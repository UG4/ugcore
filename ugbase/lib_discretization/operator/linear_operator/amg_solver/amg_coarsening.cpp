/**
 * \file amg_coarsening.cpp
 *
 * \author Martin Rupp
 *
 * \date 06.08.2010
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#include <iostream>

#include "ug.h"
#include "graph.h"
#include "amg_nodeinfo.h"
#include "maxheap.h"

using namespace std;

namespace ug
{

// CreateMeasureOfImportancePQ:
//------------------------------
/**
 * \brief Creates the "measure of importance" (=rating) priority queue.
 * calculate ratings in nodes[i].rating, build up priority queue PQ.
 * rating = unassigned neighbors + 2 * fine neighbors ( neighbors in strongT)
 * \param 	strongT 	transpose of the strong connectivity graph
 * \param	PQ			maxheap priority queue for sorting of the nodes wrt the rating
 * \param 	unassigned	nr of nodes which are now to be assigned coarse or fine
 * \param	nodes
 * \sa	CreateStrongConnectionGraph
 */
void CreateMeasureOfImportancePQ(cgraph &strong, cgraph &strongT, maxheap<amg_nodeinfo> &PQ, int &unassigned, amg_nodeinfo *nodes)
{
	PQ.create(strongT.size(), nodes);
	unassigned = 0;
	for(size_t i=0; i < strongT.size(); i++)
	{
		if(strong.is_isolated(i))
			nodes[i].setFineDirect();
		else
		{
			//UG_ASSERT(graph.iNrOfConnections[i] > 0, "node " << i << " has " << graph.iNrOfConnections[i] << " connections?");
			nodes[i].rating = strongT.num_connections(i);
			PQ.insert_item(i);
			unassigned++;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateAggressiveCoarseningGraph
//----------------------------------
/**
 * Creates the graph for aggressive coarsening from the strong connectivity graph
 * only among nodes which are marked coarse.
 * two coarse nodes a and b in graph2 are connected, if
 * - there exist at least 2 ways of length 2 from a to b (A2-Coarsening)
 * - there exist a ways of length 2 from a to b (A1-Coarsening)
 * because of the coarsening process, every coarse node has only fine neighbors in graph,
 * that means those ways are from a coarse node over a fine node to a coarse node.
 *
 * \param	graph				old graph of normal strong connectivity (from CreateGraph)
 * \param 	graph2				new graph of "distance-2-strong-connectivity"
 * \param	nodes				info about nodes if nodes are fine or coarse
 * \param	minNrOfPaths		2 or 1, depending on A2- or A1-Coarsening.
 * \param 	posInConnections	array of size graph.size for speedup of neighbor-neighbor-calculation inited with -1.
 */
void CreateAggressiveCoarseningGraph(cgraph &graph, cgraph &graph2, amg_nodeinfo *nodes,
		int minNrOfPaths, int *posInConnections)
{
	vector<int> connection(255);
	vector<int> nrOfPaths(255);

	//graph.print();
	for(size_t i=0; i < graph.size(); i++)
	{
		if(nodes[i].isFineDirect())
			continue;

		connection.clear();
		nrOfPaths.clear();
		// first calculate all nodes reachable with paths of length 2

		// ! i is coarse -> has only fine neighbors
		for(cgraph::cRowIterator conn = graph.begin_row(i); conn != graph.end_row(i); ++conn)
		{
			size_t indexN = (*conn);
			for(cgraph::cRowIterator connN = graph.begin_row(indexN); connN != graph.end_row(indexN); ++connN)
			{
				size_t indexNN = (*connN);

				if(indexNN == i || nodes[indexNN].isFineDirect())
					continue;
				int pos = posInConnections[indexNN];
				if(pos == -1)
				{
					// never reached node indexNN from i, init.
					pos = posInConnections[indexNN]= connection.size();
					connection.push_back(indexNN);
					nrOfPaths.push_back(1);
				}
				else
					nrOfPaths[pos]++;
			}
		}

		// then sort out those which were reached #aggressiveCoarseningNrOfPaths (2 or 1) times
		for(size_t j=0; j<connection.size(); j++)
		{
			if(nrOfPaths[j] >= minNrOfPaths)
			{
				// add connection i -> node
				graph2.set_connection(i, connection[j]);
				// increase rating of i
			}
			// reset posInConnections for further use
			posInConnections[connection[j]] = -1;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CreateMeasureOfImportanceAggressiveCoarseningPQ
//----------------------------------
/**
 * this function could be merged with CreateMeasureOfImportancePQ
 * \param	graphAC			graph created by CreateAggressiveCoarseningGraph
 * \param 	PQ				maxheap priority queue for sorting of the nodes wrt the rating
 * \param 	unassigned		nr of nodes which are now to be assigned coarse or fine
 * \param	iNrOfCoarse		since we need to assign coarse nodes, we need this number
 * \param	newIndex		newIndex of new coarse nodes
 * \param	nodes			fine/coarse marks (skip fine nodes)
 */
void CreateMeasureOfImportanceAggressiveCoarseningPQ(cgraph &graphAC, maxheap<amg_nodeinfo> &PQ, int &unassigned,
		int &iNrOfCoarse, int *newIndex, amg_nodeinfo *nodes)
{
	PQ.create(graphAC.size(), nodes);

	iNrOfCoarse = 0;
	for(size_t i=0; i<graphAC.size(); i++)
	{
		if(nodes[i].isFineDirect())
			continue;

		nodes[i].rating = graphAC.num_connections(i);
		PQ.insert_item(i);
		unassigned++;
	}

	//cout << endl << endl;
	//graph2.print();
}

#define AMG_PRINT_COARSEN

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coarsen:
//-------------------------
/**
 * \brief Coarsens the graph with ratings of nodes in nodes[i].rating, set up in a priority queue PQ
 * \param newIndex		store in newIndex[i] new index of node i in coarse nodes (>0, if fine < 0)
 * \param unassigned	nr of nodes to assign
 * \param bIndirect		if true, this is 2nd stage of Aggressive Coarsening, then fine nodes get marker "IndirectFine"
 *						instead of just "fine". Used later in CreateProlongation and CreateIndirectProlongation
 * \param A				matrix A (for debug)
 * \return				returns number of new coarse nodes.
 */
int Coarsen(cgraph &graph, maxheap<amg_nodeinfo> &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, amg_nodeinfo *nodes)
{
	while(unassigned > 0)
	{
		// get Node with best rating
		int best = PQ.remove_max();

#ifdef AMG_PRINT_COARSEN
		cout << endl << "set coarse: " << best << " rating " << nodes[best].rating  << ". then fine: ";
#endif

		UG_ASSERT(!nodes[best].isAssigned(), "node " << best << " is already assigned??? (rating = " << nodes[best].rating << ", unassigend = " << unassigned << ")");

		nodes[best].setCoarse();				// mark as coarse/assigned
		newIndex[best] = iNrOfCoarse++;		// new coarse index
		unassigned--;

		// unassigned neighbors will be marked as fine, so remove from PQ
		for(cgraph::cRowIterator conn = graph.begin_row(best); conn != graph.end_row(best); ++conn)
		{
			int indexN = (*conn);
			if(nodes[indexN].isAssigned()) continue;
			PQ.remove(indexN);
		}

		// now mark unassigned neighbors as fine
		for(cgraph::cRowIterator conn = graph.begin_row(best); conn != graph.end_row(best); ++conn)
		{
			int indexN = (*conn);

#ifdef AMG_PRINT_COARSEN
			cout << indexN << " ";
			if(nodes[indexN].isAssigned())
				cout << (nodes[indexN].isCoarse() ? "(c) " : "(f) ");
#endif

			if(nodes[indexN].isAssigned()) continue;

			//if(bIndirect) nodes[indexN].setFineIndirect();
			//else

			nodes[indexN].setFineDirect();
			unassigned--;

			// increase rating of neighbors of this node (= neighbors of neighbors of node "best")
			// rating = unassigned neighbors + 2 * fine neighbors
			for(cgraph::cRowIterator connN = graph.begin_row(indexN); connN != graph.end_row(indexN); ++connN)
			{
				int indexNN = (*connN);
				// TODO: perhaps we could create a the f-f candidate list here

				if(nodes[indexNN].isAssigned())
					continue;
				nodes[indexNN].rating++;
				PQ.update(indexNN);
			}
		}
		//coarse.print();
		//cout << "Ranking: " << endl;
		//PQ.print();
	}
	//cout << endl;

	UG_ASSERT(iNrOfCoarse > 0, "no coarse nodes???");

	cout << iNrOfCoarse << " Coarse Nodes. ";
	return iNrOfCoarse;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PreventFFConnections:
//-------------------------
/**
 * \brief When two fine nodes are strong connected, they need to have a common interpolating
 * node. Otherwise, we set
 * \param newIndex		store in newIndex[i] new index of node i in coarse nodes (>0, if fine < 0)
 * \param unassigned	nr of nodes to assign

 * \return				returns number of new coarse nodes.
 */
void PreventFFConnections(cgraph &graphS, cgraph &graphST, amg_nodeinfo *nodes, int *newIndex, int &iNrOfCoarse)
{
	// second pass
	//----------------
	// seems to work but doesnt help with convergence rates on complicated geometries???

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(graphS.size());

	for(size_t i=0; i< graphS.size(); i++)
	{
		if(nodes[i].isCoarse() || graphS.num_connections(i)==0)
			continue;

		// mark coarse nodes interpolating this fine node
		for(cgraph::cRowIterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
		{
			if(nodes[(*it)].isCoarse())
				marks[(*it)] = true;
		}

		// prevent strong F-F connections without common Interpolation node
		for(cgraph::cRowIterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
		{
			if(nodes[*it].isCoarse() || graphS.is_isolated(*it))
				continue;

			cgraph::cRowIterator it2 = graphST.begin_row(*it);
			cgraph::cRowIterator it2end = graphST.end_row(*it);
			for(; it2 != it2end; ++it2)
			{
				if(nodes[*it2].isCoarse() && marks[*it2])
					break;
			}

			if(it2 == it2end)
			{
				// TODO: calculate ratings, add to PQ and do coarsening again on those candidates
				// rating of a fine node = nr of f-f pairs which this node is adjacent to BOTH
				// that is: all common fine node neighbors of i and it() get rating++.
				// problem: updating

				//cout << endl << "prevent F-F-connection between (2) " << i << "[" << GetOriginalIndex(A.tolevel, i) << "] and " << it() << "[" << GetOriginalIndex(A.tolevel, it()) << "], setting " << i << " coarse.";
				nodes[i].setCoarse();
				newIndex[i] = iNrOfCoarse++;

				nrOfFFCoarseNodes++;
				break;
			}
		}

		// remove marks
		for(cgraph::cRowIterator conn = graphS.begin_row(i); conn != graphS.end_row(i); ++conn)
		{
			size_t index = (*conn);
			if(nodes[index].isCoarse())
				marks[index] = false;
		}
	}

	if(nrOfFFCoarseNodes)
		cout << "F-F prevention, now " << iNrOfCoarse << " coarse nodes." << endl;
}

} // namespace ug
