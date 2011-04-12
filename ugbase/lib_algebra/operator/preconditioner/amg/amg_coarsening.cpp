/**
 * \file amg_coarsening.cpp
 *
 * \author Martin Rupp
 *
 * \date 06.08.2010
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */
///
#include <iostream>

#include "common/assert.h"
#include "common/log.h"

using namespace std;

#include "graph.h"
#include "amg_nodeinfo.h"
#include "maxheap.h"

#include "lib_algebra/common/stl_debug.h"
#include "amg_profiling.h"

//#define AMG_PRINT_COARSEN

namespace ug
{

// CreateMeasureOfImportancePQ:
//------------------------------
/**
 * \brief Creates the "measure of importance" (=rating) priority queue.
 * calculate ratings in nodes[i].rating, build up priority queue PQ.
 * rating = unassigned neighbors + 2 * fine neighbors ( neighbors in strongT)
 * \param 	strong
 * \param 	strongT 	transpose of the strong connectivity graph
 * \param	PQ			maxheap priority queue for sorting of the nodes wrt the rating
 * \param	nodes
 * \sa	CreateStrongConnectionGraph
 */
void CreateMeasureOfImportancePQ(cgraph &strong, cgraph &strongT, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(nodes.size() == strongT.size(), "");

	PQ.create(nodes.size(), &nodes[0]);
	for(size_t i=0; i < nodes.size(); i++)
	{
		if(strong.is_isolated(i))
			nodes.set_isolated(i);
		else
		{
			//UG_ASSERT(graph.iNrOfConnections[i] > 0, "node " << i << " has " << graph.iNrOfConnections[i] << " connections?");
			nodes.set_rating(i, strongT.num_connections(i));
			PQ.insert_item(i);
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
void CreateAggressiveCoarseningGraph(cgraph &graph, cgraph &graph2, AMGNodes &nodes,
		int minNrOfPaths, int *posInConnections)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(nodes.size() == graph.size(), "");
	vector<int> connection(255);
	vector<int> nrOfPaths(255);

	//graph.print();
	for(size_t i=0; i < graph.size(); i++)
	{
		if(nodes[i].is_fine_direct())
			continue;

		connection.clear();
		nrOfPaths.clear();
		// first calculate all nodes reachable with paths of length 2

		// ! i is coarse -> has only fine neighbors
		for(cgraph::const_row_iterator conn = graph.begin_row(i); conn != graph.end_row(i); ++conn)
		{
			size_t indexN = (*conn);
			for(cgraph::const_row_iterator connN = graph.begin_row(indexN); connN != graph.end_row(indexN); ++connN)
			{
				size_t indexNN = (*connN);

				if(indexNN == i || nodes[indexNN].is_fine_direct())
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
 * \param	nodes			fine/coarse marks (skip fine nodes)
 */
void CreateMeasureOfImportanceAggressiveCoarseningPQ(cgraph &graphAC, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(nodes.size() == graphAC.size(), "");
	PQ.create(nodes.size(), &nodes[0]);

	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes[i].is_fine_direct())
			continue;

		nodes.set_rating(i, graphAC.num_connections(i));
		PQ.insert_item(i);
	}

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coarsen:
//-------------------------
/**
 * \brief Coarsens the graph with ratings of nodes in nodes[i].rating, set up in a priority queue PQ
 * \param graph			Graph
 * \param PQ
 *  bIndirect		if true, this is 2nd stage of Aggressive Coarsening, then fine nodes get marker "IndirectFine"
 *						instead of just "fine". Used later in CreateProlongation and CreateIndirectProlongation
 * \param nodes
 * \return				returns number of new coarse nodes.
 */
int Coarsen(cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(graph.size() == nodes.size() && graph.size() == PQ.arr_size(), "");
	while(nodes.get_unassigned() > 0)
	{
		// get Node with best rating
		int best = PQ.remove_max();

		UG_DLOG(LIB_ALG_AMG, 3, endl << "set coarse: " << best << " rating " << nodes[best].rating  << ". then fine: ");

		UG_ASSERT(!nodes[best].is_assigned(), "node " << best << " is already assigned??? (rating = " << nodes[best].rating << ", unassigend = " << nodes.get_unassigned() << ")");

		nodes.set_coarse(best);					// mark as coarse/assigned

		// unassigned neighbors will be marked as fine, so remove from PQ
		for(cgraph::const_row_iterator conn = graph.begin_row(best); conn != graph.end_row(best); ++conn)
		{
			size_t indexN = (*conn);
			if(nodes[indexN].is_assigned()) continue;
			PQ.remove(indexN);
		}

		// now mark unassigned neighbors as fine
		for(cgraph::const_row_iterator conn = graph.begin_row(best); conn != graph.end_row(best); ++conn)
		{
			int indexN = (*conn);

			IF_DEBUG(LIB_ALG_AMG, 3)
			{
				UG_LOG(indexN << " ");
				if(nodes[indexN].is_assigned())
					UG_LOG((nodes[indexN].is_coarse() ? "(c) " : "(f) "));
			}

			if(nodes[indexN].is_assigned()) continue;

			//if(bIndirect) nodes[indexN].setFineIndirect();
			//else

			nodes.set_fine_direct(indexN);

			// increase rating of neighbors of this node (= neighbors of neighbors of node "best")
			// rating = unassigned neighbors + 2 * fine neighbors
			for(cgraph::const_row_iterator connN = graph.begin_row(indexN); connN != graph.end_row(indexN); ++connN)
			{
				int indexNN = (*connN);
				// TODO: perhaps we could create a the f-f candidate list here

				if(nodes[indexNN].is_assigned())
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

	UG_ASSERT(nodes.get_nr_of_coarse() > 0, "no coarse nodes???");
	UG_DLOG(LIB_ALG_AMG, 1, nodes.get_nr_of_coarse() << " Coarse Nodes. ");
	return nodes.get_nr_of_coarse();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PreventFFConnections:
//-------------------------
/**
 * \brief When two fine nodes are strong connected, they need to have a common interpolating
 * node. Otherwise, we set
 * \param graphS
 * \param graphST
 * \param nodes
 * \return				returns number of new coarse nodes.
 */
void PreventFFConnections(cgraph &graphS, cgraph &graphST, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	size_t N = graphS.size();
	UG_ASSERT(N == graphST.size() && N == nodes.size(), "");
	// second pass
	//----------------
	// seems to work but doesnt help with convergence rates on complicated geometries???

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(N, false);

	for(size_t i=0; i< N; i++)
	{
		if(nodes[i].is_coarse() || graphS.num_connections(i)==0)
			continue;

		// mark coarse nodes interpolating this fine node
		for(cgraph::const_row_iterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
		{
			if(nodes[(*it)].is_coarse())
				marks[(*it)] = true;
		}

		// prevent strong F-F connections without common Interpolation node
		for(cgraph::const_row_iterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
		{
			if(nodes[*it].is_coarse() || graphS.is_isolated(*it))
				continue;

			cgraph::const_row_iterator it2 = graphST.begin_row(*it);
			cgraph::const_row_iterator it2end = graphST.end_row(*it);
			for(; it2 != it2end; ++it2)
			{
				if(nodes[*it2].is_coarse() && marks[*it2])
					break;
			}

			if(it2 == it2end)
			{
				// TODO: calculate ratings, add to PQ and do coarsening again on those candidates
				// rating of a fine node = nr of f-f pairs which this node is adjacent to BOTH
				// that is: all common fine node neighbors of i and it() get rating++.
				// problem: updating

				//cout << endl << "prevent F-F-connection between (2) " << i << "[" << GetOriginalIndex(A.tolevel, i) << "] and " << it() << "[" << GetOriginalIndex(A.tolevel, it()) << "], setting " << i << " coarse.";
				nodes.set_coarse(i);
				nrOfFFCoarseNodes++;
				break;
			}
		}

		// remove marks
		for(cgraph::const_row_iterator conn = graphS.begin_row(i); conn != graphS.end_row(i); ++conn)
		{
			size_t index = (*conn);
			if(nodes[index].is_coarse())
				marks[index] = false;
		}
	}

	if(nrOfFFCoarseNodes)
		UG_DLOG(LIB_ALG_AMG, 1, "F-F prevention, now " << nodes.get_nr_of_coarse() << " coarse nodes.\n");
}

} // namespace ug
