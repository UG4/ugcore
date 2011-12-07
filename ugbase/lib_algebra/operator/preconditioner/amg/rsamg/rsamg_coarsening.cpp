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

#include "../graph.h"
#include "rsamg_nodeinfo.h"
#include "../maxheap.h"

#include "lib_algebra/common/stl_debug.h"
#include "../amg_profiling.h"
#include "../boxsort.h"
//#define AMG_PRINT_COARSEN

namespace ug
{

// CreateMeasureOfImportancePQ:
//------------------------------
/**
 * \brief Creates the "measure of importance" (=rating) priority queue.
 * calculate ratings in nodes[i].rating, build up priority queue PQ.
 * rating = unassigned neighbors + 2 * fine neighbors ( neighbors in strongT)
 * \param 	strong		strong[i] = nodes which influence i
 * \param 	strongT 	strongT[i] = nodes which are influenced by i.
 * \param	PQ			maxheap priority queue for sorting of the nodes wrt the rating
 * \param	nodes
 * \sa	CreateStrongConnectionGraph
 */
void CreateMeasureOfImportancePQ(const cgraph &strong, const cgraph &strongT, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(nodes.size() == strongT.size(), "");

	PQ.create(nodes.size(), &nodes[0]);
	for(size_t i=0; i < nodes.size(); i++)
	{
		if(nodes.is_master_or_inner(i) == false)
		{
			nodes.set_parallel_dont_care(i);
			continue;
		}
		if(!strong.is_isolated(i))
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
void CreateAggressiveCoarseningGraph(const cgraph &graph, cgraph &graph2, const AMGNodes &nodes,
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
void CreateMeasureOfImportanceAggressiveCoarseningPQ(const cgraph &graphAC, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(nodes.size() == graphAC.size(), "");
	PQ.create(nodes.size(), &nodes[0]);

	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes.is_slave(i)) continue;
		if(nodes[i].is_fine_direct() || nodes[i].is_dirichlet())
			continue;

		nodes.set_rating(i, graphAC.num_connections(i));
		PQ.insert_item(i);
	}

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RemoveUnassignedNeighbors:
//-------------------------
void RemoveUnassignedNeighbors(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i)
{
	for(cgraph::const_row_iterator conn = graph.begin_row(i); conn != graph.end_row(i); ++conn)
	{
		size_t indexN = (*conn);
		if(nodes[indexN].is_assigned()) continue;
		PQ.remove(indexN);
	}
}

void ChangeRatingOfUnassignedNeighbors(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i, int change)
{
	for(cgraph::const_row_iterator conn = graph.begin_row(i); conn != graph.end_row(i); ++conn)
	{
		int indexN = (*conn);
		// TODO: perhaps we could create a the f-f candidate list here

		if(nodes[indexN].is_assigned())
			continue;
		nodes[indexN].rating += change;
		PQ.update(indexN);
	}

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MarkUnassignedNeighborsFine:
//-------------------------
void MarkUnassignedNeighborsFine(const cgraph &graph, nodeinfo_pq_type &PQ, AMGNodes &nodes, size_t i, bool bMarkAsFineIndirect)
{
	for(cgraph::const_row_iterator conn = graph.begin_row(i); conn != graph.end_row(i); ++conn)
	{
		int indexN = (*conn);

		IF_DEBUG(LIB_ALG_AMG, 3)
		{
			UG_LOG(indexN << " ");
			if(nodes[indexN].is_assigned())
				UG_LOG((nodes[indexN].is_coarse() ? "(c) " : "(f) "));
		}

		if(nodes[indexN].is_parallel_dont_care() == false && nodes[indexN].is_assigned()) continue;

		//if(bIndirect) nodes[indexN].setFineIndirect();
		//else

		if(bMarkAsFineIndirect)
			nodes.set_unassigned_fine_indirect(indexN);
		else
			nodes.set_fine_direct(indexN);

		ChangeRatingOfUnassignedNeighbors(graph, PQ, nodes, indexN, +1);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coarsen:
//-------------------------
/**
 * \brief Coarsens the graph with ratings of nodes in nodes[i].rating, set up in a priority queue PQ
 * \param S							S[i] nodes influencing i
 * \param ST						ST[i] = nodes beeing influenced by i
 * \param PQ
 * \param nodes
 * \param bUnsymmetric				set to true if S != ST
 * \param bAggressiveCoarsening
 * \return number of new coarse nodes.
 */
int Coarsen(const cgraph &S, const cgraph &ST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric, bool bAggressiveCoarsening)
{
	AMG_PROFILE_FUNC();
	UG_ASSERT(S.size() == nodes.size() && S.size() == PQ.arr_size() && S.size() == ST.size(), "");

	//UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 3);

#ifdef UG_DEBUG
	size_t n=0;
	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes[i].is_assigned() == false)
			UG_ASSERT(PQ.is_in(i), i);
		if(PQ.is_in(i))
			UG_ASSERT(nodes[i].is_assigned() == false, i);
		if(nodes[i].is_assigned() == false)
			n++;
	}
	UG_ASSERT(n == nodes.get_unassigned(), "");
#endif

	while(PQ.height() > 0)
	{
		UG_ASSERT(PQ.height() == nodes.get_unassigned(), PQ.height() << " != " << nodes.get_unassigned())
		// get Node with best rating
		int best = PQ.remove_max();

		UG_DLOG(LIB_ALG_AMG, 3, endl << "set coarse: " << best << " rating " << nodes[best].rating  << ". then fine: ");
		UG_ASSERT(!nodes[best].is_assigned(), "node " << best << " is already assigned??? (rating = " << nodes[best].rating << ", unassigend = " << nodes.get_unassigned() << ")");

		nodes.set_coarse(best);					// mark as coarse/assigned

		// unassigned neighbors will be marked as fine, so remove from PQ
		RemoveUnassignedNeighbors(ST, PQ, nodes, best);

		// now mark unassigned neighbors as fine
		MarkUnassignedNeighborsFine(ST, PQ, nodes, best, bAggressiveCoarsening);

		if(bUnsymmetric)
			ChangeRatingOfUnassignedNeighbors(S, PQ, nodes, best, -1);
	}
	//cout << endl;

	UG_ASSERT(nodes.get_nr_of_coarse() > 0, "no coarse nodes???");
	UG_DLOG(LIB_ALG_AMG, 1, nodes.get_nr_of_coarse() << " Coarse Nodes. ");
	return nodes.get_nr_of_coarse();
}


/**
 * \brief
 * \param graphS
 * \param graphST
 * \param nodes
 * \param i
 */
void PreventFFConnection(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes, size_t i, vector<bool> &marks, int &nrOfFFCoarseNodes)
{
	if(nodes[i].is_coarse() || graphS.num_connections(i)==0) return;

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
			nodes.set_coarse(i);
			nrOfFFCoarseNodes++;
			break;
		}
	}

	// remove marks
	for(cgraph::const_row_iterator conn = graphST.begin_row(i); conn != graphST.end_row(i); ++conn)
	{
		size_t index = (*conn);
		if(nodes[index].is_coarse())
			marks[index] = false;
	}
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
void PreventFFConnections(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	size_t N = graphS.size();
	UG_ASSERT(N == graphST.size() && N == nodes.size(), "");

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(N, false);

	for(size_t i=0; i< N; i++)
	{
		if(nodes.is_inner(i) && nodes[i].is_fine())
			PreventFFConnection(graphS, graphST, nodes, i, marks, nrOfFFCoarseNodes);
	}

	if(nrOfFFCoarseNodes)
		UG_DLOG(LIB_ALG_AMG, 0, "F-F prevention, now " << nodes.get_nr_of_coarse() << " coarse nodes.\n");
}

/*

template<typename matrix_type>
double GetRowMaxAbs(const matrix_type &A, size_t row)
{
	double m=0;
	for(typename matrix_type::const_row_iterator it = A.begin_row(row); it != A.end_row(row); ++it)
	{
		double a = BlockAbs(*it);
		if(a > m) m = a;
	}
	return m;
}

void PreventFFConnectionsMetsch(const matrix_type &mat, const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	size_t N = graphS.size();
	UG_ASSERT(N == graphST.size() && N == nodes.size(), "");

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(N, false);
	vector<double> maxk(N, -1.0);

	for(size_t r=0; i< N; r++)
		maxk[i] = GetRowMaxAbs(mat);


	for(size_t r=0; i< N; r++)
	{
		if(nodes.is_fine(r)==false)
			continue;

		for(cgraph::const_row_iterator it = graphS.begin_row(r); it != graphS.end_row(r); ++it)
		{
			int c=*it;
			if(!nodes[c].is_fine()) continue;



		}

	}

	if(nrOfFFCoarseNodes)
		UG_DLOG(LIB_ALG_AMG, 1, "F-F prevention, now " << nodes.get_nr_of_coarse() << " coarse nodes.\n");
}


*/

/**
 * \brief
 * \param graphS
 * \param graphST
 * \param nodes
 * \param i
 * \param G
 * \param ratings
 */
void GetFFConnectionGraph(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes, size_t i, vector<bool> &marks,
		cgraph &G, stdvector<size_t> &ratings)
{
	if(nodes[i].is_coarse() || graphS.num_connections(i)==0) return;

	// mark coarse nodes interpolating this fine node
	for(cgraph::const_row_iterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
	{
		if(nodes[(*it)].is_coarse())
			marks[(*it)] = true;
	}

	// prevent strong F-F connections without common Interpolation node
	for(cgraph::const_row_iterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
	{
		if(nodes.needs_assignment(*it)==false ) continue;
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
			G.set_connection(i, *it);
			G.set_connection(*it, i);
			ratings[i]++;
			ratings[*it]++;
			break;
		}
	}

	// remove marks
	for(cgraph::const_row_iterator conn = graphST.begin_row(i); conn != graphST.end_row(i); ++conn)
	{
		size_t index = (*conn);
		if(nodes[index].is_coarse())
			marks[index] = false;
	}
}

/**
 * \brief
 * \param graphS
 * \param graphST
 * \param nodes
 * \return				returns number of new coarse nodes.
 */
void PreventFFConnectionsMinimal(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	size_t N = graphS.size();
	UG_ASSERT(N == graphST.size() && N == nodes.size(), "");

	int nrOfFFCoarseNodes=0;
	vector<bool> marks(N, false);

	cgraph G(graphST.size());
	stdvector<size_t> ratings(graphST.size(), 0);


	for(size_t i=0; i< N; i++)
	{
		if(nodes.is_inner(i)==false)
			continue;
		GetFFConnectionGraph(graphS, graphST, nodes, i, marks, G, ratings);
	}

	BoxPriorityQueue<size_t> pq(ratings);

	for(size_t i=0; i<ratings.size(); i++)
		if(ratings[i]>0)
			pq.insert_item(i);

	while(pq.height() > 0)
	{
		int best = pq.remove_max();
		nodes.set_coarse(best);
		nrOfFFCoarseNodes++;
		for(cgraph::row_iterator it = G.begin_row(best); it != G.end_row(best); ++it)
		{
			ratings[*it]--;
			if(ratings[*it]==0)	pq.remove(*it);
			else				pq.update(*it);
		}
	}
	if(nrOfFFCoarseNodes)
		UG_DLOG(LIB_ALG_AMG, 1, "F-F prevention, now " << nodes.get_nr_of_coarse() << " coarse nodes.\n");
}




} // namespace ug

