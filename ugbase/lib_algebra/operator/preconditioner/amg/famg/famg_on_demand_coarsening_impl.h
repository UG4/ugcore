/**
 * \file famg_on_demand_coarsening_impl.h
 *
 * \author Martin Rupp
 *
 * \date 10.12.2010
 *
 * implementation file for famg, coarsening when calculating the possible parents only on demand
 *
 * Goethe-Center for Scientific Computing 2010-2011.
 *
 * for test purposes, functions are here in a cpp file.
 */


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_ON_DEMAND_COARSENING_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_ON_DEMAND_COARSENING_IMPL_H__

namespace ug{

template<typename neighborstruct, typename matrix_type, typename vector_type>
bool OnDemand_UpdateRating(size_t node, stdvector<neighborstruct> &PN, FAMGNodes &nodes,
		stdvector<bool> &prolongation_calculated, cgraph &SymmNeighGraph,
		FAMGInterpolationCalculator<matrix_type, vector_type> &calculator)
{
	AMG_PROFILE_FUNC();
	if(prolongation_calculated[node])
		return nodes.update_rating(node, PN);
	else
	{
#if 0
		int coarseNeighbors=0;
		for(cgraph::const_row_iterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
		{
			if(rating[(*conn)].is_coarse())
				coarseNeighbors++;
		}
		int r = min(0, 2-coarseNeighbors);
		if(r == rating[node])
			return false;

		rating[node] = r;
		if(r == 2)
			return true;
#else
		int hasCoarseNeighbors=false;
		for(cgraph::const_row_iterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
		{
			if(nodes[(*conn)].is_coarse())
			{
				hasCoarseNeighbors = true;
				break;
			}
		}
		if(!hasCoarseNeighbors)
			return false;

#endif
		UG_DLOG(LIB_ALG_AMG, 2, "node " << node << " has coarse neighbors, but rating is not calculated yet. calculating rating now\n");
		calculator.get_possible_parent_pairs(node, PN, nodes);
		prolongation_calculated[node] = true;
		nodes.update_rating(node, PN);
		return true;


	}
}

template<typename neighborstruct, typename matrix_type, typename vector_type>
void OnDemand_Update(size_t node, stdvector<stdvector<neighborstruct> > &possible_parents, FAMGNodes &nodes,
		maxheap<FAMGNode> &heap,
		stdvector<bool> &prolongation_calculated,	cgraph &SymmNeighGraph, FAMGInterpolationCalculator<matrix_type, vector_type> &calculator)
{
	AMG_PROFILE_FUNC();
	if(!nodes[node].is_valid_rating())
		return;
	if(OnDemand_UpdateRating(node, possible_parents[node], nodes, prolongation_calculated, SymmNeighGraph, calculator))
	{
		if(nodes[node].is_uninterpolateable() || nodes[node].is_fine())
		{
			if(heap.is_in(node))
			{
				UG_DLOG(LIB_ALG_AMG, 2, " remove node " << node << " from heap!\n");
				heap.remove(node);
			}
			else
				UG_DLOG(LIB_ALG_AMG, 2, " node " << node << " uninterpolateable, but not in heap anyway!\n");
		}
		else
		{

			if(heap.is_in(node))
			{
				heap.update(node);
				UG_DLOG(LIB_ALG_AMG, 2, " updated node " << node << " in heap!\n");
			}
			else
			{
				heap.insert_item(node);
				UG_DLOG(LIB_ALG_AMG, 2, " inserted node " << node << " into heap!\n");
			}
		}
	}
}

void AddUnmarkedNeighbors(cgraph &SymmNeighGraph, size_t i, stdvector<bool> &mark, stdvector<size_t> &list)
{
	for(cgraph::const_row_iterator conn = SymmNeighGraph.begin_row(i); conn != SymmNeighGraph.end_row(i); ++conn)
		if(mark[(*conn)] == false) list.push_back((*conn));
}



template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::on_demand_coarsening()
{
	AMG_PROFILE_FUNC();
	Stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "other coarsening... "); if(bTiming) SW.start();

	AMG_PROFILE_BEGIN(AMG_on_demand_Init)
	size_t N = rating.size();

	possible_parents.clear();
	possible_parents.resize(A_OL2.num_rows());

	// handle dirichlet boundaries
	for(size_t j=0; j<N; j++)
	{
		if(rating[j].is_valid_rating() &&
				rating.i_must_assign(j) && A_OL2.is_isolated(j))
			rating.set_fine(j);
	}

	prolongation_calculated.resize(N, false);

	// try to find a node which is "inside"
	size_t i;
	for(i=0; i<N; i++)
	{
		if(!rating[i].is_valid_rating() || !rating.i_must_assign(i))
			continue;

		if(IsCloseToBoundary(A_OL2, i , 2) == false)
		{
			calculator.get_possible_parent_pairs(i, possible_parents[i], rating);
			prolongation_calculated[i] = true;
			rating.update_rating(i, possible_parents[i]);
			if(!rating[i].is_uninterpolateable() && !rating[i].is_fine())
				break;
		}
	}

	UG_ASSERT(i != N, "did not found suitable node");
	UG_DLOG(LIB_ALG_AMG, 2, "\nStarting with node " << rating.get_original_index(i) << "\n");



	// set "visited" flags, so we do not update twice
	stdvector<bool> bvisited;
	bvisited.resize(N, false);
	stdvector<size_t> neighborsToUpdate;

	for(size_t j=0; j<N; j++)
	{
		if(!rating[i].is_valid_rating() || !rating.i_must_assign(j))
			bvisited[j] = true;
	}



	AMG_PROFILE_NEXT(AMG_on_demand_loop)

	while(1)
	{
		UG_ASSERT(possible_parents[i].size() > 0, "node " << i << "has no possible parents?");
		neighborstruct2 &n = possible_parents[i][0];

		UG_DLOG(LIB_ALG_AMG, 2, "\n\n\nSelect next node...\n");
		UG_DLOG(LIB_ALG_AMG, 2, "======================================\n");
		UG_DLOG(LIB_ALG_AMG, 2, "node " << rating.get_original_index(i) << " has rating " << rating[i] << ". now gets fine. parents: ");
		for(size_t j=0; j < n.parents.size(); j++)
			UG_DLOG(LIB_ALG_AMG, 2, rating.get_original_index(n.parents[j].from) << " ");
		UG_DLOG(LIB_ALG_AMG, 2,"\n");

		// node i gets fine. update neighbors.
		rating.set_fine(i);
		bvisited[i] = true;

		for(size_t j=0; j < n.parents.size(); j++)
			bvisited[i] = true;

		AddUnmarkedNeighbors(SymmNeighGraph, i, bvisited, neighborsToUpdate);

		if(fvalues.size() > i) fvalues[i] = n.F;

		UG_DLOG(LIB_ALG_AMG, 2, "Set coarse parents:\n");
		// get parent pair, set as coarse (if not already done), update neighbors.

		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;

			UG_ASSERT(!rating[node].is_fine(), "parent node " << rating.get_original_index(node) << " of node " << rating.get_original_index(i) << " is FINE?");

			if(rating[node].is_coarse()) { UG_DLOG(LIB_ALG_AMG, 2, "\nnode " << rating.get_original_index(node) << " is already coarse\n"); }
			else { UG_DLOG(LIB_ALG_AMG, 2, "\nnode " << rating.get_original_index(node) << " has rating " << rating[node] << ". now gets coarse.\nUpdate Neighbors of " << node << "\n"); }

			if(!rating[node].is_coarse())
			{
				if(heap.is_in(node)) heap.remove(node);
				rating.set_coarse(node);
				AddUnmarkedNeighbors(SymmNeighGraph, node, bvisited, neighborsToUpdate);
			}
			PoldIndices(i, node) = n.parents[j].value;
		}

		IF_DEBUG(LIB_ALG_AMG,4) print_vector(neighborsToUpdate, "neighborsToUpdate");

		// update neighbors
		for(size_t j=0; j<neighborsToUpdate.size(); j++)
		{
			size_t index = neighborsToUpdate[j];
			OnDemand_Update(index, possible_parents, rating, heap, prolongation_calculated, SymmNeighGraph, calculator);
			if(rating[index].is_valid_rating())
				bvisited[index] = false;
		}
		neighborsToUpdate.clear();


		UG_DLOG(LIB_ALG_AMG, 2, "\n\nSearching for next node...\n");
		AMG_PROFILE_BEGIN(AMG_Search_for_next_node)
		if(heap.height() == 0)
		{
			// heap empty, we need to get another start node
			for(i=0; i<N; i++)
			{
				if(!rating[i].is_valid_rating() || rating.i_must_assign(i) == false) continue;
				UG_DLOG(LIB_ALG_AMG, 2, "rating " << i << " is " << rating[i] << "\n");
				if(A.is_isolated(i) == false && rating[i].is_valid_rating())
				{
					if(prolongation_calculated[i] == false)
					{
						//UG_LOG(" calculating prologation...");
						calculator.get_possible_parent_pairs(i, possible_parents[i], rating);
						prolongation_calculated[i] = true;
					}
					rating.update_rating(i, possible_parents[i]);

					UG_DLOG(LIB_ALG_AMG, 2, " new rating is " << rating[i] << "\n");
					if(rating[i].is_valid_rating()) break;
				}
			}
			if(i == N) { UG_DLOG(LIB_ALG_AMG, 2, "\nno more new start nodes found\n"); break; }

			UG_DLOG(LIB_ALG_AMG, 2, "\n\nRESTARTING WITH NODE " << rating.get_original_index(i) << "!!!\n\n");

			calculator.get_possible_parent_pairs(i, possible_parents[i], rating);
			prolongation_calculated[i] = true;
			rating.update_rating(i, possible_parents[i]);
		}
		else
			while(1)
			{
				do
				{
					i = heap.get_max();
					// update rating since because we made some nodes fine that rating could be wrong
					if(rating.update_rating(i, possible_parents[i]) == false)
						break;
					if(rating[i].is_uninterpolateable() || rating[i].is_fine())
						heap.remove(i);
					else
						heap.update(i);
				} while(true);
				UG_DLOG(LIB_ALG_AMG, 2, "node " << i << " has best rating");
				UG_ASSERT(rating.i_must_assign(i), "node " << i << " is not master!");

				if(prolongation_calculated[i] == false)
				{
					UG_DLOG(LIB_ALG_AMG, 2, ", but prolongation not calculated. update.\n");
					calculator.get_possible_parent_pairs(i, possible_parents[i], rating);
					prolongation_calculated[i] = true;
					rating.update_rating(i, possible_parents[i]);
					heap.update(i);
				}
				else
				{
					UG_DLOG(LIB_ALG_AMG, 2, ", and prolongation calculated, take this node.\n");
					i = heap.remove_max();
					break;
				}
			}
	}

	for(size_t i=0; i<N; i++)
	{
		UG_ASSERT(rating.i_must_assign(i) == false || rating[i].rating != 0.0, i);
	}

	UG_DLOG(LIB_ALG_AMG, 2, "\nDone!\n");

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
}

}

#endif
