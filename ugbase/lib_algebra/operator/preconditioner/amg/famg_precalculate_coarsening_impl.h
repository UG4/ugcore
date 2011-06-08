/**
 * \file famg_precalculate_coarsening_impl.h
 *
 * \author Martin Rupp
 *
 * \date 14.02.2011
 *
 * implementation file for famg precalculated coarsening
 *
 * Goethe-Center for Scientific Computing 2011.
 *
 */


#ifndef __H__LIB_ALGEBRA__AMG__FAMG_PRECALCULATE_COARSENING_IMPL_H__
#define __H__LIB_ALGEBRA__AMG__FAMG_PRECALCULATE_COARSENING_IMPL_H__

namespace ug
{

template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::calculate_all_possible_parent_pairs()
{
	stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 1, "\ncreating possible parent list... "); if(bTiming) SW.start();
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, iDebugLevelCalculateParentPairs);
	AMG_PROFILE_FUNC();

	possible_parents.clear();
	possible_parents.resize(A.num_rows());
	prolongation_calculated.resize(A.num_rows(), false);

	UG_LOG(std::scientific);

	for(size_t i=0; i<A.num_rows(); i++)
	{
		if(rating.i_must_assign(i) == false) continue;
		if(rating[i].is_valid_rating() == false) continue; // could be set coarse from others

		calculator.get_possible_parent_pairs(i, possible_parents[i], rating);
		prolongation_calculated[i] = true;
	}


	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms");
}



template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::precalculate_coarsening()
{
	AMG_PROFILE_FUNC();
	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, iDebugLevelPrecalculateCoarsening);
	stopwatch SW;
	UG_DLOG(LIB_ALG_AMG, 1, std::endl << "coarsening... "); if(bTiming) SW.start();

	size_t N = rating.size();
	for(size_t j=0; j<N; j++)
	{
		if(rating[j].is_valid_rating() && rating.i_must_assign(j) && A_OL2.is_isolated(j))
			rating.set_fine(j);
	}

	while(heap.height() != 0)
	{
		// get node i with best rating
		size_t i = heap.remove_max();
		neighborstruct2 &n = possible_parents[i][0];

		UG_DLOG(LIB_ALG_AMG, 4, "\n\n\nSelect next node...\n");
		UG_DLOG(LIB_ALG_AMG, 4, "node " << i << " has rating " << rating[i] << ". now gets fine. parents: ");
		for(size_t j=0; j < n.parents.size(); j++)	UG_DLOG(LIB_ALG_AMG, 2, n.parents[j].from << " ");
		UG_DLOG(LIB_ALG_AMG, 4, "\nUpdate Neighbors of " << i << "\n");

		// node i then gets fine, parent nodes get coarse, and neighbors of those updated.

		// node i gets fine. update neighbors.
		rating.set_fine(i);
		UpdateNeighbors(SymmNeighGraph, i, possible_parents, rating, heap);

		UG_DLOG(LIB_ALG_AMG, 4, "Set coarse parents:\n");
		// get parent pair, set as coarse (if not already done), update neighbors.

		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;

			if(rating[node].is_coarse()) { UG_DLOG(LIB_ALG_AMG, 4, "\nnode " << node << " is already coarse\n"); }
			else { UG_DLOG(LIB_ALG_AMG, 4, "\nnode " << node << " has rating " << rating[node] << ". now gets coarse.\nUpdate Neighbors of " << node << "\n"); }

			if(!rating[node].is_coarse())
			{
				heap.remove(node);
				rating.set_coarse(node);
				UpdateNeighbors(SymmNeighGraph, node, possible_parents, rating, heap);
			}
			PoldIndices(i, node) = n.parents[j].value;
		}
	}

	if(bTiming) UG_DLOG(LIB_ALG_AMG, 1, "took " << SW.ms() << " ms.");
}

}

#endif
