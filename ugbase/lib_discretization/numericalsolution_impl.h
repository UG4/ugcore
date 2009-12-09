/*
 * numericalsolution_impl.h
 *
 *  Created on: 03.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__NUMERICALSOLUTION_IMPL__
#define __H__LIBDISCRETIZATION__NUMERICALSOLUTION_IMPL__

namespace ug{

template<typename TElem>
bool NumericalSolutionPattern::get_indices(TElem* elem, int nr_solution, int* ind)
{
	assert(nr_solution < _SingleSolutionInfoVec.size());
	int subsetIndex = _sh->get_subset_index(elem);
	int i;
	for(i = 0; i < _SingleSolutionInfoVec[nr_solution]->num_SubsetIndex; ++i)
	{
		if(subsetIndex == _SingleSolutionInfoVec[nr_solution]->SubsetIndex[i]) break;
	}
	if(i ==  _SingleSolutionInfoVec[nr_solution]->num_SubsetIndex) assert(0 && "ERROR in get_index. Solution not defined for SubsetIndex.");

	int ncomp = _SingleSolutionInfoVec[nr_solution]->group_Comp[subsetIndex];

	for(int i=0; i< geometry_traits<TElem>::Descriptor.num_vertices(); i++)
	{
		VertexBase* vert = elem->vertex(i);
		ind[i] = _aaDoFVRT[vert].nr + ncomp;
	}
}


template <class TElem>
ug::TrialSpace<TElem>& NumericalSolution::TrialSpace(int nr_func, int SubsetIndex)
{
	return TrialSpaces<TElem>::TrialSpace(_pattern->get_trial_space_type(nr_func));
}

template <class TElem>
bool NumericalSolution::get_local_DoFValues(TElem* elem, int nr_func, number* DoFValues)
{
	typename geometry_traits<TElem>::Descriptor TDesc;

	int nvalues = TDesc.num_vertices();

	int* indices = new int[nvalues];

	for(int i=0; i< nvalues; i++)
	{
		VertexBase* vert = elem->vertex(i);
		indices[i] = (int) (_pattern->get_index(vert, nr_func));
	}

	int level = 0;
	if(_pattern->num_levels() >= 2)
	{
		MultiGrid* mg = dynamic_cast<MultiGrid*>(_pattern->get_assigned_subset()->get_assigned_grid());
		assert(mg != NULL && "ERROR in get_local_DoFValues: num_levels > 0 but not a MultiGrid. Aborting.");
		level = mg->get_level(elem);
	}

	double *doubleDoFValues = new double[nvalues];
	if(_GridVector[level]->get_values(nvalues, indices, doubleDoFValues) == false) return false;

	for(int i=0; i< nvalues; i++)
	{
		DoFValues[i] = (number) doubleDoFValues[i];
	}


	return true;
}


}


#endif /* __H__LIBDISCRETIZATION__NUMERICALSOLUTION_IMPL__ */
