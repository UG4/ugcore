/*
 * ass_tuner_impl.h
 *
 *  Created on: 01.03.2013
 *      Author:	raphaelprohl, Andreas Vogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER_IMPL__

#include "ass_tuner.h"

namespace ug{

template <typename TAlgebra>
void AssemblingTuner<TAlgebra>::resize(ConstSmartPtr<DoFDistribution> dd,
                                  vector_type& vec)	const
{
	if (single_dof_index_assembling_enabled()){ vec.resize(1);}
	else{
		const size_t numIndex = dd->num_indices();
		vec.resize(numIndex);
	}
	vec.set(0.0);
}

template <typename TAlgebra>
void AssemblingTuner<TAlgebra>::resize(ConstSmartPtr<DoFDistribution> dd,
								  matrix_type& mat) const
{
	if (single_dof_index_assembling_enabled()){ mat.resize_and_clear(1, 1);
	}
	else{
		const size_t numIndex = dd->num_indices();
		mat.resize_and_clear(numIndex, numIndex);
	}
}

template <typename TAlgebra>
template <typename TElem>
bool AssemblingTuner<TAlgebra>::element_used(TElem* elem) const
{
	if(m_pBoolMarker)
		if(!m_pBoolMarker->is_marked(elem)) return false;

	return true;
}


template <typename TAlgebra>
template <typename TElem>
void AssemblingTuner<TAlgebra>::collect_selected_elements(std::vector<TElem*>& vElem,
                                                     ConstSmartPtr<DoFDistribution> dd, int si) const
{
	if (!m_pSelector)
		UG_THROW("Selector-iterator not set!")

	Selector* sel = m_pSelector;
	const ISubsetHandler& sh = *dd->subset_handler();

	for(typename Selector::traits<TElem>::iterator iter = sel->begin<TElem>();
		iter != sel->end<TElem>(); ++iter)
	{
		if(sh.get_subset_index(*iter) == si)
			vElem.push_back(*iter);
	}
}

template <typename TAlgebra>
void AssemblingTuner<TAlgebra>::set_dirichlet_row(matrix_type& mat, const DoFIndex& ind) const
{
	// 	check if assembling has been carried out with respect to one index only.
	//	For that case assembling-matrices have been resized to a block-matrix at one DoF only.
	if(single_dof_index_assembling_enabled())
	{
		UG_ASSERT(mat.num_rows() == 1, "#rows needs to be 1 for setting Dirichlet "
								"in an index-wise manner.");
		UG_ASSERT(mat.num_cols() == 1, "#cols needs to be 1 for setting Dirichlet "
				"in an index-wise manner.");

		if (ind == m_SingleAssDoFIndex)
			SetDirichletRow(mat, ind);
	}
	else{
		SetDirichletRow(mat, ind);
	}
}

template <typename TAlgebra>
void AssemblingTuner<TAlgebra>::set_dirichlet_val(vector_type& vec, const DoFIndex& ind, const double val) const
{
	//	check if assembling has been carried out with respect to one index only.
	//	For that case assembling-vectors have been resized to a block-vector at one DoF only.
	if(single_dof_index_assembling_enabled())
	{
		UG_ASSERT(vec.size() == 1, "vector-size needs to be 1 for setting Dirichlet "
			"in an index-wise manner.");

		UG_LOG("single_dof_index_enabled().set_dirichlet_val m_SingleAssDoFIndex: " << m_SingleAssDoFIndex << "\n");
		UG_LOG("single_dof_index_enabled().set_dirichlet_val ind: " << ind << "\n");

		if(ind == m_SingleAssDoFIndex)
			DoFRef(vec, ind) = val;
	}
	else{
		DoFRef(vec, ind) = val;
	}
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER_IMPL__ */
