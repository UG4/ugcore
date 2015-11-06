/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Raphael Prohl, Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER_IMPL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER_IMPL__

#include "ass_tuner.h"

namespace ug{

template <typename TAlgebra>
void AssemblingTuner<TAlgebra>::resize(ConstSmartPtr<DoFDistribution> dd,
                                  vector_type& vec)	const
{
	if (single_index_assembling_enabled()){ vec.resize(1);}
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
	if (single_index_assembling_enabled()){ mat.resize_and_clear(1, 1);
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
	if(single_index_assembling_enabled())
	{
		if (mat.num_rows() != 1 || mat.num_cols() != 1)
			UG_THROW("#rows and #cols need to be 1 for setting dirichlet rows"
					" in an index-wise manner.")

		const size_t index = ind[0];
		if (index == m_SingleAssIndex)
			SetDirichletRow(mat, 0, ind[1]);
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
	if(single_index_assembling_enabled())
	{
		if(vec.size() != 1)
			UG_THROW("vector-size needs to be 1 for setting dirichlet values"
					" in an index-wise manner.");

		const size_t index = ind[0];
		if(index == m_SingleAssIndex)
			BlockRef(vec[0], ind[1]) = val;
	}
	else{
		DoFRef(vec, ind) = val;
	}
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ASS_TUNER_IMPL__ */
