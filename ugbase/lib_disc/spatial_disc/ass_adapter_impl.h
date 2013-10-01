/*
 * ass_adapter_impl.h
 *
 *  Created on: 01.03.2013
 *      Author:	raphaelprohl
 */

#ifndef ASS_ADAPTER_IMPL_H_
#define ASS_ADAPTER_IMPL_H_

#include "ass_adapter.h"

namespace ug{

template <typename TAlgebra>
void AssAdapter<TAlgebra>::resize(ConstSmartPtr<DoFDistribution> dd,
		vector_type& vec)
{
	if (m_assIndex.index_set){ vec.resize(1);}
	else{
		const size_t numIndex = dd->num_indices();
		vec.resize(numIndex);
	}
	vec.set(0.0);
}

template <typename TAlgebra>
void AssAdapter<TAlgebra>::resize(ConstSmartPtr<DoFDistribution> dd,
		matrix_type& mat)
{
	if (m_assIndex.index_set){ mat.resize_and_clear(1, 1);
	}
	else{
		const size_t numIndex = dd->num_indices();
		mat.resize_and_clear(numIndex, numIndex);
	}
}

template <typename TAlgebra>
template <typename TElem>
void AssAdapter<TAlgebra>::elemIter_fromSel(ConstSmartPtr<DoFDistribution> dd,
    	int si, std::vector<TElem*>& elems)
{
	if (!m_pSelector)
		UG_THROW("Selector-iterator not set!")

	Selector* sel = m_pSelector;
	const ISubsetHandler& sh = *dd->subset_handler();

	for(typename Selector::traits<TElem>::iterator iter = sel->begin<TElem>();
		iter != sel->end<TElem>(); ++iter)
	{
		if(sh.get_subset_index(*iter) == si)
			elems.push_back(*iter);
	}
}

template <typename TAlgebra>
void AssAdapter<TAlgebra>::adjust_matrix(matrix_type& mat, const DoFIndex& ind)
{
	UG_ASSERT(mat.num_rows() == 1, "#rows needs to be 1 for setting Dirichlet "
			"in an index-wise manner.");
	UG_ASSERT(mat.num_cols() == 1, "#cols needs to be 1 for setting Dirichlet "
			"in an index-wise manner.");

	const size_t index = ind[0];
	const size_t alpha = ind[1];
	if (index == m_assIndex.index)
	{
		typename matrix_type::value_type& block = mat(0,0);

		BlockRef(block, alpha, alpha) = 1.0;
		for(size_t beta = 0; beta < (size_t) GetCols(block); ++beta){
			if(beta != alpha) BlockRef(block, alpha, beta) = 0.0;
		}
	}
}

template <typename TAlgebra>
void AssAdapter<TAlgebra>::adjust_vector(vector_type& vec, const DoFIndex& ind, const double val)
{
	UG_ASSERT(vec.size() == 1, "vector-size needs to be 1 for setting Dirichlet "
			"in an index-wise manner.");

	const size_t index = ind[0];
	const size_t alpha = ind[1];
	if (index == m_assIndex.index){
		BlockRef(vec[0], alpha) = val;
	}
}


} // end namespace ug

#endif /* ASS_ADAPTER_IMPL_H_ */
