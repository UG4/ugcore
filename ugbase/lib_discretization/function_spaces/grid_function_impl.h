/*
 * grid_function_impl.h
 *
 *  Created on: 13.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_IMPL__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_IMPL__

#include "grid_function.h"

namespace ug{

template <typename TDoFDistribution>
void
IGridFunction<TDoFDistribution>::
assign_dof_distribution(typename IGridFunction<TDoFDistribution>::dof_distribution_type& DoFDistr, bool adapt)
{
//	unregister from dof distribution iff already dd set
	if(m_pDD != NULL)
		m_pDD->unmanage_grid_function(*this);

//	remember new dof distribution
	m_pDD = &DoFDistr;

//	schedule for adaption
	if(adapt)
		m_pDD->manage_grid_function(*this);

//	resize the vector
	resize_values(num_indices());
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void
GridFunction<TDomain, TDoFDistribution, TAlgebra>::
clone_pattern(const this_type& v)
{
// 	copy approximation space
	assign_approximation_space(*v.m_pApproxSpace);

//	assign dof distribution (resizes vector)
	assign_dof_distribution(*v.m_pDD);
};

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void
GridFunction<TDomain, TDoFDistribution, TAlgebra>::
resize_values(size_t s, number defaultValue)
{
//	remember old values
	const size_t oldSize = vector_type::size();

//	resize vector
	vector_type::resize(s);

//	set vector to zero-values
	for(size_t i = oldSize; i < s; ++i)
		this->operator[](i) = defaultValue;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
GridFunction<TDomain, TDoFDistribution, TAlgebra>::
permute_values(const std::vector<size_t>& vIndNew)
{
//	check sizes
	if(vIndNew.size() != this->size())
	{
		UG_LOG("ERROR in GridFunction::swap_values: For a permutation the"
				" index set must have same cardinality as vector.\n");
		return false;
	}

// \todo: avoid tmp vector, only copy values into new vector and use that one
//	create tmp vector
	vector_type vecTmp; vecTmp.resize(this->size());

//	loop indices and copy values
	for(size_t i = 0; i < vIndNew.size(); ++i)
		vecTmp[vIndNew[i]] = this->operator[](i);

//	copy tmp vector into this vector
	if(!this->assign(vecTmp)) return false;

//	we're done
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
GridFunction<TDomain, TDoFDistribution, TAlgebra>::
copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,bool bDisjunct)
{
//	disjunct case
	if(bDisjunct)
		for(size_t i = 0; i < vIndexMap.size(); ++i)
			this->operator[](vIndexMap[i].second)
				= this->operator[](vIndexMap[i].first);
//	other case not implemented
	else return false;

//	we're done
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
GridFunction<TDomain, TDoFDistribution, TAlgebra>::
assign(const vector_type& v)
{
//	check size
	if(v.size() != vector_type::size())
	{
		UG_LOG("ERROR in GridFunction::assign:"
				"Assigned vector has incorrect size.\n");
		return false;
	}

//	assign vector
	*(dynamic_cast<vector_type*>(this)) = v;

//	we're done
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
GridFunction<TDomain, TDoFDistribution, TAlgebra>::
assign(const this_type& v)
{
// 	copy approximation space
	assign_approximation_space(*v.m_pApproxSpace);

//	assign dof distribution (resizes vector)
	assign_dof_distribution(*v.m_pDD);

//  copy values
	*(dynamic_cast<vector_type*>(this)) = *dynamic_cast<const vector_type*>(&v);

//	we're done
	return true;
}

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_IMPL__ */
