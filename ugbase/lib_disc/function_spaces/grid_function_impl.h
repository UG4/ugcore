/*
 * grid_function_impl.h
 *
 *  Created on: 13.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__

#include "grid_function.h"

#include "lib_algebra/algebra_type.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// GridFunction
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>::
GridFunction(SmartPtr<approximation_space_type> approxSpace,
             SmartPtr<DoFDistribution> spDoFDistr, bool bManage)
 : m_spDD(spDoFDistr), m_spApproxSpace(approxSpace)
{
	if(!m_spDD.valid()) UG_THROW("DoF Distribution is null.");

	if(bManage)
		m_spDD->manage_grid_function(*this);

	check_algebra();
	resize_values(num_indices());
#ifdef UG_PARALLEL
//	set layouts
	this->set_layouts(m_spDD->layouts());

//	set storage type
	this->set_storage_type(PST_UNDEFINED);
#endif
};

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>::
GridFunction(SmartPtr<approximation_space_type> approxSpace, bool bManage)
	: m_spDD(approxSpace->surface_dof_distribution()),
	  m_spApproxSpace(approxSpace)
{
	if(!m_spDD.valid()) UG_THROW("DoF Distribution is null.");

	if(bManage)
		m_spDD->manage_grid_function(*this);

	m_bManaged = bManage;

	check_algebra();
	resize_values(num_indices());
#ifdef UG_PARALLEL
//	set layouts
	this->set_layouts(m_spDD->layouts());

//	set storage type
	this->set_storage_type(PST_UNDEFINED);
#endif
};

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>::
GridFunction(SmartPtr<approximation_space_type> approxSpace, int level, bool bManage)
	: m_spDD(approxSpace->surface_dof_distribution(level)),
	  m_spApproxSpace(approxSpace)
{
	if(!m_spDD.valid()) UG_THROW("DoF Distribution is null.");

	if(bManage)
		m_spDD->manage_grid_function(*this);

	check_algebra();
	resize_values(num_indices());
#ifdef UG_PARALLEL
//	set layouts
	this->set_layouts(m_spDD->layouts());

//	set storage type
	this->set_storage_type(PST_UNDEFINED);
#endif
};

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>::
GridFunction(SmartPtr<approximation_space_type> approxSpace, const GridLevel& gl, bool bManage)
	: m_spDD(approxSpace->dof_distribution(gl)),
	  m_spApproxSpace(approxSpace)
{
	if(!m_spDD.valid()) UG_THROW("DoF Distribution is null.");

	if(bManage)
		m_spDD->manage_grid_function(*this);

	check_algebra();
	resize_values(num_indices());
#ifdef UG_PARALLEL
//	set layouts
	this->set_layouts(m_spDD->layouts());

//	set storage type
	this->set_storage_type(PST_UNDEFINED);
#endif
};

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::check_algebra()
{
//	get blocksize of algebra
	const int blockSize = algebra_type::blockSize;

//	a)	If blocksize fixed and > 1, we need grouping.
	if(blockSize > 1 && !this->m_spDD->grouped())
	{
		UG_THROW("Fixed block algebra needs grouped dofs.");
	}
//	b) 	If blocksize flexible, we group
	else if (blockSize == AlgebraType::VariableBlockSize
			&& !this->m_spDD->grouped())
	{
		UG_THROW("Variable block algebra needs grouped dofs.");
	}
//	c)	If blocksize == 1, we do not group. This will allow us to handle
//		this case for any problem.
	else if (blockSize == 1 && this->m_spDD->grouped())
	{
		UG_THROW("block 1x1 algebra needs non-grouped dofs.");
	}
}

template <typename TDomain, typename TAlgebra>
template <typename TElem>
bool
GridFunction<TDomain, TAlgebra>::
dof_positions(TElem* elem, size_t fct, std::vector<MathVector<dim> >& vPos) const
{
	return DoFPosition(vPos, elem, *domain(),
	                   this->local_finite_element_id(fct),
	                   this->dim(fct));
};

template <typename TDomain,typename TAlgebra>
template <typename TElem>
bool
GridFunction<TDomain, TAlgebra>::
inner_dof_positions(TElem* elem, size_t fct, std::vector<MathVector<dim> >& vPos) const
{
	return InnerDoFPosition(vPos, elem, *domain(),
	                        this->local_finite_element_id(fct),
	                        this->dim(fct));
};

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
clone_pattern(const this_type& v)
{
// 	copy approximation space
	m_spApproxSpace = v.m_spApproxSpace;

//	assign dof distribution (resizes vector)
	this->m_spDD = v.m_spDD;

//	resize the vector
	resize_values(v.size());

#ifdef UG_PARALLEL
//	set layouts
	this->set_layouts(v.layouts());

//	copy storage type
	this->set_storage_type(v.get_storage_mask());
#endif
};

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
resize_values(size_t s, number defaultValue)
{
//	remember old values
	const size_t oldSize = vector_type::size();

//	resize vector
	vector_type::resize_sloppy(s);

//	set vector to zero-values
	for(size_t i = oldSize; i < s; ++i)
		this->operator[](i) = defaultValue;
}

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
permute_values(const std::vector<size_t>& vIndNew)
{
//	check sizes
	if(vIndNew.size() != this->size())
		UG_THROW("GridFunction::permute_values: For a permutation the"
				" index set must have same cardinality as vector.");

// \todo: avoid tmp vector, only copy values into new vector and use that one
//	create tmp vector
	vector_type vecTmp; vecTmp.resize(this->size());
#ifdef UG_PARALLEL
//	copy storage type
	vecTmp.set_storage_type(this->get_storage_mask());
#endif

//	loop indices and copy values
	for(size_t i = 0; i < vIndNew.size(); ++i)
		vecTmp[vIndNew[i]] = this->operator[](i);

//	copy tmp vector into this vector
	this->assign(vecTmp);
}

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,bool bDisjunct)
{
//	disjunct case
	if(bDisjunct)
		for(size_t i = 0; i < vIndexMap.size(); ++i)
			this->operator[](vIndexMap[i].second)
				= this->operator[](vIndexMap[i].first);
	else {
		typedef typename vector_type::value_type value_type;
		std::vector<value_type> values;
		values.resize(vIndexMap[vIndexMap.size()-1].first);
		for(size_t i = 0; i < vIndexMap.size(); ++i){ 
			const size_t index = vIndexMap[i].first;
			if (index>=values.size()) values.resize(index+1);
			values[index] = this->operator[](index);
		}
		for(size_t i = 0; i < vIndexMap.size(); ++i)
			this->operator[](vIndexMap[i].second)
				= values[vIndexMap[i].first];
	}
}

template <typename TDomain, typename TAlgebra>
void GridFunction<TDomain, TAlgebra>::assign(const vector_type& v)
{
//	check size
	if(v.size() != vector_type::size())
		UG_THROW("GridFunction: Assigned vector has incorrect size.");

//	assign vector
	*(dynamic_cast<vector_type*>(this)) = v;

#ifdef UG_PARALLEL
//	set layouts
	this->set_layouts(v.layouts());

//	copy storage type
	this->set_storage_type(v.get_storage_mask());
#endif
}

template <typename TDomain, typename TAlgebra>
void GridFunction<TDomain, TAlgebra>::assign(const this_type& v)
{
// 	copy approximation space
	m_spApproxSpace = v.m_spApproxSpace;

//	assign dof distribution (resizes vector)
	this->m_spDD = v.m_spDD;

	if(v.m_bManaged==true)
		this->m_spDD->manage_grid_function(*this);

//	resize the vector
	resize_values(v.size());

//  copy values
	*(dynamic_cast<vector_type*>(this)) = *dynamic_cast<const vector_type*>(&v);

#ifdef UG_PARALLEL
//	set layouts
	this->set_layouts(v.layouts());

//	copy storage type
	this->set_storage_type(v.get_storage_mask());
#endif
}

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>*
GridFunction<TDomain, TAlgebra>::virtual_clone_without_values() const
{
	GridFunction<TDomain, TAlgebra>* v
		= new GridFunction<TDomain, TAlgebra>(m_spApproxSpace, this->m_spDD);
	v->resize_values(this->size());
#ifdef UG_PARALLEL
	v->set_layouts(this->layouts());
	v->set_storage_type(PST_UNDEFINED);
#endif
	return v;
}


template <typename TDomain, typename TAlgebra>
void GridFunction<TDomain, TAlgebra>::
add_transfer(SmartPtr<ILocalTransferAlgebra<TAlgebra> > spTransfer)
{
	spTransfer->set_vector(this);
	m_spDD->add_transfer(spTransfer);
}

template <typename TDomain, typename TAlgebra>
void GridFunction<TDomain, TAlgebra>::add_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	m_spDD->add_transfer(spTransfer);
}

template <typename TDomain, typename TAlgebra>
void GridFunction<TDomain, TAlgebra>::remove_transfer(SmartPtr<ILocalTransfer> spTransfer)
{
	m_spDD->remove_transfer(spTransfer);
}

template <typename TDomain, typename TAlgebra>
void GridFunction<TDomain, TAlgebra>::clear_transfers()
{
	m_spDD->clear_transfers();
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__ */
