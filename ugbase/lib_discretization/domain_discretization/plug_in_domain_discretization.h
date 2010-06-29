/*
 * plug_in_spacial_discretization.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__PLUG_IN_SPACIAL_DISCRETIZATION__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__PLUG_IN_SPACIAL_DISCRETIZATION__

// extern includes
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/domain_discretization/domain_discretization_interface.h"
#include "lib_discretization/domain_discretization/disc_helper/fvgeom.h"
#include "lib_discretization/domain_discretization/dirichlet_bnd_values.h"

namespace ug {

template <	typename TDiscreteFunction,
			int ref_dim,
			template<typename TDomain, typename TAlgebra > class TElemDisc>
class PlugInDomainDiscretization : public IDimensionDomainDiscretization<TDiscreteFunction> {
	protected:
		// type of discrete function
		typedef TDiscreteFunction discrete_function_type;

		// type of domain
		typedef typename TDiscreteFunction::domain_type domain_type;

		// type of algebra
		typedef typename TDiscreteFunction::algebra_type algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of local matrix
		typedef typename matrix_type::local_matrix_type local_matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

		// type of algebra vector
		typedef typename vector_type::local_vector_type local_vector_type;

		// type of multi_index used in algebra
		typedef typename matrix_type::index_type index_type;

		// type of local index container
		typedef typename matrix_type::local_index_type local_index_type;

	public:
		PlugInDomainDiscretization(TElemDisc<domain_type, algebra_type>& elemDisc)
			: m_elemDisc(elemDisc)
		  {}

		// Assemble routines for time independent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, int si);
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, int si);
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, int si);
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si);

		// Assemble routines for time dependent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, int si, number time, number s_m, number s_a);
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, int si, number time, number s_m, number s_a);
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, int si, number time, number s_m, number s_a);
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si, number time, number s_m, number s_a);

		virtual size_t num_fct() const{ return m_elemDisc.num_fct();}

	protected:
		// Element assembling routines
		template <typename TElem>
		bool assemble_jacobian_defect	(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											matrix_type& J, vector_type& d, const discrete_function_type& u,
											number time, number s_m, number s_a);
		template <typename TElem>
		bool assemble_jacobian			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											matrix_type& J, const discrete_function_type& u,
											number time, number s_m, number s_a);
		template <typename TElem>
		bool assemble_defect			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											vector_type& d, const discrete_function_type& u,
											number time, number s_m, number s_a);

		template <typename TElem>
		bool assemble_linear			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											matrix_type& mat, vector_type& rhs, const discrete_function_type& u);

		TElemDisc<domain_type, algebra_type>& m_elemDisc;
};


// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_jacobian_defect(	typename geometry_traits<TElem>::iterator iterBegin,
							typename geometry_traits<TElem>::iterator iterEnd,
							matrix_type& J, vector_type& d, const discrete_function_type& u,
							number time, number s_m, number s_a)
{
	TElem* elem = NULL;

	// get total number of DoF's on this element type 'TElem'
	const size_t num_sh = m_elemDisc.num_sh(elem);

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_d(num_sh);
	local_vector_type loc_d_temp(num_sh);
	local_vector_type loc_u(num_sh);
	local_matrix_type loc_J(num_sh, num_sh);
	local_matrix_type loc_J_temp(num_sh, num_sh);

	m_elemDisc.prepare_element_loop(elem);

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// reset index offset
		size_t offset = 0;

		// loop over all functions
		for(size_t i = 0; i < m_elemDisc.num_fct(); i++)
		{
			offset += u.get_multi_indices(elem, m_elemDisc.fct(i), glob_ind, offset);
		}
		UG_ASSERT(offset == num_sh, offset << " indices are read in, but we have " << num_sh << " dofs on this element.\n");

		// read dof values
		u.get_dof_values(loc_u, glob_ind);

		// reset local matrix and rhs
		loc_d.set(0.0);
		loc_J.set(0.0);

		// prepare element
		m_elemDisc.prepare_element(elem, loc_u, glob_ind);

		// Assemble JA
		loc_J_temp.set(0.0);
		m_elemDisc.assemble_element_JA(elem, loc_J_temp, loc_u, time);
		loc_J += loc_J_temp * s_a;

		// Assemble JM
		loc_J_temp.set(0.0);
		m_elemDisc.assemble_element_JM(elem, loc_J_temp, loc_u, time);
		loc_J += loc_J_temp * s_m;

		// Assemble A
		loc_d_temp.set(0.0);
		m_elemDisc.assemble_element_A(elem, loc_d_temp, loc_u, time);
		loc_d += loc_d_temp * s_a;

		// Assemble M
		loc_d_temp.set(0.0);
		m_elemDisc.assemble_element_M(elem, loc_d_temp, loc_u, time);
		loc_d += loc_d_temp * s_m;

		// Assemble f
		loc_d_temp.set(0.0);
		m_elemDisc.assemble_element_f(elem, loc_d_temp, time);
		loc_d -= loc_d_temp * s_a;

		// send local to global matrix
		J.add(loc_J, glob_ind, glob_ind);
		d.add(loc_d, glob_ind);
	}

	// finish element loop
	m_elemDisc.finish_element_loop(elem);

	return true;
}

// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_jacobian(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					matrix_type& J, const discrete_function_type& u,
					number time, number s_m, number s_a)
{
	TElem* elem = NULL;

	// get total number of DoF's on this element type 'TElem'
	const size_t num_sh = m_elemDisc.num_sh(elem);

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_u(num_sh);
	local_matrix_type loc_J(num_sh, num_sh);
	local_matrix_type loc_J_temp(num_sh, num_sh);

	m_elemDisc.prepare_element_loop(elem);

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// reset index offset
		size_t offset = 0;

		// loop over all functions
		for(size_t i = 0; i < m_elemDisc.num_fct(); i++)
		{
			offset += u.get_multi_indices(elem, m_elemDisc.fct(i), glob_ind, offset);
		}
		UG_ASSERT(offset == num_sh, offset << " indices are read in, but we have " << num_sh << " dofs on this element.\n");

		// read local values of current solution
		u.get_dof_values(loc_u, glob_ind);

		// reset local matrix and rhs
		loc_J.set(0.0);

		// prepare element
		m_elemDisc.prepare_element(elem, loc_u, glob_ind);

		// Assemble JA
		loc_J_temp.set(0.0);
		m_elemDisc.assemble_element_JA(elem, loc_J_temp, loc_u, time);
		loc_J += loc_J_temp * s_a;

		// Assemble JM
		loc_J_temp.set(0.0);
		m_elemDisc.assemble_element_JM(elem, loc_J_temp, loc_u, time);
		loc_J += loc_J_temp * s_m;

		// send local to global matrix
		J.add(loc_J, glob_ind, glob_ind);

		UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Adding local Matrix: \n" << loc_J << " for indices " << glob_ind << " x " << glob_ind << ".\n");
	}

	// finish element loop
	m_elemDisc.finish_element_loop(elem);

	return true;
}

// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_defect(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					vector_type& d, const discrete_function_type& u,
					number time, number s_m, number s_a)
{
	TElem* elem = NULL;

	// get total number of DoF's on this element type 'TElem'
	const size_t num_sh = m_elemDisc.num_sh(elem);

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_d(num_sh);
	local_vector_type loc_d_temp(num_sh);
	local_vector_type loc_u(num_sh);

	m_elemDisc.prepare_element_loop(elem);

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// reset index offset
		size_t offset = 0;

		// loop over all functions
		for(size_t i = 0; i < m_elemDisc.num_fct(); i++)
		{
			offset += u.get_multi_indices(elem, m_elemDisc.fct(i), glob_ind, offset);
		}
		UG_ASSERT(offset == num_sh, offset << " indices are read in, but we have " << num_sh << " dofs on this element.\n");

		// read values
		u.get_dof_values(loc_u, glob_ind);

		// reset local matrix and rhs
		loc_d.set(0.0);

		m_elemDisc.prepare_element(elem, loc_u, glob_ind);

		// Assemble A
		loc_d_temp.set(0.0);
		m_elemDisc.assemble_element_A(elem, loc_d_temp, loc_u, time);
		loc_d += loc_d_temp * s_a;

		// Assemble M
		loc_d_temp.set(0.0);
		m_elemDisc.assemble_element_M(elem, loc_d_temp, loc_u, time);
		loc_d += loc_d_temp * s_m;

		// Assemble f
		loc_d_temp.set(0.0);
		m_elemDisc.assemble_element_f(elem, loc_d_temp, time);
		loc_d -= loc_d_temp * s_a;

		// send local to global matrix
		d.add(loc_d, glob_ind);
	}

	// finish element loop
	m_elemDisc.finish_element_loop(elem);

	return true;
}



// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_linear(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
{
	TElem* elem = NULL;

	// get total number of DoF's on this element type 'TElem'
	const size_t num_sh = m_elemDisc.num_sh(elem);

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_rhs(num_sh);
	local_vector_type loc_u(num_sh);
	local_matrix_type loc_mat(num_sh, num_sh);

	// prepare for loop
	if(m_elemDisc.prepare_element_loop(elem) != IPlugInReturn_OK) return false;;

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// reset index offset
		size_t offset = 0;

		// loop over all functions
		for(size_t i = 0; i < m_elemDisc.num_fct(); i++)
		{
			offset += u.get_multi_indices(elem, m_elemDisc.fct(i), glob_ind, offset);
		}
		UG_ASSERT(offset == num_sh, offset << " indices are read in, but we have " << num_sh << " dofs on this element.\n");

		u.get_dof_values(loc_u, glob_ind);

		// reset local matrix and rhs
		loc_mat.set(0.0);
		loc_rhs.set(0.0);

		// prepare element data
		if(m_elemDisc.prepare_element(elem, loc_u, glob_ind) != IPlugInReturn_OK) return false;

		// assemble stiffness matrix for inner elements
		if(m_elemDisc.assemble_element_JA(elem, loc_mat, loc_u) != IPlugInReturn_OK) return false;

		// assemble rhs for inner elements
		if(m_elemDisc.assemble_element_f(elem, loc_rhs) != IPlugInReturn_OK) return false;

		// send local to global (matrix and rhs) [this is a virtual call]
		mat.add(loc_mat, glob_ind, glob_ind);
		rhs.add(loc_rhs, glob_ind);
	}

	// finish element loop
	if(m_elemDisc.finish_element_loop(elem) != IPlugInReturn_OK) return false;

	return true;
}


template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, int si)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u, int si)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	//return IAssemble_NOT_IMPLEMENTED;

	return this->assemble_jacobian(J, u, si, 0.0, 0.0, 1.0);
}

template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_defect(vector_type& d, const discrete_function_type& u, int si)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(size_t i = 0; i < m_elemDisc.num_fct(); i++){
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
			{UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");return IAssemble_ERROR;}
		}

	if(assemble_linear<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), mat, rhs, u)==false)
		{UG_LOG("Error in assemble_linear, aborting.\n");return IAssemble_ERROR;}
	if(assemble_linear<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), mat, rhs, u)==false)
		{UG_LOG("Error in assemble_linear, aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, int si, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(size_t i = 0; i < m_elemDisc.num_fct(); i++){
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
			{UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");return IAssemble_ERROR;}
		}

	if(assemble_jacobian_defect<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), J, d, u, time, s_m, s_a)==false)
		{UG_LOG("Error in assemble_linear, aborting.\n");return IAssemble_ERROR;}
	if(assemble_jacobian_defect<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), J, d, u, time, s_m, s_a)==false)
		{UG_LOG("Error in assemble_linear, aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u, int si, number time, number s_m, number s_a)
{

	// check that Solution number 'nr' matches trial space required by Discretization
	for(size_t i = 0; i < m_elemDisc.num_fct(); i++){
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
			{UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");return IAssemble_ERROR;}
		}

	if(assemble_jacobian<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), J, u, time, s_m, s_a)==false)
		{UG_LOG("Error in 'assemble_jacobian' while calling 'assemble_jacobian<Triangle>', aborting.\n");return IAssemble_ERROR;}
	if(assemble_jacobian<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), J, u, time, s_m, s_a)==false)
		{UG_LOG("Error in 'assemble_jacobian' while calling 'assemble_jacobian<Quadrilateral>', aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_defect(vector_type& d, const discrete_function_type& u, int si, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(size_t i = 0; i < m_elemDisc.num_fct(); i++){
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
		{UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");return IAssemble_ERROR;}
		}

	if(assemble_defect<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), d, u, time, s_m, s_a)==false)
		{UG_LOG("Error in assemble_defect, aborting.\n");return IAssemble_ERROR;}
	if(assemble_defect<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), d, u, time, s_m, s_a)==false)
		{UG_LOG("Error in assemble_defect, aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, int ref_dim, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, ref_dim, TElemDisc>::
assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si, number time, number s_m, number s_a)
{
	return IAssemble_NOT_IMPLEMENTED;
}



}
#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__PLUG_IN_SPACIAL_DISCRETIZATION__ */
