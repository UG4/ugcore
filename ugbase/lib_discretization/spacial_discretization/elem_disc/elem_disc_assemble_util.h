/*
 * assemble_elem_disc.h
 *
 *  Created on: 08.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// intern headers
#include "lib_discretization/reference_element/reference_elements.h"
#include "elem_disc_interface.h"
#include "../function_group.h"

namespace ug {

//////////////////////////////////
// Assemble Jacobian
//////////////////////////////////

// assemble elements of type TElem in d dimensions
template <	typename TElem,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleJacobian(	IElemDisc<TAlgebra>& elemDisc,
					typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					typename TAlgebra::matrix_type& J,
					const TDiscreteFunction& u,
					const FunctionGroup& fcts,
					number time, number s_m, number s_a)
{
	typedef typename TAlgebra::matrix_type::local_matrix_type local_matrix_type;
	typedef typename TAlgebra::vector_type::local_vector_type local_vector_type;
	typedef typename TAlgebra::matrix_type::local_index_type local_index_type;
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// set elem type
	if(!elemDisc.set_geometric_object_type(reference_element_type::REFERENCE_ELEMENT_TYPE, IEDN_LINEAR))
		{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

	// get total number of DoF's on this element type 'TElem'
	const size_t num_sh = elemDisc.num_total_sh();

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_u(num_sh);
	local_matrix_type loc_J(num_sh, num_sh);
	local_matrix_type loc_J_temp(num_sh, num_sh);

	// prepare loop
	elemDisc.prepare_element_loop();

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem* elem = *iter;

		// reset index offset
		size_t offset = 0;

		// loop over all functions
		for(size_t i = 0; i < fcts.num_functions(); i++)
		{
			offset += u.get_multi_indices(elem, fcts[i], glob_ind, offset);
		}
		UG_ASSERT(offset == num_sh, offset << " indices are read in, but we have " << num_sh << " dofs on this element.\n");

		// read local values of current solution
		u.get_dof_values(loc_u, glob_ind);

		// reset local matrix and rhs
		loc_J.set(0.0);

		// prepare element
		elemDisc.prepare_element(elem, loc_u, glob_ind);

		// Assemble JA
		loc_J_temp.set(0.0);
		elemDisc.assemble_JA(loc_J_temp, loc_u, time);
		loc_J += loc_J_temp * s_a;

		// Assemble JM
		loc_J_temp.set(0.0);
		elemDisc.assemble_JM(loc_J_temp, loc_u, time);
		loc_J += loc_J_temp * s_m;

		// send local to global matrix
		J.add(loc_J, glob_ind, glob_ind);

		UG_DLOG(LIB_DISC_ASSEMBLE, 3, "Adding local Matrix: \n" << loc_J << " for indices " << glob_ind << " x " << glob_ind << ".\n");
	}

	// finish element loop
	elemDisc.finish_element_loop();

	return true;
}

//////////////////////////////////
// Assemble Defect
//////////////////////////////////

template <	typename TElem,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleDefect(	IElemDisc<TAlgebra>& elemDisc,
				typename geometry_traits<TElem>::iterator iterBegin,
				typename geometry_traits<TElem>::iterator iterEnd,
				typename TAlgebra::vector_type& d,
				const TDiscreteFunction& u,
				const FunctionGroup& fcts,
				number time, number s_m, number s_a)
{
	typedef typename TAlgebra::matrix_type::local_matrix_type local_matrix_type;
	typedef typename TAlgebra::vector_type::local_vector_type local_vector_type;
	typedef typename TAlgebra::matrix_type::local_index_type local_index_type;
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// set elem type
	if(!elemDisc.set_geometric_object_type(reference_element_type::REFERENCE_ELEMENT_TYPE, IEDN_LINEAR))
		{UG_LOG("ERROR in AssembleDefect: Cannot set geometric object type.\n"); return false;}

	// get total number of DoF's on this element type 'TElem'
	const size_t num_sh = elemDisc.num_total_sh();

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_d(num_sh);
	local_vector_type loc_d_temp(num_sh);
	local_vector_type loc_u(num_sh);

	// prepare loop
	elemDisc.prepare_element_loop();

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem* elem = *iter;

		// reset index offset
		size_t offset = 0;

		// loop over all functions
		for(size_t i = 0; i < fcts.num_functions(); i++)
		{
			offset += u.get_multi_indices(elem, fcts[i], glob_ind, offset);
		}
		UG_ASSERT(offset == num_sh, offset << " indices are read in, but we have " << num_sh << " dofs on this element.\n");

		// read values
		u.get_dof_values(loc_u, glob_ind);

		// reset local matrix and rhs
		loc_d.set(0.0);

		elemDisc.prepare_element(elem, loc_u, glob_ind);

		// Assemble A
		loc_d_temp.set(0.0);
		elemDisc.assemble_A(loc_d_temp, loc_u, time);
		loc_d += loc_d_temp * s_a;

		// Assemble M
		loc_d_temp.set(0.0);
		elemDisc.assemble_M(loc_d_temp, loc_u, time);
		loc_d += loc_d_temp * s_m;

		// Assemble f
		loc_d_temp.set(0.0);
		elemDisc.assemble_f(loc_d_temp, time);
		loc_d -= loc_d_temp * s_a;

		// send local to global matrix
		d.add(loc_d, glob_ind);
	}

	// finish element loop
	elemDisc.finish_element_loop();

	return true;
}

//////////////////////////////////
// Assemble linear
//////////////////////////////////

template <	typename TElem,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleLinear(	IElemDisc<TAlgebra>& elemDisc,
				typename geometry_traits<TElem>::iterator iterBegin,
				typename geometry_traits<TElem>::iterator iterEnd,
				typename TAlgebra::matrix_type& mat,
				typename TAlgebra::vector_type& rhs,
				const TDiscreteFunction& u,
				const FunctionGroup& fcts)
{
	typedef typename TAlgebra::matrix_type::local_matrix_type local_matrix_type;
	typedef typename TAlgebra::vector_type::local_vector_type local_vector_type;
	typedef typename TAlgebra::matrix_type::local_index_type local_index_type;
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// set elem type
	if(!elemDisc.set_geometric_object_type(reference_element_type::REFERENCE_ELEMENT_TYPE, IEDN_LINEAR))
		{UG_LOG("ERROR in AssembleLinear: Cannot set geometric object type.\n"); return false;}

	// get total number of DoF's on this element type 'TElem'
	const size_t num_sh = elemDisc.num_total_sh();

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_rhs(num_sh);
	local_vector_type loc_u(num_sh);
	local_matrix_type loc_mat(num_sh, num_sh);

	// prepare loop
	elemDisc.prepare_element_loop();

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem* elem = *iter;

		// reset index offset
		size_t offset = 0;

		// loop over all functions
		for(size_t i = 0; i < fcts.num_functions(); i++)
		{
			offset += u.get_multi_indices(elem, fcts[i], glob_ind, offset);
		}
		UG_ASSERT(offset == num_sh, offset << " indices are read in, but we have " << num_sh << " dofs on this element.\n");

		u.get_dof_values(loc_u, glob_ind);

		// reset local matrix and rhs
		loc_mat.set(0.0);
		loc_rhs.set(0.0);

		// prepare element data
		if(!elemDisc.prepare_element(elem, loc_u, glob_ind)) return false;

		// assemble stiffness matrix for inner elements
		if(!elemDisc.assemble_JA(loc_mat, loc_u)) return false;

		// assemble rhs for inner elements
		if(!elemDisc.assemble_f(loc_rhs)) return false;

		// send local to global (matrix and rhs)
		mat.add(loc_mat, glob_ind, glob_ind);
		rhs.add(loc_rhs, glob_ind);
	}

	// finish element loop
	if(!elemDisc.finish_element_loop()) return false;

	return true;
};

} // end namespace ug


#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__ */
