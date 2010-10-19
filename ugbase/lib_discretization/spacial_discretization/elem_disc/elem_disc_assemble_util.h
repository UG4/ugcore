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
#include "../../reference_element/reference_element.h"
#include "./elem_disc_interface.h"
#include "lib_discretization/common/common.h"

namespace ug {

//////////////////////////////////
// Assemble Stiffness Matrix
//////////////////////////////////

// assemble elements of type TElem in d dimensions
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleStiffnessMatrix(	IElemDisc<TAlgebra>& elemDisc,
						typename geometry_traits<TElem>::const_iterator iterBegin,
						typename geometry_traits<TElem>::const_iterator iterEnd,
						int si,
						typename TAlgebra::matrix_type& J,
						const typename TAlgebra::vector_type& u,
						const TDoFDistribution& dofDistr)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// get function group
	const FunctionGroup& fcts = elemDisc.get_function_group();

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u;
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_J;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!dofDistr.prepare_indices(refID, si, ind))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IEDN_STIFFNESS))
		{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

	// adjust local algebra
	loc_u.set_indices(ind);
	loc_J.set_indices(ind, ind);

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element loop.\n"); return false;}

	// loop over all elements
	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		dofDistr.update_indices(elem, ind);

		// read local values of u
		loc_u.read_values(u);

		// prepare element
		if(!elemDisc.prepare_element(elem, loc_u, ind))
			{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element.\n"); return false;}

		// reset local matrix and rhs
		loc_J.set(0.0);

		// Assemble JA
		if(!elemDisc.assemble_JA(loc_J, loc_u))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix.\n"); return false;}

		// send local to global matrix
		J.add(loc_J);
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n"); return false;}

	return true;
}


//////////////////////////////////
// Assemble Mass Matrix
//////////////////////////////////

// assemble elements of type TElem in d dimensions
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleMassMatrix(		IElemDisc<TAlgebra>& elemDisc,
						typename geometry_traits<TElem>::const_iterator iterBegin,
						typename geometry_traits<TElem>::const_iterator iterEnd,
						int si,
						typename TAlgebra::matrix_type& J,
						const typename TAlgebra::vector_type& u,
						const TDoFDistribution& dofDistr)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// get function group
	const FunctionGroup& fcts = elemDisc.get_function_group();

	// local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u;
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_J;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!dofDistr.prepare_indices(refID, si, ind))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IEDN_MASS))
		{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

	// adjust local algebra
	loc_u.set_indices(ind);
	loc_J.set_indices(ind, ind);

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element loop.\n"); return false;}

	// loop over all elements
	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		dofDistr.update_indices(elem, ind);

		// read local values of u
		loc_u.read_values(u);

		// prepare element
		if(!elemDisc.prepare_element(elem, loc_u, ind))
			{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element.\n"); return false;}

		// reset local matrix and rhs
		loc_J.set(0.0);

		// Assemble JA
		if(!elemDisc.assemble_JM(loc_J, loc_u))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix.\n"); return false;}

		// send local to global matrix
		J.add(loc_J);
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n"); return false;}

	return true;
}

//////////////////////////////////
// Assemble Jacobian
//////////////////////////////////

// assemble elements of type TElem in d dimensions
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleJacobian(	IElemDisc<TAlgebra>& elemDisc,
					typename geometry_traits<TElem>::const_iterator iterBegin,
					typename geometry_traits<TElem>::const_iterator iterEnd,
					int si,
					typename TAlgebra::matrix_type& J,
					const typename TAlgebra::vector_type& u,
					const TDoFDistribution& dofDistr,
					number time, number s_m, number s_a)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// get function group
	const FunctionGroup& fcts = elemDisc.get_function_group();

	// flag, wheather to use haning nodes as well
	bool useHanging = elemDisc.use_hanging();

	// local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u;
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_J;
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_J_temp;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!dofDistr.prepare_indices(refID, si, ind, useHanging))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IEDN_JACOBIAN))
		{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

	// adjust local algebra
	if(!useHanging)
	{
		loc_u.set_indices(ind);
		loc_J.set_indices(ind, ind);
		loc_J_temp.set_indices(ind, ind);
	}

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element loop.\n"); return false;}

	// loop over all elements
	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		dofDistr.update_indices(elem, ind, useHanging);

		// adapt local algebra
		if(useHanging)
		{
			loc_u.set_indices(ind);
			loc_J.set_indices(ind, ind);
			loc_J_temp.set_indices(ind, ind);
		}

		// read local values of u
		loc_u.read_values(u);

		// reset local matrix and rhs
		loc_J.set(0.0);

		// prepare element
		if(!elemDisc.prepare_element(elem, loc_u, ind))
			{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element.\n"); return false;}

		// Assemble JA
		loc_J_temp.set(0.0);
		if(!elemDisc.assemble_JA(loc_J_temp, loc_u, time))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix.\n"); return false;}
		loc_J += loc_J_temp * s_a;

		// Assemble JM
		loc_J_temp.set(0.0);
		if(!elemDisc.assemble_JM(loc_J_temp, loc_u, time))
			{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Mass Matrix.\n"); return false;}
		loc_J += loc_J_temp * s_m;

		// send local to global matrix
		J.add(loc_J);
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop.\n"); return false;}

	return true;
}

//////////////////////////////////
// Assemble Defect
//////////////////////////////////

template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleDefect(	IElemDisc<TAlgebra>& elemDisc,
				typename geometry_traits<TElem>::const_iterator iterBegin,
				typename geometry_traits<TElem>::const_iterator iterEnd,
				int si,
				typename TAlgebra::vector_type& d,
				const typename TAlgebra::vector_type& u,
				const TDoFDistribution& dofDistr,
				number time, number s_m, number s_a)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// get function group
	const FunctionGroup& fcts = elemDisc.get_function_group();

	// flag, wheather to use haning nodes as well
	bool useHanging = elemDisc.use_hanging();

	// local indices
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u;
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_d;
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_d_temp;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!dofDistr.prepare_indices(refID, si, ind, useHanging))
		{UG_LOG("ERROR in AssembleDefect: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IEDN_DEFECT))
		{UG_LOG("ERROR in AssembleDefect: Cannot set geometric object type.\n"); return false;}

	// adjust local algebra
	if(!useHanging)
	{
		loc_u.set_indices(ind);
		loc_d.set_indices(ind);
		loc_d_temp.set_indices(ind);
	}

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleDefect: Cannot prepare element loop.\n"); return false;}

	// loop over all elements
	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		dofDistr.update_indices(elem, ind, useHanging);

		// adjust local algebra
		if(useHanging)
		{
			loc_u.set_indices(ind);
			loc_d.set_indices(ind);
			loc_d_temp.set_indices(ind);
		}

		// read values
		loc_u.read_values(u);

		// reset local matrix and rhs
		loc_d.set(0.0);

		// prepare element
		if(!elemDisc.prepare_element(elem, loc_u, ind))
			{UG_LOG("ERROR in AssembleDefect: Cannot prepare element.\n"); return false;}

		// Assemble A
		loc_d_temp.set(0.0);
		if(!elemDisc.assemble_A(loc_d_temp, loc_u, time))
			{UG_LOG("ERROR in AssembleDefect: Cannot assemble local Stiffness Defect.\n"); return false;}
		loc_d += loc_d_temp * s_a;

		// Assemble M
		loc_d_temp.set(0.0);
		if(!elemDisc.assemble_M(loc_d_temp, loc_u, time))
			{UG_LOG("ERROR in AssembleDefect: Cannot assemble local Mass Defect.\n"); return false;}
		loc_d += loc_d_temp * s_m;

		// Assemble f
		loc_d_temp.set(0.0);
		if(!elemDisc.assemble_f(loc_d_temp, time))
			{UG_LOG("ERROR in AssembleDefect: Cannot assemble local Right-Hand Side.\n"); return false;}
		loc_d -= loc_d_temp * s_a;

		// send local to global matrix
		d.add(loc_d);
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleDefect: Cannot finish element loop.\n"); return false;}

	return true;
}

//////////////////////////////////
// Assemble linear
//////////////////////////////////

template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleLinear(	IElemDisc<TAlgebra>& elemDisc,
				typename geometry_traits<TElem>::const_iterator iterBegin,
				typename geometry_traits<TElem>::const_iterator iterEnd,
				int si,
				typename TAlgebra::matrix_type& mat,
				typename TAlgebra::vector_type& rhs,
				const typename TAlgebra::vector_type& u,
				const TDoFDistribution& dofDistr)
{
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// get function group
	const FunctionGroup& fcts = elemDisc.get_function_group();

	// flag, wheather to use haning nodes as well
	bool useHanging = elemDisc.use_hanging();

	// local indices and local algebra
	LocalIndices ind;
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u;
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_rhs;
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_mat;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices for elem type
	if(!dofDistr.prepare_indices(refID, si, ind, useHanging))
		{UG_LOG("ERROR in AssembleLinear: Cannot prepare indices.\n"); return false;}

	// set elem type in elem disc
	if(!elemDisc.set_geometric_object_type(refID, IEDN_LINEAR))
		{UG_LOG("ERROR in AssembleLinear: Cannot set geometric object type.\n"); return false;}

	// adjust local algebra
	if(!useHanging)
	{
		loc_u.set_indices(ind);
		loc_rhs.set_indices(ind);
		loc_mat.set_indices(ind, ind);
	}

	// prepare loop
	if(!elemDisc.prepare_element_loop())
		{UG_LOG("ERROR in AssembleLinear: Cannot prepare element loop.\n"); return false;}

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		dofDistr.update_indices(elem, ind, useHanging);

		// adjust local algebra
		if(useHanging)
		{
			loc_u.set_indices(ind);
			loc_rhs.set_indices(ind);
			loc_mat.set_indices(ind, ind);
		}

		// read local values of u
		loc_u.read_values(u);

		// reset local matrix and rhs
		loc_mat.set(0.0);
		loc_rhs.set(0.0);

		// prepare element data
		if(!elemDisc.prepare_element(elem, loc_u, ind))
			{UG_LOG("ERROR in AssembleLinear: Cannot prepare element.\n"); return false;}

		// assemble stiffness matrix for inner elements
		if(!elemDisc.assemble_JA(loc_mat, loc_u))
			{UG_LOG("ERROR in AssembleLinear: Cannot assemble local Stiffness Matrix.\n"); return false;}

		// assemble rhs for inner elements
		if(!elemDisc.assemble_f(loc_rhs))
			{UG_LOG("ERROR in AssembleLinear: Cannot assemble local Right-Hand Side.\n"); return false;}

		// send local to global (matrix and rhs)
		mat.add(loc_mat);
		rhs.add(loc_rhs);
	}

	// finish element loop
	if(!elemDisc.finish_element_loop())
		{UG_LOG("ERROR in AssembleLinear: Cannot finish element loop.\n"); return false;}

	return true;
};

} // end namespace ug


#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_ASSEMBLE_UTIL__ */
