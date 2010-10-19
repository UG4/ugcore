/*
 * coupled_elem_disc_assemble_util.h
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_ASSEMBLE_UTIL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_algebra/lib_algebra.h"

// intern headers
#include "../../reference_element/reference_element.h"
#include "./coupled_elem_disc_interface.h"
#include "./coupled_system.h"
#include "lib_discretization/common/function_group.h"


namespace ug {

//////////////////////////////////
// Assemble Jacobian
//////////////////////////////////

// assemble elements of type TElem in d dimensions
template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleJacobian(	CoupledSystem<TAlgebra>& cplElemDisc,
					typename geometry_traits<TElem>::const_iterator iterBegin,
					typename geometry_traits<TElem>::const_iterator iterEnd,
					int si,
					typename TAlgebra::matrix_type& J,
					const typename TAlgebra::vector_type& u,
					const TDoFDistribution& dofDistr,
					number time, number s_m, number s_a)
{
	using std::vector;
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// get function group
	const FunctionGroup& fcts = cplElemDisc.get_function_group();

	// clear identification ( may be skipped if identification is the same for all GeomObject types)
	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	//ElemDataContainer.clear_identification();

	// local indices
	LocalIndices rowInd;
	LocalIndices colInd;

	// set function group
	rowInd.set_function_group(fcts);
	colInd.set_function_group(fcts);

	// prepare local indices
	if(!dofDistr.prepare_indices(refID, si, rowInd))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}
	// prepare local indices
	if(!dofDistr.prepare_indices(refID, si, colInd))
		{UG_LOG("ERROR in AssembleJacobian: Cannot prepare indices.\n"); return false;}

	for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
	{
		// get system
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.system(sys);

		// set elem type
		if(!system.set_geometric_object_type(refID, IEDN_JACOBIAN))
			{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type for system " << sys << ".\n"); return false;}

		// prepare element loop:
		// Imports: set_positions and num_eq
		// Exports: set new eval function + set num_shapes + sys
		if(!system.prepare_element_loop())
			{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element loop for system " << sys << ".\n"); return false;}
	}

	// identify again the exports to avoid double computation
	//ElemDataContainer.identify_exports();

	// local algebra
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u(colInd);
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_J(rowInd, colInd);
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_J_temp(rowInd, colInd);

	// loop over all elements
	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		dofDistr.update_indices(elem, rowInd);
		dofDistr.update_indices(elem, colInd);

		// read local values of u
		loc_u.read_values(u);

		// loop all systems
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			rowInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));
			colInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			// prepare export
			// Exports: set_local_solution(elem, local_vector_type& u, local_index_type& glob_ind)
			if(!cplElemDisc.system(sys).prepare_element(elem, loc_u, colInd))
				{UG_LOG("ERROR in AssembleJacobian: Cannot prepare element for system " << sys << ".\n"); return false;}
		}

		ElemDataContainer.compute(true);

		// reset local matrix and rhs
		loc_J.set(0.0);

		// Assemble JA
		loc_J_temp.set(0.0);
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			rowInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));
			colInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			if(!cplElemDisc.system(sys).assemble_JA(loc_J_temp, loc_u, time))
				{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix for system " << sys << ".\n"); return false;}
		}
		loc_J += loc_J_temp * s_a;

		// Assemble JM
		loc_J_temp.set(0.0);
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			rowInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));
			colInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			if(!cplElemDisc.system(sys).assemble_JM(loc_J_temp, loc_u, time))
				{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Mass Matrix for system " << sys << ".\n"); return false;}
		}
		loc_J += loc_J_temp * s_m;

		// compute matrix entries induced by coupling
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			for(size_t i = 0; i < cplElemDisc.system(sys).num_imports(); ++i)
			{
				DataImportItem* Imp = cplElemDisc.system(sys).import(i);

				for(size_t r = 0; r < Imp->num_sys(); ++r)
				{
					rowInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));
					colInd.access_sub_function_group(cplElemDisc.sub_function_map(r));

					// TODO: currently we assume, that imports are only in the stiffness part
					if(!Imp->add_offdiagonal(loc_J, r, s_a))
						{UG_LOG("Offdiagonal coupling went wrong.\n"); return false;}
				}
			}
		}

		// send local to global matrix
		rowInd.access_all();
		colInd.access_all();
		J.add(loc_J);
	}

	// finish element loop for each system
	for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
	{
		if(!cplElemDisc.system(sys).finish_element_loop())
			{UG_LOG("ERROR in AssembleJacobian: Cannot finish element loop for system " << sys << ".\n"); return false;}
	}

	return true;
}

//////////////////////////////////
// Assemble Defect
//////////////////////////////////

template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleDefect(	CoupledSystem<TAlgebra>& cplElemDisc,
				typename geometry_traits<TElem>::const_iterator iterBegin,
				typename geometry_traits<TElem>::const_iterator iterEnd,
				int si,
				typename TAlgebra::vector_type& d,
				const typename TAlgebra::vector_type& u,
				const TDoFDistribution& dofDistr,
				number time, number s_m, number s_a)
{
	using std::vector;
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// get function group
	const FunctionGroup& fcts = cplElemDisc.get_function_group();

	// TODO : Handle identification
	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	//ElemDataContainer.clear_identification();

	// local indices
	LocalIndices ind;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices
	if(!dofDistr.prepare_indices(refID, si, ind))
		{UG_LOG("ERROR in AssembleDefect: Cannot prepare indices.\n"); return false;}

	// local algebra
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u(ind);
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_d(ind);
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_d_temp(ind);

	// prepare each systems
	for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.system(sys);

		// set elem type
		if(!system.set_geometric_object_type(refID, IEDN_DEFECT))
			{UG_LOG("ERROR in AssembleDefect: Cannot set geometric object type for system " << sys << ".\n"); return false;}

		// prepare element loop
		if(!system.prepare_element_loop())
			{UG_LOG("ERROR in AssembleDefect: Cannot prepare element loop for system " << sys << ".\n"); return false;}
	}

	// TODO : Handle identification
	// identify again the exports to avoid double computation
	//ElemDataContainer.identify_exports();

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
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			if(!cplElemDisc.system(sys).prepare_element(elem, loc_u, ind))
				{UG_LOG("ERROR in AssembleDefect: Cannot prepare element for system " << sys << ".\n"); return false;}
		}

		ElemDataContainer.compute(false);

		// reset local matrix and rhs
		loc_d.set(0.0);

		// Assemble A
		loc_d_temp.set(0.0);
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			if(!cplElemDisc.system(sys).assemble_A(loc_d_temp, loc_u, time))
				{UG_LOG("ERROR in AssembleDefect: Cannot assemble local Stiffness Matrix for system " << sys << ".\n"); return false;}
		}
		loc_d += loc_d_temp * s_a;

		// Assemble M
		loc_d_temp.set(0.0);
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			if(!cplElemDisc.system(sys).assemble_M(loc_d_temp, loc_u, time))
				{UG_LOG("ERROR in AssembleDefect: Cannot assemble local Mass Matrix for system " << sys << ".\n"); return false;}
		}
		loc_d += loc_d_temp * s_m;

		// Assemble f
		loc_d_temp.set(0.0);
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			if(!cplElemDisc.system(sys).assemble_f(loc_d_temp, time))
				{UG_LOG("ERROR in AssembleDefect: Cannot assemble local Right-Hand Side for system " << sys << ".\n"); return false;}
		}
		loc_d -= loc_d_temp * s_a;

		// send local to global matrix
		ind.access_all();
		d.add(loc_d);
	}

	// finish element loop
	for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
	{
		if(!cplElemDisc.system(sys).finish_element_loop())
			{UG_LOG("ERROR in AssembleDefect: Cannot finish element loop for system " << sys << ".\n"); return false;}
	}

	return true;
}



//////////////////////////////////
// Assemble linear
//////////////////////////////////

template <	typename TElem,
			typename TDoFDistribution,
			typename TAlgebra>
bool
AssembleLinear(	CoupledSystem<TAlgebra>& cplElemDisc,
				typename geometry_traits<TElem>::const_iterator iterBegin,
				typename geometry_traits<TElem>::const_iterator iterEnd,
				int si,
				typename TAlgebra::matrix_type& mat,
				typename TAlgebra::vector_type& rhs,
				const typename TAlgebra::vector_type& u,
				const TDoFDistribution& dofDistr)
{
	using std::vector;
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// get function group
	const FunctionGroup& fcts = cplElemDisc.get_function_group();

	// clear identification ( may be skipped if identification is the same for all GeomObject types)
	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	//ElemDataContainer.clear_identification();

	// local indices
	LocalIndices rowInd;
	LocalIndices colInd;

	// set function group
	rowInd.set_function_group(fcts);
	colInd.set_function_group(fcts);

	// prepare local indices
	if(!dofDistr.prepare_indices(refID, si, rowInd))
		{UG_LOG("ERROR in AssembleLinear: Cannot prepare indices.\n"); return false;}
	// prepare local indices
	if(!dofDistr.prepare_indices(refID, si, colInd))
		{UG_LOG("ERROR in AssembleLinear: Cannot prepare indices.\n"); return false;}

	for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
	{
		// get system
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.system(sys);

		// set elem type
		if(!system.set_geometric_object_type(refID, IEDN_LINEAR))
			{UG_LOG("ERROR in AssembleLinear: Cannot set geometric object type for system " << sys << ".\n"); return false;}

		// prepare element loop:
		// Imports: set_positions and num_eq
		// Exports: set new eval function + set num_shapes + sys
		if(!system.prepare_element_loop())
			{UG_LOG("ERROR in AssembleLinear: Cannot prepare element loop for system " << sys << ".\n"); return false;}
	}

	// identify again the exports to avoid double computation
	//ElemDataContainer.identify_exports();

	// local algebra
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u(colInd);
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_rhs(colInd);
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_mat(rowInd, colInd);

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		dofDistr.update_indices(elem, rowInd);
		dofDistr.update_indices(elem, colInd);

		// read local values of u
		loc_u.read_values(u);

		// loop all systems
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			rowInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));
			colInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			// prepare export
			// Exports: set_local_solution(elem, local_vector_type& u, local_index_type& glob_ind)
			if(!cplElemDisc.system(sys).prepare_element(elem, loc_u, colInd))
				{UG_LOG("ERROR in AssembleLinear: Cannot prepare element for system " << sys << ".\n"); return false;}
		}

		ElemDataContainer.compute(true);

		// reset local matrix and rhs
		loc_mat.set(0.0);
		loc_rhs.set(0.0);

		// assemble stiffness matrix
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			rowInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));
			colInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			// assemble stiffness matrix
			if(!cplElemDisc.system(sys).assemble_JA(loc_mat, loc_u))
				{UG_LOG("ERROR in AssembleLinear: Cannot assemble local Stiffness Matrix for system " << sys << ".\n"); return false;}

			// assemble rhs for inner elements
			if(!cplElemDisc.system(sys).assemble_f(loc_rhs))
				{UG_LOG("ERROR in AssembleLinear: Cannot assemble local Right-Hand side for system " << sys << ".\n"); return false;}
		}

		// compute matrix entries induced by coupling
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			for(size_t i = 0; i < cplElemDisc.system(sys).num_imports(); ++i)
			{
				DataImportItem* Imp = cplElemDisc.system(sys).import(i);

				for(size_t r = 0; r < Imp->num_sys(); ++r)
				{
					rowInd.access_sub_function_group(cplElemDisc.sub_function_map(sys));
					colInd.access_sub_function_group(cplElemDisc.sub_function_map(r));

					if(!Imp->add_offdiagonal(loc_mat, r, 1.0))
						{UG_LOG("Offdiagonal coupling went wrong.\n"); return false;}
				}
			}
		}

		// send local to global (matrix and rhs)
		rowInd.access_all();
		colInd.access_all();
		mat.add(loc_mat);
		rhs.add(loc_rhs);
	}

	// finish element loop
	for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
	{
		if(!cplElemDisc.system(sys).finish_element_loop())
			{UG_LOG("ERROR in AssembleLinear: Cannot finish element loop for system " << sys << ".\n"); return false;}
	}

	return true;
}

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_ASSEMBLE_UTIL__*/
