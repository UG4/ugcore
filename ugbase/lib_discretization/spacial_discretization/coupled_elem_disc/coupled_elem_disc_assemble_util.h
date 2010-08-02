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
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleJacobian(	CoupledSystem<TDiscreteFunction, TAlgebra>& cplElemDisc,
					typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					int si,
					typename TAlgebra::matrix_type& J,
					const TDiscreteFunction& u,
					const FunctionGroup& fcts,
					number time, number s_m, number s_a)
{
	using std::vector;
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// clear identification ( may be skipped if identification is the same for all GeomObject types)
	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	ElemDataContainer.clear_identification();

	// local indices
	LocalIndices ind;

	// set function group
	ind.set_function_group(fcts);

	// prepare local indices
	if(!u.prepare_indices(refID, si, ind))
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
	ElemDataContainer.identify_exports();

	// local algebra
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u(ind);
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_J(ind, ind);
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_J_temp(ind, ind);

	// TODO: I think this is nor needed anymore
	// allocate matrix for off-diagonal part (induced by coupling)
/*	for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.system(sys);

		loc_J_coupl[sys].resize( system.num_imports());
		for(size_t i = 0; i < system.num_imports(); ++i)
		{
			DataImportItem* Imp = system.import(i);

			for(size_t s = 0; s < Imp->num_sys(); ++s)
			{
				loc_J_coupl[sys][i].resize(Imp->num_sys());

				UG_LOG("Coupling system " << sys << " ("<< num_sh[sys]<< ") with "<< s << " (" << num_sh[Imp->sys_id(s)] <<").\n");
				loc_J_coupl[sys][i][s].resize(num_sh[sys], num_sh[Imp->sys_id(s)]);
			}
		}
	}*/

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		u.update_indices(elem, ind);

		// read local values of u
		const typename TAlgebra::vector_type& u_vec = u.get_vector();
		loc_u.read_values(u_vec);

		// loop all systems
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			// prepare export
			// Exports: set_local_solution(elem, local_vector_type& u, local_index_type& glob_ind)
			if(!cplElemDisc.system(sys).prepare_element(elem, loc_u, ind))
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
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			if(!cplElemDisc.system(sys).assemble_JA(loc_J_temp, loc_u, time))
				{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Stiffness Matrix for system " << sys << ".\n"); return false;}
		}
		loc_J += loc_J_temp * s_a;

		// Assemble JM
		loc_J_temp.set(0.0);
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			if(!cplElemDisc.system(sys).assemble_JM(loc_J_temp, loc_u, time))
				{UG_LOG("ERROR in AssembleJacobian: Cannot assemble local Mass Matrix for system " << sys << ".\n"); return false;}
		}
		loc_J += loc_J_temp * s_m;

		// send local to global matrix
		J.add(loc_J);

		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			// compute matrix entries induced by coupling
/*			for(size_t i = 0; i < system.num_imports(); ++i)
			{
				DataImportItem* Imp = system.import(i);

				for(size_t r = 0; r < Imp->num_sys(); ++r)
				{
					loc_J_coupl[sys][i][r].set(0.0);

					// TODO: currently we assume, that imports are only in the stiffness part
					if(!Imp->add_offdiagonal(loc_J_coupl[sys][i][r], r, s_a))
						{UG_LOG("Offdiagonal coupling went wrong.\n"); return false;}

					//UG_LOG("Adding Local Couling: "<< loc_J_coupl[sys][i][r] << "\n");

					// add coupling
					J.add(loc_J_coupl[sys][i][r], glob_ind[sys], glob_ind[Imp->sys_id(r)]);
				}
			}
*/		}
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
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleDefect(	CoupledSystem<TDiscreteFunction, TAlgebra>& cplElemDisc,
				typename geometry_traits<TElem>::iterator iterBegin,
				typename geometry_traits<TElem>::iterator iterEnd,
				int si,
				typename TAlgebra::vector_type& d,
				const TDiscreteFunction& u,
				const FunctionGroup& fcts,
				number time, number s_m, number s_a)
{
	using std::vector;
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// TODO : Handle identification
	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	//ElemDataContainer.clear_identification();

	// local indices
	LocalIndices ind;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices
	if(!u.prepare_indices(refID, si, ind))
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
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		u.update_indices(elem, ind);

		// read local values of u
		const typename TAlgebra::vector_type& u_vec = u.get_vector();
		loc_u.read_values(u_vec);

		// prepare element
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			const SubFunctionMap& map = cplElemDisc.sub_function_map(sys);
			for(size_t i = 0; i < map.num_fct(); ++i)
			{
				UG_LOG("Map: " << i << " -> " << map[i] << ".\n");
			}

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
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleLinear(	CoupledSystem<TDiscreteFunction, TAlgebra>& cplElemDisc,
				typename geometry_traits<TElem>::iterator iterBegin,
				typename geometry_traits<TElem>::iterator iterEnd,
				int si,
				typename TAlgebra::matrix_type& mat,
				typename TAlgebra::vector_type& rhs,
				const TDiscreteFunction& u,
				const FunctionGroup& fcts)
{
	using std::vector;
	const ReferenceObjectID refID = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	ElemDataContainer.clear_identification();

	// local indices
	LocalIndices ind;

	// set functions
	ind.set_function_group(fcts);

	// prepare local indices
	if(!u.prepare_indices(refID, si, ind))
		{UG_LOG("ERROR in AssembleLinear: Cannot prepare indices.\n"); return false;}

	for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.system(sys);

		// set elem type
		if(!system.set_geometric_object_type(refID, IEDN_LINEAR))
			{UG_LOG("ERROR in AssembleLinear: Cannot set geometric object type for system " << sys << ".\n"); return false;}

		// prepare for loop
		if(!system.prepare_element_loop()) return false;;
			{UG_LOG("ERROR in AssembleLinear: Cannot prepare element loop for system " << sys << ".\n"); return false;}
	}

	// local algebra
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_u(ind);
	LocalVector<typename TAlgebra::vector_type::entry_type> loc_rhs(ind);
	LocalMatrix<typename TAlgebra::matrix_type::entry_type> loc_mat(ind, ind);

	// identify again the exports to avoid double computation
	ElemDataContainer.identify_exports();

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get Element
		TElem* elem = *iter;

		// get global indices
		u.update_indices(elem, ind);

		// read local values of u
		const typename TAlgebra::vector_type& u_vec = u.get_vector();
		loc_u.read_values(u_vec);

		// loop all systems
		for(size_t sys = 0; sys < cplElemDisc.num_system(); ++sys)
		{
			//set SubFunctionMap
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			// prepare element
			cplElemDisc.system(sys).prepare_element(elem, loc_u, ind);
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
			ind.access_sub_function_group(cplElemDisc.sub_function_map(sys));

			// assemble stiffness matrix
			if(!cplElemDisc.system(sys).assemble_JA(loc_mat, loc_u))
				{UG_LOG("ERROR in AssembleLinear: Cannot assemble local Stiffness Matrix for system " << sys << ".\n"); return false;}

			// assemble rhs for inner elements
			if(!cplElemDisc.system(sys).assemble_f(loc_rhs))
				{UG_LOG("ERROR in AssembleLinear: Cannot assemble local Right-Hand side for system " << sys << ".\n"); return false;}
		}

		// send local to global (matrix and rhs)
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
