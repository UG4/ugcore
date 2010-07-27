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
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	// clear identification ( may be skipped if identification is the same for all GeomObject types)
	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	ElemDataContainer.clear_identification();

	// number of systems
	const size_t num_sys = cplElemDisc.num_sys();

	// local indices
	vector<LocalIndices> ind(num_sys);

	for(size_t sys = 0; sys < num_sys; ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

		// set elem type
		if(!system.set_geometric_object_type(refID, IEDN_LINEAR))
			{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

		// set functions
		// TODO: set correct function group
		ind[sys].set_function_group(fcts);

		// prepare local indices
		if(!u.prepare_indices(refID, si, ind[sys]))
			{UG_LOG("Cannot prepare indices.\n"); return false;}

		// prepare element loop:
		// Imports: set_positions and num_eq
		// Exports: set new eval function + set num_shapes + sys
		system.prepare_element_loop();
	}

	// identify again the exports to avoid double computation
	ElemDataContainer.identify_exports();

	// local algebra
	vector<LocalVector<typename TAlgebra::vector_type::entry_type> > loc_u(num_sys);
	vector<LocalMatrix<typename TAlgebra::matrix_type::entry_type> > loc_J(num_sys);
	vector<LocalMatrix<typename TAlgebra::matrix_type::entry_type> > loc_J_temp(num_sys);

	// TODO: set indices in local algebra

	// TODO: set couple matrices

	// allocating memory
	for(size_t sys = 0; sys < num_sys; ++sys)
	{
		// local values of unknowns of system
		loc_u[sys].set_indices(ind[sys]);

		// local matrix of unknowns (diagonal part)
		loc_J[sys].set_indices(ind[sys], ind[sys]);
		loc_J_temp[sys].set_indices(ind[sys], ind[sys]);
	}

	// allocate matrix for off-diagonal part (induced by coupling)
	/*
	for(size_t sys = 0; sys < num_sys; ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

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
	}
	*/

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem* elem = *iter;

		// loop all systems
		for(size_t sys = 0; sys < cplElemDisc.num_sys(); ++sys)
		{
			ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

			// get global indices
			u.update_indices(elem, ind[sys]);

			// read local values of u
			const typename TAlgebra::vector_type& u_vec = u.get_vector();
			loc_u[sys].read_values(u_vec);

			// prepare export
			// Exports: set_local_solution(elem, local_vector_type& u, local_index_type& glob_ind)
			system.prepare_element(elem, loc_u[sys], ind[sys]);
		}

		ElemDataContainer.compute(true);

		for(size_t sys = 0; sys < num_sys; ++sys)
		{
			ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

			// reset local matrix and rhs
			loc_J[sys].set(0.0);

			// Assemble JA
			loc_J_temp[sys].set(0.0);
			system.assemble_JA(loc_J_temp[sys], loc_u[sys], time);
			loc_J[sys] += loc_J_temp[sys] * s_a;

			// Assemble JM
			loc_J_temp[sys].set(0.0);
			system.assemble_JM(loc_J_temp[sys], loc_u[sys], time);
			loc_J[sys] += loc_J_temp[sys] * s_m;

			// send local to global matrix
			J.add(loc_J[sys]);

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
	for(size_t sys = 0; sys < cplElemDisc.num_sys(); ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);
		system.finish_element_loop();
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
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	ElemDataContainer.clear_identification();

	// number of systems
	const size_t num_sys = cplElemDisc.num_sys();

	// local indices
	vector<LocalIndices> ind(num_sys);

	// local algebra
	vector<LocalVector<typename TAlgebra::vector_type::entry_type> > loc_u(num_sys);
	vector<LocalVector<typename TAlgebra::vector_type::entry_type> > loc_d(num_sys);
	vector<LocalVector<typename TAlgebra::vector_type::entry_type> > loc_d_temp(num_sys);

	// prepare each systems
	for(size_t sys = 0; sys < num_sys; ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

		// set elem type
		if(!system.set_geometric_object_type(refID, IEDN_LINEAR))
			{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

		// prepare element loop
		system.prepare_element_loop();

		// set functions
		// TODO: set correct function group
		ind[sys].set_function_group(fcts);

		// prepare local indices
		if(!u.prepare_indices(refID, si, ind[sys]))
			{UG_LOG("Cannot prepare indices.\n"); return false;}
	}

	// identify again the exports to avoid double computation
	ElemDataContainer.identify_exports();

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem* elem = *iter;

		// loop all systems
		for(size_t sys = 0; sys < num_sys; ++sys)
		{
			ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

			// get global indices
			u.update_indices(elem, ind[sys]);

			// read local values of u
			const typename TAlgebra::vector_type& u_vec = u.get_vector();
			loc_u[sys].read_values(u_vec);

			// prepare element
			system.prepare_element(elem, loc_u[sys], ind[sys]);
		}

		ElemDataContainer.compute(false);

		for(size_t sys = 0; sys < num_sys; ++sys)
		{
			ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

			// reset local matrix and rhs
			loc_d[sys].set(0.0);

			// Assemble A
			loc_d_temp[sys].set(0.0);
			system.assemble_A(loc_d_temp[sys], loc_u[sys], time);
			loc_d[sys] += loc_d_temp[sys] * s_a;

			// Assemble M
			loc_d_temp[sys].set(0.0);
			system.assemble_M(loc_d_temp[sys], loc_u[sys], time);
			loc_d[sys] += loc_d_temp[sys] * s_m;

			// Assemble f
			loc_d_temp[sys].set(0.0);
			system.assemble_f(loc_d_temp[sys], time);
			loc_d[sys] -= loc_d_temp[sys] * s_a;

			// send local to global matrix
			d.add(loc_d[sys]);
		}
	}

	// finish element loop
	for(size_t sys = 0; sys < cplElemDisc.num_sys(); ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);
		system.finish_element_loop();
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
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	const ReferenceObjectID refID = reference_element_type::REFERENCE_OBJECT_ID;

	// check if at least on element exist, else return
	if(iterBegin == iterEnd) return true;

	DataContainer& ElemDataContainer = cplElemDisc.get_elem_data_container();
	ElemDataContainer.clear_identification();

	// number of systems
	const size_t num_sys = cplElemDisc.num_sys();

	// local indices
	vector<LocalIndices> ind(num_sys);

	for(size_t sys = 0; sys < num_sys; ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

		// set elem type
		if(!system.set_geometric_object_type(refID, IEDN_LINEAR))
			{UG_LOG("ERROR in AssembleJacobian: Cannot set geometric object type.\n"); return false;}

		// set functions
		// TODO: set correct function group
		ind[sys].set_function_group(fcts);

		// prepare local indices
		if(!u.prepare_indices(refID, si, ind[sys]))
			{UG_LOG("Cannot prepare indices.\n"); return false;}

		// prepare for loop
		if(system.prepare_element_loop() != true) return false;;
	}

	// local algebra
	vector<LocalVector<typename TAlgebra::vector_type::entry_type> > loc_u(num_sys);
	vector<LocalVector<typename TAlgebra::vector_type::entry_type> > loc_rhs(num_sys);
	vector<LocalMatrix<typename TAlgebra::matrix_type::entry_type> > loc_mat(num_sys);

	// identify again the exports to avoid double computation
	ElemDataContainer.identify_exports();

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem* elem = *iter;

		// loop all systems
		for(size_t sys = 0; sys < num_sys; ++sys)
		{
			ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

			// get global indices
			u.update_indices(elem, ind[sys]);

			// read local values of u
			const typename TAlgebra::vector_type& u_vec = u.get_vector();
			loc_u[sys].read_values(u_vec);

			// prepare element
			system.prepare_element(elem, loc_u[sys], ind[sys]);
		}

		ElemDataContainer.compute(true);

		for(size_t sys = 0; sys < num_sys; ++sys)
		{
			ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);

			// reset local matrix and rhs
			loc_mat[sys].set(0.0);
			loc_rhs[sys].set(0.0);

			// assemble stiffness matrix for inner elements
			if(system.assemble_JA(loc_mat[sys], loc_u[sys]) != true) return false;

			// assemble rhs for inner elements
			if(system.assemble_f(loc_rhs[sys]) != true) return false;

			// send local to global (matrix and rhs) [this is a virtual call]
			mat.add(loc_mat[sys]);
			rhs.add(loc_rhs[sys]);
		}
	}

	// finish element loop
	for(size_t sys = 0; sys < cplElemDisc.num_sys(); ++sys)
	{
		ICoupledElemDisc<TAlgebra>& system = cplElemDisc.sys(sys);
		if(system.finish_element_loop() != true) return false;
	}

	return true;
}

} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_ASSEMBLE_UTIL__*/
