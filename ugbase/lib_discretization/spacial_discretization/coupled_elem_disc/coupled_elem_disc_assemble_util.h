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
#include "lib_discretization/reference_element/reference_elements.h"
#include "elem_disc_interface.h"

namespace ug {

//////////////////////////////////
// Assemble Jacobian
//////////////////////////////////

// assemble elements of type TElem in d dimensions
template <	typename TElem,
			typename TDiscreteFunction,
			typename TAlgebra>
bool
AssembleJacobian(	ICoupledElemDisc<TAlgebra>& elemDisc,
					typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					typename TAlgebra::matrix_type& J,
					const TDiscreteFunction& u,
					const std::vector<size_t>& u_comp,
					number time, number s_m, number s_a)
{
	// element
	TElem* elem = NULL;

	// clear identification ( may be skipped if identification is the same for all GeomObject types)
	m_ElemDataContainer.clear_identification();

	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		// prepare element loop:
		// Imports: set_positions and num_eq
		// Exports: set new eval function + set num_shapes + sys
		m_systems[sys]->prepare_element_loop(elem);
	}

	// identify again the exports to avoid double computation
	m_ElemDataContainer.identify_exports();


	size_t num_sys = m_systems.size();
	std::vector<size_t> num_sh(num_sys);

	std::vector<local_index_type> glob_ind(num_sys);
	std::vector<local_vector_type> loc_u(num_sys);
	std::vector<local_matrix_type> loc_J(num_sys);
	std::vector<local_matrix_type> loc_J_temp(num_sys);

	std::vector<std::vector<std::vector<local_matrix_type> > > loc_J_coupl(num_sys);

	// allocating memory
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		// get total number of DoF's on this element type 'TElem'
		num_sh[sys] = m_systems[sys]->num_sh(elem);

		// global indices of system
		glob_ind[sys].resize(num_sh[sys]);

		// local values of unknowns of system
		loc_u[sys].resize(num_sh[sys]);

		// local matrix of unknowns (diagonal part)
		loc_J[sys].resize(num_sh[sys], num_sh[sys]);
		loc_J_temp[sys].resize(num_sh[sys], num_sh[sys]);
	}

	// allocate matrix for off-diagonal part (induced by coupling)
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		loc_J_coupl[sys].resize(m_systems[sys]->num_imports());
		for(size_t i = 0; i < m_systems[sys]->num_imports(); ++i)
		{
			DataImportItem* Imp = m_systems[sys]->import(i);

			for(size_t s = 0; s < Imp->num_sys(); ++s)
			{
				loc_J_coupl[sys][i].resize(Imp->num_sys());
				loc_J_coupl[sys][i][s].resize(num_sh[sys], num_sh[Imp->sys(s)]);
			}
		}
	}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// loop all systems
		for(size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset index offset
			size_t offset = 0;

			// get local indices and fill local matrix pattern
			for(size_t i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				offset += u.get_multi_indices(elem, m_systems[sys]->fct(i), glob_ind[sys], offset);
			}
			UG_ASSERT(offset == num_sh[sys], offset << " indices are read in, but we have " << num_sh[sys] << " dofs on this element.\n");

			// get local dof values of all unknowns
			u.get_dof_values(loc_u[sys], glob_ind[sys]);

/*			UG_LOG("sys = " << sys << "\n");
			for(size_t i = 0; i < loc_u[sys].size(); ++i)
			{
				UG_LOG("i = " << i << ", u= " << loc_u[sys][i] << "ind = "<<  glob_ind[sys][i][0] << "\n");
			}
*/
			// prepare export
			// Exports: set_local_solution(elem, local_vector_type& u, local_index_type& glob_ind)
			m_systems[sys]->prepare_element(elem, loc_u[sys], glob_ind[sys]);
		}

		m_ElemDataContainer.compute(true);

		for(size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset local matrix and rhs
			loc_J[sys].set(0.0);

			// Assemble JA
			loc_J_temp[sys].set(0.0);
			m_systems[sys]->assemble_element_JA(elem, loc_J_temp[sys], loc_u[sys], time);
			loc_J[sys] += loc_J_temp[sys] * s_a;

			// Assemble JM
			loc_J_temp[sys].set(0.0);
			m_systems[sys]->assemble_element_JM(elem, loc_J_temp[sys], loc_u[sys], time);
			loc_J[sys] += loc_J_temp[sys] * s_m;

			// send local to global matrix
			J.add(loc_J[sys], glob_ind[sys], glob_ind[sys]);

			// compute matrix entries induced by coupling
			for(size_t i = 0; i < m_systems[sys]->num_imports(); ++i)
			{
				DataImportItem* Imp = m_systems[sys]->import(i);

				for(size_t r = 0; r < Imp->num_sys(); ++r)
				{
					loc_J_coupl[sys][i][r].set(0.0);

					// TODO: currently we assume, that imports are only in the stiffness part
					if(!Imp->add_offdiagonal(loc_J_coupl[sys][i][r], r, s_a))
						{UG_LOG("Offdiagonal coupling went wrong.\n"); return false;}

					//UG_LOG("Adding Local Couling: "<< loc_J_coupl[sys][i][r] << "\n");

					// add coupling
					J.add(loc_J_coupl[sys][i][r], glob_ind[sys], glob_ind[Imp->sys(r)]);
				}
			}
		}
	}

	// finish element loop for each system
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		m_systems[sys]->finish_element_loop(elem);
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
AssembleDefect(	ICoupledElemDisc<TAlgebra>& elemDisc,
				typename geometry_traits<TElem>::iterator iterBegin,
				typename geometry_traits<TElem>::iterator iterEnd,
				typename TAlgebra::vector_type& d,
				const TDiscreteFunction& u,
				const std::vector<size_t>& u_comp,
				number time, number s_m, number s_a)
{
	// element
	TElem* elem = NULL;

	size_t num_systems = m_systems.size();
	std::vector<size_t> num_sh(num_systems);

	std::vector<local_index_type> glob_ind(num_systems);
	std::vector<local_vector_type> loc_u(num_systems);
	std::vector<local_vector_type> loc_d(num_systems);
	std::vector<local_vector_type> loc_d_temp(num_systems);

	// prepare each systems
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		// prepare element loop
		m_systems[sys]->prepare_element_loop(elem);

		// get total number of DoF's on this element type 'TElem'
		num_sh[sys] = m_systems[sys]->num_sh(elem);

		// allocate memory for local rhs and local Stiffness matrix_type
		glob_ind[sys].resize(num_sh[sys]);
		loc_u[sys].resize(num_sh[sys]);
		loc_d[sys].resize(num_sh[sys]);
		loc_d_temp[sys].resize(num_sh[sys]);
	}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// loop all systems
		for(size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset index offset
			size_t offset = 0;

			// get local indices and fill local matrix pattern
			for(size_t i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				offset += u.get_multi_indices(elem, m_systems[sys]->fct(i), glob_ind[sys], offset);
			}
			UG_ASSERT(offset == num_sh[sys], offset << " indices are read in, but we have " << num_sh[sys] << " dofs on this element.\n");

			// get local dof values of all unknowns
			u.get_dof_values(loc_u[sys], glob_ind[sys]);

			// prepare element
			m_systems[sys]->prepare_element(elem, loc_u[sys], glob_ind[sys]);
		}

		m_ElemDataContainer.compute(false);

		for(size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset local matrix and rhs
			loc_d[sys].set(0.0);

			// Assemble A
			loc_d_temp[sys].set(0.0);
			m_systems[sys]->assemble_element_A(elem, loc_d_temp[sys], loc_u[sys], time);
			loc_d[sys] += loc_d_temp[sys] * s_a;

			// Assemble M
			loc_d_temp[sys].set(0.0);
			m_systems[sys]->assemble_element_M(elem, loc_d_temp[sys], loc_u[sys], time);
			loc_d[sys] += loc_d_temp[sys] * s_m;

			// Assemble f
			loc_d_temp[sys].set(0.0);
			m_systems[sys]->assemble_element_f(elem, loc_d_temp[sys], time);
			loc_d[sys] -= loc_d_temp[sys] * s_a;

			// send local to global matrix
			d.add(loc_d[sys], glob_ind[sys]);
		}
	}

	// finish element loop
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		m_systems[sys]->finish_element_loop(elem);
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
AssembleLinear(	ICoupledElemDisc<TAlgebra>& elemDisc,
				typename geometry_traits<TElem>::iterator iterBegin,
				typename geometry_traits<TElem>::iterator iterEnd,
				typename TAlgebra::matrix_type& mat,
				typename TAlgebra::vector_type& rhs,
				const TDiscreteFunction& u,
				const std::vector<size_t>& u_comp)
{
	// element
	TElem* elem = NULL;

	size_t num_systems = m_systems.size();
	std::vector<size_t> num_sh(num_systems);

	std::vector<local_index_type> glob_ind(num_systems);
	std::vector<local_vector_type> loc_u(num_systems);
	std::vector<local_vector_type> loc_rhs(num_systems);
	std::vector<local_matrix_type> loc_mat(num_systems);

	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		// get total number of DoF's on this element type 'TElem'
		num_sh[sys] = m_systems[sys]->num_sh(elem);

		// allocate memory for local rhs and local Stiffness matrix_type
		glob_ind[sys].resize(num_sh[sys]);
		loc_u[sys].resize(num_sh[sys]);
		loc_rhs[sys].resize(num_sh[sys]);
		loc_mat[sys].resize(num_sh[sys], num_sh[sys]);

		// prepare for loop
		if(m_systems[sys]->prepare_element_loop(elem) != true) return false;;
	}

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// loop all systems
		for(size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset index offset
			size_t offset = 0;

			// get local indices and fill local matrix pattern
			for(size_t i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				offset += u.get_multi_indices(elem, m_systems[sys]->fct(i), glob_ind[sys], offset);
			}
			UG_ASSERT(offset == num_sh[sys], offset << " indices are read in, but we have " << num_sh[sys] << " dofs on this element.\n");

			// get local dof values of all unknowns
			u.get_dof_values(loc_u[sys], glob_ind[sys]);  //<-- not needed, since linear, but maybe for linker

			// prepare element
			m_systems[sys]->prepare_element(elem, loc_u[sys], glob_ind[sys]);
		}

		m_ElemDataContainer.compute(true);

		for(size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset local matrix and rhs
			loc_mat[sys].set(0.0);
			loc_rhs[sys].set(0.0);

			// assemble stiffness matrix for inner elements
			if(m_systems[sys]->assemble_element_JA(elem, loc_mat[sys], loc_u[sys]) != true) return false;

			// assemble rhs for inner elements
			if(m_systems[sys]->assemble_element_f(elem, loc_rhs[sys]) != true) return false;

			// send local to global (matrix and rhs) [this is a virtual call]
			mat.add(loc_mat[sys], glob_ind[sys], glob_ind[sys]);
			rhs.add(loc_rhs[sys], glob_ind[sys]);
		}
	}

	// finish element loop
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		if(m_systems[sys]->finish_element_loop(elem) != true) return false;
	}

	return true;
}


//////////////////////////////////
//////////////////////////////////
// Subset assembling
//////////////////////////////////
//////////////////////////////////

//////////////////////////////////
// Jacobian subset assembling
//////////////////////////////////


} // end namespace ug

#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__COUPLED_ELEM_DISC__COUPLED_ELEM_DISC_ASSEMBLE_UTIL__*/
