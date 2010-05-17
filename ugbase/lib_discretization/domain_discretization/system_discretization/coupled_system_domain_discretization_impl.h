/*
 * problem.cpp
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION_IMPL__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION_IMPL__

namespace ug {


template <typename TDiscreteFunction>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
clear_dirichlet_jacobian_defect(	geometry_traits<Vertex>::iterator iterBegin,
									geometry_traits<Vertex>::iterator iterEnd,
									matrix_type& J, vector_type& d,
									const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();
	grid_type& grid = *const_cast<grid_type*>(&u.get_domain().get_grid());

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	number val;
	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(IsBoundaryVertex2D(grid, vert))
		{
			for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
			{
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				if(u.get_multi_indices_of_geom_obj(vert, m_systems[sys]->fct(i), ind) != 1) assert(0);
				if(m_systems[sys]->boundary_value(val, corner, i, time))
				{
					glob_ind.push_back(ind[0]);
					dirichlet_vals.push_back(0.0);
				}
			}}
		}
	}

	if(d.set(dirichlet_vals, glob_ind) != true)
		return false;

	if(J.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
clear_dirichlet_defect(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						vector_type& d,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();
	grid_type& grid = *const_cast<grid_type*>(&u.get_domain().get_grid());

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	number val;
	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(IsBoundaryVertex2D(grid, vert))
		{
			for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
			{
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				if(u.get_multi_indices_of_geom_obj(vert, m_systems[sys]->fct(i), ind) != 1) assert(0);
				if(m_systems[sys]->boundary_value(val, corner, i, time))
				{
					glob_ind.push_back(ind[0]);
					dirichlet_vals.push_back(0.0);
				}
			}}
		}
	}

	if(d.set(dirichlet_vals, glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
clear_dirichlet_jacobian(	geometry_traits<Vertex>::iterator iterBegin,
							geometry_traits<Vertex>::iterator iterEnd,
							matrix_type& J,
							const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();
	grid_type& grid = *const_cast<grid_type*>(&u.get_domain().get_grid());

	local_index_type ind(1);
	local_index_type glob_ind;

	number val;
	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(IsBoundaryVertex2D(grid, vert))
		{
			for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
			{
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				if(u.get_multi_indices_of_geom_obj(vert, m_systems[sys]->fct(i), ind) != 1) assert(0);
				if(m_systems[sys]->boundary_value(val, corner, i, time))
				{
					glob_ind.push_back(ind[0]);
				}
			}}
		}
	}

	if(J.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
set_dirichlet_solution(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						vector_type& x,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();
	grid_type& grid = *const_cast<grid_type*>(&u.get_domain().get_grid());

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	number val;
	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(IsBoundaryVertex2D(grid, vert))
		{
			for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
			{
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				if(u.get_multi_indices_of_geom_obj(vert, m_systems[sys]->fct(i), ind) != 1) assert(0);
				if(m_systems[sys]->boundary_value(val, corner, i, time))
				{
					glob_ind.push_back(ind[0]);
					dirichlet_vals.push_back(val);
				}
			}}
		}
	}

	if(x.set(dirichlet_vals, glob_ind) != true)
		return false;

	return true;
}


template <typename TDiscreteFunction>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
set_dirichlet_linear(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						matrix_type& mat, vector_type& rhs,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();
	grid_type& grid = *const_cast<grid_type*>(&u.get_domain().get_grid());

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	number val;
	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(IsBoundaryVertex2D(grid, vert))
		{
			for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
			{
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				if(u.get_multi_indices_of_geom_obj(vert, m_systems[sys]->fct(i), ind) != 1) assert(0);
				if(m_systems[sys]->boundary_value(val, corner, i, time))
				{
					dirichlet_vals.push_back(val);
					glob_ind.push_back(ind[0]);

				}
			}}
		}
	}

	if(rhs.set(dirichlet_vals, glob_ind) != true)
		return false;

	if(mat.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction>
template <typename TElem>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian_defect(	typename geometry_traits<TElem>::iterator iterBegin,
							typename geometry_traits<TElem>::iterator iterEnd,
							matrix_type& J, vector_type& d, const discrete_function_type& u,
							number time, number s_m, number s_a)
{
	// element
	TElem* elem = NULL;

	std::size_t num_systems = m_systems.size();
	std::vector<uint> num_sh(num_systems);

	std::vector<local_index_type> glob_ind(num_systems);
	std::vector<local_vector_type> loc_u(num_systems);
	std::vector<local_vector_type> loc_d(num_systems);
	std::vector<local_vector_type> loc_d_temp(num_systems);
	std::vector<local_matrix_type> loc_J(num_systems);
	std::vector<local_matrix_type> loc_J_temp(num_systems);

	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		// get total number of DoF's on this element type 'TElem'
		num_sh[sys] = m_systems[sys]->num_sh(elem);

		// allocate memory for local rhs and local Stiffness matrix_type
		glob_ind[sys].resize(num_sh[sys]);
		loc_u[sys].resize(num_sh[sys]);
		loc_d[sys].resize(num_sh[sys]);
		loc_d_temp[sys].resize(num_sh[sys]);
		loc_J[sys].resize(num_sh[sys], num_sh[sys]);
		loc_J_temp[sys].resize(num_sh[sys], num_sh[sys]);

		m_systems[sys]->prepare_element_loop(elem);
	}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// loop all systems
		for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset index offset
			uint offset = 0;

			// get local indices and fill local matrix pattern
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				offset += u.get_multi_indices(elem, m_systems[sys]->fct(i), glob_ind[sys], offset);
			}
			UG_ASSERT(offset == num_sh[sys], offset << " indices are read in, but we have " << num_sh[sys] << " dofs on this element.\n");

			// get local dof values of all unknowns
			u.get_dof_values(loc_u[sys], glob_ind[sys]);

			// prepare export
			m_systems[sys]->prepare_element(elem, loc_u[sys], glob_ind[sys]);
		}

		m_ElemDataContainer.compute(true);

		// loop all systems
		for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset local matrix and rhs
			loc_d[sys].set(0.0);
			loc_J[sys].set(0.0);

			// Assemble JA
			loc_J_temp[sys].set(0.0);
			m_systems[sys]->assemble_element_JA(elem, loc_J_temp[sys], loc_u[sys], time);
			loc_J[sys] += loc_J_temp[sys] * s_a;

			// Assemble JM
			loc_J_temp[sys].set(0.0);
			m_systems[sys]->assemble_element_JM(elem, loc_J_temp[sys], loc_u[sys], time);
			loc_J[sys] += loc_J_temp[sys] * s_m;

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
			J.add(loc_J[sys], glob_ind[sys], glob_ind[sys]);
			d.add(loc_d[sys], glob_ind[sys]);
		}
	}

	// finish element loop for each system
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		m_systems[sys]->finish_element_loop(elem);
	}

	return true;
}

// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction>
template <typename TElem>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					matrix_type& J, const discrete_function_type& u,
					number time, number s_m, number s_a)
{
	// element
	TElem* elem = NULL;

	// clear identification ( may be skipped if identification is the same for all GeomObject types)
	m_ElemDataContainer.clear_identification();

	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		// prepare element loop:
		// Imports: set_positions and num_eq
		// Exports: set new eval function + set num_shapes + sys
		m_systems[sys]->prepare_element_loop(elem);
	}

	// identify again the exports to avoid double computation
	m_ElemDataContainer.identify_exports();


	std::size_t num_sys = m_systems.size();
	std::vector<uint> num_sh(num_sys);

	std::vector<local_index_type> glob_ind(num_sys);
	std::vector<local_vector_type> loc_u(num_sys);
	std::vector<local_matrix_type> loc_J(num_sys);
	std::vector<local_matrix_type> loc_J_temp(num_sys);

	//std::vector<std::vector<local_matrix_type> > loc_J_coupl(num_sys);

	// allocating memory
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
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

		// allocate matrix for off-diagonal part (induced by coupling)
/*		loc_J_coupl[sys].resize(m_systems[sys]->num_imports());
		for(std::size_t i = 0; i < m_systems[sys]->num_imports(); ++i)
		{
			DataImport* Imp = m_systems[sys]->import(i);

			for(std::size_t s = 0; s < Imp->num_sys(); ++s)
			{
				loc_J_coupl[sys][i][s].resize(num_sh[sys], num_sh[Imp->sys(s)]);
			}
		}
*/	}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// loop all systems
		for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset index offset
			uint offset = 0;

			// get local indices and fill local matrix pattern
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
			{
				offset += u.get_multi_indices(elem, m_systems[sys]->fct(i), glob_ind[sys], offset);
			}
			UG_ASSERT(offset == num_sh[sys], offset << " indices are read in, but we have " << num_sh[sys] << " dofs on this element.\n");

			// get local dof values of all unknowns
			u.get_dof_values(loc_u[sys], glob_ind[sys]);

			// prepare export
			// Exports: set_local_solution(elem, local_vector_type& u, local_index_type& glob_ind)
			m_systems[sys]->prepare_element(elem, loc_u[sys], glob_ind[sys]);
		}

		m_ElemDataContainer.compute(true);

		for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
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

/*
			// compute matrix entries induced by coupling
			loc_J_coupl[sys].set(0.0);

			for(std::size_t i = 0; i < m_systems[sys]->num_imports(); ++i)
			{
				DataImport* Imp = m_systems[sys]->import(i);

				for(std::size_t r = 0; r < Imp->num_sys(); ++r)
				{
					s = Imp->sys(r);
					for(std::size_t k = 0; k < Imp->num_sh(r); ++k)
					{
						for(std::size_t j = 0; j < num_sh[sys]; ++j)
						{
							for(std::size_t ip = 0; ip < Imp->num_ip(); ++ip)
							{

								loc_J_coupl[sys][i][s](j, k) += Imp->defect(j, ip) * (*Imp)(s, ip, k);
							}
						}
					}
				}
				if(Imp->in_mass_part() != true) loc_J_coupl *= s_a;

				J.add(loc_J_coupl[sys][i][s], glob_ind[sys], glob_ind[s]);
				}
			}
*/
		}
	}

	// finish element loop for each system
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		m_systems[sys]->finish_element_loop(elem);
	}

	return true;
}

// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction>
template <typename TElem>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_defect(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					vector_type& d, const discrete_function_type& u,
					number time, number s_m, number s_a)
{
	// element
	TElem* elem = NULL;

	std::size_t num_systems = m_systems.size();
	std::vector<uint> num_sh(num_systems);

	std::vector<local_index_type> glob_ind(num_systems);
	std::vector<local_vector_type> loc_u(num_systems);
	std::vector<local_vector_type> loc_d(num_systems);
	std::vector<local_vector_type> loc_d_temp(num_systems);

	// prepare each systems
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		// get total number of DoF's on this element type 'TElem'
		num_sh[sys] = m_systems[sys]->num_sh(elem);

		// allocate memory for local rhs and local Stiffness matrix_type
		glob_ind[sys].resize(num_sh[sys]);
		loc_u[sys].resize(num_sh[sys]);
		loc_d[sys].resize(num_sh[sys]);
		loc_d_temp[sys].resize(num_sh[sys]);

		// prepare element loop
		m_systems[sys]->prepare_element_loop(elem);
	}

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// loop all systems
		for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset index offset
			uint offset = 0;

			// get local indices and fill local matrix pattern
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
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

		for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
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
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		m_systems[sys]->finish_element_loop(elem);
	}

	return true;
}



// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction>
template <typename TElem>
bool
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_linear(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
{
	// element
	TElem* elem = NULL;

	std::size_t num_systems = m_systems.size();
	std::vector<uint> num_sh(num_systems);

	std::vector<local_index_type> glob_ind(num_systems);
	std::vector<local_vector_type> loc_u(num_systems);
	std::vector<local_vector_type> loc_rhs(num_systems);
	std::vector<local_matrix_type> loc_mat(num_systems);

	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		// get total number of DoF's on this element type 'TElem'
		num_sh[sys] = m_systems[sys]->num_sh(elem);

		// allocate memory for local rhs and local Stiffness matrix_type
		glob_ind[sys].resize(num_sh[sys]);
		loc_u[sys].resize(num_sh[sys]);
		loc_rhs[sys].resize(num_sh[sys]);
		loc_mat[sys].resize(num_sh[sys], num_sh[sys]);

		// prepare for loop
		if(m_systems[sys]->prepare_element_loop(elem) != IPlugInReturn_OK) return false;;
	}

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		elem = *iter;

		// loop all systems
		for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset index offset
			uint offset = 0;

			// get local indices and fill local matrix pattern
			for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
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

		for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
		{
			// reset local matrix and rhs
			loc_mat[sys].set(0.0);
			loc_rhs[sys].set(0.0);

			// assemble stiffness matrix for inner elements
			if(m_systems[sys]->assemble_element_JA(elem, loc_mat[sys], loc_u[sys]) != IPlugInReturn_OK) return false;

			// assemble rhs for inner elements
			if(m_systems[sys]->assemble_element_f(elem, loc_rhs[sys]) != IPlugInReturn_OK) return false;

			// send local to global (matrix and rhs) [this is a virtual call]
			mat.add(loc_mat[sys], glob_ind[sys], glob_ind[sys]);
			rhs.add(loc_rhs[sys], glob_ind[sys]);
		}
	}

	// finish element loop
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		if(m_systems[sys]->finish_element_loop(elem) != IPlugInReturn_OK) return false;
	}

	return true;
}


template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	//return IAssemble_NOT_IMPLEMENTED;

	return this->assemble_jacobian(J, u, 0.0, 0.0, 1.0);
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_defect(vector_type& d, const discrete_function_type& u)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_solution(discrete_function_type& u)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
		{
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
			{
				UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
				return IAssemble_ERROR;
			}
		}
	}

	uint si = 0;
	UG_ASSERT(u.num_subsets() == 1, "Currently only one subset allowed.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes in solution vector ...");
	if(set_dirichlet_solution(u.template begin<Vertex>(si), u.template end<Vertex>(si), u.get_vector(), u) == false)
	{
		std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
{
	UG_DLOG(LIB_DISC_ASSEMBLE, 0, " ---- START: 'assemble_linear' (level = " << u.get_level() << ") ----\n");

	// check that Solution number 'nr' matches trial space required by Discretization
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
		{
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
			{
				UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
				return IAssemble_ERROR;
			}
		}
	}

	uint si = 0;
	UG_ASSERT(u.num_subsets() == 1, "Currently only one subset allowed.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Triangle>(si) << " Triangle(s) on Level " << u.get_level() << " ... ");
	if(assemble_linear<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), mat, rhs, u)==false)
	{
		UG_LOG("Error in assemble_linear, aborting.\n");
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "done.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Quadrilateral>(si) << " Quadrilateral(s) on Level " << u.get_level() << " ... ");
	if(assemble_linear<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), mat, rhs, u)==false)
	{
		UG_LOG("Error in assemble_linear, aborting.\n");
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "done.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes ..." );
	if(set_dirichlet_linear(u.template begin<Vertex>(si), u.template end<Vertex>(si), mat, rhs, u) == false)
	{
		UG_LOG("Error in assemble_linear, aborting.\n");
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "done.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 0, " ---- END: 'assemble_linear' ----\n");
	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
		{
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
			{
				UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
				return IAssemble_ERROR;
			}
		}
	}

	uint si = 0;
	UG_ASSERT(u.num_subsets() == 1, "Currently only one subset allowed.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Triangle>(si) << " Triangle(s) on Level " << u.get_level() << " ... ");
	if(assemble_jacobian_defect<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), J, d, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Quadrilateral>(si) << " Quadrilateral(s) on Level " << u.get_level() << " ... ");
	if(assemble_jacobian_defect<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), J, d, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes ...");
	if(clear_dirichlet_jacobian_defect(u.template begin<Vertex>(si),u.template end<Vertex>(si), J, d, u, time) == false)
	{
		std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u, number time, number s_m, number s_a)
{
	UG_DLOG(LIB_DISC_ASSEMBLE, 0, " ---- START: 'assembling_jacobian' (time = " << time << ", s_m = " << s_m << ", s_a = " << s_a << ", level = " << 	u.get_level() << ") ----\n");

	// check that Solution number 'nr' matches trial space required by Discretization
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
		{
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
			{
				UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
				return IAssemble_ERROR;
			}
		}
	}

	uint si = 0;
	UG_ASSERT(u.num_subsets() == 1, "Currently only one subset allowed.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Triangle>(si) << " Triangle(s) on Level " << u.get_level() << " ... ");
	if(assemble_jacobian<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), J, u, time, s_m, s_a)==false)
	{
		UG_LOG("Error in 'assemble_jacobian' while calling 'assemble_jacobian<Triangle>', aborting." << std::endl);
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Quadrilateral>(si) << " Quadrilateral(s) on Level " << u.get_level() << " ... ");
	if(assemble_jacobian<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), J, u, time, s_m, s_a)==false)
	{
		UG_LOG("Error in 'assemble_jacobian' while calling 'assemble_jacobian<Quadrilateral>', aborting." << std::endl);
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes ...");
	if(clear_dirichlet_jacobian(u.template begin<Vertex>(si),u.template end<Vertex>(si), J, u, time) == false)
	{
		UG_LOG("Error in 'assemble_jacobian' while calling 'clear_dirichlet_jacobian', aborting." << std::endl);
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	UG_DLOG(LIB_DISC_ASSEMBLE, 0, " ---- END: 'assembling_jacobian' ----\n");

	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_defect(vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
		{
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
			{
				UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
				return IAssemble_ERROR;
			}
		}
	}

	uint si = 0;
	UG_ASSERT(u.num_subsets() == 1, "Currently only one subset allowed.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Triangle>(si) << " Triangle(s) on Level " << u.get_level() << " ... ");
	if(assemble_defect<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), d, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Quadrilateral>(si) << " Quadrilateral(s) on Level " << u.get_level() << " ... ");
	if(assemble_defect<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), d, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes ...");
	if(clear_dirichlet_defect(u.template begin<Vertex>(si),u.template end<Vertex>(si), d, u, time) == false)
	{
		std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time, number s_m, number s_a)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_solution(discrete_function_type& u, number time)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(std::size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		for(uint i = 0; i < m_systems[sys]->num_fct(); i++)
		{
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
			{
				UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
				return IAssemble_ERROR;
			}
		}
	}

	uint si = 0;
	UG_ASSERT(u.num_subsets() == 1, "Currently only one subset allowed.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes in solution vector ...");
	if(set_dirichlet_solution(u.template begin<Vertex>(si), u.template end<Vertex>(si), u.get_vector(), u, time) == false)
	{
		std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);

	return IAssemble_OK;
}


}
#endif /*__H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION_IMPL__*/
