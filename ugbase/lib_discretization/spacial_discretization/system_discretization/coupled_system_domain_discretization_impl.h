/*
 * coupled_system_domain_discretization_impl.h
 *
 *  Created on: 25.06.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION_IMPL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION_IMPL__

namespace ug {


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

	size_t num_systems = m_systems.size();
	std::vector<size_t> num_sh(num_systems);

	std::vector<local_index_type> glob_ind(num_systems);
	std::vector<local_vector_type> loc_u(num_systems);
	std::vector<local_vector_type> loc_d(num_systems);
	std::vector<local_vector_type> loc_d_temp(num_systems);
	std::vector<local_matrix_type> loc_J(num_systems);
	std::vector<local_matrix_type> loc_J_temp(num_systems);

	for(size_t sys = 0; sys < m_systems.size(); ++sys)
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

			// prepare export
			m_systems[sys]->prepare_element(elem, loc_u[sys], glob_ind[sys]);
		}

		m_ElemDataContainer.compute(true);

		// loop all systems
		for(size_t sys = 0; sys < m_systems.size(); ++sys)
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
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
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
					size_t s = Imp->sys(r);

					for(size_t k = 0; k < Imp->num_sh(r); ++k)
					{
						for(size_t j = 0; j < num_sh[sys]; ++j)
						{
						//	for(size_t ip = 0; ip < Imp->num_ip(); ++ip)
							{
								//loc_J_coupl[sys][i][s](j, k) += Imp->lin_defect(j, ip) * (*Imp)(s, ip, k);
							}
						}
					}
				}
				//if(Imp->in_mass_part() != true) loc_J_coupl *= s_a;

				//J.add(loc_J_coupl[sys][i][s], glob_ind[sys], glob_ind[s]);
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
		if(m_systems[sys]->prepare_element_loop(elem) != IPlugInReturn_OK) return false;;
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
			if(m_systems[sys]->assemble_element_JA(elem, loc_mat[sys], loc_u[sys]) != IPlugInReturn_OK) return false;

			// assemble rhs for inner elements
			if(m_systems[sys]->assemble_element_f(elem, loc_rhs[sys]) != IPlugInReturn_OK) return false;

			// send local to global (matrix and rhs) [this is a virtual call]
			mat.add(loc_mat[sys], glob_ind[sys], glob_ind[sys]);
			rhs.add(loc_rhs[sys], glob_ind[sys]);
		}
	}

	// finish element loop
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
	{
		if(m_systems[sys]->finish_element_loop(elem) != IPlugInReturn_OK) return false;
	}

	return true;
}


template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, int si)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u, int si)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	//return IAssemble_NOT_IMPLEMENTED;

	return this->assemble_jacobian(J, u, si, 0.0, 0.0, 1.0);
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_defect(vector_type& d, const discrete_function_type& u, int si)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si)
{
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
		for(size_t i = 0; i < m_systems[sys]->num_fct(); i++)
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
				{UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");return IAssemble_ERROR;}

	if(assemble_linear<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), mat, rhs, u)==false)
		{UG_LOG("Error in assemble_linear, aborting.\n");return IAssemble_ERROR;}
	if(assemble_linear<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), mat, rhs, u)==false)
		{UG_LOG("Error in assemble_linear, aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, int si, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
		for(size_t i = 0; i < m_systems[sys]->num_fct(); i++)
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
				{UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n"); return IAssemble_ERROR;}

	if(assemble_jacobian_defect<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), J, d, u, time, s_m, s_a)==false)
		{UG_LOG("Error in assemble_linear, aborting.\n");return IAssemble_ERROR;}
	if(assemble_jacobian_defect<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), J, d, u, time, s_m, s_a)==false)
		{UG_LOG("Error in assemble_linear, aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u, int si, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
		for(size_t i = 0; i < m_systems[sys]->num_fct(); i++)
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
				{UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n"); return IAssemble_ERROR;}

	if(assemble_jacobian<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), J, u, time, s_m, s_a)==false)
		{UG_LOG("Error in 'assemble_jacobian' while calling 'assemble_jacobian<Triangle>', aborting.\n");return IAssemble_ERROR;}
	if(assemble_jacobian<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), J, u, time, s_m, s_a)==false)
		{UG_LOG("Error in 'assemble_jacobian' while calling 'assemble_jacobian<Quadrilateral>', aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_defect(vector_type& d, const discrete_function_type& u, int si, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(size_t sys = 0; sys < m_systems.size(); ++sys)
		for(size_t i = 0; i < m_systems[sys]->num_fct(); i++)
			if(u.get_local_shape_function_set_id(m_systems[sys]->fct(i)) != m_systems[sys]->local_shape_function_set(i))
				{UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n"); return IAssemble_ERROR;}

	if(assemble_defect<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), d, u, time, s_m, s_a)==false)
		{UG_LOG("Error in assemble_defect, aborting.\n");return IAssemble_ERROR;}
	if(assemble_defect<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), d, u, time, s_m, s_a)==false)
		{UG_LOG("Error in assemble_defect, aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction>
IAssembleReturn
CoupledSystemDomainDiscretization<TDiscreteFunction>::
assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si, number time, number s_m, number s_a)
{
	return IAssemble_NOT_IMPLEMENTED;
}


}
#endif /*__H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SYSTEM_DISCREZIZATION__COUPLED_SYSTEM_DOMAIN_DISCRETIZATION_IMPL__*/
