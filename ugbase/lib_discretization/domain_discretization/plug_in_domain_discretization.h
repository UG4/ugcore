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
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/domain_discretization/domain_discretization_interface.h"
#include "lib_discretization/domain_discretization/disc_helper/fvgeom.h"

namespace ug {

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
class PlugInDomainDiscretization : public IDomainDiscretization<typename TDiscreteFunction::algebra_type, TDiscreteFunction> {
	protected:
		// forward constants and typenames

		// type of discrete function
		typedef TDiscreteFunction discrete_function_type;

		// type of domain
		typedef typename TDiscreteFunction::domain_type domain_type;

		// type of grid used
		typedef typename domain_type::grid_type grid_type;

		// type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

		// type of dof manager
		typedef typename TDiscreteFunction::dof_manager_type dof_manager_type;

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
		PlugInDomainDiscretization(TElemDisc<domain_type, algebra_type>& elemDisc) :
		  m_elemDisc(elemDisc)
		{};

		// Assemble routines for time independent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u);
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u);
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u);
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u);
		IAssembleReturn assemble_solution(discrete_function_type& u);

		// Assemble routines for time dependent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a);
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u, number time, number s_m, number s_a);
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a);
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time, number s_m, number s_a);
		IAssembleReturn assemble_solution(discrete_function_type& u, number time);

		virtual std::size_t num_fct() const
		{
			return m_elemDisc.num_fct();
		}

		virtual bool is_dirichlet(int s, uint fct)
		{
			//UG_LOG("Checking Dirichlet value");
			return m_elemDisc.is_dirichlet(s, fct);
		}

	protected:
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

		bool clear_dirichlet_jacobian_defect(	geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, uint fct, int s,matrix_type& J,vector_type& d, const discrete_function_type& u, number time = 0.0);
		bool clear_dirichlet_jacobian(			geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, uint fct, int s,matrix_type& J, const discrete_function_type& u, number time = 0.0);
		bool clear_dirichlet_defect(			geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, uint fct, int s,vector_type& d, const discrete_function_type& u, number time = 0.0);
		bool set_dirichlet_solution( 			geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, uint fct, int s,vector_type& x, const discrete_function_type& u, number time = 0.0);
		bool set_dirichlet_linear(				geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, uint fct, int s,matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time = 0.0);

		TElemDisc<domain_type, algebra_type>& m_elemDisc;
};



template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
clear_dirichlet_jacobian_defect(	geometry_traits<Vertex>::iterator iterBegin,
									geometry_traits<Vertex>::iterator iterEnd,
									uint fct, int s,
									matrix_type& J, vector_type& d,
									const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(u.get_multi_indices_of_geom_obj(vert, m_elemDisc.fct(fct), ind) != 1) assert(0);

		glob_ind.push_back(ind[0]);
		dirichlet_vals.push_back(0.0);
	}

	if(d.set(dirichlet_vals, glob_ind) != true)
		return false;

	if(J.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
clear_dirichlet_defect(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						uint fct, int s,
						vector_type& d,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(u.get_multi_indices_of_geom_obj(vert, m_elemDisc.fct(fct), ind) != 1) assert(0);

		glob_ind.push_back(ind[0]);
		dirichlet_vals.push_back(0.0);
	}

	if(d.set(dirichlet_vals, glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
clear_dirichlet_jacobian(	geometry_traits<Vertex>::iterator iterBegin,
							geometry_traits<Vertex>::iterator iterEnd,
							uint fct, int s,
							matrix_type& J,
							const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();

	local_index_type ind(1);
	local_index_type glob_ind;

	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(u.get_multi_indices_of_geom_obj(vert, m_elemDisc.fct(fct), ind) != 1) assert(0);

		glob_ind.push_back(ind[0]);
	}

	if(J.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
set_dirichlet_solution(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						uint fct, int s,
						vector_type& x,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_domain().get_position_accessor();

	local_index_type ind(1);
	local_index_type glob_ind;
	local_vector_type dirichlet_vals;

	number val;
	position_type corner;

	for(geometry_traits<Vertex>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		VertexBase *vert = *iter;
		corner = aaPos[vert];

		if(u.get_multi_indices_of_geom_obj(vert, m_elemDisc.fct(fct), ind) != 1) assert(0);

		if(m_elemDisc.boundary_value(val, corner, time, s, fct))
		{
			glob_ind.push_back(ind[0]);
			dirichlet_vals.push_back(val);
		}
	}

	if(x.set(dirichlet_vals, glob_ind) != true)
		return false;

	return true;
}


template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
set_dirichlet_linear(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						uint fct, int s,
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
			for(uint i = 0; i < m_elemDisc.num_fct(); i++)
			{
				if(u.get_multi_indices_of_geom_obj(vert, m_elemDisc.fct(fct), ind) != 1) assert(0);
				if(m_elemDisc.boundary_value(val, corner, time, s, fct))
				{
					dirichlet_vals.push_back(val);
					glob_ind.push_back(ind[0]);

				}
			}
		}
	}

	if(rhs.set(dirichlet_vals, glob_ind) != true)
		return false;

	if(mat.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

// assemble elements of type TElem in d dimensions
template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_jacobian_defect(	typename geometry_traits<TElem>::iterator iterBegin,
							typename geometry_traits<TElem>::iterator iterEnd,
							matrix_type& J, vector_type& d, const discrete_function_type& u,
							number time, number s_m, number s_a)
{
	TElem* elem = NULL;

	// get total number of DoF's on this element type 'TElem'
	const uint num_sh = m_elemDisc.num_sh(elem);

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
		uint offset = 0;

		// loop over all functions
		for(uint i = 0; i < m_elemDisc.num_fct(); i++)
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
template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_jacobian(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					matrix_type& J, const discrete_function_type& u,
					number time, number s_m, number s_a)
{
	TElem* elem = NULL;

	// get total number of DoF's on this element type 'TElem'
	const uint num_sh = m_elemDisc.num_sh(elem);

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
		uint offset = 0;

		// loop over all functions
		for(uint i = 0; i < m_elemDisc.num_fct(); i++)
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
template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_defect(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					vector_type& d, const discrete_function_type& u,
					number time, number s_m, number s_a)
{
	TElem* elem = NULL;

	// get total number of DoF's on this element type 'TElem'
	const uint num_sh = m_elemDisc.num_sh(elem);

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
		uint offset = 0;

		// loop over all functions
		for(uint i = 0; i < m_elemDisc.num_fct(); i++)
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
template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_linear(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
{
	TElem* elem = NULL;

	// get total number of DoF's on this element type 'TElem'
	const uint num_sh = m_elemDisc.num_sh(elem);

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
		uint offset = 0;

		// loop over all functions
		for(uint i = 0; i < m_elemDisc.num_fct(); i++)
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


template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u)
{
	// TODO: This is a costly quick hack, compute matrices directly (without time assembling) !
	//return IAssemble_NOT_IMPLEMENTED;

	return this->assemble_jacobian(J, u, 0.0, 0.0, 1.0);
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_defect(vector_type& d, const discrete_function_type& u)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_solution(discrete_function_type& u)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(uint i = 0; i < m_elemDisc.num_fct(); i++)
	{
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
		{
			UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
			return IAssemble_ERROR;
		}
	}

	for(uint i = 0; i < m_elemDisc.num_bnd_subsets(0); ++i)
	{
		int si = m_elemDisc.bnd_subset(0, i);

		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes in solution vector ...");
		for(uint fct = 0; fct < m_elemDisc.num_fct(); ++fct)
		{
			if(!m_elemDisc.is_dirichlet(si, fct)) continue;
			if(set_dirichlet_solution(u.template begin<Vertex>(si), u.template end<Vertex>(si), fct, si, u.get_vector(), u) == false)
			{
				std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
				return IAssemble_ERROR;
			}
		}
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);
	}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
{
	UG_DLOG(LIB_DISC_ASSEMBLE, 0, " ---- START: 'assemble_linear' (level = " << u.get_level() << ") ----\n");

	// check that Solution number 'nr' matches trial space required by Discretization
	for(uint i = 0; i < m_elemDisc.num_fct(); i++)
	{
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
		{
			UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
			return IAssemble_ERROR;
		}
	}

	for(uint i = 0; i < m_elemDisc.num_elem_subsets(2); ++i)
	{
		int si = m_elemDisc.elem_subset(2, i);

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
	}

	for(uint i = 0; i < m_elemDisc.num_bnd_subsets(0); ++i)
	{
		int si = m_elemDisc.bnd_subset(0, i);
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes in solution vector ...");
		for(uint fct = 0; fct < m_elemDisc.num_fct(); ++fct)
		{
			if(!m_elemDisc.is_dirichlet(si, fct)) continue;
			UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes ..." );
			if(set_dirichlet_linear(u.template begin<Vertex>(si), u.template end<Vertex>(si), fct, si, mat, rhs, u) == false)
			{
				UG_LOG("Error in assemble_linear, aborting.\n");
				return IAssemble_ERROR;
			}
			UG_DLOG(LIB_DISC_ASSEMBLE, 1, "done.\n");
		}
	}

	UG_DLOG(LIB_DISC_ASSEMBLE, 0, " ---- END: 'assemble_linear' ----\n");
	return IAssemble_OK;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(uint i = 0; i < m_elemDisc.num_fct(); i++)
	{
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
		{
			UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
			return IAssemble_ERROR;
		}
	}

	for(uint i = 0; i < m_elemDisc.num_elem_subsets(2); ++i)
	{
		int si = m_elemDisc.elem_subset(2, i);

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
	}

	for(uint i = 0; i < m_elemDisc.num_bnd_subsets(0); ++i)
	{
		int si = m_elemDisc.bnd_subset(0, i);
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes in solution vector ...");
		for(uint fct = 0; fct < m_elemDisc.num_fct(); ++fct)
		{
			if(!m_elemDisc.is_dirichlet(si, fct)) continue;

			UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes ...");
			if(clear_dirichlet_jacobian_defect(u.template begin<Vertex>(si),u.template end<Vertex>(si), fct, si, J, d, u, time) == false)
			{
				std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
				return IAssemble_ERROR;
			}
			UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);
		}
	}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u, number time, number s_m, number s_a)
{
	UG_DLOG(LIB_DISC_ASSEMBLE, 0, " ---- START: 'assembling_jacobian' (time = " << time << ", s_m = " << s_m << ", s_a = " << s_a << ", level = " << 	u.get_level() << ") ----\n");

	// check that Solution number 'nr' matches trial space required by Discretization
	for(uint i = 0; i < m_elemDisc.num_fct(); i++)
	{
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
		{
			UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
			return IAssemble_ERROR;
		}
	}

	for(uint i = 0; i < m_elemDisc.num_elem_subsets(2); ++i)
	{
		int si = m_elemDisc.elem_subset(2, i);

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
	}

	for(uint i = 0; i < m_elemDisc.num_bnd_subsets(0); ++i)
	{
		int si = m_elemDisc.bnd_subset(0, i);

		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes in solution vector ...");
		for(uint fct = 0; fct < m_elemDisc.num_fct(); ++fct)
		{
			if(!m_elemDisc.is_dirichlet(si, fct)) continue;
			UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes ...");
			if(clear_dirichlet_jacobian(u.template begin<Vertex>(si),u.template end<Vertex>(si), fct, si, J, u, time) == false)
			{
				UG_LOG("Error in 'assemble_jacobian' while calling 'clear_dirichlet_jacobian', aborting." << std::endl);
				return IAssemble_ERROR;
			}
			UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done." << std::endl);
		}
	}

	UG_DLOG(LIB_DISC_ASSEMBLE, 0, " ---- END: 'assembling_jacobian' ----\n");

	return IAssemble_OK;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_defect(vector_type& d, const discrete_function_type& u, number time, number s_m, number s_a)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(uint i = 0; i < m_elemDisc.num_fct(); i++)
	{
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
		{
			UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
			return IAssemble_ERROR;
		}
	}

	for(uint i = 0; i < m_elemDisc.num_elem_subsets(2); ++i)
	{
		int si = m_elemDisc.elem_subset(2, i);

		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Triangle>(si) << " Triangle(s) on Level " << u.get_level() << " ... ");
		if(assemble_defect<Triangle>(u.template begin<Triangle>(si), u.template end<Triangle>(si), d, u, time, s_m, s_a)==false)
		{
			std::cout << "Error in assemble_linear, aborting." << std::endl;
			return IAssemble_ERROR;
		}
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done. \n");

		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << u.template num<Quadrilateral>(si) << " Quadrilateral(s) on Level " << u.get_level() << " ... ");
		if(assemble_defect<Quadrilateral>(u.template begin<Quadrilateral>(si), u.template end<Quadrilateral>(si), d, u, time, s_m, s_a)==false)
		{
			std::cout << "Error in assemble_linear, aborting." << std::endl;
			return IAssemble_ERROR;
		}
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done. \n");
	}

	for(uint i = 0; i < m_elemDisc.num_bnd_subsets(0); ++i)
	{
		int si = m_elemDisc.bnd_subset(0, i);
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting " << u.template num<Vertex>(si) << " dirichlet nodes in Subset " << si <<" ...");
		for(uint fct = 0; fct < m_elemDisc.num_fct(); ++fct)
		{
			if(!m_elemDisc.is_dirichlet(si, fct)) continue;

			if(clear_dirichlet_defect(u.template begin<Vertex>(si),u.template end<Vertex>(si), fct, si, d, u, time) == false)
			{
				std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
				return IAssemble_ERROR;
			}
		}
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done.\n");
	}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, number time, number s_m, number s_a)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TDiscreteFunction, template<typename TDomain, typename TAlgebra > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TDiscreteFunction, TElemDisc>::
assemble_solution(discrete_function_type& u, number time)
{
	// check that Solution number 'nr' matches trial space required by Discretization
	for(uint i = 0; i < m_elemDisc.num_fct(); i++)
	{
		if(u.get_local_shape_function_set_id(m_elemDisc.fct(i)) != m_elemDisc.local_shape_function_set(i))
		{
			UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
			return IAssemble_ERROR;
		}
	}

	for(uint i = 0; i < m_elemDisc.num_bnd_subsets(0); ++i)
	{
		int si = m_elemDisc.bnd_subset(0, i);
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting " << u.template num<Vertex>(si) << " dirichlet values in solution vector ...");

		for(uint fct = 0; fct < m_elemDisc.num_fct(); ++fct)
		{
			if(!m_elemDisc.is_dirichlet(si, fct)) continue;
			if(set_dirichlet_solution(u.template begin<Vertex>(si), u.template end<Vertex>(si), fct, si, u.get_vector(), u, time) == false)
			{
				std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
				return IAssemble_ERROR;
			}
		}
		UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Done.\n");
	}

	return IAssemble_OK;
}




}
#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__PLUG_IN_SPACIAL_DISCRETIZATION__ */
