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

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
class PlugInDomainDiscretization : public IDomainDiscretization<TAlgebra, TDiscreteFunction> {
	protected:
		// forward constants and typenames

		// type of discrete function
		typedef TDiscreteFunction discrete_function_type;

		// type of domain
		typedef typename TDiscreteFunction::domain_type domain_type;

		// type of grid used
		typedef typename domain_type::grid_type grid_type;

		// type of dof manager
		typedef typename TDiscreteFunction::dof_manager_type dof_manager_type;

		// type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

		// type of algebra
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename TAlgebra::matrix_type matrix_type;

		// type of local matrix
		typedef typename matrix_type::local_matrix_type local_matrix_type;

		// type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

		// type of algebra vector
		typedef typename vector_type::local_vector_type local_vector_type;

		// type of multi_index used in algebra
		typedef typename matrix_type::index_type index_type;

		// type of local index container
		typedef typename matrix_type::local_index_type local_index_type;

	public:
		PlugInDomainDiscretization(domain_type& domain, dof_manager_type& dofmanager, uint num_func, bool (*bndfct)(number&, position_type, uint, number),
									TElemDisc<domain_type, Triangle>& ElemDiscTriangle, TElemDisc<domain_type, Quadrilateral>& ElemDiscQuadrilateral
		)
		: m_bndfct(bndfct), m_domain(domain), m_pattern(dofmanager), _num_func(num_func),
		  m_ElemDiscTriangle(ElemDiscTriangle), m_ElemDiscQuadrilateral(ElemDiscQuadrilateral)
		{};

		// Assemble routines for time independent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, discrete_function_type& u, uint level = 0);
		IAssembleReturn assemble_jacobian(matrix_type& J, discrete_function_type& u, uint level = 0);
		IAssembleReturn assemble_defect(vector_type& d, discrete_function_type& u, uint level = 0);
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, uint level = 0);
		IAssembleReturn assemble_solution(discrete_function_type& u, uint level = 0);

		// Assemble routines for time dependent problems
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, discrete_function_type& u, number time, number s_m, number s_a, uint level = 0);
		IAssembleReturn assemble_jacobian(matrix_type& J, discrete_function_type& u, number time, number s_m, number s_a, uint level = 0);
		IAssembleReturn assemble_defect(vector_type& d, discrete_function_type& u, number time, number s_m, number s_a, uint level = 0);
		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs,  number time, number s_m, number s_a, uint level = 0);
		IAssembleReturn assemble_solution(discrete_function_type& u, number time, uint level = 0);


	protected:
		template <typename TElem>
		bool assemble_jacobian_defect	(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											TElemDisc<domain_type, TElem>& elemDisc,
											matrix_type& J, vector_type& d, discrete_function_type& u,
											number time, number s_m, number s_a);
		template <typename TElem>
		bool assemble_jacobian			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											TElemDisc<domain_type, TElem>& elemDisc,
											matrix_type& J, discrete_function_type& u,
											number time, number s_m, number s_a);
		template <typename TElem>
		bool assemble_defect			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											TElemDisc<domain_type, TElem>& elemDisc,
											vector_type& d, discrete_function_type& u,
											number time, number s_m, number s_a);

		template <typename TElem>
		bool assemble_linear			(	typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd,
											TElemDisc<domain_type, TElem>& elemDisc,
											matrix_type& mat, vector_type& rhs);

		inline bool clear_dirichlet_jacobian_defect(	geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, matrix_type& J, vector_type& d);
		inline bool clear_dirichlet_jacobian(			geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, matrix_type& J);
		inline bool clear_dirichlet_defect(				geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, vector_type& d);
		inline bool set_dirichlet_solution( 			geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, vector_type& x, number time = 0.0);
		inline bool set_dirichlet_linear(				geometry_traits<Vertex>::iterator iterBegin, geometry_traits<Vertex>::iterator iterEnd, matrix_type& mat, vector_type& rhs, number time = 0.0);

	protected:
		bool (*m_bndfct)(number&, position_type, uint, number);
		domain_type& m_domain;
		dof_manager_type& m_pattern;
		uint _num_func;

		TElemDisc<domain_type, Triangle>& m_ElemDiscTriangle;
		TElemDisc<domain_type, Quadrilateral>& m_ElemDiscQuadrilateral;
};



template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
inline
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
clear_dirichlet_jacobian_defect(	geometry_traits<Vertex>::iterator iterBegin,
									geometry_traits<Vertex>::iterator iterEnd,
									matrix_type& J, vector_type& d)
{
	typename domain_type::position_accessor_type aaPos = m_domain.get_position_accessor();
	grid_type& grid = m_domain.get_grid();

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
			for(uint fct = 0; fct < _num_func; fct++)
			{
				if(m_pattern.get_multi_indices_of_geom_obj(vert, fct, ind) != 1) assert(0);
				if(m_bndfct(val, corner, fct, 0.0))
				{
					glob_ind.push_back(ind[0]);
					dirichlet_vals.push_back(0.0);
				}
			}
		}
	}

	if(d.set(dirichlet_vals, glob_ind) != true)
		return false;

	if(J.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
inline
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
clear_dirichlet_defect(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						vector_type& d)
{
	typename domain_type::position_accessor_type aaPos = m_domain.get_position_accessor();
	grid_type& grid = m_domain.get_grid();

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
			for(uint fct = 0; fct < _num_func; fct++)
			{
				if(m_pattern.get_multi_indices_of_geom_obj(vert, fct, ind) != 1) assert(0);
				if(m_bndfct(val, corner, fct, 0.0))
				{
					glob_ind.push_back(ind[0]);
					dirichlet_vals.push_back(0.0);
				}
			}
		}
	}

	if(d.set(dirichlet_vals, glob_ind) != true)
		return false;

	return true;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
inline
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
clear_dirichlet_jacobian(	geometry_traits<Vertex>::iterator iterBegin,
							geometry_traits<Vertex>::iterator iterEnd,
							matrix_type& J)
{
	typename domain_type::position_accessor_type aaPos = m_domain.get_position_accessor();
	grid_type& grid = m_domain.get_grid();

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
			for(uint fct = 0; fct < _num_func; fct++)
			{
				if(m_pattern.get_multi_indices_of_geom_obj(vert, fct, ind) != 1) assert(0);
				if(m_bndfct(val, corner, fct, 0.0))
				{
					glob_ind.push_back(ind[0]);
				}
			}
		}
	}

	if(J.set_dirichlet_rows(glob_ind) != true)
		return false;

	return true;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
inline
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
set_dirichlet_solution(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						vector_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = m_domain.get_position_accessor();
	grid_type& grid = m_domain.get_grid();

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
			for(uint fct = 0; fct < _num_func; fct++)
			{
				if(m_pattern.get_multi_indices_of_geom_obj(vert, fct, ind) != 1) assert(0);
				if(m_bndfct(val, corner, fct, 0.0))
				{
					glob_ind.push_back(ind[0]);
					dirichlet_vals.push_back(val);
				}
			}
		}
	}

	if(u.set(dirichlet_vals, glob_ind) != true)
		return false;

	return true;
}


template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
inline
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
set_dirichlet_linear(	geometry_traits<Vertex>::iterator iterBegin,
						geometry_traits<Vertex>::iterator iterEnd,
						matrix_type& mat, vector_type& rhs, number time)
{
	typename domain_type::position_accessor_type aaPos = m_domain.get_position_accessor();
	grid_type& grid = m_domain.get_grid();

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
			for(uint fct = 0; fct < _num_func; fct++)
			{
				if(m_pattern.get_multi_indices_of_geom_obj(vert, fct, ind) != 1) assert(0);
				if(m_bndfct(val, corner, fct, time))

				dirichlet_vals.push_back(val);
				glob_ind.push_back(ind[0]);
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
template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_jacobian_defect(	typename geometry_traits<TElem>::iterator iterBegin,
							typename geometry_traits<TElem>::iterator iterEnd,
							TElemDisc<domain_type, TElem>& elemDisc,
							matrix_type& J, vector_type& d, discrete_function_type& u,
							number time, number s_m, number s_a)
{
	// get total number of DoF's on this element type 'TElem'
	const uint num_sh = elemDisc.num_sh();

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_d(num_sh);
	local_vector_type loc_d_temp(num_sh);
	local_vector_type loc_u(num_sh);
	local_matrix_type loc_J(num_sh, num_sh);
	local_matrix_type loc_J_temp(num_sh, num_sh);

	elemDisc.prepare_element_loop();

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem *elem = *iter;

		// get local indices and fill local matrix pattern
		for(uint fct = 0; fct < _num_func; fct++)
		{
			// reset local matrix and rhs
			loc_d.set(0.0);
			loc_J.set(0.0);

			m_pattern.get_multi_indices(elem, fct, glob_ind);
			uint level = m_pattern.get_grid().get_level(elem);
			u.get_dof_values(level, loc_u, glob_ind);

			elemDisc.prepare_element(elem);

			// Assemble JA
			loc_J_temp.set(0.0);
			elemDisc.assemble_element_JA(loc_J_temp, loc_u, time);
			loc_J += loc_J_temp * s_a;

			// Assemble JM
			loc_J_temp.set(0.0);
			elemDisc.assemble_element_JM(loc_J_temp, loc_u, time);
			loc_J += loc_J_temp * s_m;

			// Assemble A
			loc_d_temp.set(0.0);
			elemDisc.assemble_element_A(loc_d_temp, loc_u, time);
			loc_d += loc_d_temp * s_a;

			// Assemble M
			loc_d_temp.set(0.0);
			elemDisc.assemble_element_M(loc_d_temp, loc_u, time);
			loc_d += loc_d_temp * s_m;

			// Assemble f
			loc_d_temp.set(0.0);
			elemDisc.assemble_element_f(loc_d_temp, time);
			loc_d -= loc_d_temp * s_a;

			// send local to global matrix
			J.add(loc_J, glob_ind, glob_ind);
			d.add(loc_d, glob_ind);
		}
	}

	// finish element loop
	elemDisc.finish_element_loop();

	return true;
}

// assemble elements of type TElem in d dimensions
template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_jacobian(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					TElemDisc<domain_type, TElem>& elemDisc,
					matrix_type& J, discrete_function_type& u,
					number time, number s_m, number s_a)
{
	// get total number of DoF's on this element type 'TElem'
	const uint num_sh = elemDisc.num_sh();

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_u(num_sh);
	local_matrix_type loc_J(num_sh, num_sh);
	local_matrix_type loc_J_temp(num_sh, num_sh);

	elemDisc.prepare_element_loop();

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem *elem = *iter;

		// get local indices and fill local matrix pattern
		for(uint fct = 0; fct < _num_func; fct++)
		{
			// reset local matrix and rhs
			loc_J.set(0.0);

			m_pattern.get_multi_indices(elem, fct, glob_ind);
			uint level = m_pattern.get_grid().get_level(elem);
			u.get_dof_values(level, loc_u, glob_ind);

			elemDisc.prepare_element(elem);

			// Assemble JA
			loc_J_temp.set(0.0);
			elemDisc.assemble_element_JA(loc_J_temp, loc_u, time);
			loc_J += loc_J_temp * s_a;

			// Assemble JM
			loc_J_temp.set(0.0);
			elemDisc.assemble_element_JM(loc_J_temp, loc_u, time);
			loc_J += loc_J_temp * s_m;

			// send local to global matrix
			J.add(loc_J, glob_ind, glob_ind);
		}
	}

	// finish element loop
	elemDisc.finish_element_loop();

	return true;
}

// assemble elements of type TElem in d dimensions
template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_defect(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					TElemDisc<domain_type, TElem>& elemDisc,
					vector_type& d, discrete_function_type& u,
					number time, number s_m, number s_a)
{
	// get total number of DoF's on this element type 'TElem'
	const uint num_sh = elemDisc.num_sh();

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_d(num_sh);
	local_vector_type loc_d_temp(num_sh);
	local_vector_type loc_u(num_sh);

	elemDisc.prepare_element_loop();

	// loop over all elements
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem *elem = *iter;

		// loop over all functions
		for(uint fct = 0; fct < _num_func; fct++)
		{
			// reset local matrix and rhs
			loc_d.set(0.0);

			m_pattern.get_multi_indices(elem, fct, glob_ind);
			uint level = m_pattern.get_grid().get_level(elem);
			u.get_dof_values(level, loc_u, glob_ind);

			elemDisc.prepare_element(elem);

			// Assemble A
			loc_d_temp.set(0.0);
			elemDisc.assemble_element_A(loc_d_temp, loc_u, time);
			loc_d += loc_d_temp * s_a;

			// Assemble M
			loc_d_temp.set(0.0);
			elemDisc.assemble_element_M(loc_d_temp, loc_u, time);
			loc_d += loc_d_temp * s_m;

			// Assemble f
			loc_d_temp.set(0.0);
			elemDisc.assemble_element_f(loc_d_temp, time);
			loc_d -= loc_d_temp * s_a;

			// send local to global matrix
			d.add(loc_d, glob_ind);
		}
	}

	// finish element loop
	elemDisc.finish_element_loop();

	return true;
}



// assemble elements of type TElem in d dimensions
template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
template <typename TElem>
bool
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_linear(	typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					TElemDisc<domain_type, TElem>& elemDisc,
					matrix_type& mat, vector_type& rhs)
{
	// get total number of DoF's on this element type 'TElem'
	const uint num_sh = elemDisc.num_sh();

	// allocate memory for local rhs and local Stiffness matrix_type
	local_index_type glob_ind(num_sh);
	local_vector_type loc_rhs(num_sh);
	local_vector_type loc_u(num_sh);
	local_matrix_type loc_mat(num_sh, num_sh);

	// prepare for loop
	elemDisc.prepare_element_loop();

	// loop over all elements of type TElem
	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		// get Element
		TElem *elem = *iter;

		// loop over all functions
		for(uint fct = 0; fct < _num_func; fct++)
		{
			// reset local matrix and rhs
			loc_mat.set(0.0);
			loc_rhs.set(0.0);

			// get multiindices
			m_pattern.get_multi_indices(elem, fct, glob_ind);
			// ** u.get_dof_values(loc_u, glob_ind); ** <- not needed, since linear

			// prepare element data
			elemDisc.prepare_element(elem);

			// assemble stiffness matrix for inner elements
			elemDisc.assemble_element_JA(loc_mat, loc_u);

			// assemble stiffness matrix for bnd elements (bnd part only)
			elemDisc.assemble_element_JA_bnd(loc_mat, loc_u);

			// assemble rhs for inner elements
			elemDisc.assemble_element_f(loc_rhs);

			// assemble rhs for bnd elements (bnd part only)
			elemDisc.assemble_element_f_bnd(loc_rhs);

			// send local to global (matrix and rhs) [this is a virtual call]
			mat.add(loc_mat, glob_ind, glob_ind);
			rhs.add(loc_rhs, glob_ind);
		}
	}

	// finish element loop
	elemDisc.finish_element_loop();

	return true;
}


template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, discrete_function_type& u, uint level)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_jacobian(matrix_type& J, discrete_function_type& u, uint level)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_defect(vector_type& d, discrete_function_type& u, uint level)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_solution(discrete_function_type& u, uint level)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_linear(matrix_type& mat, vector_type& rhs, uint level)
{
	// check if Solution number 'nr' has linear trial space
	for(uint fct = 0; fct < _num_func; fct++)
	{
		if(m_pattern.get_local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
		{
			UG_LOG("Choosen Trial Space does not match this Discretization scheme. Abort discretization.\n");
			return IAssemble_ERROR;
		}
	}

	grid_type& mg = m_domain.get_grid();

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << mg.template num<Triangle>(level) << " Triangle(s) on Level " << level << " ... ");
	if(assemble_linear<Triangle>(mg.template begin<Triangle>(level), mg.template end<Triangle>(level), m_ElemDiscTriangle, mat, rhs)==false)
	{
		UG_LOG("Error in assemble_linear, aborting.\n");
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "done.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Assembling " << mg.template num<Quadrilateral>(level) << " Quadrilateral(s) on Level " << level << " ... ");
	if(assemble_linear<Quadrilateral>(mg.template begin<Quadrilateral>(level), mg.template end<Quadrilateral>(level), m_ElemDiscQuadrilateral, mat, rhs)==false)
	{
		UG_LOG("Error in assemble_linear, aborting.\n");
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "done.\n");

	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "Setting dirichlet nodes ..." );
	if(set_dirichlet_linear(mg.template begin<Vertex>(level),mg.template end<Vertex>(level), mat, rhs) == false)
	{
		UG_LOG("Error in assemble_linear, aborting.\n");
		return IAssemble_ERROR;
	}
	UG_DLOG(LIB_DISC_ASSEMBLE, 1, "done.\n");

	return IAssemble_OK;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, discrete_function_type& u, number time, number s_m, number s_a, uint level)
{
	// check if Solution number 'nr' has linear trial space
	for(uint fct = 0; fct < _num_func; fct++)
	{
		if(m_pattern.get_local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
		{
			std::cout << "Choosen Trial Space does not match this Discretization scheme. Abort discretization." << std::endl;
			return IAssemble_ERROR;
		}
	}

	grid_type& mg = m_domain.get_grid();

	//std::cout << "Assembling " << mg.num<Triangle>(level) << " Triangle(s) on Level " << level << " ... " << std::flush;
	if(assemble_jacobian_defect<Triangle>(mg.template begin<Triangle>(level), mg.template end<Triangle>(level), m_ElemDiscTriangle, J, d, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	//std::cout << "Assembling " << mg.num<Quadrilateral>(level) << " Quadrilateral(s) on Level " << level << " ... " << std::flush;
	if(assemble_jacobian_defect<Quadrilateral>(mg.template begin<Quadrilateral>(level), mg.template end<Quadrilateral>(level), m_ElemDiscQuadrilateral, J, d, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	//std::cout << "Setting dirichlet nodes ..." << std::flush;
	if(clear_dirichlet_jacobian_defect(mg.template begin<Vertex>(level),mg.template end<Vertex>(level), J, d) == false)
	{
		std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	return IAssemble_OK;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_jacobian(matrix_type& J, discrete_function_type& u, number time, number s_m, number s_a, uint level)
{
	// check if Solution number 'nr' has linear trial space
	for(uint fct = 0; fct < _num_func; fct++)
	{
		if(m_pattern.get_local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
		{
			std::cout << "Choosen Trial Space does not match this Discretization scheme. Abort discretization." << std::endl;
			return IAssemble_ERROR;
		}
	}

	grid_type& mg = m_domain.get_grid();

	//std::cout << "Assembling " << mg.num<Triangle>(level) << " Triangle(s) on Level " << level << " ... " << std::flush;
	if(assemble_jacobian<Triangle>(mg.template begin<Triangle>(level), mg.template end<Triangle>(level), m_ElemDiscTriangle, J, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	//std::cout << "Assembling " << mg.num<Quadrilateral>(level) << " Quadrilateral(s) on Level " << level << " ... " << std::flush;
	if(assemble_jacobian<Quadrilateral>(mg.template begin<Quadrilateral>(level), mg.template end<Quadrilateral>(level), m_ElemDiscQuadrilateral, J, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	//std::cout << "Setting dirichlet nodes ..." << std::flush;
	if(clear_dirichlet_jacobian(mg.template begin<Vertex>(level),mg.template end<Vertex>(level), J) == false)
	{
		std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	return IAssemble_OK;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_defect(vector_type& d, discrete_function_type& u, number time, number s_m, number s_a, uint level)
{
	// check if Solution number 'nr' has linear trial space
	for(uint fct = 0; fct < _num_func; fct++)
	{
		if(m_pattern.get_local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
		{
			std::cout << "Choosen Trial Space does not match this Discretization scheme. Abort discretization." << std::endl;
			return IAssemble_ERROR;
		}
	}

	grid_type& mg = m_domain.get_grid();

	//std::cout << "Assembling " << mg.num<Triangle>(level) << " Triangle(s) on Level " << level << " ... " << std::flush;
	if(assemble_defect<Triangle>(mg.template begin<Triangle>(level), mg.template end<Triangle>(level),m_ElemDiscTriangle, d, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	//std::cout << "Assembling " << mg.num<Quadrilateral>(level) << " Quadrilateral(s) on Level " << level << " ... " << std::flush;
	if(assemble_defect<Quadrilateral>(mg.template begin<Quadrilateral>(level), mg.template end<Quadrilateral>(level), m_ElemDiscQuadrilateral, d, u, time, s_m, s_a)==false)
	{
		std::cout << "Error in assemble_linear, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	//std::cout << "Setting dirichlet nodes ..." << std::flush;
	if(clear_dirichlet_defect(mg.template begin<Vertex>(level),mg.template end<Vertex>(level), d) == false)
	{
		std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	return IAssemble_OK;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_linear(matrix_type& mat, vector_type& rhs, number time, number s_m, number s_a, uint level)
{
	return IAssemble_NOT_IMPLEMENTED;
}

template <typename TAlgebra, typename TDiscreteFunction, template<typename TDomain, typename TElem > class TElemDisc>
IAssembleReturn
PlugInDomainDiscretization<TAlgebra, TDiscreteFunction, TElemDisc>::
assemble_solution(discrete_function_type& u, number time, uint level)
{
	// check if Solution number 'nr' has linear trial space
	for(uint fct = 0; fct < _num_func; fct++)
	{
		if(m_pattern.get_local_shape_function_set_id(fct) != LSFS_LAGRANGEP1)
		{
			std::cout << "Choosen Trial Space does not match this Discretization scheme. Abort discretization." << std::endl;
			return IAssemble_ERROR;
		}
	}

	grid_type& mg = m_domain.get_grid();

	//std::cout << "Setting dirichlet nodes in solution vector ..." << std::flush;
	if(set_dirichlet_solution(mg.template begin<Vertex>(level),mg.template end<Vertex>(level), u.get_vector(level), time) == false)
	{
		std::cout << "Error in set_dirichlet_nodes, aborting." << std::endl;
		return IAssemble_ERROR;
	}
	//std::cout << "Done." << std::endl;

	return IAssemble_OK;
}




}
#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__PLUG_IN_SPACIAL_DISCRETIZATION__ */
