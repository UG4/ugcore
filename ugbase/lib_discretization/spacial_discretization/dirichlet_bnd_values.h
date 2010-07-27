/*
 * dirichlet_bnd_values.h
 *
 *  Created on: 08.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DIRICHLET_BND_VALUES__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DIRICHLET_BND_VALUES__

#include "lib_discretization/spacial_discretization/domain_discretization_interface.h"

namespace ug{

enum BND_TYPE {
	BND_TYPE_NONE = 0,
	BND_TYPE_DIRICHLET,
	BND_TYPE_NEUMANN
};

template <	typename TDiscreteFunction,
			int ref_dim>
class DirichletBNDValues : public IDirichletBoundaryValues<TDiscreteFunction> {
	public:
		// type of discrete function
		typedef TDiscreteFunction discrete_function_type;

		// type of domain
		typedef typename TDiscreteFunction::domain_type domain_type;

		// type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

		// type of algebra
		typedef typename TDiscreteFunction::algebra_type algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename TDiscreteFunction::vector_type vector_type;

		// type of local index container
		typedef typename TDiscreteFunction::local_index_type local_index_type;

		// type of multi index vector
		typedef typename TDiscreteFunction::multi_index_vector_type multi_index_vector_type;

	public:
		typedef bool (*Boundary_fct)(number&, const position_type&, number);

		DirichletBNDValues(size_t fct, Boundary_fct bnd_fct, domain_type& domain) :
		  m_bndfct(bnd_fct), m_fct(fct), m_domain(domain)
		  {}

		IAssembleReturn clear_dirichlet_jacobian(matrix_type& J, const discrete_function_type& u, int si, number time = 0.0);
		IAssembleReturn clear_dirichlet_defect(vector_type& d, const discrete_function_type& u, int si,number time = 0.0);
		IAssembleReturn set_dirichlet_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si, number time = 0.0);

		IAssembleReturn set_dirichlet_solution(discrete_function_type& u, int si, number time = 0.0);

	protected:
		template <typename TElem>
		bool clear_dirichlet_jacobian(			typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, matrix_type& J,
												const discrete_function_type& u, number time = 0.0);

		template <typename TElem>
		bool clear_dirichlet_defect(			typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, vector_type& d,
												const discrete_function_type& u, number time = 0.0);

		template <typename TElem>
		bool set_dirichlet_solution( 			typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, vector_type& x,
												const discrete_function_type& u, number time = 0.0);

		template <typename TElem>
		bool set_dirichlet_linear(				typename geometry_traits<TElem>::iterator iterBegin,
												typename geometry_traits<TElem>::iterator iterEnd,
												size_t fct, int si, matrix_type& mat, vector_type& rhs,
												const discrete_function_type& u, number time = 0.0);

		bool is_dirichlet(size_t fct) {return fct == m_fct;}

	protected:
		Boundary_fct m_bndfct;
		size_t m_fct;
		domain_type& m_domain;
};


template <typename TDiscreteFunction, int ref_dim>
IAssembleReturn
DirichletBNDValues<TDiscreteFunction, ref_dim>::
clear_dirichlet_jacobian(matrix_type& J, const discrete_function_type& u, int si, number time)
{
	if(clear_dirichlet_jacobian<VertexBase>(u.template begin<VertexBase>(si),u.template end<VertexBase>(si), m_fct, si, J, u, time) == false)
		{UG_LOG("Error in 'assemble_jacobian' while calling 'clear_dirichlet_jacobian', aborting.\n"); return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, int ref_dim>
IAssembleReturn
DirichletBNDValues<TDiscreteFunction, ref_dim>::
clear_dirichlet_defect(vector_type& d, const discrete_function_type& u, int si, number time)
{
	if(clear_dirichlet_defect<VertexBase>(u.template begin<VertexBase>(si),u.template end<VertexBase>(si), m_fct, si, d, u, time) == false)
		{UG_LOG("Error in set_dirichlet_nodes, aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, int ref_dim>
IAssembleReturn
DirichletBNDValues<TDiscreteFunction, ref_dim>::
set_dirichlet_solution(discrete_function_type& u, int si, number time)
{
	if(set_dirichlet_solution<VertexBase>(u.template begin<VertexBase>(si), u.template end<VertexBase>(si), m_fct, si, u.get_vector(), u, time) == false)
		{UG_LOG("Error in set_dirichlet_nodes, aborting.\n"); return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDiscreteFunction, int ref_dim>
IAssembleReturn
DirichletBNDValues<TDiscreteFunction, ref_dim>::
set_dirichlet_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u, int si, number time)
{
	if(set_dirichlet_linear<VertexBase>(u.template begin<VertexBase>(si), u.template end<VertexBase>(si), m_fct, si, mat, rhs, u) == false)
		{UG_LOG("Error in assemble_linear, aborting.\n"); return IAssemble_ERROR;}

	return IAssemble_OK;
}


template <typename TDiscreteFunction, int ref_dim>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction, ref_dim>::
clear_dirichlet_defect(	typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						size_t fct, int si,
						vector_type& d,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_approximation_space().get_domain().get_position_accessor();
	multi_index_vector_type multInd;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		corner = aaPos[elem];

		if(u.template get_multi_indices_of_geom_obj<TElem>(elem, fct, multInd) != 1)
			return false;

		if(m_bndfct(val, corner, time))
		{
			BlockRef(d[multInd[0][0]], multInd[0][1]) = 0.0;
		}
	}

	return true;
}

template <typename TDiscreteFunction, int ref_dim>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction, ref_dim>::
clear_dirichlet_jacobian(	typename geometry_traits<TElem>::iterator iterBegin,
							typename geometry_traits<TElem>::iterator iterEnd,
							size_t fct, int si,
							matrix_type& J,
							const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_approximation_space().get_domain().get_position_accessor();
	multi_index_vector_type multInd;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		corner = aaPos[elem];

		if(u.template get_multi_indices_of_geom_obj<TElem>(elem, fct, multInd) != 1)
			return false;

		if(m_bndfct(val, corner, time))
		{
			SetDirichletRow(J, multInd[0][0], multInd[0][1]);
		}

	}
	return true;
}

template <typename TDiscreteFunction, int ref_dim>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction, ref_dim>::
set_dirichlet_solution(	typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						size_t fct, int si,
						vector_type& x,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_approximation_space().get_domain().get_position_accessor();

	multi_index_vector_type multInd;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		// TODO: if TElem != Vertex we have to do something else
		corner = aaPos[elem];

		if(u.template get_multi_indices_of_geom_obj<TElem>(elem, fct, multInd) != 1)
			return false;

		if(m_bndfct(val, corner, time))
		{
			BlockRef(x[multInd[0][0]], multInd[0][1]) = val;
		}
	}
	return true;
}


template <typename TDiscreteFunction, int ref_dim>
template <typename TElem>
bool
DirichletBNDValues<TDiscreteFunction, ref_dim>::
set_dirichlet_linear(	typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						size_t fct, int si,
						matrix_type& mat, vector_type& rhs,
						const discrete_function_type& u, number time)
{
	typename domain_type::position_accessor_type aaPos = u.get_approximation_space().get_domain().get_position_accessor();

	multi_index_vector_type multInd;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		// TODO: if TElem != Vertex we have to do something else
		corner = aaPos[elem];

		if(u.template get_multi_indices_of_geom_obj<TElem>(elem, fct, multInd) != 1)
			return false;
		if(m_bndfct(val, corner, time))
		{
			BlockRef(rhs[multInd[0][0]], multInd[0][1]) = val;
			SetDirichletRow(mat, multInd[0][0], multInd[0][1]);
		}
	}

	return true;
}




} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DIRICHLET_BND_VALUES__ */
