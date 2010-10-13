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

template <int dim>
class DirichletBoundaryFunction
{
	public:
	//	Function Type
		typedef bool (*Boundary_fct)(number&, const MathVector<dim>&, number);

	public:
		virtual Boundary_fct get_bnd_function() const = 0;

		virtual ~DirichletBoundaryFunction() {};
};

template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class DirichletBNDValues : public IDirichletBoundaryValues<TDoFDistribution, TAlgebra> {
	public:
	// 	Type of discrete function
		typedef TDoFDistribution dof_distribution_type;

	// 	Type of domain
		typedef TDomain domain_type;

	// 	Type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

	// 	Type of algebra
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	// 	Type of multi index vector
		typedef typename dof_distribution_type::multi_index_vector_type multi_index_vector_type;

	public:
		typedef bool (*Boundary_fct)(number&, const position_type&, number);

		DirichletBNDValues() : m_bndfct(NULL), m_fct(-1), m_pDomain(NULL) {}
		DirichletBNDValues(size_t fct, Boundary_fct bnd_fct, domain_type& domain) :
		  m_bndfct(bnd_fct), m_fct(fct), m_pDomain(&domain)
		  {}

		void set_function(size_t fct) {m_fct = fct;}
		void set_dirichlet_function(DirichletBoundaryFunction<domain_type::dim>& f) {m_bndfct = f.get_bnd_function();}
		void set_user_function(Boundary_fct bnd_fct) {m_bndfct = bnd_fct;}
		void set_domain(domain_type& domain) {m_pDomain = &domain;}

		IAssembleReturn clear_dirichlet_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr, int si, number time = 0.0);
		IAssembleReturn clear_dirichlet_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr, int si,number time = 0.0);
		IAssembleReturn set_dirichlet_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, int si, number time = 0.0);

		IAssembleReturn set_dirichlet_solution(vector_type& u, const dof_distribution_type& dofDistr, int si, number time = 0.0);

	protected:
		template <typename TElem>
		bool clear_dirichlet_jacobian(			typename geometry_traits<TElem>::const_iterator iterBegin,
												typename geometry_traits<TElem>::const_iterator iterEnd,
												size_t fct, int si, matrix_type& J,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		template <typename TElem>
		bool clear_dirichlet_defect(			typename geometry_traits<TElem>::const_iterator iterBegin,
												typename geometry_traits<TElem>::const_iterator iterEnd,
												size_t fct, int si, vector_type& d,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		template <typename TElem>
		bool set_dirichlet_solution( 			typename geometry_traits<TElem>::const_iterator iterBegin,
												typename geometry_traits<TElem>::const_iterator iterEnd,
												size_t fct, int si, vector_type& x,
												const dof_distribution_type& dofDistr, number time = 0.0);

		template <typename TElem>
		bool set_dirichlet_linear(				typename geometry_traits<TElem>::const_iterator iterBegin,
												typename geometry_traits<TElem>::const_iterator iterEnd,
												size_t fct, int si, matrix_type& mat, vector_type& rhs,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		bool is_dirichlet(size_t fct) {return fct == m_fct;}

	protected:
		Boundary_fct m_bndfct;
		size_t m_fct;
		domain_type* m_pDomain;
};


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DirichletBNDValues<TDomain, TDoFDistribution, TAlgebra>::
clear_dirichlet_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr, int si, number time)
{
	if(clear_dirichlet_jacobian<VertexBase>(dofDistr.template begin<VertexBase>(si), dofDistr.template end<VertexBase>(si), m_fct, si, J, u, dofDistr, time) == false)
		{UG_LOG("Error in 'assemble_jacobian' while calling 'clear_dirichlet_jacobian', aborting.\n"); return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DirichletBNDValues<TDomain, TDoFDistribution, TAlgebra>::
clear_dirichlet_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr, int si, number time)
{
	if(clear_dirichlet_defect<VertexBase>(dofDistr.template begin<VertexBase>(si),dofDistr.template end<VertexBase>(si), m_fct, si, d, u, dofDistr, time) == false)
		{UG_LOG("Error in set_dirichlet_nodes, aborting.\n");return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DirichletBNDValues<TDomain, TDoFDistribution, TAlgebra>::
set_dirichlet_solution(vector_type& u, const dof_distribution_type& dofDistr, int si, number time)
{
	if(set_dirichlet_solution<VertexBase>(dofDistr.template begin<VertexBase>(si), dofDistr.template end<VertexBase>(si), m_fct, si, u, dofDistr, time) == false)
		{UG_LOG("Error in set_dirichlet_nodes, aborting.\n"); return IAssemble_ERROR;}

	return IAssemble_OK;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
DirichletBNDValues<TDomain, TDoFDistribution, TAlgebra>::
set_dirichlet_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, int si, number time)
{
	if(set_dirichlet_linear<VertexBase>(dofDistr.template begin<VertexBase>(si), dofDistr.template end<VertexBase>(si), m_fct, si, mat, rhs, u, dofDistr) == false)
		{UG_LOG("Error in assemble_linear, aborting.\n"); return IAssemble_ERROR;}

	return IAssemble_OK;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TElem>
bool
DirichletBNDValues<TDomain, TDoFDistribution, TAlgebra>::
clear_dirichlet_defect(	typename geometry_traits<TElem>::const_iterator iterBegin,
						typename geometry_traits<TElem>::const_iterator iterEnd,
						size_t fct, int si,
						vector_type& d,
						const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
	typename domain_type::position_accessor_type aaPos = m_pDomain->get_position_accessor();
	multi_index_vector_type multInd;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		corner = aaPos[elem];

		if(dofDistr.template get_inner_multi_indices<TElem>(elem, fct, multInd) != 1)
			return false;

		if(m_bndfct(val, corner, time))
		{
			BlockRef(d[multInd[0][0]], multInd[0][1]) = 0.0;
		}
	}

	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TElem>
bool
DirichletBNDValues<TDomain, TDoFDistribution, TAlgebra>::
clear_dirichlet_jacobian(	typename geometry_traits<TElem>::const_iterator iterBegin,
							typename geometry_traits<TElem>::const_iterator iterEnd,
							size_t fct, int si,
							matrix_type& J,
							const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
	typename domain_type::position_accessor_type aaPos = m_pDomain->get_position_accessor();
	multi_index_vector_type multInd;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		corner = aaPos[elem];

		if(dofDistr.template get_inner_multi_indices<TElem>(elem, fct, multInd) != 1)
			return false;

		if(m_bndfct(val, corner, time))
		{
			SetDirichletRow(J, multInd[0][0], multInd[0][1]);
		}

	}
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TElem>
bool
DirichletBNDValues<TDomain, TDoFDistribution, TAlgebra>::
set_dirichlet_solution(	typename geometry_traits<TElem>::const_iterator iterBegin,
						typename geometry_traits<TElem>::const_iterator iterEnd,
						size_t fct, int si,
						vector_type& x,
						const dof_distribution_type& dofDistr, number time)
{
	typename domain_type::position_accessor_type aaPos = m_pDomain->get_position_accessor();

	multi_index_vector_type multInd;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		// TODO: if TElem != Vertex we have to do something else
		corner = aaPos[elem];

		if(dofDistr.template get_inner_multi_indices<TElem>(elem, fct, multInd) != 1)
			return false;

		if(m_bndfct(val, corner, time))
		{
			BlockRef(x[multInd[0][0]], multInd[0][1]) = val;
		}
	}
	return true;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TElem>
bool
DirichletBNDValues<TDomain, TDoFDistribution, TAlgebra>::
set_dirichlet_linear(	typename geometry_traits<TElem>::const_iterator iterBegin,
						typename geometry_traits<TElem>::const_iterator iterEnd,
						size_t fct, int si,
						matrix_type& mat, vector_type& rhs,
						const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
	typename domain_type::position_accessor_type aaPos = m_pDomain->get_position_accessor();

	multi_index_vector_type multInd;

	number val;
	position_type corner;

	for(typename geometry_traits<TElem>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
	{
		TElem* elem = *iter;
		// TODO: if TElem != Vertex we have to do something else
		corner = aaPos[elem];

		if(m_bndfct(val, corner, time))
		{
			if(dofDistr.template get_inner_multi_indices<TElem>(elem, fct, multInd) != 1)
				return false;
			const size_t index = multInd[0][0];
			const size_t alpha = multInd[0][1];
			BlockRef(rhs[index], alpha) = val;
			SetDirichletRow(mat, multInd[0][0], multInd[0][1]);
		}
	}

	return true;
}




} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DIRICHLET_BND_VALUES__ */
