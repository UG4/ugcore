/*
 * lagrange_dirichlet_boundary.h
 *
 *  Created on: 08.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__LAGRANGE_DIFFERENCE_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__LAGRANGE_DIFFERENCE_BOUNDARY__

#include "lib_discretization/common/function_group.h"
#include "lib_discretization/spatial_discretization/domain_discretization_interface.h"
#include "lib_discretization/function_spaces/approximation_space.h"
#include "lib_discretization/spatial_discretization/ip_data/const_user_data.h"
#include "lib_grid/tools/subset_handler_interface.h"

#include <map>
#include <vector>

namespace ug{

template <	typename TDomain, typename TDoFDistribution, typename TAlgebra>
class LagrangeDifferenceBoundary
	: public IConstraint<TDoFDistribution, TAlgebra>
{
	public:
	///	Type of dof distribution
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Type of domain
		typedef TDomain domain_type;

	///	world Dimension
		static const int dim = domain_type::dim;

	///	Type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of multi index vector
		typedef typename dof_distribution_type::multi_index_vector_type multi_index_vector_type;

	public:
	///	constructor
		LagrangeDifferenceBoundary() :
			m_pDomain(NULL), m_pPattern(NULL),
			m_differenceCmp(-1), m_minuendCmp(-1), m_subtrahendCmp(-1)
		{}

	///	destructor
		~LagrangeDifferenceBoundary() {}

	///	adds a constant value as dirichlet condition for a function on subsets
		bool set(const char* subsets, const char* difference,
				 const char* minuend, const char* subtrahend)
		{
			if(m_pPattern == NULL)
			{
				UG_LOG("ERROR in 'LagrangeDifferenceBoundary::set':"
						" No Approximation Space given, please set it.\n");
				return false;
			}

			if(!ConvertStringToSubsetGroup(m_ssGrp, *m_pPattern, subsets))
			{
				UG_LOG("ERROR in 'LagrangeDifferenceBoundary:set':"
						" Subsets '"<<subsets<<"' not"
						" all contained in ApproximationSpace.\n");
				return false;
			}

			FunctionGroup fctGrp;

			if(!ConvertStringToFunctionGroup(fctGrp, *m_pPattern, difference))
			{
				UG_LOG("ERROR in 'LagrangeDifferenceBoundary::set':"
						" Functions '"<<difference<<"' not"
						" all contained in ApproximationSpace.\n");
				return false;
			}
			if(fctGrp.num_fct() != 1)
			{
				UG_LOG("ERROR in 'LagrangeDifferenceBoundary::set':"
						" Only one function allowed, '"<<difference<<"' given.\n");
				return false;
			}
			m_differenceCmp = fctGrp[0];

			if(!ConvertStringToFunctionGroup(fctGrp, *m_pPattern, minuend))
			{
				UG_LOG("ERROR in 'LagrangeDifferenceBoundary::set':"
						" Functions '"<<minuend<<"' not"
						" all contained in ApproximationSpace.\n");
				return false;
			}
			if(fctGrp.num_fct() != 1)
			{
				UG_LOG("ERROR in 'LagrangeDifferenceBoundary::set':"
						" Only one function allowed, '"<<minuend<<"' given.\n");
				return false;
			}
			m_minuendCmp = fctGrp[0];

			if(!ConvertStringToFunctionGroup(fctGrp, *m_pPattern, subtrahend))
			{
				UG_LOG("ERROR in 'LagrangeDifferenceBoundary::set':"
						" Functions '"<<subtrahend<<"' not"
						" all contained in ApproximationSpace.\n");
				return false;
			}
			if(fctGrp.num_fct() != 1)
			{
				UG_LOG("ERROR in 'LagrangeDifferenceBoundary::set':"
						" Only one function allowed, '"<<subtrahend<<"' given.\n");
				return false;
			}
			m_subtrahendCmp = fctGrp[0];

			return true;
		}

	///	sets the approximation space to work on
		void set_approximation_space(IApproximationSpace<domain_type>& approxSpace)
		{
			m_pDomain = &approxSpace.get_domain();
			m_aaPos = m_pDomain->get_position_accessor();
			m_pPattern = &approxSpace;
		}

	public:
	///////////////////////////////
	// 	Implement Interface
	///////////////////////////////

	/// sets a unity row for all dirichlet indices
		bool adjust_jacobian(matrix_type& J, const vector_type& u,
		                     const dof_distribution_type& dd, number time = 0.0);

	/// sets a zero value in the defect for all dirichlet indices
		bool adjust_defect(vector_type& d, const vector_type& u,
		                   const dof_distribution_type& dd, number time = 0.0);

	/// sets the dirichlet value in the solution for all dirichlet indices
		bool adjust_solution(vector_type& u,
		                     const dof_distribution_type& dd, number time = 0.0);

	///	sets unity rows in A and dirichlet values in right-hand side b
		bool adjust_linear(matrix_type& A, vector_type& b, const vector_type& u,
		                   const dof_distribution_type& dd, number time = 0.0);

	///	sets the dirichlet value in the right-hand side
		bool adjust_rhs(vector_type& b, const vector_type& u,
		                const dof_distribution_type& dd, number time = 0.0);

		template <typename TBaseElem>
		bool adjust_jacobian(int si, matrix_type& J, const vector_type& u,
							 const dof_distribution_type& dd, number time);

		template <typename TBaseElem>
		bool adjust_defect(int si, vector_type& d, const vector_type& u,
							const dof_distribution_type& dd, number time);

	///	returns the type of the constraints
		virtual int type()	{return CT_DIRICHLET;}

	protected:
	///	current domain
		domain_type* m_pDomain;

	///	current position accessor
		typename domain_type::position_accessor_type m_aaPos;

	///	current function pattern
		const FunctionPattern* m_pPattern;

	///	result cmp
		int m_differenceCmp;

	///	minuend cmp
		int m_minuendCmp;

	///	subtrahend cmp
		int m_subtrahendCmp;

	///	subsets
		SubsetGroup m_ssGrp;
};

////////////////////////////////////////////////////////////////////////////////
//	adjust JACOBIAN
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDifferenceBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dd, number time)
{
//	loop boundary subsets
	for(size_t i = 0; i < m_ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = m_ssGrp[i];

	//	adapt jacobian for dofs in each base element type
		bool bRes = true;
		if(dd.has_indices_on(VERTEX))
			bRes &= adjust_jacobian<VertexBase>(si, J, u, dd, time);
		if(dd.has_indices_on(EDGE))
			bRes &= adjust_jacobian<EdgeBase>(si, J, u, dd, time);
		if(dd.has_indices_on(FACE))
			bRes &= adjust_jacobian<Face>(si, J, u, dd, time);
		if(dd.has_indices_on(VOLUME))
			bRes &= adjust_jacobian<Volume>(si, J, u, dd, time);

	//	check success
		if(!bRes)
		{
			UG_LOG("ERROR in 'LagrangeDifferenceBoundary::adjust_jacobian':"
					" While calling 'adapt_jacobian', aborting.\n");
			return false;
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TBaseElem>
bool LagrangeDifferenceBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_jacobian(int si, matrix_type& J, const vector_type& u,
           	    const dof_distribution_type& dd, number time)
{
//	create Multiindex
	multi_index_vector_type DiffMultInd, MinuendMultInd, SubtrahendMultInd;

//	iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd.template begin<TBaseElem>(si);
	iterEnd = dd.template end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	get multi indices
		dd.inner_multi_indices(elem, m_differenceCmp, DiffMultInd);
		dd.inner_multi_indices(elem, m_minuendCmp, MinuendMultInd);
		dd.inner_multi_indices(elem, m_subtrahendCmp, SubtrahendMultInd);

		UG_ASSERT(DiffMultInd.size() == MinuendMultInd.size(), "wrong size");
		UG_ASSERT(DiffMultInd.size() == SubtrahendMultInd.size(), "wrong size");

	//	loop dofs on element
		for(size_t j = 0; j < DiffMultInd.size(); ++j)
		{
			SetDirichletRow(J, DiffMultInd[j][0], DiffMultInd[j][1]);

			BlockRef(J(DiffMultInd[j][0], DiffMultInd[j][0]),
					   DiffMultInd[j][1], DiffMultInd[j][1]) = 1.0;

			BlockRef(J(DiffMultInd[j][0], MinuendMultInd[j][0]),
					   DiffMultInd[j][1], MinuendMultInd[j][1]) = -1.0;

			BlockRef(J(DiffMultInd[j][0], SubtrahendMultInd[j][0]),
					   DiffMultInd[j][1], SubtrahendMultInd[j][1]) = -1.0;
		}
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	adjust DEFECT
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDifferenceBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              const dof_distribution_type& dd, number time)
{
//	loop boundary subsets
	for(size_t i = 0; i < m_ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = m_ssGrp[i];

	//	adapt jacobian for dofs in each base element type
		bool bRes = true;
		if(dd.has_indices_on(VERTEX))
			bRes &= adjust_defect<VertexBase>(si, d, u, dd, time);
		if(dd.has_indices_on(EDGE))
			bRes &= adjust_defect<EdgeBase>(si, d, u, dd, time);
		if(dd.has_indices_on(FACE))
			bRes &= adjust_defect<Face>(si, d, u, dd, time);
		if(dd.has_indices_on(VOLUME))
			bRes &= adjust_defect<Volume>(si, d, u, dd, time);

	//	check success
		if(!bRes)
		{
			UG_LOG("ERROR in 'LagrangeDifferenceBoundary::adjust_defect':"
					" While calling 'adjust_defect', aborting.\n");
			return false;
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TBaseElem>
bool LagrangeDifferenceBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_defect(int si, vector_type& d, const vector_type& u,
              const dof_distribution_type& dd, number time)
{
//	create Multiindex
	multi_index_vector_type DiffMultInd, MinuendMultInd, SubtrahendMultInd;

//	iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd.template begin<TBaseElem>(si);
	iterEnd = dd.template end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	get multi indices
		dd.inner_multi_indices(elem, m_differenceCmp, DiffMultInd);
		dd.inner_multi_indices(elem, m_minuendCmp, MinuendMultInd);
		dd.inner_multi_indices(elem, m_subtrahendCmp, SubtrahendMultInd);

		UG_ASSERT(DiffMultInd.size() == MinuendMultInd.size(), "wrong size");
		UG_ASSERT(DiffMultInd.size() == SubtrahendMultInd.size(), "wrong size");

	//	loop dofs on element
		for(size_t j = 0; j < DiffMultInd.size(); ++j)
		{
			BlockRef(d[DiffMultInd[j][0]], DiffMultInd[j][1])
				=
				BlockRef(u[DiffMultInd[j][0]], DiffMultInd[j][1])
				-BlockRef(u[MinuendMultInd[j][0]], MinuendMultInd[j][1])
				-BlockRef(u[SubtrahendMultInd[j][0]], SubtrahendMultInd[j][1]);
		}
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	adjust SOLUTION
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDifferenceBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_solution(vector_type& u, const dof_distribution_type& dd, number time)
{

	return true;
}


////////////////////////////////////////////////////////////////////////////////
//	adjust LINEAR
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDifferenceBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_linear(matrix_type& A, vector_type& b,
              const vector_type& u, const dof_distribution_type& dd, number time)
{
	return false;
}

////////////////////////////////////////////////////////////////////////////////
//	adjust RHS
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDifferenceBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_rhs(vector_type& b, const vector_type& u,
           const dof_distribution_type& dd, number time)
{
	return false;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__LAGRANGE_DIFFERENCE_BOUNDARY__ */
