/*
 * p1_dirichlet_boundary.h
 *
 *  Created on: 08.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS_DIRICHLET_BOUNDARY__P1_DIRICHLET_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS_DIRICHLET_BOUNDARY__P1_DIRICHLET_BOUNDARY__

#include "lib_discretization/common/common.h"
#include "lib_discretization/spatial_discretization/domain_discretization_interface.h"
#include "lib_grid/tools/subset_handler_interface.h"

#include <map>
#include <vector>

namespace ug{

template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class P1DirichletBoundary : public IPostProcess<TDoFDistribution, TAlgebra> {
	public:
	// 	Type of discrete function
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	// 	Type of domain
		typedef TDomain domain_type;

	//	Dimension
		static const int dim = domain_type::dim;

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
		typedef typename IBoundaryNumberProvider<dim>::functor_type BNDNumberFunctor;

		P1DirichletBoundary() :
			m_pDomain(NULL), m_pPattern(NULL) {	m_mBoundarySegment.clear();}

		~P1DirichletBoundary()
		{
			for(size_t i = 0; i < m_vConstBoundaryNumber.size(); ++i)
				if(m_vConstBoundaryNumber[i] != NULL)
					delete m_vConstBoundaryNumber[i];
		}

		bool add_boundary_value(typename IBoundaryNumberProvider<dim>::functor_type func,
								const char* function, const char* subsets)
		{
		//	check that function pattern exists
			if(m_pPattern == NULL)
			{
				UG_LOG("P1DirichletBoundary:add_boundary_value: Function Pattern not set.\n");
				return false;
			}

		//	create Function Group and Subset Group
			FunctionGroup functionGroup;
			SubsetGroup subsetGroup;

		//	convert strings
			if(!ConvertStringToSubsetGroup(subsetGroup, *m_pPattern, subsets))
			{
				UG_LOG("ERROR while parsing Subsets.\n");
				return false;
			}
			if(!ConvertStringToFunctionGroup(functionGroup, *m_pPattern, function))
			{
				UG_LOG("ERROR while parsing Functions.\n");
				return false;
			}

		//	only one function allowed
			if(functionGroup.num_fct() != 1)
			{
				UG_LOG("P1DirichletBoundary:add_boundary_value: Exactly one function needed, but given '"<<function<<"' as functions.\n");
				return false;
			}

			return add_boundary_value(func, functionGroup.unique_id(0), subsetGroup);
		}

		bool add_boundary_value(IBoundaryNumberProvider<dim>& user,
								const char* function, const char* subsets)
		{
		//	forward request
			return add_boundary_value(user.get_functor(), function, subsets);
		}

		bool add_boundary_value(typename IBoundaryNumberProvider<dim>::functor_type func,
								size_t fct, SubsetGroup subsetGroup)
		{
		//	check that function pattern exists
			if(m_pPattern == NULL)
			{
				UG_LOG("P1DirichletBoundary:add_boundary_value: Function Pattern not set.\n");
				return false;
			}

		//	get subsethandler
			const ISubsetHandler* pSH = m_pPattern->get_subset_handler();

		// 	loop subsets
			for(size_t si = 0; si < subsetGroup.num_subsets(); ++si)
			{
			//	get subset index
				const int subsetIndex = subsetGroup[si];

			//	check that subsetIndex is valid
				if(subsetIndex < 0 || subsetIndex >= pSH->num_subsets())
				{
					UG_LOG("P1DirichletBoundary:add_boundary_value: Invalid subset Index "
							<< subsetIndex << ". (Valid is 0, .. , " << pSH->num_subsets() <<").\n");
					return false;
				}

			// 	check if function exist
				if(fct >= m_pPattern->num_fct())
				{
					UG_LOG("P1DirichletBoundary:add_boundary_value: Function "
							<< fct << " does not exist in pattern.\n");
					return false;
				}

			// 	check that function is defined for segment
				if(!m_pPattern->is_def_in_subset(fct, subsetIndex))
				{
					UG_LOG("P1DirichletBoundary:add_boundary_value: Function "
							<< fct << " not defined in subset " << subsetIndex << ".\n");
					return false;
				}

			//	get Boundary segment from map
				std::vector<UserDataFunction>& vSegmentFunction = m_mBoundarySegment[subsetIndex];

			//	remember functor and function
				vSegmentFunction.push_back(UserDataFunction(fct, func));
			}

		//	we're done
			return true;
		}

		bool add_boundary_value(IBoundaryNumberProvider<dim>& user,
								size_t fct, SubsetGroup subsetGroup)
		{
		//	forward request
			return add_boundary_value(user.get_functor(), fct, subsetGroup);
		}

		bool add_constant_boundary_value(number value,  const char* function, const char* subsets)
		{
			ConstBoundaryNumber<dim>* valProvider = new ConstBoundaryNumber<dim>;
			valProvider->set(value);
			m_vConstBoundaryNumber.push_back(valProvider);
			return add_boundary_value(*valProvider, function, subsets);
		}

		void set_approximation_space(IApproximationSpace<domain_type>& approxSpace)
		{
			set_domain(approxSpace.get_domain());
			set_pattern(approxSpace.get_function_pattern());
		}

		void set_domain(domain_type& domain)
		{
			m_pDomain = &domain;
			m_aaPos = m_pDomain->get_position_accessor();
		}

		void set_pattern(const FunctionPattern& pattern)
		{
			m_pPattern = &pattern;
			m_mBoundarySegment.clear();
		}

	public:
	// 	Implement Interface
		IAssembleReturn post_process_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);
		IAssembleReturn post_process_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);
		IAssembleReturn post_process_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);
		IAssembleReturn post_process_solution(vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		virtual int type()	{return PPT_DIRICHLET;}

	protected:
		struct UserDataFunction
		{
			UserDataFunction(size_t fct_, BNDNumberFunctor functor_)
				: fct(fct_), functor(functor_) {}

			size_t fct;
			BNDNumberFunctor functor;
		};

	protected:
		bool clear_dirichlet_jacobian(			geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<UserDataFunction>& userData, int si, matrix_type& J,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		bool clear_dirichlet_defect(			geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<UserDataFunction>& userData, int si, vector_type& d,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		bool set_dirichlet_solution( 			geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<UserDataFunction>& userData, int si, vector_type& x,
												const dof_distribution_type& dofDistr, number time = 0.0);

		bool set_dirichlet_linear(				geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<UserDataFunction>& userData, int si, matrix_type& mat, vector_type& rhs,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

	protected:
		domain_type* m_pDomain;

		// remember created ConstNumbers
		std::vector<ConstBoundaryNumber<dim>*> m_vConstBoundaryNumber;

		const FunctionPattern* m_pPattern;

		std::map<int, std::vector<UserDataFunction> > m_mBoundarySegment;

		typename domain_type::position_accessor_type m_aaPos;
};


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
post_process_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<UserDataFunction> >::const_iterator iter;
	for(iter = m_mBoundarySegment.begin(); iter != m_mBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<UserDataFunction>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!clear_dirichlet_jacobian(iterBegin, iterEnd, userData, si, J, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::post_process_jacobian':"
					" while calling 'clear_dirichlet_jacobian', aborting.\n");
			return IAssemble_ERROR;
		}
	}

	return IAssemble_OK;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
post_process_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<UserDataFunction> >::const_iterator iter;
	for(iter = m_mBoundarySegment.begin(); iter != m_mBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<UserDataFunction>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!clear_dirichlet_defect(iterBegin, iterEnd, userData, si, d, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::post_process_jacobian':"
					" while calling 'clear_dirichlet_jacobian', aborting.\n");
			return IAssemble_ERROR;
		}
	}

	return IAssemble_OK;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
post_process_solution(vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<UserDataFunction> >::const_iterator iter;
	for(iter = m_mBoundarySegment.begin(); iter != m_mBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<UserDataFunction>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!set_dirichlet_solution(iterBegin, iterEnd, userData, si, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::post_process_jacobian':"
					" while calling 'clear_dirichlet_jacobian', aborting.\n");
			return IAssemble_ERROR;
		}
	}
	return IAssemble_OK;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
IAssembleReturn
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
post_process_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<UserDataFunction> >::const_iterator iter;
	for(iter = m_mBoundarySegment.begin(); iter != m_mBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<UserDataFunction>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!set_dirichlet_linear(iterBegin, iterEnd, userData, si, mat, rhs, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::post_process_jacobian':"
					" while calling 'clear_dirichlet_jacobian', aborting.\n");
			return IAssemble_ERROR;
		}
	}
	return IAssemble_OK;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
clear_dirichlet_defect(	geometry_traits<VertexBase>::const_iterator iterBegin,
						geometry_traits<VertexBase>::const_iterator iterEnd,
						const std::vector<UserDataFunction>& userData, int si,
						vector_type& d,
						const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	for readin
	number val;
	position_type corner;

//	loop vertices
	for(geometry_traits<VertexBase>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
	{
	//	get vertex
		VertexBase* vertex = *iter;

	//	get corner position
		corner = m_aaPos[vertex];

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < userData.size(); ++i)
		{
		// 	check if function is dirichlet
			if(!userData[i].functor(val, corner, time)) continue;

		//	get function index
			const size_t fct = userData[i].fct;

		//	get multi indices
			if(dofDistr.template get_inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
				return false;

		//	set dirichlet value
			BlockRef(d[multInd[0][0]], multInd[0][1]) = 0.0;
		}
	}
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
clear_dirichlet_jacobian(	geometry_traits<VertexBase>::const_iterator iterBegin,
							geometry_traits<VertexBase>::const_iterator iterEnd,
							const std::vector<UserDataFunction>& userData, int si,
							matrix_type& J,
							const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	for readin
	number val;
	position_type corner;

//	loop vertices
	for(geometry_traits<VertexBase>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
	{
	//	get vertex
		VertexBase* vertex = *iter;

	//	get corner position
		corner = m_aaPos[vertex];

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < userData.size(); ++i)
		{
		// 	check if function is dirichlet
			if(!userData[i].functor(val, corner, time)) continue;

		//	get function index
			const size_t fct = userData[i].fct;

		//	get multi indices
			if(dofDistr.template get_inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
				return false;

		//	set dirichlet row
			SetDirichletRow(J, multInd[0][0], multInd[0][1]);
		}
	}
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
set_dirichlet_solution(	geometry_traits<VertexBase>::const_iterator iterBegin,
						geometry_traits<VertexBase>::const_iterator iterEnd,
						const std::vector<UserDataFunction>& userData, int si,
						vector_type& x,
						const dof_distribution_type& dofDistr, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	for readin
	number val;
	position_type corner;

//	loop vertices
	for(geometry_traits<VertexBase>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
	{
	//	get vertex
		VertexBase* vertex = *iter;

	//	get corner position
		corner = m_aaPos[vertex];

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < userData.size(); ++i)
		{
		// 	check if function is dirichlet
			if(!userData[i].functor(val, corner, time)) continue;

		//	get function index
			const size_t fct = userData[i].fct;

		//	get multi indices
			if(dofDistr.template get_inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
				return false;

		//	set dirichlet value
			BlockRef(x[multInd[0][0]], multInd[0][1]) = val;
		}
	}
	return true;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
set_dirichlet_linear(	geometry_traits<VertexBase>::const_iterator iterBegin,
						geometry_traits<VertexBase>::const_iterator iterEnd,
						const std::vector<UserDataFunction>& userData, int si,
						matrix_type& mat, vector_type& rhs,
						const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	for readin
	number val;
	position_type corner;

//	loop vertices
	for(geometry_traits<VertexBase>::const_iterator iter = iterBegin; iter != iterEnd; iter++)
	{
	//	get vertex
		VertexBase* vertex = *iter;

	//	get corner position
		corner = m_aaPos[vertex];

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < userData.size(); ++i)
		{
		// 	check if function is dirichlet
			if(!userData[i].functor(val, corner, time)) continue;

		//	get function index
			const size_t fct = userData[i].fct;

		//	get multi indices
			if(dofDistr.template get_inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
				return false;

			const size_t index = multInd[0][0];
			const size_t alpha = multInd[0][1];

		//	set dirichlet value
			BlockRef(rhs[index], alpha) = val;

		//	set dirichlet row
			SetDirichletRow(mat, index, alpha);
		}
	}
	return true;
}




} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__POST_PROCESS_DIRICHLET_BOUNDARY__P1_DIRICHLET_BOUNDARY__ */
