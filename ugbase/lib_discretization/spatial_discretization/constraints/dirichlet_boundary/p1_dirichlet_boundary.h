/*
 * p1_dirichlet_boundary.h
 *
 *  Created on: 08.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__P1_DIRICHLET_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__P1_DIRICHLET_BOUNDARY__

#include "lib_discretization/common/function_group.h"
#include "lib_discretization/spatial_discretization/domain_discretization_interface.h"
#include "lib_discretization/function_spaces/approximation_space.h"
#include "lib_discretization/spatial_discretization/ip_data/const_user_data.h"
#include "lib_grid/tools/subset_handler_interface.h"

#include <map>
#include <vector>

namespace ug{

template <	typename TDomain, typename TDoFDistribution, typename TAlgebra>
class P1DirichletBoundary
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

	///	Type of conditional boundary functor for a number (bool return value)
		typedef boost::function<bool (number& value, const MathVector<dim>& x, number time)> BNDNumberFunctor;

	///	Type of non-conditional boundary functor for a number (always true)
		typedef boost::function<void (number& value, const MathVector<dim>& x, number time)> NumberFunctor;

	public:
	///	constructor
		P1DirichletBoundary() :
			m_pDomain(NULL), m_pPattern(NULL) {	m_mConditionalBoundarySegment.clear();}

	///	destructor
		~P1DirichletBoundary()
		{
			for(size_t i = 0; i < m_vConstBoundaryNumber.size(); ++i)
				if(m_vConstBoundaryNumber[i] != NULL)
					delete m_vConstBoundaryNumber[i];
		}

		bool add_boundary_value(BNDNumberFunctor& func,
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

		bool add_boundary_value(BNDNumberFunctor func,
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
				std::vector<BNDNumberData>& vSegmentFunction = m_mConditionalBoundarySegment[subsetIndex];

			//	remember functor and function
				vSegmentFunction.push_back(BNDNumberData(fct, func));
			}

		//	we're done
			return true;
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
			set_pattern(approxSpace);
		}

		void set_domain(domain_type& domain)
		{
			m_pDomain = &domain;
			m_aaPos = m_pDomain->get_position_accessor();
		}

		void set_pattern(const FunctionPattern& pattern)
		{
			m_pPattern = &pattern;
			m_mConditionalBoundarySegment.clear();
		}

		void clear()
		{
			m_mConditionalBoundarySegment.clear();
		}

	///	Sets dirichlet rows for all registered dirichlet values
	/**	(implemented by Mr. Xylouris and Mr. Reiter)
	 *
	 * This method is just an attempt to allow to set dirichlet rows in a matrix.
	 * It should probably be a virtual method derived from IDirichletPostProcess.
	 *
	 * Note that adjust_jacobian does the same (...!!!)
	 * You should thus use adjust_jacobian.
	 *
	 * This method is probably removed in the near future!
	 *
	 * It could make sense to keep it but implement it as an overload of
	 * adjust_jacobian...
	 *
	 * If Mr. Vogel decides that this is nonsense, he may of course remove it!!!
	 */
		void assemble_dirichlet_rows(matrix_type& mat, const dof_distribution_type& dofDistr, number time = 0.0);

	public:
	// 	Implement Interface
		bool adjust_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);
		bool adjust_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);
		bool adjust_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);
		bool adjust_rhs(vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);
		bool adjust_solution(vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		virtual int type()	{return CT_DIRICHLET;}

	protected:
	///	grouping for subset and conditional data
		struct BNDNumberData
		{
			BNDNumberData(size_t fct_, BNDNumberFunctor functor_)
				: fct(fct_), functor(functor_) {}
			size_t fct;
			BNDNumberFunctor functor;
		};

	///	grouping for subset and non-conditional data
		struct NumberData
		{
			NumberData(size_t fct_, NumberFunctor functor_)
				: fct(fct_), functor(functor_) {}
			size_t fct;
			NumberFunctor functor;
		};

	protected:
		bool clear_dirichlet_jacobian(			geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<BNDNumberData>& userData, int si, matrix_type& J,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		bool clear_dirichlet_defect(			geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<BNDNumberData>& userData, int si, vector_type& d,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		bool set_dirichlet_solution( 			geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<BNDNumberData>& userData, int si, vector_type& x,
												const dof_distribution_type& dofDistr, number time = 0.0);

		bool set_dirichlet_linear(				geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<BNDNumberData>& userData, int si, matrix_type& mat, vector_type& rhs,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

		bool set_dirichlet_rhs(				geometry_traits<VertexBase>::const_iterator iterBegin,
												geometry_traits<VertexBase>::const_iterator iterEnd,
												const std::vector<BNDNumberData>& userData, int si, vector_type& rhs,
												const vector_type& u, const dof_distribution_type& dofDistr, number time = 0.0);

	protected:
	///	current domain
		domain_type* m_pDomain;

	///	current position accessor
		typename domain_type::position_accessor_type m_aaPos;

	///	current function pattern
		const FunctionPattern* m_pPattern;

	///	remember created ConstNumbers
		std::vector<ConstBoundaryNumber<dim>*> m_vConstBoundaryNumber;

	///	conditional boundary values for all subsets
		std::map<int, std::vector<BNDNumberData> > m_mConditionalBoundarySegment;

	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<NumberData> > m_mBoundarySegment;
};


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
assemble_dirichlet_rows(matrix_type& mat, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator iter;
	for(iter = m_mConditionalBoundarySegment.begin(); iter != m_mConditionalBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<BNDNumberData>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

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
				if(dofDistr.template inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
					return;

				const size_t index = multInd[0][0];
				const size_t alpha = multInd[0][1];

			//	set dirichlet row
				SetDirichletRow(mat, index, alpha);
			}
		}
	}
}



template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator iter;
	for(iter = m_mConditionalBoundarySegment.begin(); iter != m_mConditionalBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<BNDNumberData>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!clear_dirichlet_jacobian(iterBegin, iterEnd, userData, si, J, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::adjust_jacobian':"
					" while calling 'clear_dirichlet_jacobian', aborting.\n");
			return false;
		}
	}

	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator iter;
	for(iter = m_mConditionalBoundarySegment.begin(); iter != m_mConditionalBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<BNDNumberData>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!clear_dirichlet_defect(iterBegin, iterEnd, userData, si, d, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::adjust_jacobian':"
					" while calling 'clear_dirichlet_jacobian', aborting.\n");
			return false;
		}
	}

	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_solution(vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator iter;
	for(iter = m_mConditionalBoundarySegment.begin(); iter != m_mConditionalBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<BNDNumberData>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!set_dirichlet_solution(iterBegin, iterEnd, userData, si, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::adjust_jacobian':"
					" while calling 'clear_dirichlet_jacobian', aborting.\n");
			return false;
		}
	}
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator iter;
	for(iter = m_mConditionalBoundarySegment.begin(); iter != m_mConditionalBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<BNDNumberData>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!set_dirichlet_linear(iterBegin, iterEnd, userData, si, mat, rhs, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::adjust_jacobian':"
					" while calling 'clear_dirichlet_jacobian', aborting.\n");
			return false;
		}
	}
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_rhs(vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator iter;
	for(iter = m_mConditionalBoundarySegment.begin(); iter != m_mConditionalBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<BNDNumberData>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dofDistr.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dofDistr.template end<VertexBase>(si);

		if(!set_dirichlet_rhs(iterBegin, iterEnd, userData, si, rhs, u, dofDistr, time))
		{
			UG_LOG("ERROR in 'P1DirichletBoundary::adjust_rhs':"
					" while calling 'clear_dirichlet_rhs', aborting.\n");
			return false;
		}
	}
	return true;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
clear_dirichlet_defect(	geometry_traits<VertexBase>::const_iterator iterBegin,
						geometry_traits<VertexBase>::const_iterator iterEnd,
						const std::vector<BNDNumberData>& userData, int si,
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
			if(dofDistr.template inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
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
							const std::vector<BNDNumberData>& userData, int si,
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
			if(dofDistr.template inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
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
						const std::vector<BNDNumberData>& userData, int si,
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
			if(dofDistr.template inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
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
						const std::vector<BNDNumberData>& userData, int si,
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
			if(dofDistr.template inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
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


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool
P1DirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
set_dirichlet_rhs(	geometry_traits<VertexBase>::const_iterator iterBegin,
					geometry_traits<VertexBase>::const_iterator iterEnd,
					const std::vector<BNDNumberData>& userData, int si,
					vector_type& rhs,
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
			if(dofDistr.template inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
				return false;

			const size_t index = multInd[0][0];
			const size_t alpha = multInd[0][1];

		//	set dirichlet value
			BlockRef(rhs[index], alpha) = val;
		}
	}
	return true;
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__P1_DIRICHLET_BOUNDARY__ */
