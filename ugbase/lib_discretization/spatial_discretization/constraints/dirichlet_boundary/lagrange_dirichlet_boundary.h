/*
 * lagrange_dirichlet_boundary.h
 *
 *  Created on: 08.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__

#include "lib_discretization/common/function_group.h"
#include "lib_discretization/spatial_discretization/domain_discretization_interface.h"
#include "lib_discretization/function_spaces/approximation_space.h"
#include "lib_discretization/spatial_discretization/ip_data/const_user_data.h"
#include "lib_grid/tools/subset_handler_interface.h"

#include <map>
#include <vector>

namespace ug{

template <	typename TDomain, typename TDoFDistribution, typename TAlgebra>
class LagrangeDirichletBoundary
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
		LagrangeDirichletBoundary() :
			m_pDomain(NULL), m_pPattern(NULL) {	m_mConditionalBoundarySegment.clear();}

	///	destructor
		~LagrangeDirichletBoundary()
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
				UG_LOG("LagrangeDirichletBoundary:add_boundary_value: Function Pattern not set.\n");
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
				UG_LOG("LagrangeDirichletBoundary:add_boundary_value: Exactly one function needed, but given '"<<function<<"' as functions.\n");
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
				UG_LOG("LagrangeDirichletBoundary:add_boundary_value: Function Pattern not set.\n");
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
					UG_LOG("LagrangeDirichletBoundary:add_boundary_value: Invalid subset Index "
							<< subsetIndex << ". (Valid is 0, .. , " << pSH->num_subsets() <<").\n");
					return false;
				}

			// 	check if function exist
				if(fct >= m_pPattern->num_fct())
				{
					UG_LOG("LagrangeDirichletBoundary:add_boundary_value: Function "
							<< fct << " does not exist in pattern.\n");
					return false;
				}

			// 	check that function is defined for segment
				if(!m_pPattern->is_def_in_subset(fct, subsetIndex))
				{
					UG_LOG("LagrangeDirichletBoundary:add_boundary_value: Function "
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
		void assemble_dirichlet_rows(matrix_type& mat, const dof_distribution_type& dd, number time = 0.0);

	public:
	// 	Implement Interface
		bool adjust_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dd, number time = 0.0);
		bool adjust_defect(vector_type& d, const vector_type& d, const dof_distribution_type& dd, number time = 0.0);
		bool adjust_solution(vector_type& u, const dof_distribution_type& dd, number time = 0.0);

		bool adjust_linear(matrix_type& A, vector_type& b, const vector_type& u, const dof_distribution_type& dd, number time = 0.0);
		bool adjust_rhs(vector_type& b, const vector_type& u, const dof_distribution_type& dd, number time = 0.0);

		virtual int type()	{return CT_DIRICHLET;}

	protected:
	///	grouping for subset and conditional data
		struct BNDNumberData
		{
			const static bool isConditional = true;
			BNDNumberData(size_t fct_, BNDNumberFunctor functor_)
				: fct(fct_), functor(functor_) {}
			bool operator()(number& val, const MathVector<dim> x, number time) const
			{
				return functor(val, x, time);
			}
			size_t fct;
			BNDNumberFunctor functor;
		};

	///	grouping for subset and non-conditional data
		struct NumberData
		{
			const static bool isConditional = false;
			NumberData(size_t fct_, NumberFunctor functor_)
				: fct(fct_), functor(functor_) {}
			bool operator()(number& val, const MathVector<dim> x, number time) const
			{
				functor(val, x, time); return true;
			}
			size_t fct;
			NumberFunctor functor;
		};

	protected:
		template <typename TUserData>
		bool adjust_jacobian(const std::map<int, std::vector<TUserData> >& mvUserData,
		                     matrix_type& J, const vector_type& u,
		                     const dof_distribution_type& dd, number time);

		template <typename TBaseElem, typename TUserData>
		bool adjust_jacobian(const std::vector<TUserData>& vUserData, int si,
		                     matrix_type& J, const vector_type& u,
		                     const dof_distribution_type& dd, number time);

		template <typename TUserData>
		bool adjust_defect(const std::map<int, std::vector<TUserData> >& mvUserData,
		                   vector_type& d, const vector_type& u,
		                   const dof_distribution_type& dd, number time);

		template <typename TBaseElem, typename TUserData>
		bool adjust_defect(const std::vector<TUserData>& vUserData, int si,
		                   vector_type& d, const vector_type& u,
		                   const dof_distribution_type& dd, number time);

		template <typename TUserData>
		bool adjust_solution(const std::map<int, std::vector<TUserData> >& mvUserData,
		                     vector_type& u, const dof_distribution_type& dd, number time);

		template <typename TBaseElem, typename TUserData>
		bool adjust_solution(const std::vector<TUserData>& vUserData, int si,
		                     vector_type& u, const dof_distribution_type& dd, number time);

		template <typename TUserData>
		bool adjust_linear(const std::map<int, std::vector<TUserData> >& mvUserData,
		                   matrix_type& A, vector_type& b, const vector_type& u,
		                   const dof_distribution_type& dd, number time);

		template <typename TBaseElem, typename TUserData>
		bool adjust_linear(const std::vector<TUserData>& vUserData, int si,
		                   matrix_type& A, vector_type& b, const vector_type& u,
		                   const dof_distribution_type& dd, number time);

		template <typename TUserData>
		bool adjust_rhs(const std::map<int, std::vector<TUserData> >& mvUserData,
		                vector_type& b, const vector_type& u,
		                const dof_distribution_type& dd, number time);

		template <typename TBaseElem, typename TUserData>
		bool adjust_rhs(const std::vector<TUserData>& vUserData, int si,
		                vector_type& b, const vector_type& u,
		                const dof_distribution_type& dd, number time);

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

////////////////////////////////////////////////////////////////////////////////
//	assemble_dirichlet_rows
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
assemble_dirichlet_rows(matrix_type& mat, const dof_distribution_type& dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator iter;
	for(iter = m_mConditionalBoundarySegment.begin(); iter != m_mConditionalBoundarySegment.end(); ++iter)
	{
		const int si = (*iter).first;
		const std::vector<BNDNumberData>& userData = (*iter).second;

		typename geometry_traits<VertexBase>::const_iterator iterBegin 	= dd.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd 	= dd.template end<VertexBase>(si);

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
				if(dd.template inner_multi_indices<VertexBase>(vertex, fct, multInd) != 1)
					return;

				const size_t index = multInd[0][0];
				const size_t alpha = multInd[0][1];

			//	set dirichlet row
				SetDirichletRow(mat, index, alpha);
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	adjust JACOBIAN
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dd, number time)
{
	bool bRet = true;
	bRet &= adjust_jacobian<BNDNumberData>(m_mConditionalBoundarySegment, J, u, dd, time);
	bRet &= adjust_jacobian<NumberData>(m_mBoundarySegment, J, u, dd, time);
	return bRet;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_jacobian(const std::map<int, std::vector<TUserData> >& mvUserData,
                matrix_type& J, const vector_type& u,
           	    const dof_distribution_type& dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		bool bRes = true;
		if(dd.has_indices_on(VERTEX))
			bRes &= adjust_jacobian<VertexBase, TUserData>(vUserData, si, J, u, dd, time);
		if(dd.has_indices_on(EDGE))
			bRes &= adjust_jacobian<EdgeBase, TUserData>(vUserData, si, J, u, dd, time);
		if(dd.has_indices_on(FACE))
			bRes &= adjust_jacobian<Face, TUserData>(vUserData, si, J, u, dd, time);
		if(dd.has_indices_on(VOLUME))
			bRes &= adjust_jacobian<Volume, TUserData>(vUserData, si, J, u, dd, time);

	//	check success
		if(!bRes)
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary::adjust_jacobian':"
					" While calling 'adapt_jacobian' for TUserData, aborting.\n");
			return false;
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_jacobian(const std::vector<TUserData>& vUserData, int si,
                matrix_type& J, const vector_type& u,
           	    const dof_distribution_type& dd, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	dummy for readin
	number val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd.template begin<TBaseElem>(si);
	iterEnd = dd.template end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
		//	get function index
			const size_t fct = vUserData[i].fct;

		//	get local finite element id
			const LFEID& lfeID = dd.local_finite_element_id(fct);

		//	get dof position
			if(TUserData::isConditional)
				InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);

		//	get multi indices
			dd.inner_multi_indices(elem, fct, multInd);

			UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

		//	loop dofs on element
			for(size_t j = 0; j < vPos.size(); ++j)
			{
			// 	check if function is dirichlet
				if(TUserData::isConditional){
					if(!vUserData[i](val, vPos[j], time)) continue;
				}

			//	set dirichlet row
				SetDirichletRow(J, multInd[j][0], multInd[j][1]);
			}
		}
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	adjust DEFECT
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_defect(vector_type& d, const vector_type& u,
              const dof_distribution_type& dd, number time)
{
	bool bRet = true;
	bRet &= adjust_defect<BNDNumberData>(m_mConditionalBoundarySegment, d, u, dd, time);
	bRet &= adjust_defect<NumberData>(m_mBoundarySegment, d, u, dd, time);
	return bRet;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_defect(const std::map<int, std::vector<TUserData> >& mvUserData,
               vector_type& d, const vector_type& u,
               const dof_distribution_type& dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		bool bRes = true;
		if(dd.has_indices_on(VERTEX))
			bRes &= adjust_defect<VertexBase, TUserData>(vUserData, si, d, u, dd, time);
		if(dd.has_indices_on(EDGE))
			bRes &= adjust_defect<EdgeBase, TUserData>(vUserData, si, d, u, dd, time);
		if(dd.has_indices_on(FACE))
			bRes &= adjust_defect<Face, TUserData>(vUserData, si, d, u, dd, time);
		if(dd.has_indices_on(VOLUME))
			bRes &= adjust_defect<Volume, TUserData>(vUserData, si, d, u, dd, time);

	//	check success
		if(!bRes)
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary::adjust_defect':"
					" While calling 'adjust_defect' for TUserData, aborting.\n");
			return false;
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_defect(const std::vector<TUserData>& vUserData, int si,
              vector_type& d, const vector_type& u,
              const dof_distribution_type& dd, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	dummy for readin
	number val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd.template begin<TBaseElem>(si);
	iterEnd = dd.template end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
		//	get function index
			const size_t fct = vUserData[i].fct;

		//	get local finite element id
			const LFEID& lfeID = dd.local_finite_element_id(fct);

		//	get dof position
			if(TUserData::isConditional)
				InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);

		//	get multi indices
			dd.inner_multi_indices(elem, fct, multInd);

			UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

		//	loop dofs on element
			for(size_t j = 0; j < vPos.size(); ++j)
			{
			// 	check if function is dirichlet
				if(TUserData::isConditional){
					if(!vUserData[i](val, vPos[j], time)) continue;
				}

			//	set zero for dirichlet values
				BlockRef(d[multInd[j][0]], multInd[j][1]) = 0.0;
			}
		}
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	adjust SOLUTION
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_solution(vector_type& u, const dof_distribution_type& dd, number time)
{
	bool bRet = true;
	bRet &= adjust_solution<BNDNumberData>(m_mConditionalBoundarySegment, u, dd, time);
	bRet &= adjust_solution<NumberData>(m_mBoundarySegment, u, dd, time);
	return bRet;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_solution(const std::map<int, std::vector<TUserData> >& mvUserData,
                vector_type& u, const dof_distribution_type& dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		bool bRes = true;
		if(dd.has_indices_on(VERTEX))
			bRes &= adjust_solution<VertexBase, TUserData>(vUserData, si, u, dd, time);
		if(dd.has_indices_on(EDGE))
			bRes &= adjust_solution<EdgeBase, TUserData>(vUserData, si, u, dd, time);
		if(dd.has_indices_on(FACE))
			bRes &= adjust_solution<Face, TUserData>(vUserData, si, u, dd, time);
		if(dd.has_indices_on(VOLUME))
			bRes &= adjust_solution<Volume, TUserData>(vUserData, si, u, dd, time);

	//	check success
		if(!bRes)
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary::adjust_solution':"
					" While calling 'adjust_solution' for TUserData, aborting.\n");
			return false;
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_solution(const std::vector<TUserData>& vUserData, int si,
                vector_type& u, const dof_distribution_type& dd, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	value readin
	number val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd.template begin<TBaseElem>(si);
	iterEnd = dd.template end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
		//	get function index
			const size_t fct = vUserData[i].fct;

		//	get local finite element id
			const LFEID& lfeID = dd.local_finite_element_id(fct);

		//	get dof position
			InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);

		//	get multi indices
			dd.inner_multi_indices(elem, fct, multInd);

			UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

		//	loop dofs on element
			for(size_t j = 0; j < vPos.size(); ++j)
			{
			//  get dirichlet value
				if(!vUserData[i](val, vPos[j], time)) continue;

			//	set zero for dirichlet values
				BlockRef(u[multInd[j][0]], multInd[j][1]) = val;
			}
		}
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	adjust LINEAR
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_linear(matrix_type& A, vector_type& b,
              const vector_type& u, const dof_distribution_type& dd, number time)
{
	bool bRet = true;
	bRet &= adjust_linear<BNDNumberData>(m_mConditionalBoundarySegment, A, b, u, dd, time);
	bRet &= adjust_linear<NumberData>(m_mBoundarySegment, A, b, u, dd, time);
	return bRet;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_linear(const std::map<int, std::vector<TUserData> >& mvUserData,
              matrix_type& A, vector_type& b, const vector_type& u,
           	  const dof_distribution_type& dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		bool bRes = true;
		if(dd.has_indices_on(VERTEX))
			bRes &= adjust_linear<VertexBase, TUserData>(vUserData, si, A, b, u, dd, time);
		if(dd.has_indices_on(EDGE))
			bRes &= adjust_linear<EdgeBase, TUserData>(vUserData, si, A, b, u, dd, time);
		if(dd.has_indices_on(FACE))
			bRes &= adjust_linear<Face, TUserData>(vUserData, si, A, b, u, dd, time);
		if(dd.has_indices_on(VOLUME))
			bRes &= adjust_linear<Volume, TUserData>(vUserData, si, A, b, u, dd, time);

	//	check success
		if(!bRes)
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary::adjust_linear':"
					" While calling 'adjust_linear' for TUserData, aborting.\n");
			return false;
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_linear(const std::vector<TUserData>& vUserData, int si,
              matrix_type& A, vector_type& b, const vector_type& u,
              const dof_distribution_type& dd, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	readin value
	number val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd.template begin<TBaseElem>(si);
	iterEnd = dd.template end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
		//	get function index
			const size_t fct = vUserData[i].fct;

		//	get local finite element id
			const LFEID& lfeID = dd.local_finite_element_id(fct);

		//	get dof position
			InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);

		//	get multi indices
			dd.inner_multi_indices(elem, fct, multInd);

			UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

		//	loop dofs on element
			for(size_t j = 0; j < vPos.size(); ++j)
			{
			// 	check if function is dirichlet and read value
				if(!vUserData[i](val, vPos[j], time)) continue;

				const size_t index = multInd[j][0];
				const size_t alpha = multInd[j][1];

			//	set dirichlet row
				SetDirichletRow(A, index, alpha);

			//	set dirichlet value in rhs
				BlockRef(b[index], alpha) = val;
			}
		}
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//	adjust RHS
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_rhs(vector_type& b, const vector_type& u,
           const dof_distribution_type& dd, number time)
{
	bool bRet = true;
	bRet &= adjust_rhs<BNDNumberData>(m_mConditionalBoundarySegment, b, u, dd, time);
	bRet &= adjust_rhs<NumberData>(m_mBoundarySegment, b, u, dd, time);
	return bRet;
}


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_rhs(const std::map<int, std::vector<TUserData> >& mvUserData,
           vector_type& b, const vector_type& u,
           const dof_distribution_type& dd, number time)
{
//	loop boundary subsets
	typename std::map<int, std::vector<TUserData> >::const_iterator iter;
	for(iter = mvUserData.begin(); iter != mvUserData.end(); ++iter)
	{
	//	get subset index
		const int si = (*iter).first;

	//	get vector of scheduled dirichlet data on this subset
		const std::vector<TUserData>& vUserData = (*iter).second;

	//	adapt jacobian for dofs in each base element type
		bool bRes = true;
		if(dd.has_indices_on(VERTEX))
			bRes &= adjust_rhs<VertexBase, TUserData>(vUserData, si, b, u, dd, time);
		if(dd.has_indices_on(EDGE))
			bRes &= adjust_rhs<EdgeBase, TUserData>(vUserData, si, b, u, dd, time);
		if(dd.has_indices_on(FACE))
			bRes &= adjust_rhs<Face, TUserData>(vUserData, si, b, u, dd, time);
		if(dd.has_indices_on(VOLUME))
			bRes &= adjust_rhs<Volume, TUserData>(vUserData, si, b, u, dd, time);

	//	check success
		if(!bRes)
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary::adjust_rhs':"
					" While calling 'adjust_rhs' for TUserData, aborting.\n");
			return false;
		}
	}

//	ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TBaseElem, typename TUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
adjust_rhs(const std::vector<TUserData>& vUserData, int si,
           vector_type& b, const vector_type& u,
           const dof_distribution_type& dd, number time)
{
//	create Multiindex
	multi_index_vector_type multInd;

//	readin value
	number val;

//	position of dofs
	std::vector<position_type> vPos;

//	iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;
	iter = dd.template begin<TBaseElem>(si);
	iterEnd = dd.template end<TBaseElem>(si);

//	loop elements
	for( ; iter != iterEnd; iter++)
	{
	//	get vertex
		TBaseElem* elem = *iter;

	//	loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size(); ++i)
		{
		//	get function index
			const size_t fct = vUserData[i].fct;

		//	get local finite element id
			const LFEID& lfeID = dd.local_finite_element_id(fct);

		//	get dof position
			InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);

		//	get multi indices
			dd.inner_multi_indices(elem, fct, multInd);

			UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

		//	loop dofs on element
			for(size_t j = 0; j < vPos.size(); ++j)
			{
			// 	check if function is dirichlet and read value
				if(!vUserData[i](val, vPos[j], time)) continue;

				const size_t index = multInd[j][0];
				const size_t alpha = multInd[j][1];

			//	set dirichlet value in rhs
				BlockRef(b[index], alpha) = val;
			}
		}
	}

//	done
	return true;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__ */
