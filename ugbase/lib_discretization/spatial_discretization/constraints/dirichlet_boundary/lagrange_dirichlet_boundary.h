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
			m_pDomain(NULL), m_pPattern(NULL) {clear();}

	///	destructor
		~LagrangeDirichletBoundary() {}

	///	adds a conditional user-defined value as dirichlet condition for a function on subsets
		void add(BNDNumberFunctor& func, const char* function, const char* subsets);

	///	adds a user-defined value as dirichlet condition for a function on subsets
		void add(NumberFunctor& func, const char* function, const char* subsets);

	///	adds a constant value as dirichlet condition for a function on subsets
		void add(number value, const char* function, const char* subsets);

	///	sets the approximation space to work on
		void set_approximation_space(IApproximationSpace<domain_type>& approxSpace)
		{
			m_pDomain = &approxSpace.get_domain();
			m_aaPos = m_pDomain->get_position_accessor();
			m_pPattern = &approxSpace;
			clear();
		}

	///	removes all scheduled dirichlet data.
		void clear()
		{
			m_vScheduledBNDNumberData.clear();
			m_vScheduledNumberData.clear();
			m_vScheduledConstNumberData.clear();
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

	///	returns the type of the constraints
		virtual int type()	{return CT_DIRICHLET;}

	protected:
		bool check_functions_and_subsets(FunctionGroup& functionGroup,
		                                 SubsetGroup& subsetGroup) const;

		bool extract_scheduled_data();

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

	///	grouping for subset and conditional data
		struct ConstNumberData
		{
			const static bool isConditional = false;
			ConstNumberData(size_t fct_, number value_)
				: fct(fct_), value(value_) {}
			inline bool operator()(number& val, const MathVector<dim> x, number time) const
			{
				val = value; return true;
			}
			size_t fct;
			number value;
		};

	///	to remember the scheduled data
		struct ScheduledBNDNumberData
		{
			ScheduledBNDNumberData(BNDNumberFunctor functor_,
								   std::string fctName_, std::string ssName_)
				: functor(functor_), fctName(fctName_), ssName(ssName_)
			{}

			BNDNumberFunctor functor;
			std::string fctName;
			std::string ssName;
		};

	///	to remember the scheduled data
		struct ScheduledNumberData
		{
			ScheduledNumberData(NumberFunctor functor_,
			                    std::string fctName_, std::string ssName_)
				: functor(functor_), fctName(fctName_), ssName(ssName_)
			{}

			NumberFunctor functor;
			std::string fctName;
			std::string ssName;
		};

	///	to remember the scheduled data
		struct ScheduledConstNumberData
		{
			ScheduledConstNumberData(number value_,
			                         std::string fctName_, std::string ssName_)
				: value(value_), fctName(fctName_), ssName(ssName_)
			{}

			number value;
			std::string fctName;
			std::string ssName;
		};

		std::vector<ScheduledBNDNumberData> m_vScheduledBNDNumberData;
		std::vector<ScheduledNumberData> m_vScheduledNumberData;
		std::vector<ScheduledConstNumberData> m_vScheduledConstNumberData;

	protected:
	///	current domain
		domain_type* m_pDomain;

	///	current position accessor
		typename domain_type::position_accessor_type m_aaPos;

	///	current function pattern
		const FunctionPattern* m_pPattern;

	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<NumberData> > m_mBoundarySegment;

	///	constant boundary values for all subsets
		std::map<int, std::vector<ConstNumberData> > m_mConstBoundarySegment;

	///	conditional boundary values for all subsets
		std::map<int, std::vector<BNDNumberData> > m_mConditionalBoundarySegment;
};


////////////////////////////////////////////////////////////////////////////////
//	setup
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
add(BNDNumberFunctor& func,
						const char* function, const char* subsets)
{
	m_vScheduledBNDNumberData.push_back(ScheduledBNDNumberData(func, function, subsets));
	extract_scheduled_data();
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
add(NumberFunctor& func,
						const char* function, const char* subsets)
{
	m_vScheduledNumberData.push_back(ScheduledNumberData(func, function, subsets));
	extract_scheduled_data();
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
add(number value,  const char* function, const char* subsets)
{
	m_vScheduledConstNumberData.push_back(ScheduledConstNumberData(value, function, subsets));
	extract_scheduled_data();
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
check_functions_and_subsets(FunctionGroup& functionGroup, SubsetGroup& subsetGroup) const
{
//	only one function allowed
	if(functionGroup.num_fct() != 1)
	{
		UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
				" Only one function allowed in specification of each"
				" Dirichlet Value, but the following functions given:"
				<<functionGroup<<"\n");
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
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Invalid Subset Index " << subsetIndex << ". (Valid is"
					" 0, .. , " << pSH->num_subsets() <<").\n");
			return false;
		}

		const size_t fct = functionGroup[0];

	// 	check if function exist
		if(fct >= m_pPattern->num_fct())
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Function "<< fct << " does not exist in pattern.\n");
			return false;
		}

	// 	check that function is defined for segment
		if(!m_pPattern->is_def_in_subset(fct, subsetIndex))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
				" Function "<<fct<<" not defined on subset "<<subsetIndex<<".\n");
			return false;
		}
	}

//	everything ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
extract_scheduled_data()
{
//	check that function pattern exists
	if(m_pPattern == NULL)
	{
		UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
				" Approximation Space not set.\n");
		return false;
	}

//	clear data
	m_mConditionalBoundarySegment.clear();
	m_mBoundarySegment.clear();
	m_mConstBoundarySegment.clear();

	for(size_t i = 0; i < m_vScheduledBNDNumberData.size(); ++i)
	{
	//	create Function Group and Subset Group
		FunctionGroup fctGrp;
		SubsetGroup ssGrp;

	//	convert strings
		if(!ConvertStringToSubsetGroup(ssGrp, *m_pPattern,
									   m_vScheduledBNDNumberData[i].ssName.c_str()))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Subsets '"<<m_vScheduledBNDNumberData[i].ssName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}
		if(!ConvertStringToFunctionGroup(fctGrp, *m_pPattern,
										 m_vScheduledBNDNumberData[i].fctName.c_str()))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Functions '"<<m_vScheduledBNDNumberData[i].fctName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}

	//	check functions and subsets
		if(!check_functions_and_subsets(fctGrp, ssGrp)) return false;

	// 	loop subsets
		for(size_t si = 0; si < ssGrp.num_subsets(); ++si)
		{
		//	get subset index and function
			const int subsetIndex = ssGrp[si];
			const size_t fct = fctGrp[0];

		//	get Boundary segment from map
			std::vector<BNDNumberData>& vSegmentFunction
									= m_mConditionalBoundarySegment[subsetIndex];

		//	remember functor and function
			vSegmentFunction.push_back(BNDNumberData(fct, m_vScheduledBNDNumberData[i].functor));
		}
	}


	for(size_t i = 0; i < m_vScheduledNumberData.size(); ++i)
	{
	//	create Function Group and Subset Group
		FunctionGroup fctGrp;
		SubsetGroup ssGrp;

	//	convert strings
		if(!ConvertStringToSubsetGroup(ssGrp, *m_pPattern,
		                               m_vScheduledNumberData[i].ssName.c_str()))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Subsets '"<<m_vScheduledNumberData[i].ssName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}
		if(!ConvertStringToFunctionGroup(fctGrp, *m_pPattern,
		                                 m_vScheduledNumberData[i].fctName.c_str()))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Functions '"<<m_vScheduledNumberData[i].fctName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}

	//	check functions and subsets
		if(!check_functions_and_subsets(fctGrp, ssGrp)) return false;

	// 	loop subsets
		for(size_t si = 0; si < ssGrp.num_subsets(); ++si)
		{
		//	get subset index and function
			const int subsetIndex = ssGrp[si];
			const size_t fct = fctGrp[0];

		//	get Boundary segment from map
			std::vector<NumberData>& vSegmentFunction
									= m_mBoundarySegment[subsetIndex];

		//	remember functor and function
			vSegmentFunction.push_back(NumberData(fct, m_vScheduledNumberData[i].functor));
		}
	}

	for(size_t i = 0; i < m_vScheduledConstNumberData.size(); ++i)
	{
	//	create Function Group and Subset Group
		FunctionGroup fctGrp;
		SubsetGroup ssGrp;

	//	convert strings
		if(!ConvertStringToSubsetGroup(ssGrp, *m_pPattern,
		                               m_vScheduledConstNumberData[i].ssName.c_str()))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Subsets '"<<m_vScheduledConstNumberData[i].ssName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}
		if(!ConvertStringToFunctionGroup(fctGrp, *m_pPattern,
		                                 m_vScheduledConstNumberData[i].fctName.c_str()))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Functions '"<<m_vScheduledConstNumberData[i].fctName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}

	//	check functions and subsets
		if(!check_functions_and_subsets(fctGrp, ssGrp)) return false;

	// 	loop subsets
		for(size_t si = 0; si < ssGrp.num_subsets(); ++si)
		{
		//	get subset index and function
			const int subsetIndex = ssGrp[si];
			const size_t fct = fctGrp[0];

		//	get Boundary segment from map
			std::vector<ConstNumberData>& vSegmentFunction
									= m_mConstBoundarySegment[subsetIndex];

		//	remember functor and function
			vSegmentFunction.push_back(ConstNumberData(fct, m_vScheduledConstNumberData[i].value));
		}
	}

//	we're done
	return true;
}

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
	bRet &= adjust_jacobian<ConstNumberData>(m_mConstBoundarySegment, J, u, dd, time);
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
	bRet &= adjust_defect<ConstNumberData>(m_mConstBoundarySegment, d, u, dd, time);
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
	bRet &= adjust_solution<ConstNumberData>(m_mConstBoundarySegment, u, dd, time);
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
	bRet &= adjust_linear<ConstNumberData>(m_mConstBoundarySegment, A, b, u, dd, time);
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

			UG_ASSERT(multInd.size() == vPos.size(),
			          "Mismatch: numInd="<<multInd.size()<<", numPos="<<vPos.size());

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
	bRet &= adjust_rhs<ConstNumberData>(m_mConstBoundarySegment, b, u, dd, time);
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
