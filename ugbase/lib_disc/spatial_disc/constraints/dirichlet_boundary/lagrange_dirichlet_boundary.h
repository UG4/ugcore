/*
 * lagrange_dirichlet_boundary.h
 *
 *  Created on: 08.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/spatial_disc/ip_data/const_user_data.h"
#include "lib_grid/tools/subset_handler_interface.h"

#include <map>
#include <vector>

namespace ug{

template <	typename TDomain, typename TDoFDistribution, typename TAlgebra>
class LagrangeDirichletBoundary
	: public IDomainConstraint<TDomain, TDoFDistribution, TAlgebra>
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

	///	Type of non-conditional boundary functor for a Vector (always true)
		typedef boost::function<void (MathVector<dim>& value, const MathVector<dim>& x, number time)> VectorFunctor;

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

	///	adds a user-defined vector as dirichlet condition for a vector-function on subsets
		void add(VectorFunctor& func, const char* functions, const char* subsets);

	///	sets the approximation space to work on
		void set_approximation_space(IApproximationSpace<TDomain>& approxSpace)
		{
			m_pDomain = &approxSpace.get_domain();
			m_aaPos = m_pDomain->get_position_accessor();
			m_pPattern = &approxSpace;
		}

	///	removes all scheduled dirichlet data.
		void clear()
		{
			m_vScheduledBNDNumberData.clear();
			m_vScheduledNumberData.clear();
			m_vScheduledConstNumberData.clear();
			m_vScheduledVectorData.clear();
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
		                                 SubsetGroup& subsetGroup, size_t numFct) const;

		bool extract_scheduled_data();

		template <typename TUserData, typename TScheduledUserData>
		bool extract_scheduled_data(std::map<int, std::vector<TUserData> >& mvUserDataBndSegment,
		                            const std::vector<TScheduledUserData>& vScheduledUserData);

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
	///	grouping for subset and non-conditional data
		struct NumberData
		{
			const static bool isConditional = false;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			NumberData(const FunctionGroup& fcts_, NumberFunctor functor_)
				: functor(functor_)
			{
				UG_ASSERT(fcts_.num_fct()==numFct, "Must be exactly numFct functions.");
				for(size_t i=0;i<numFct; ++i) fct[i]=fcts_[i];
			}
			bool operator()(MathVector<1>& val, const MathVector<dim> x, number time) const
			{
				functor(val[0], x, time); return true;
			}
			size_t fct[numFct];
			NumberFunctor functor;
		};

	///	grouping for subset and conditional data
		struct BNDNumberData
		{
			const static bool isConditional = true;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			BNDNumberData(const FunctionGroup& fcts_, BNDNumberFunctor functor_)
				: functor(functor_)
			{
				UG_ASSERT(fcts_.num_fct()==numFct, "Must be exactly numFct functions.");
				for(size_t i=0;i<numFct; ++i) fct[i]=fcts_[i];
			}
			bool operator()(MathVector<1>& val, const MathVector<dim> x, number time) const
			{
				return functor(val[0], x, time);
			}
			size_t fct[numFct];
			BNDNumberFunctor functor;
		};

	///	grouping for subset and conditional data
		struct ConstNumberData
		{
			const static bool isConditional = false;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			ConstNumberData(const FunctionGroup& fcts_, number value_)
				: functor(value_)
			{
				UG_ASSERT(fcts_.num_fct()==numFct, "Must be exactly numFct functions.");
				for(size_t i=0;i<numFct; ++i) fct[i]=fcts_[i];
			}
			inline bool operator()(MathVector<1>& val, const MathVector<dim> x, number time) const
			{
				val[0] = functor; return true;
			}
			size_t fct[numFct];
			number functor;
		};

	///	grouping for subset and non-conditional data
		struct VectorData
		{
			const static bool isConditional = false;
			const static size_t numFct = dim;
			typedef MathVector<dim> value_type;
			VectorData(const FunctionGroup& fcts_, VectorFunctor functor_)
				: functor(functor_)
			{
				UG_ASSERT(fcts_.num_fct()==numFct, "Must be exactly numFct functions.");
				for(size_t i=0;i<numFct; ++i) fct[i]=fcts_[i];
			}
			bool operator()(MathVector<dim>& val, const MathVector<dim> x, number time) const
			{
				functor(val, x, time); return true;
			}
			size_t fct[numFct];
			VectorFunctor functor;
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
				: functor(value_), fctName(fctName_), ssName(ssName_)
			{}

			number functor;
			std::string fctName;
			std::string ssName;
		};

	///	to remember the scheduled data
		struct ScheduledVectorData
		{
			ScheduledVectorData(VectorFunctor functor_,
								std::string fctName_, std::string ssName_)
				: functor(functor_), fctName(fctName_), ssName(ssName_)
			{}

			VectorFunctor functor;
			std::string fctName;
			std::string ssName;
		};

		std::vector<ScheduledBNDNumberData> m_vScheduledBNDNumberData;
		std::vector<ScheduledNumberData> m_vScheduledNumberData;
		std::vector<ScheduledConstNumberData> m_vScheduledConstNumberData;

		std::vector<ScheduledVectorData> m_vScheduledVectorData;

	protected:
	///	current domain
		domain_type* m_pDomain;

	///	current position accessor
		typename domain_type::position_accessor_type m_aaPos;

	///	current function pattern
		const FunctionPattern* m_pPattern;

	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<NumberData> > m_mNumberBndSegment;

	///	constant boundary values for all subsets
		std::map<int, std::vector<ConstNumberData> > m_mConstNumberBndSegment;

	///	conditional boundary values for all subsets
		std::map<int, std::vector<BNDNumberData> > m_mBNDNumberBndSegment;

	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<VectorData> > m_mVectorBndSegment;
};


////////////////////////////////////////////////////////////////////////////////
//	setup
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
add(BNDNumberFunctor& func, const char* function, const char* subsets)
{
	m_vScheduledBNDNumberData.push_back(ScheduledBNDNumberData(func, function, subsets));
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
add(NumberFunctor& func, const char* function, const char* subsets)
{
	m_vScheduledNumberData.push_back(ScheduledNumberData(func, function, subsets));
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
add(number value, const char* function, const char* subsets)
{
	m_vScheduledConstNumberData.push_back(ScheduledConstNumberData(value, function, subsets));
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
add(VectorFunctor& func, const char* functions, const char* subsets)
{
	m_vScheduledVectorData.push_back(ScheduledVectorData(func, functions, subsets));
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
check_functions_and_subsets(FunctionGroup& functionGroup, SubsetGroup& subsetGroup, size_t numFct) const
{
//	only number of functions allowed
	if(functionGroup.num_fct() != numFct)
	{
		UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
				" Only "<<numFct<<" function(s) allowed in specification of a"
				" Dirichlet Value, but the following functions given:"
				<<functionGroup<<"\n");
		return false;
	}

//	get subsethandler
	const ISubsetHandler& rSH = m_pPattern->get_subset_handler();

// 	loop subsets
	for(size_t si = 0; si < subsetGroup.num_subsets(); ++si)
	{
	//	get subset index
		const int subsetIndex = subsetGroup[si];

	//	check that subsetIndex is valid
		if(subsetIndex < 0 || subsetIndex >= rSH.num_subsets())
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Invalid Subset Index " << subsetIndex << ". (Valid is"
					" 0, .. , " << rSH.num_subsets() <<").\n");
			return false;
		}

	//	check all functions
		for(size_t i=0; i < functionGroup.num_fct(); ++i)
		{
			const size_t fct = functionGroup[i];

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
	}

//	everything ok
	return true;
}

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
template <typename TUserData, typename TScheduledUserData>
bool LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
extract_scheduled_data(std::map<int, std::vector<TUserData> >& mvUserDataBndSegment,
                       const std::vector<TScheduledUserData>& vScheduledUserData)
{
//	clear the extracted data
	mvUserDataBndSegment.clear();

	for(size_t i = 0; i < vScheduledUserData.size(); ++i)
	{
	//	create Function Group and Subset Group
		FunctionGroup fctGrp;
		SubsetGroup ssGrp;

	//	convert strings
		if(!ConvertStringToSubsetGroup(ssGrp, *m_pPattern,
		                               vScheduledUserData[i].ssName.c_str()))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Subsets '"<<vScheduledUserData[i].ssName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}
		if(!ConvertStringToFunctionGroup(fctGrp, *m_pPattern,
		                                 vScheduledUserData[i].fctName.c_str()))
		{
			UG_LOG("ERROR in 'LagrangeDirichletBoundary:extract_scheduled_data':"
					" Functions '"<<vScheduledUserData[i].fctName<<"' not"
					" all contained in ApproximationSpace.\n");
			return false;
		}

	//	check functions and subsets
		if(!check_functions_and_subsets(fctGrp, ssGrp, TUserData::numFct)) return false;

	// 	loop subsets
		for(size_t si = 0; si < ssGrp.num_subsets(); ++si)
		{
		//	get subset index and function
			const int subsetIndex = ssGrp[si];

		//	get Boundary segment from map
			std::vector<TUserData>& vSegmentFunction = mvUserDataBndSegment[subsetIndex];

		//	remember functor and function
			vSegmentFunction.push_back(TUserData(fctGrp, vScheduledUserData[i].functor));
		}
	}

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

	bool bRet = true;
	bRet &= extract_scheduled_data(m_mNumberBndSegment, m_vScheduledNumberData);
	bRet &= extract_scheduled_data(m_mBNDNumberBndSegment, m_vScheduledBNDNumberData);
	bRet &= extract_scheduled_data(m_mConstNumberBndSegment, m_vScheduledConstNumberData);
	bRet &= extract_scheduled_data(m_mVectorBndSegment, m_vScheduledVectorData);

//	we're done
	return bRet;
}

////////////////////////////////////////////////////////////////////////////////
//	assemble_dirichlet_rows
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
void LagrangeDirichletBoundary<TDomain, TDoFDistribution, TAlgebra>::
assemble_dirichlet_rows(matrix_type& mat, const dof_distribution_type& dd, number time)
{
	extract_scheduled_data();

//	loop boundary subsets
	typename std::map<int, std::vector<BNDNumberData> >::const_iterator iter;
	for(iter = m_mBNDNumberBndSegment.begin(); iter != m_mBNDNumberBndSegment.end(); ++iter)
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
				const size_t fct = userData[i].fct[0];

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
	extract_scheduled_data();

	bool bRet = true;
	bRet &= adjust_jacobian<BNDNumberData>(m_mBNDNumberBndSegment, J, u, dd, time);
	bRet &= adjust_jacobian<NumberData>(m_mNumberBndSegment, J, u, dd, time);
	bRet &= adjust_jacobian<ConstNumberData>(m_mConstNumberBndSegment, J, u, dd, time);

	bRet &= adjust_jacobian<VectorData>(m_mVectorBndSegment, J, u, dd, time);

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
	typename TUserData::value_type val;

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
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i].fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd.local_finite_element_id(fct);

			//	get multi indices
				dd.inner_multi_indices(elem, fct, multInd);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);
					UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");
				}

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
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
	extract_scheduled_data();

	bool bRet = true;
	bRet &= adjust_defect<BNDNumberData>(m_mBNDNumberBndSegment, d, u, dd, time);
	bRet &= adjust_defect<NumberData>(m_mNumberBndSegment, d, u, dd, time);
	bRet &= adjust_defect<ConstNumberData>(m_mConstNumberBndSegment, d, u, dd, time);

	bRet &= adjust_defect<VectorData>(m_mVectorBndSegment, d, u, dd, time);

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
	typename TUserData::value_type val;

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
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i].fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd.local_finite_element_id(fct);

			//	get multi indices
				dd.inner_multi_indices(elem, fct, multInd);

			//	get dof position
				if(TUserData::isConditional){
					InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);
					UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch. (multInd.size()="<<
					          multInd.size()<<", vPos.size()="<<vPos.size()<<")");
				}

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
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
	extract_scheduled_data();

	bool bRet = true;
	bRet &= adjust_solution<BNDNumberData>(m_mBNDNumberBndSegment, u, dd, time);
	bRet &= adjust_solution<NumberData>(m_mNumberBndSegment, u, dd, time);
	bRet &= adjust_solution<ConstNumberData>(m_mConstNumberBndSegment, u, dd, time);

	bRet &= adjust_solution<VectorData>(m_mVectorBndSegment, u, dd, time);

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
	typename TUserData::value_type val;

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
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i].fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd.local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);

			//	get multi indices
				dd.inner_multi_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				//  get dirichlet value
					if(!vUserData[i](val, vPos[j], time)) continue;

				//	set zero for dirichlet values
					BlockRef(u[multInd[j][0]], multInd[j][1]) = val[f];
				}
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
	extract_scheduled_data();

	bool bRet = true;
	bRet &= adjust_linear<BNDNumberData>(m_mBNDNumberBndSegment, A, b, u, dd, time);
	bRet &= adjust_linear<NumberData>(m_mNumberBndSegment, A, b, u, dd, time);
	bRet &= adjust_linear<ConstNumberData>(m_mConstNumberBndSegment, A, b, u, dd, time);

	bRet &= adjust_linear<VectorData>(m_mVectorBndSegment, A, b, u, dd, time);

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
	typename TUserData::value_type val;

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
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i].fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd.local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);

			//	get multi indices
				dd.inner_multi_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(),
						  "Mismatch: numInd="<<multInd.size()<<", numPos="<<vPos.size());

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet and read value
					if(!vUserData[i](val, vPos[j], time)) continue;

					const size_t index = multInd[j][0];
					const size_t alpha = multInd[j][1];

				//	set dirichlet row
					SetDirichletRow(A, index, alpha);

				//	set dirichlet value in rhs
					BlockRef(b[index], alpha) = val[f];
				}
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
	extract_scheduled_data();

	bool bRet = true;
	bRet &= adjust_rhs<BNDNumberData>(m_mBNDNumberBndSegment, b, u, dd, time);
	bRet &= adjust_rhs<NumberData>(m_mNumberBndSegment, b, u, dd, time);
	bRet &= adjust_rhs<ConstNumberData>(m_mConstNumberBndSegment, b, u, dd, time);

	bRet &= adjust_rhs<VectorData>(m_mVectorBndSegment, b, u, dd, time);

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
	typename TUserData::value_type val;

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
			for(size_t f = 0; f < TUserData::numFct; ++f)
			{
			//	get function index
				const size_t fct = vUserData[i].fct[f];

			//	get local finite element id
				const LFEID& lfeID = dd.local_finite_element_id(fct);

			//	get dof position
				InnerDoFPosition(vPos, elem, *m_pDomain, lfeID);

			//	get multi indices
				dd.inner_multi_indices(elem, fct, multInd);

				UG_ASSERT(multInd.size() == vPos.size(), "Size mismatch");

			//	loop dofs on element
				for(size_t j = 0; j < multInd.size(); ++j)
				{
				// 	check if function is dirichlet and read value
					if(!vUserData[i](val, vPos[j], time)) continue;

					const size_t index = multInd[j][0];
					const size_t alpha = multInd[j][1];

				//	set dirichlet value in rhs
					BlockRef(b[index], alpha) = val[f];
				}
			}
		}
	}

//	done
	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__ */
