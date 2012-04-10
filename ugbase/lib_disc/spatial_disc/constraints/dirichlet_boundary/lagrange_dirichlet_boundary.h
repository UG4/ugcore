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

#include "lib_disc/spatial_disc/constraints/constraint_base.h"

#include <map>
#include <vector>

namespace ug{

template <	typename TDomain, typename TAlgebra>
class LagrangeDirichletBoundary
	: public ConstraintBase<TDomain, TAlgebra,
	  	  	  	  	  	  	LagrangeDirichletBoundary<TDomain, TAlgebra> >
{
	public:
	///	Base Type
		typedef ConstraintBase<TDomain, TAlgebra,
	  	  	  	  	LagrangeDirichletBoundary<TDomain, TAlgebra> > base_type;

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

	///	Type of conditional boundary functor for a number (bool return value)
		typedef boost::function<bool (number& value, const MathVector<dim>& x, number time)> BNDNumberFunctor;

	///	Type of non-conditional boundary functor for a number (always true)
		typedef boost::function<void (number& value, const MathVector<dim>& x, number time)> NumberFunctor;

	///	Type of non-conditional boundary functor for a Vector (always true)
		typedef boost::function<void (MathVector<dim>& value, const MathVector<dim>& x, number time)> VectorFunctor;

	public:
	///	constructor
		LagrangeDirichletBoundary() {clear();}

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
		void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	///	removes all scheduled dirichlet data.
		void clear();

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
		template <typename TDD>
		void assemble_dirichlet_rows(matrix_type& mat, ConstSmartPtr<TDD> dd, number time = 0.0);

	public:
	///////////////////////////////
	// 	Implement Interface
	///////////////////////////////

	/// sets a unity row for all dirichlet indices
		template <typename TDD>
		void adjust_jacobian(matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<TDD> dd, number time = 0.0);

	/// sets a zero value in the defect for all dirichlet indices
		template <typename TDD>
		void adjust_defect(vector_type& d, const vector_type& u,
		                   ConstSmartPtr<TDD> dd, number time = 0.0);

	/// sets the dirichlet value in the solution for all dirichlet indices
		template <typename TDD>
		void adjust_solution(vector_type& u,
		                     ConstSmartPtr<TDD> dd, number time = 0.0);

	///	sets unity rows in A and dirichlet values in right-hand side b
		template <typename TDD>
		void adjust_linear(matrix_type& A, vector_type& b,
		                   ConstSmartPtr<TDD> dd, number time = 0.0);

	///	sets the dirichlet value in the right-hand side
		template <typename TDD>
		void adjust_rhs(vector_type& b, const vector_type& u,
		                ConstSmartPtr<TDD> dd, number time = 0.0);

	///	returns the type of the constraints
		virtual int type()	{return CT_DIRICHLET;}

	protected:
		void check_functions_and_subsets(FunctionGroup& functionGroup,
		                                 SubsetGroup& subsetGroup, size_t numFct) const;

		void extract_data();

		template <typename TUserData, typename TScheduledUserData>
		void extract_data(std::map<int, std::vector<TUserData*> >& mvUserDataBndSegment,
		                  std::vector<TScheduledUserData>& vUserData);

		template <typename TUserData, typename TDD>
		void adjust_jacobian(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                     matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<TDD> dd, number time);

		template <typename TBaseElem, typename TUserData, typename TDD>
		void adjust_jacobian(const std::vector<TUserData*>& vUserData, int si,
		                     matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<TDD> dd, number time);

		template <typename TUserData, typename TDD>
		void adjust_defect(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                   vector_type& d, const vector_type& u,
		                   ConstSmartPtr<TDD> dd, number time);

		template <typename TBaseElem, typename TUserData, typename TDD>
		void adjust_defect(const std::vector<TUserData*>& vUserData, int si,
		                   vector_type& d, const vector_type& u,
		                   ConstSmartPtr<TDD> dd, number time);

		template <typename TUserData, typename TDD>
		void adjust_solution(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                     vector_type& u, ConstSmartPtr<TDD> dd, number time);

		template <typename TBaseElem, typename TUserData, typename TDD>
		void adjust_solution(const std::vector<TUserData*>& vUserData, int si,
		                     vector_type& u, ConstSmartPtr<TDD> dd, number time);

		template <typename TUserData, typename TDD>
		void adjust_linear(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                   matrix_type& A, vector_type& b,
		                   ConstSmartPtr<TDD> dd, number time);

		template <typename TBaseElem, typename TUserData, typename TDD>
		void adjust_linear(const std::vector<TUserData*>& vUserData, int si,
		                   matrix_type& A, vector_type& b,
		                   ConstSmartPtr<TDD> dd, number time);

		template <typename TUserData, typename TDD>
		void adjust_rhs(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                vector_type& b, const vector_type& u,
		                ConstSmartPtr<TDD> dd, number time);

		template <typename TBaseElem, typename TUserData, typename TDD>
		void adjust_rhs(const std::vector<TUserData*>& vUserData, int si,
		                vector_type& b, const vector_type& u,
		                ConstSmartPtr<TDD> dd, number time);

	protected:
	///	grouping for subset and non-conditional data
		struct NumberData
		{
			const static bool isConditional = false;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			NumberData(NumberFunctor functor_,
			           std::string fctName_, std::string ssName_)
				: functor(functor_), fctName(fctName_), ssName(ssName_)
			{}

			bool operator()(MathVector<1>& val, const MathVector<dim> x, number time) const
			{
				functor(val[0], x, time); return true;
			}

			NumberFunctor functor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

	///	grouping for subset and conditional data
		struct BNDNumberData
		{
			const static bool isConditional = true;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			BNDNumberData(BNDNumberFunctor functor_,
			              std::string fctName_, std::string ssName_)
				: functor(functor_), fctName(fctName_), ssName(ssName_)
			{}
			bool operator()(MathVector<1>& val, const MathVector<dim> x, number time) const
			{
				return functor(val[0], x, time);
			}

			BNDNumberFunctor functor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

	///	grouping for subset and conditional data
		struct ConstNumberData
		{
			const static bool isConditional = false;
			const static size_t numFct = 1;
			typedef MathVector<1> value_type;
			ConstNumberData(number value_,
			              std::string fctName_, std::string ssName_)
				: functor(value_), fctName(fctName_), ssName(ssName_)
			{}
			inline bool operator()(MathVector<1>& val, const MathVector<dim> x, number time) const
			{
				val[0] = functor; return true;
			}

			number functor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

	///	grouping for subset and non-conditional data
		struct VectorData
		{
			const static bool isConditional = false;
			const static size_t numFct = dim;
			typedef MathVector<dim> value_type;
			VectorData(VectorFunctor value_,
			           std::string fctName_, std::string ssName_)
				: functor(value_), fctName(fctName_), ssName(ssName_)
			{}
			bool operator()(MathVector<dim>& val, const MathVector<dim> x, number time) const
			{
				functor(val, x, time); return true;
			}

			VectorFunctor functor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

		std::vector<BNDNumberData> m_vBNDNumberData;
		std::vector<NumberData> m_vNumberData;
		std::vector<ConstNumberData> m_vConstNumberData;

		std::vector<VectorData> m_vVectorData;

	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<NumberData*> > m_mNumberBndSegment;

	///	constant boundary values for all subsets
		std::map<int, std::vector<ConstNumberData*> > m_mConstNumberBndSegment;

	///	conditional boundary values for all subsets
		std::map<int, std::vector<BNDNumberData*> > m_mBNDNumberBndSegment;

	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<VectorData*> > m_mVectorBndSegment;

	protected:
	///	current ApproxSpace
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	current domain
		SmartPtr<TDomain> m_spDomain;

	///	current position accessor
		typename domain_type::position_accessor_type m_aaPos;
};

} // end namespace ug

#include "lagrange_dirichlet_boundary_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__ */
