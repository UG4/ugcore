/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__CONSTRAINTS__LAGRANGE_DIRICHLET_BOUNDARY__

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_grid/tools/subset_handler_interface.h"

#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

#include <map>
#include <vector>


// #define LAGRANGE_DIRICHLET_ADJ_TRANSFER_FIX 

namespace ug{

template <typename TDomain, typename TAlgebra>
class DirichletBoundary
	: public IDomainConstraint<TDomain, TAlgebra>
{
	public:
	///	Base Type
		using base_type = IDomainConstraint<TDomain, TAlgebra>;

	///	Type of domain
		using domain_type = TDomain;

	///	world Dimension
		static constexpr int dim = domain_type::dim;

	///	Type of position coordinates (e.g. position_type)
		using position_type = typename domain_type::position_type;

	///	Type of algebra
		using algebra_type = TAlgebra;

	///	Type of algebra matrix
		using matrix_type = typename algebra_type::matrix_type;

	///	Type of algebra vector
		using vector_type = typename algebra_type::vector_type;

	/// Type of value type
		using value_type = typename matrix_type::value_type;

	/// error estimator type
		using err_est_type = MultipleSideAndElemErrEstData<TDomain>;


	public:
	/// If you want to have Dirichlet columns than use the constructor with the flag
	/// The default ist without Dirichlet columns.

	///	constructor
		DirichletBoundary()
			:	m_bInvertSubsetSelection(false),
				m_bDirichletColumns(false),
				m_A(nullptr)
#ifdef LAGRANGE_DIRICHLET_ADJ_TRANSFER_FIX
		 , m_bAdjustTransfers(true)
#endif
			{clear();}

	/// constructor with flag for Dirichlet-Columns.
		DirichletBoundary(bool DirichletColumns)
			:	m_bInvertSubsetSelection(false),
				m_bDirichletColumns(DirichletColumns),
				m_A(nullptr)
#ifdef LAGRANGE_DIRICHLET_ADJ_TRANSFER_FIX
		 , m_bAdjustTransfers(true)
#endif
			{clear();}

#ifdef LAGRANGE_DIRICHLET_ADJ_TRANSFER_FIX
		/// constructor with flag for Dirichlet-Columns.
		DirichletBoundary(bool DirichletColumns, bool bAdjustTransfers)
			:	m_bInvertSubsetSelection(false),
				m_bDirichletColumns(DirichletColumns),
				m_A(nullptr),
				m_bAdjustTransfers(bAdjustTransfers)

			{clear();}
#endif

	///	destructor
		~DirichletBoundary() {}

	///	adds a lua callback (cond and non-cond)
#ifdef UG_FOR_LUA
		void add(const char* name, const char* function, const char* subsets);
		void add(const char* name, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets);
		void add(LuaFunctionHandle fct, const char* function, const char* subsets);
		void add(LuaFunctionHandle fct, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets);
#endif

	///	adds a conditional user-defined value as dirichlet condition for a function on subsets
		void add(SmartPtr<UserData<number, dim, bool> > func, const char* function, const char* subsets);
		void add(SmartPtr<UserData<number, dim, bool> > func, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets);

	///	adds a user-defined value as dirichlet condition for a function on subsets
		void add(SmartPtr<UserData<number, dim> > func, const char* function, const char* subsets);
		void add(SmartPtr<UserData<number, dim> > func, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets);

	///	adds a constant value as dirichlet condition for a function on subsets
		void add(number value, const char* function, const char* subsets);
		void add(number value, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets);

	///	adds a user-defined vector as dirichlet condition for a vector-function on subsets
		void add(SmartPtr<UserData<MathVector<dim>, dim> > func, const char* functions, const char* subsets);
		void add(SmartPtr<UserData<MathVector<dim>, dim> > func, const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets);

	///	adds a number Dirichlet condition whose value is taken from the solution vector itself
		void add(const char* functions, const char* subsets);
		void add(const std::vector<std::string>& Fcts, const std::vector<std::string>& Subsets);
		
	///	inverts the subset selection making the conditions be imposed on the rest of the domain
		void invert_subset_selection() {m_bInvertSubsetSelection = true;};
	
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
		void assemble_dirichlet_rows(matrix_type& mat, ConstSmartPtr<DoFDistribution> dd, number time = 0.0);

	public:
	///////////////////////////////
	// 	Implement Interface
	///////////////////////////////

	/// sets a unity row for all dirichlet indices
		void adjust_jacobian(matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
                             ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = nullptr,
							 const number s_a0 = 1.0);

	/// sets a zero value in the defect for all dirichlet indices
		void adjust_defect(vector_type& d, const vector_type& u,
		                   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0,
                           ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = nullptr,
						   const std::vector<number>* vScaleMass = nullptr,
						   const std::vector<number>* vScaleStiff = nullptr);

	/// sets the dirichlet value in the solution for all dirichlet indices
		void adjust_solution(vector_type& u,
		                     ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

	/// sets zero to correction
		virtual void adjust_correction(vector_type& c, ConstSmartPtr<DoFDistribution> dd, int type,
									   number time = 0.0);

	///	sets unity rows in A and dirichlet values in right-hand side b
		void adjust_linear(matrix_type& A, vector_type& b,
		                   ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

	///	sets the dirichlet value in the right-hand side
		void adjust_rhs(vector_type& b, const vector_type& u,
		                ConstSmartPtr<DoFDistribution> dd, int type, number time = 0.0);

	/// @copydoc IConstraint::adjust_error()
		virtual void adjust_error(const vector_type& u, ConstSmartPtr<DoFDistribution> dd, int type,
								  number time = 0.0,
								  ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = nullptr,
								  const std::vector<number>* vScaleMass = nullptr,
								  const std::vector<number>* vScaleStiff = nullptr);

	///	sets constraints in prolongation
		virtual void adjust_prolongation(matrix_type& P,
										 ConstSmartPtr<DoFDistribution> ddFine,
										 ConstSmartPtr<DoFDistribution> ddCoarse,
										 int type,
										 number time = 0.0);

	///	sets constraints in restriction
		virtual void adjust_restriction(matrix_type& R,
										ConstSmartPtr<DoFDistribution> ddCoarse,
										ConstSmartPtr<DoFDistribution> ddFine,
										int type,
										number time = 0.0);

	///	returns the type of the constraints
		virtual int type() const {return CT_DIRICHLET;}

	protected:
		void check_functions_and_subsets(FunctionGroup& functionGroup,
		                                 SubsetGroup& subsetGroup, size_t numFct) const;

		void extract_data();

		template <typename TUserData, typename TScheduledUserData>
		void extract_data(std::map<int, std::vector<TUserData*> >& mvUserDataBndSegment,
		                  std::vector<TScheduledUserData>& vUserData);

		template <typename TUserData>
		void adjust_jacobian(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                     matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TBaseElem, typename TUserData>
		void adjust_jacobian(const std::vector<TUserData*>& vUserData, int si,
		                     matrix_type& J, const vector_type& u,
		                     ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TUserData>
		void adjust_defect(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                   vector_type& d, const vector_type& u,
		                   ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TBaseElem, typename TUserData>
		void adjust_defect(const std::vector<TUserData*>& vUserData, int si,
		                   vector_type& d, const vector_type& u,
		                   ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TUserData>
		void adjust_solution(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                     vector_type& u, ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TUserData>
		void adjust_correction(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                     vector_type& c, ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TBaseElem, typename TUserData>
		void adjust_solution(const std::vector<TUserData*>& vUserData, int si,
		                     vector_type& u, ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TBaseElem, typename TUserData>
		void adjust_correction(const std::vector<TUserData*>& vUserData, int si,
		                     vector_type& c, ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TUserData>
		void adjust_linear(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                   matrix_type& A, vector_type& b,
		                   ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TBaseElem, typename TUserData>
		void adjust_linear(const std::vector<TUserData*>& vUserData, int si,
		                   matrix_type& A, vector_type& b,
		                   ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TUserData>
		void adjust_rhs(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                vector_type& b, const vector_type& u,
		                ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TBaseElem, typename TUserData>
		void adjust_rhs(const std::vector<TUserData*>& vUserData, int si,
		                vector_type& b, const vector_type& u,
		                ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TUserData>
		void adjust_error(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                  const vector_type& u, ConstSmartPtr<DoFDistribution> dd, number time);

		template <typename TUserData>
		void adjust_prolongation(const std::map<int, std::vector<TUserData*> >& mvUserData,
		                         matrix_type& P,
		                         ConstSmartPtr<DoFDistribution> ddFine,
		                         ConstSmartPtr<DoFDistribution> ddCoarse,
		                         number time);

		template <typename TBaseElem, typename TUserData>
		void adjust_prolongation(const std::vector<TUserData*>& vUserData, int si,
								matrix_type& P,
								ConstSmartPtr<DoFDistribution> ddFine,
								ConstSmartPtr<DoFDistribution> ddCoarse,
								number time);

		template <typename TUserData>
		void adjust_restriction(const std::map<int, std::vector<TUserData*> >& mvUserData,
							   matrix_type& R,
							   ConstSmartPtr<DoFDistribution> ddCoarse,
							   ConstSmartPtr<DoFDistribution> ddFine,
							   number time);

		template <typename TBaseElem, typename TUserData>
		void adjust_restriction(const std::vector<TUserData*>& vUserData, int si,
							   matrix_type& R,
							   ConstSmartPtr<DoFDistribution> ddCoarse,
							   ConstSmartPtr<DoFDistribution> ddFine,
							   number time);

	protected:
	///	grouping for subset and non-conditional data
		struct NumberData
		{
			static constexpr bool isConditional = false;
			static constexpr bool setSolValue = true;
			static constexpr size_t numFct = 1;
			using value_type = MathVector<1>;
			NumberData(SmartPtr<UserData<number, dim> > functor_,
			           std::string fctName_, std::string ssName_)
				: spFunctor(functor_), fctName(fctName_), ssName(ssName_)
			{}

			bool operator()(MathVector<1>& val, const MathVector<dim> x,
			                number time, int si) const
			{
				(*spFunctor)(val[0], x, time, si); return true;
			}

			SmartPtr<UserData<number, dim> > spFunctor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

	///	grouping for subset and conditional data
		struct CondNumberData
		{
			static constexpr bool isConditional = true;
			static constexpr bool setSolValue = true;
			static constexpr size_t numFct = 1;
			using value_type = MathVector<1>;
			CondNumberData(SmartPtr<UserData<number, dim, bool> > functor_,
			              std::string fctName_, std::string ssName_)
				: spFunctor(functor_), fctName(fctName_), ssName(ssName_)
			{}
			bool operator()(MathVector<1>& val, const MathVector<dim> x,
			                number time, int si) const
			{
				return (*spFunctor)(val[0], x, time, si);
			}

			SmartPtr<UserData<number, dim, bool> > spFunctor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

	///	grouping for subset and constant data
		struct ConstNumberData
		{
			static constexpr bool isConditional = false;
			static constexpr bool setSolValue = true;
			static constexpr size_t numFct = 1;
			using value_type = MathVector<1>;
			ConstNumberData(number value_,
			              std::string fctName_, std::string ssName_)
				: functor(value_), fctName(fctName_), ssName(ssName_)
			{}
			inline bool operator()(MathVector<1>& val, const MathVector<dim> x,
			                       number time, int si) const
			{
				val[0] = functor; return true;
			}

			number functor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

	///	grouping for subset and non-conditional vector data
		struct VectorData
		{
			static constexpr bool isConditional = false;
			static constexpr bool setSolValue = true;
			static const size_t numFct = dim;
			using value_type = MathVector<dim>;
			VectorData(SmartPtr<UserData<MathVector<dim>, dim> > value_,
			           std::string fctName_, std::string ssName_)
				: spFunctor(value_), fctName(fctName_), ssName(ssName_)
			{}
			bool operator()(MathVector<dim>& val, const MathVector<dim> x,
			                number time, int si) const
			{
				(*spFunctor)(val, x, time, si); return true;
			}

			SmartPtr<UserData<MathVector<dim>, dim> > spFunctor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

	///	grouping for subset and the data already stored in the solution
		struct OldNumberData
		{
			static constexpr bool isConditional = false;
			static constexpr bool setSolValue = false;
			static constexpr size_t numFct = 1;
			using value_type = MathVector<1>;
			OldNumberData(std::string fctName_, std::string ssName_)
				: fctName(fctName_), ssName(ssName_)
			{}
			inline bool operator()(MathVector<1>& val, const MathVector<dim> x,
			                       number time, int si) const
			{
				return true; // note that we do not set val because setSolValue == false
			}

			number functor;
			std::string fctName;
			std::string ssName;
			size_t fct[numFct];
			SubsetGroup ssGrp;
		};

		std::vector<CondNumberData> m_vBNDNumberData;
		std::vector<NumberData> m_vNumberData;
		std::vector<ConstNumberData> m_vConstNumberData;

		std::vector<VectorData> m_vVectorData;

		std::vector<OldNumberData> m_vOldNumberData;
		
	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<NumberData*> > m_mNumberBndSegment;

	///	constant boundary values for all subsets
		std::map<int, std::vector<ConstNumberData*> > m_mConstNumberBndSegment;

	///	conditional boundary values for all subsets
		std::map<int, std::vector<CondNumberData*> > m_mBNDNumberBndSegment;

	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<VectorData*> > m_mVectorBndSegment;

	///	non-conditional boundary values for all subsets
		std::map<int, std::vector<OldNumberData*> > m_mOldNumberBndSegment;

	protected:
	///	flag for inverting the subset selection: use Dirichlet throughout except for the given subsets
		bool m_bInvertSubsetSelection;
	
	/// flag for setting dirichlet columns
		bool m_bDirichletColumns;

	/// maps a column dirichlet index to the
	/// row and its corresponding matrix entry.
		std::map<int, std::map<int, std::map<int, value_type> > > m_dirichletMap;
	///
		matrix_type* m_A;

	///	current ApproxSpace
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	current domain
		SmartPtr<TDomain> m_spDomain;

	///	current position accessor
		typename domain_type::position_accessor_type m_aaPos;
#ifdef LAGRANGE_DIRICHLET_ADJ_TRANSFER_FIX
		/// flag for setting dirichlet columns
		bool m_bAdjustTransfers;
#endif
};

} // end namespace ug

#include "lagrange_dirichlet_boundary_impl.h"

#endif