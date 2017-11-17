/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NI_SOLVER__NI__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NI_SOLVER__NI__

#include <cmath>

#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/debug_writer.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/assemble_interface.h"

#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/function_spaces/error_elem_marking_strategy.h"

#include "lib_algebra/operator/debug_writer.h"

namespace ug {

//! Nested iteration solver (e.g. for full multigrid)
/*! If error estimators are available and the convergence rate of the FE/FV method is known,
 * 	these methods allow to construct an optimal order method.
 * */
template <typename TDomain, typename TAlgebra>
class NestedIterationSolver
	: 	public IOperatorInverse<typename TAlgebra::vector_type>,
		public DebugWritingObject<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	GridFunction type
		typedef GridFunction<TDomain, TAlgebra> grid_function_type;

	public:
	///	default constructor
		NestedIterationSolver();

	///	constructor setting operator
		NestedIterationSolver(SmartPtr<IOperator<vector_type> > N);

	///	constructor using assembling
		NestedIterationSolver(SmartPtr<IAssemble<TAlgebra> > spAss);

	///	constructor using assembling
		NestedIterationSolver(SmartPtr<IAssemble<TAlgebra> > spAss, SmartPtr<IAssemble<TAlgebra> > spDomErr);

	///	constructor
		NestedIterationSolver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver);

	///	sets the linear solver (this should be a fixed number of multi-grid cycles )
		void set_linear_solver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver) {m_spLinearSolver = LinearSolver;}

	/** @name IOperatorInverse interface*/
	///@{
	/// This operator inverts the Operator N: Y -> X
		virtual bool init(SmartPtr<IOperator<vector_type> > N);

	/// prepare Operator
		virtual bool prepare(vector_type& u);

	/// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(vector_type& u);
	///@}


	///	returns information about configuration parameters
		/**
		 * this should return necessary information about parameters and possibly
		 * calling config_string of subcomponents.
		 *
		 * \returns std::string	necessary information about configuration parameters
		 */

		virtual std::string config_string() const;




	/** @name Adaptive refinement
	 *  for adaptive nested iteration
	 */
	///@{

	//! getter/setter for top level
		int top_level() const			{return m_topLevel;}
		void set_top_level(int lev) 	{ m_topLevel = lev;}

	//! getter/setter for base level
		int base_level() const			{return m_baseLevel;}
		void set_base_level(int lev) 	{ m_baseLevel = lev;}

	//! set grid refiner
		void set_refiner(SmartPtr<IRefiner> r) {m_spRefiner =r; }

	//! set marking strategy
		void set_refinement_marking(SmartPtr<IElementMarkingStrategy<TDomain> > m) {m_spRefinementMarking =m; }
		void set_coarsening_marking(SmartPtr<IElementMarkingStrategy<TDomain> > m) {m_spCoarseningMarking =m; }
		void set_tolerance(number tol) {m_TOL = tol;}

	//! indicates, if grids should be refined further (if desired accuracy has not been achieved)
		bool use_adaptive_refinement() const {return m_bUseAdaptiveRefinement;}
	//! disable grid refinement
		void disable_adaptive_refinement() {m_bUseAdaptiveRefinement= false;}
	//! enable grid refinement
		void enable_adaptive_refinement()  {m_bUseAdaptiveRefinement= true;}


		void set_max_steps(int steps) {m_maxSteps = steps;}
		number last_error() const {return m_lastError;}
	///@}

	/// for debug output
		using DebugWritingObject<TAlgebra>::write_debug;
		using DebugWritingObject<TAlgebra>::vector_debug_writer;
		void set_debug_elem_error(SmartPtr<grid_function_type> spErrEta)
		{m_spElemError = spErrEta;}

	protected:
		//number estimate_error(const grid_function_type& u);
		number refine_domain(const grid_function_type& u);
		number coarsen_domain(const grid_function_type& u);


	private:
	///	linear solver
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spLinearSolver;

	///	assembling routine
		SmartPtr<AssembledOperator<algebra_type> > m_N;
	///	jacobi operator
		SmartPtr<AssembledLinearOperator<algebra_type> > m_J;
	///	assembling
		SmartPtr<IAssemble<TAlgebra> > m_spAss;
		SmartPtr<IAssemble<TAlgebra> > m_spDomErr;

	///	call counter
		int m_dgbCall;
		int m_lastNumSteps;


		int m_baseLevel;
		int m_topLevel;

	/// (adaptive) refinement
		SmartPtr<IRefiner> m_spRefiner;

		SmartPtr<IElementMarkingStrategy<TDomain> > m_spRefinementMarking;
		SmartPtr<IElementMarkingStrategy<TDomain> > m_spCoarseningMarking;
		bool m_bUseAdaptiveRefinement;

		number m_TOL;
		int m_maxSteps;
		number m_lastError;

	/// (optional) debug info for adaptive refinement
		SmartPtr<grid_function_type> m_spElemError;

};

}

#include "nested_iteration_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NI_SOLVER__NI__ */
