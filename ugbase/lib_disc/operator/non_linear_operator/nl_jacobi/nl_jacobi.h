/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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

/*
 *  (main parts are based on the structure of
 *  	newton.h by Andreas Vogel)
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_JACOBI__NL_JACOBIL_H_
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_JACOBI__NL_JACOBIL_H_

#include "lib_algebra/operator/interface/operator_inverse.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"

namespace ug {

/// Nonlinear Jacobi-method
/**
* 	Let L(u) denote a nonlinear functional of n components (l_1,...,l_n).
* 	Then the basic step of the nonlinear Jacobi method is to solve the
* 	i-th equation
*
* 	l_i(u_1^{k},...,u_{i-1}^{k},u_i,u_{i+1}^{k},...,u_{n}^{k}) = 0
*
* 	for u_i and to set u_i^{k+1} = u_i. Here k denotes the iteration-index.
* 	Thus, in order to obtain u^{k+1} from u^k, we solve successively the n
* 	dimensional nonlinear equations for i = 1,...,n. Here this is done
* 	by a scalar newton step for every i. But every other scalar nonlinear
* 	method could be applied as well.
*
* 	Using a damped version of the nonlinear jacobi method results in the
* 	following update of the variables
*
* 	u_i^{k+1} = u_i^k + damp * (u_i -u_i^k).
*
* References:
* <ul>
* <li> J. M. Ortega and W. C. Rheinbolt. Iterative Solution of nonlinear equations in several variables.(1970)
* </ul>
*
* \tparam 	TAlgebra	Algebra type
*/
template <typename TAlgebra>
class NLJacobiSolver
	: 	public IOperatorInverse<typename TAlgebra::vector_type>,
		public DebugWritingObject<TAlgebra>
{
	public:
	///	Algebra type
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	protected:
		using base_writer_type = DebugWritingObject<TAlgebra>;

	public:
	///	default constructor
		NLJacobiSolver();

	///	constructor
		NLJacobiSolver(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

	    ~NLJacobiSolver() override = default;

		void set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_damp(number damp) {m_damp = damp;}

	///	returns information about configuration parameters
		std::string config_string() const override {
			std::stringstream ss;
			ss << "NonlinearJacobiSolver( damp = " << m_damp << ")\n";
			ss << " ConvergenceCheck: ";
			if(m_spConvCheck.valid())	ss << ConfigShift(m_spConvCheck->config_string()) << "\n";
			else						ss << " NOT SET!\n";

			return ss.str();

		}

	///////////////////////////////////////////////////////////////////////////
	//	OperatorInverse interface methods
	///////////////////////////////////////////////////////////////////////////

	/// This operator inverts the Operator op: Y -> X
		bool init(SmartPtr<IOperator<vector_type> > op) override;

	/// prepare Operator
		bool prepare(vector_type& u) override;

	/// apply Operator, i.e. op^{-1}(0) = u
		bool apply(vector_type& u) override;

	private:
	///	help functions for debug output
	///	\{
		void write_debug(const vector_type& vec, const char* filename);
		void write_debug(const matrix_type& mat, const char* filename);
	/// \}

	private:
		SmartPtr<IConvergenceCheck<vector_type> > m_spConvCheck;

		///	damping factor
		number m_damp;

		SmartPtr<AssembledOperator<algebra_type> > m_spAssOp;
		SmartPtr<AssembledLinearOperator<algebra_type> > m_spJ;
		SmartPtr<IAssemble<TAlgebra> > m_spAss;

		///	call counter
		int m_dgbCall;
};

}

#include "nl_jacobi_impl.h"

#endif