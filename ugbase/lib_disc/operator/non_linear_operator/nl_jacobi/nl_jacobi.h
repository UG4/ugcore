/*
 * nl_jacobi.h
 *
 *  Created on: 07.01.2013
 *  (main parts are based on the structure of
 *  	newton.h by Andreas Vogel)
 *
 *      Author: raphaelprohl
 */

#ifndef NL_JACOBI_H_
#define NL_JACOBI_H_

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
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	protected:
		typedef DebugWritingObject<TAlgebra> base_writer_type;

	public:
	///	default constructor
		NLJacobiSolver();

	///	constructor
		NLJacobiSolver(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_damp(number damp) {m_damp = damp;}

	///////////////////////////////////////////////////////////////////////////
	//	OperatorInverse interface methods
	///////////////////////////////////////////////////////////////////////////

	/// This operator inverts the Operator op: Y -> X
		virtual bool init(SmartPtr<IOperator<vector_type> > op);

	/// prepare Operator
		virtual bool prepare(vector_type& u);

	/// apply Operator, i.e. op^{-1}(0) = u
		virtual bool apply(vector_type& u);

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

#endif /* NL_JACOBI_H_ */
