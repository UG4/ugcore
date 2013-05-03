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
		using base_writer_type::write_debug;

	public:
	///	default constructor
		NLJacobiSolver();

	///	constructor
		NLJacobiSolver(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);
		void set_damp(number damp) {m_damp = damp;}

		/// This operator inverts the Operator N: Y -> X
		virtual bool init(SmartPtr<IOperator<vector_type> > N);

		/// prepare Operator
		virtual bool prepare(vector_type& u);

		/// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(vector_type& u);

	private:
	///	help functions for debug output
	///	\{
		void write_debug(const vector_type& vec, const char* filename);
		void write_debug(const matrix_type& mat, const char* filename);
	/// \}

	private:
		SmartPtr<IConvergenceCheck<vector_type> > m_spConvCheck;
		number m_damp;

		SmartPtr<AssembledOperator<algebra_type> > m_N;
		SmartPtr<AssembledLinearOperator<algebra_type> > m_J;
		SmartPtr<IAssemble<TAlgebra> > m_spAss;

		///	call counter
		int m_dgbCall;
};

}

#include "nl_jacobi_impl.h"

#endif /* NL_JACOBI_H_ */
