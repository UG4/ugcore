/*
 * nl_gauss_seidel.h
 *
 *  Created on: 07.01.2013
 *  (main parts are based on the structure of
 *  	newton.h and some ideas of Sebastian Reiter & Andreas Vogel)
 *
 *      Author: raphaelprohl
 */

#ifndef NL_GAUSS_SEIDEL_H_
#define NL_GAUSS_SEIDEL_H_


#include "lib_algebra/operator/interface/operator_inverse.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"

namespace ug {

template <typename TDomain, typename TAlgebra>
class NLGaussSeidelSolver
	: public IOperatorInverse<typename TAlgebra::vector_type>,
	  public DebugWritingObject<TAlgebra>
{
	private:
	///	own type
		typedef NLGaussSeidelSolver<TDomain, TAlgebra> this_type;

	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Domain type
		typedef TDomain domain_type;

	///	Type of approximation space
		typedef ApproximationSpace<domain_type>	approx_space_type;

	protected:
		typedef DebugWritingObject<TAlgebra> base_writer_type;
		using base_writer_type::write_debug;

	public:
	///	default constructor
		NLGaussSeidelSolver();

	///	constructor
		NLGaussSeidelSolver(SmartPtr<approx_space_type> spApproxSpace,
					SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_approximation_space(SmartPtr<approx_space_type> spApproxSpace)
		{m_spApproxSpace = spApproxSpace;}
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
		///	Approximation Space
		SmartPtr<approx_space_type> m_spApproxSpace;

		///	DoF distribution pointer
		ConstSmartPtr<LevelDoFDistribution> m_spLevDD;
		ConstSmartPtr<SurfaceDoFDistribution> m_spSurfDD;

		/// DoF Distribution used
		GridLevel m_gridLevel;

		SmartPtr<IConvergenceCheck<vector_type> > m_spConvCheck;

		vector_type m_d;
		vector_type m_c;

		SmartPtr<AssembledOperator<algebra_type> > m_N;
		SmartPtr<AssembledLinearOperator<algebra_type> > m_J;
		IAssemble<algebra_type>* m_pAss;

		number m_damp;

		///	call counter
		int m_dgbCall;

		///	selector of elements with contributions to a specific DoF
		Selector m_vElemSelector;
};

}

#include "nl_gauss_seidel_impl.h"

#endif /* NL_GAUSS_SEIDEL_H_ */
