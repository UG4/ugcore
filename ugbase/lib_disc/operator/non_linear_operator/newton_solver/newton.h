/*
 * newton.h
 *
 *  Created on: 26.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__
#define __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__

#include <cmath>

#include "lib_algebra/operator/operator_inverse_interface.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "../line_search.h"
#include "lib_algebra/operator/debug_writer.h"

namespace ug {

template <typename TDoFDistribution, typename TAlgebra>
class NewtonSolver : public IOperatorInverse<	typename TAlgebra::vector_type,
												typename TAlgebra::vector_type>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	//	DoFDistribution Type
		typedef TDoFDistribution dof_distribution_type;

	public:
		NewtonSolver(ILinearOperatorInverse<vector_type, vector_type>& LinearSolver,
					IConvergenceCheck& ConvCheck,
					ILineSearch<vector_type>* LineSearch, bool reallocate) :
					m_pLinearSolver(&LinearSolver),
					m_pConvCheck(&ConvCheck),
					m_pLineSearch(LineSearch),
					m_reallocate(reallocate), m_allocated(false),
					m_pDebugWriter(NULL), m_dgbCall(0)
			{};

		NewtonSolver() :
			m_pLinearSolver(NULL), m_pConvCheck(NULL), m_pLineSearch(NULL),
			m_reallocate(false), m_allocated(false), m_pDebugWriter(NULL),
			m_dgbCall(0)
			{};

		void set_linear_solver(ILinearOperatorInverse<vector_type, vector_type>& LinearSolver) {m_pLinearSolver = &LinearSolver;}
		void set_convergence_check(IConvergenceCheck& ConvCheck)
		{
			m_pConvCheck = &ConvCheck;
			m_pConvCheck->set_offset(3);
			m_pConvCheck->set_symbol('#');
			m_pConvCheck->set_name("Newton Solver");
		}
		void set_line_search(ILineSearch<vector_type>& LineSearch) {m_pLineSearch = &LineSearch;}

	///	set debug output
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

		// init: This operator inverts the Operator N: Y -> X
		virtual bool init(IOperator<vector_type, vector_type>& N);

		// prepare Operator
		virtual bool prepare(vector_type& u);

		// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(vector_type& u);

		~NewtonSolver();

	private:
		bool allocate_memory(const vector_type& u);
		bool deallocate_memory();

		bool write_debug(const vector_type& vec, const char* filename)
		{
		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return true;

		//	add iter count to name
			std::string name(filename);
			char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
			name.append(ext);

		//	write
			return m_pDebugWriter->write_vector(vec, name.c_str());
		}

		bool write_debug(const matrix_type& mat, const char* filename)
		{
		//	if no debug writer set, we're done
			if(m_pDebugWriter == NULL) return true;

		//	add iter count to name
			std::string name(filename);
			char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
			name.append(ext);

		//	write
			return m_pDebugWriter->write_matrix(mat, name.c_str());
		}

	private:
		ILinearOperatorInverse<vector_type, vector_type>* m_pLinearSolver;

		// Convergence Check
		IConvergenceCheck* m_pConvCheck;

		// LineSearch
		ILineSearch<vector_type>* m_pLineSearch;

		vector_type m_d;
		vector_type m_c;

		AssembledOperator<dof_distribution_type, algebra_type>* m_N;
		AssembledLinearOperator<dof_distribution_type, algebra_type>* m_J;
		IAssemble<dof_distribution_type, algebra_type>* m_pAss;

		// line search parameters
		int m_maxLineSearch;
		number m_lambda_start;
		number m_lambda_reduce;

		bool m_reallocate;
		bool m_allocated;

	//	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;

		int m_dgbCall;
};

}

#include "newton_impl.h"

#endif /* __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__ */
