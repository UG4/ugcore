/*
 * lib_algebra.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__
#define __H__LIB_ALGEBRA__

#include <iomanip>
// other ug4 modules
#include "common/common.h"

#include "local_matrix_vector/flex_local_matrix_vector.h"

// library intern includes
#include "lib_algebra/multi_index/multi_indices.h"



/////////////////////////////////////////////
/////////////////////////////////////////////
//   Linear Solver Interface
/////////////////////////////////////////////
/////////////////////////////////////////////

namespace ug{

/// Linear Solver Interface
/**
 * Interface for linear solvers
 */
template <typename TAlgebra>
class ILinearSolver {

	// algebra type
	typedef TAlgebra algebra_type;

	// matrix type used
	typedef typename TAlgebra::matrix_type matrix_type;

	// vector type used
	typedef typename TAlgebra::vector_type vector_type;

	public:
	/// this function prepares the execution of solve
	/**
	 * This function prepares the execution of solve.
	 * It must be called before solve.
	 *
	 * Typically, here memory is allocated and submatrices for
	 * MultiGrid solvers are constructed.
	 */
	virtual bool prepare() = 0;

	/// solves A*x = b
	/**
	 * Solves the linear equation system A*x = b
	 *
	 * The current value of x will be used as starting iterative
	 * in iterative solvers
	 *
	 * \param[in]  A Matrix
	 * \param[in]  b Right-Hand side
	 * \param[out] x Solution (to be computed)
	 *
	 */
	virtual bool solve(matrix_type& A, vector_type& x, vector_type &b) = 0;

	/// this function finishes the execution of solve.
	/**
	 * This function finishes the execution of solve.
	 * It must be called after 'solve'.
	 *
	 * Typically, here memory is deallocated.
	 */
	virtual bool finish() = 0;


	/// virtual destructor
	virtual ~ILinearSolver()
	{}
};

/// Iterative Linear Solver Interface
/**
 * Interface for iterative linear solvers
 *
 * Can be used as a linear solver using method 'solve'.
 * Also only on iterative step can be performed using method 'step'.
 */
template <typename TAlgebra>
class IIterativeStep{

	// algebra type
	typedef TAlgebra algebra_type;

	// matrix type used
	typedef typename TAlgebra::matrix_type matrix_type;

	// vector type used
	typedef typename TAlgebra::vector_type vector_type;

	public:

	/// prepares the execution of step.
	/**
	 * This function prepares the execution of step.
	 * It must be called before step.
	 *
	 * Typically, here memory is allocated and submatrices for
	 * MultiGrid solvers are constructed.
	 */
	virtual bool prepare() = 0;

	/// performs an iterative step to solve A*c = d
	/**
	 * Performs a step of an iterative solver: c = B*d with B \approx A^{-1}
	 *
	 * If d = b - A*x, the step computes a correction c such that
	 * (hopefully) x+c is closer to the solution of A*x = b than x.
	 *
	 * c is not initialized on input, i.e. it may have any value. Therefore, if the implementation
	 * needs a vector c=0 as an input vector, the implementation itself must reset the vector c
	 * inside of the function step.
	 *
	 * The new defect is returned in d, i.e. d_{k+1} = b - A*(x_k + c) = d_k - A*c.
	 *
	 * Usually, then the solution is updated x_{k+1} = x_k + c,
	 *
	 * \param[in]  A Matrix
	 * \param[in]  d current defect (i.e. d := b - A*x if the linear system A*x = b has to be solved)
	 * \param[out] d new defect (i.e. d := d - A*c = b - A*(x+c) )
	 * \param[out] c corretion
	 */
	virtual bool step(matrix_type& A, vector_type& c, vector_type &d) = 0;

	/// finishes the execution of step.
	/**
	 * This function finishes the execution of step.
	 * It must be called after 'step'.
	 *
	 * Typically, here memory is deallocated.
	 */
	virtual bool finish() = 0;


	/// virtual destructor
	virtual ~IIterativeStep()
	{}
};


// Standard implementation for iterative solvers using a Iterative Step
template <typename TAlgebra>
class IterativeLinearSolver : public ILinearSolver<TAlgebra>{
		// algebra type
		typedef TAlgebra algebra_type;

		// matrix type used
		typedef typename TAlgebra::matrix_type matrix_type;

		// vector type used
		typedef typename TAlgebra::vector_type vector_type;

	public:
		IterativeLinearSolver(IIterativeStep<TAlgebra>& step, int maxIter, number absTol, number relTol) :
			m_step(step), m_maxIter(maxIter), m_absTol(absTol), m_relTol(relTol)
			{};

		bool prepare()
		{
			if(m_step.prepare() != true)
			{
				UG_LOG("ERROR on 'IterativeLinearSolver::prepare': Cannot prepare step routine.\n");
				return false;
			}

			return true;
		}

		bool solve(matrix_type& A, vector_type& x, vector_type &b)
		{
			// rename b as d (b will be overwritten)
			vector_type& d = b;

			// build defect:  d := b - A*x
			A.matmul_minus(d, x);

			// create correction, that has same memory pattern as x
			vector_type c; c.create(x);

			// compute start norm ||d||_2
			number norm, norm_old, norm_start;
			norm = norm_old = norm_start = d.two_norm();

			// Print Start information
			UG_LOG("\n   ######### Iterative Linear Solver #########\n");
			UG_LOG("  Iter     Defect         Rate \n");
			UG_LOG(std::setw(4) << 0 << ":  " << std::scientific << norm_old <<  "      -------" << std::endl);

			// Iteration loop
			for(int i = 1; i <= m_maxIter; ++i)
			{
				// check if defect is small enough (absolute)
				if(norm < m_absTol)
				{
					UG_LOG("\n ##### Absolute defect " << m_absTol << " reached. Linear Solver converged. #####\n");
					return true;
				}

				// check if defect is small enough (absolute)
				if(norm/norm_start < m_relTol)
				{
					UG_LOG("\n ##### Relative defect " << m_relTol << " reached. Linear Solver converged. #####\n");
					return true;
				}

				// Compute a correction c using one the iterative step
				// Internally the defect is updated d := d - A*c = d - A*(x+c)
				if(m_step.step(A, c, d) == false)
				{
					UG_LOG("Iterative Linear Solver: Error in Step Routine. Aborting.\n");
					return false;
				}

				// add correction to solution
				UG_DLOG(LIB_ALG_LINEAR_SOLVER, 3, "   ||c||_2 = " << c.two_norm() << "\n");
				UG_DLOG(LIB_ALG_LINEAR_SOLVER, 10, "  c = \n" << c << "\n");
				UG_DLOG(LIB_ALG_LINEAR_SOLVER, 3, "   ||x||_2 = " << x.two_norm() << "\n");
				UG_DLOG(LIB_ALG_LINEAR_SOLVER, 10, "  x = \n" << x << "\n");
				x += c;
				UG_DLOG(LIB_ALG_LINEAR_SOLVER, 3, "   ||x+c||_2 = " << x.two_norm() << "\n");
				UG_DLOG(LIB_ALG_LINEAR_SOLVER, 10, "  x+c = \n" << x << "\n");

				// compute new defect norm
				norm = d.two_norm();

				// print convergence rate
				UG_LOG(std::setw(4) << i << ":  " << std::scientific << norm << "    " << norm/norm_old << std::endl);

				// remember current norm
				norm_old = norm;
			}

			UG_LOG("\n ##### Absolute defect " << m_absTol << " and relative defect " << m_relTol << " NOT reached after " << m_maxIter << " Iterations. Iterative Linear Solver did NOT CONVERGE. ##### \n");
			return false;
		}

		bool finish()
		{
			return m_step.finish();
		}

	protected:
		IIterativeStep<TAlgebra>&	m_step;
		int m_maxIter;
		number m_absTol;
		number m_relTol;
};


/////////////////////////////////////////////
/////////////////////////////////////////////
//   Jacobi
/////////////////////////////////////////////
/////////////////////////////////////////////

template <typename TAlgebra>
class JacobiStep : public IIterativeStep<TAlgebra>{

	// algebra type
	typedef TAlgebra algebra_type;

	// matrix type used
	typedef typename TAlgebra::matrix_type matrix_type;

	// vector type used
	typedef typename TAlgebra::vector_type vector_type;

	public:

	JacobiStep(number damp) : m_damp(damp) {};


	/// prepares the execution of step.
	/**
	 * This function prepares the execution of step.
	 * It must be called before step.
	 *
	 * Typically, here memory is allocated and submatrices for
	 * MultiGrid solvers are constructed.
	 */
	bool prepare()
	{
		return true;
	}

	/// performs an iterative step to solve A*c = d
	/**
	 * Performs a step of an iterative solver: c = B*d with B \approx A^{-1}
	 *
	 * If d = b - A*x, the step computes a correction c such that
	 * (hopefully) x+c is closer to the solution of A*x = b than x.
	 *
	 * c is not initialized on input, i.e. it may have any value. Therefore, if the implementation
	 * needs a vector c=0 as an input vector, the implementation itself must reset the vector c
	 * inside of the function step.
	 *
	 * The new defect is returned in d, i.e. d_{k+1} = b - A*(x_k + c) = d_k - A*c.
	 *
	 * Usually, then the solution is updated x_{k+1} = x_k + c,
	 *
	 * \param[in]  A Matrix
	 * \param[in]  d current defect (i.e. d := b - A*x if the linear system A*x = b has to be solved)
	 * \param[out] d new defect (i.e. d := d - A*c = b - A*(x+c) )
	 * \param[out] c corretion
	 */
	bool step(matrix_type& A, vector_type& c, vector_type &d)
	{
		diag_step(A, c, d, m_damp);

		return true;
	}

	/// finishes the execution of step.
	/**
	 * This function finishes the execution of step.
	 * It must be called after 'step'.
	 *
	 * Typically, here memory is deallocated.
	 */
	bool finish()
	{
		return true;
	}


	/// virtual destructor
	~JacobiStep()
	{}

	private:
		number m_damp;
};


}


/////////////////////////////////////////////
/////////////////////////////////////////////
//   Arne Algebra
/////////////////////////////////////////////
/////////////////////////////////////////////

#include "arne_algebra/arnematrix.h"
#include "arne_algebra/arnevector.h"
#include "arne_algebra/arnelinearsolver.h"

namespace ug {

/** Define different algebra types.
 *  An Algebra should export the following typedef:
 *  - matrix_type
 *  - vector_type
 *  - index_type
 */
class ArneAlgebra{
	public:
		// matrix type
		typedef ArneMatrix matrix_type;

		// vector type
		typedef ArneVector vector_type;

		// index_type
		typedef MultiIndex<1> index_type;
};


} // namespace ug

/////////////////////////////////////////////
/////////////////////////////////////////////
//   Hypre Algebra
/////////////////////////////////////////////
/////////////////////////////////////////////


#if 0
//def HYPRELIB_LIB_DIR

#include "hypre_algebra/hyprematrix.h"
#include "hypre_algebra/hyprevector.h"
#include "hypre_algebra/hyprelinearsolver.h"

namespace ug{
class HypreAlgebra{
	public:
		// matrix type
		typedef HypreMatrix matrix_type;

		// vector type
		typedef HypreVector vector_type;

		// index_type
		typedef MultiIndex<1> index_type;

		typedef HYPREboomerAMG linear_solver_type;

};
}

#endif

/////////////////////////////////////////////
/////////////////////////////////////////////
//   Martin Algebra
/////////////////////////////////////////////
/////////////////////////////////////////////
#ifdef USE_MARTIN_ALGEBRA
#include "martin_algebra/vector.h"

namespace ug
{
class MartinAlgebra
	{
	public:
		// matrix type
		typedef SparseMatrix<double> matrix_type;
		
		// vector type
		typedef Vector<double> vector_type;
		
		// index_type
		typedef MultiIndex<1> index_type;
		
		//	typedef HYPREboomerAMG linear_solver_type;
	};

  // this will soon be moved
bool diag_step(const SparseMatrix<number>& A, Vector<number>& x, Vector<number>& b, number damp)
{
	UG_ASSERT(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");

	for(int j=0; j < A.getLength(); j++)
		x[j] += b[j] / A.getDiag(j);

	// update defect
	// b -= A*x;

	return true;
}
}

#endif

#endif /* __H__LIB_ALGEBRA__ */
