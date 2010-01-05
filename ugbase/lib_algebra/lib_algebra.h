/*
 * lib_algebra.h
 *
 *  Created on: 02.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__
#define __H__LIB_ALGEBRA__

/*
#include "hypre_algebra/hyprematrix.h"
#include "hypre_algebra/hyprevector.h"
#include "hypre_algebra/hyprelinearsolver.h"
*/

#include "arne_algebra/arnematrix.h"
#include "arne_algebra/arnevector.h"
#include "arne_algebra/arnelinearsolver.h"

namespace ug {

/*
typedef HYPREboomerAMG LinearSolver;
typedef HypreMatrix Matrix;
typedef HypreVector Vector;
*/


typedef ArneJacobi LinearSolver;
typedef ArneMatrix Matrix;
typedef ArneVector Vector;

/// Vector Interface
/**
 * Interface for Vectors (using virtual functions)
 */
class IVector {
	public:
		/// creates a vector of length nentries
		bool create_vector(int nentries);

		/// deletes memory
		bool delete_vector();

		/// set selected entries
		bool set_values(int nvalues, int* indices, double* values);

		/// add to selected entries
		bool add_values(int nvalues, int* indices, double* values);

		/// get selected entries
		bool get_values(int nvalues, int* indices, double* values) const;

		/// finalize
		bool finalize();

		/// print to ascii - file
		bool printToFile(const char* filename);

		/// virtual destructor
		virtual ~ArneVector()
		{}

		virtual ArneVector& operator+= (const ArneVector& v);

		ArneVector& operator-= (const ArneVector& v);

		number norm2();

		bool set(number w);

		int length();

		ArneVector::ScalarVector* getStorage();


};


/// Matrix Interface
/**
 * Interface for Matrices (using virtual functions)
 */
class IMatrix {
	public:
	/// creates a (nrow x ncol) - matrix
	virtual bool create_matrix(int nrow, int ncol) = 0;

	/// deletes the matrix memory
	virtual bool delete_matrix() = 0;

	/// set selected entries
	virtual bool set_values(int nrows, int* ncols, int* rows, int* cols, double* values) = 0;

	/// add to selected entries
	virtual bool add_values(int nrows, int* ncols, int* rows, int* cols, double* values) = 0;

	/// normalize dirichlet rows
	virtual bool set_dirichletrows(int nrows, int* rows) = 0;

	/// set all entries to value
	virtual bool set(number w) = 0;

	/// finalize
	virtual bool finalize() = 0;

	/// b := A*x (A = this Object)
	virtual	bool apply(IVector&b, IVector& x) = 0;

	/// b := A^T * x (A^T = transposed of this object)
	virtual	bool applyTransposed(IVector&b, IVector& x) = 0;

	/// b := b - A * x (A = this object)
	virtual bool matmul_minus(IVector&b, IVector& x) = 0;

	/// export to ascii - file
	virtual bool printToFile(const char* filename) = 0;

	/// virtual destructor
	virtual ~IMatrix()
	{}
};



/// Linear Solver Interface
/**
 * Interface for linear solvers
 */
class ILinearSolver {
	public:

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
	virtual bool solve(Matrix& A, Vector& x, Vector &b) = 0;


	/// virtual destructor
	virtual ~ILinearSolver()
	{}
};


}

#endif /* __H__LIB_ALGEBRA__ */
