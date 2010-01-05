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

/// Vector Block Interface
class IVectorBlock{
	public:
		/// access component i of this block
		virtual number& operator() (uint block_i) = 0;

};


/// SubVector Inteface
/**
 * Interface for SubVectors (using virtual functions)
 */
class ISubVector{
	public:
		/// constructor
		/**
		 * @param[in] indices Array of required indices
		 * @param[in] unknowns Array of number of unknowns for each index
		 * @param[in] num_indices number of indices required (size of arrays)
		 */
		ISubVector(int* indices, int* unknowns, int num_indices);

		/// returns entry i (which may be a small block)
		virtual IVectorBlock& operator() (uint local_i) = 0;

		///
		virtual uint get_global_index(uint local_i) = 0;

};


/// Vector Interface
/**
 * Interface for Vectors (using virtual functions)
 */
class IVector {
	public:
		/// creates a vector of length nentries (allocates memory, etc.)
		/**
		 * The length of the vector is the number of blocks.
		 */
		virtual bool create(uint nentries) = 0;

		/// creates a vector with the entry-pattern of vector v
		virtual bool create(const IVector& v) = 0;

		/// deletes memory
		virtual bool destroy() = 0;

		/// adds subvector to vector
		virtual bool add(const ISubVector& SubV) = 0;

		/// set subvector entries in this vector
		virtual bool set(const ISubVector& SubV) = 0;

		/// get values of a subvector
		virtual bool get(ISubVector& SubV) const = 0;

		/// finalize - internal memory rearrangement
		virtual bool finalize()
		{}

		/// print to ascii - file
		virtual bool printToFile(const char* filename) = 0;

		/// virtual destructor
		virtual ~ArneVector()
		{}

		virtual IVector& operator+= (const IVector& v) = 0;
		virtual IVector& operator-= (const IVector& v) = 0;
		virtual IVector& operator= (const IVector& v) = 0;

		/// set all entries of the vector to the value
		virtual bool set(number w) = 0;
		virtual bool operator= (number w) = 0;
		virtual bool operator*= (number w) = 0;

		/// scalar product
		virtual number operator *(const IVector& v) = 0;

		virtual number one_norm() = 0;
		virtual number two_norm() = 0;

		virtual uint length() = 0;

	private:
		// do not allow copy constructor
		IVector(const IVector& v);

};

/// Matrix Block Interface
class IMatrixBlock{
	public:
		/// access component i of this block
		virtual number& operator() (uint block_i, uint block_j) = 0;

};


/// SubMatrix Inteface
/**
 * Interface for SubMatrices (using virtual functions)
 */
class ISubMatrix{
	public:
		/// constructor
		/**
		 * creates a quadratic submatrix for requested indices,
		 * entries may be small blocks
		 *
		 * @param[in] indices Array of required indices
		 * @param[in] unknowns Array of number of unknowns for each index
		 * @param[in] num_indices number of indices required (size of arrays)
		 */
		ISubMatrix(int* indices, int* unknowns, int num_indices);

		/// returns entry (i,j) (which may be a small block)
		virtual IMatrixBlock& operator() (uint local_i, uint local_j) = 0;

		///
		virtual uint get_global_index(uint local_i) = 0;

};


/// Matrix Interface
/**
 * Interface for Matrices (using virtual functions)
 */
class IMatrix {
	public:
	/// creates a (num_rows x num_ncols) - matrix
	virtual bool create(int num_rows, int num_cols) = 0;

	/// creates a matrix with the same structure as matrix m
	virtual bool create(IMatrix& m) = 0;

	/// deletes the matrix memory
	virtual bool destroy() = 0;

	/// add submatrix entries
	virtual bool add(const ISubMatrix& SubM) = 0;

	/// set submatrix entries
	virtual bool set(const ISubMatrix& SubM) = 0;

	/// get submatrix entries
	virtual bool get(ISubMatrix& SubM) = 0;

	/// normalize dirichlet rows
	/**
	 * sets the nrows indicated in the rows-array to identity matrix
	 * (i.e. coupling 1 to itself, coupling 0 to all others)
	 * In the case of a block matrix, the diagonal coupling is an identity matrix
	 *
	 * Attention: the resulting matrix is no longer symmetric.
	 */
	virtual bool set_dirichlet_rows(int nrows, int* rows) = 0;

	// TODO: add eliminate_dirichlet_cols

	/// set all entries to value
	virtual bool operator= (number w) = 0;
	virtual bool set(number w) = 0;

	/// finalize
	bool finalize()
	{};

	// TODO: Add matrix norms here.

	/// b := A*x (A = this Object)
	virtual	bool apply(IVector&b, const IVector& x) const = 0;

	/// b := A^T * x (A^T = transposed of this object)
	virtual	bool applyTransposed(IVector&b, const IVector& x) const = 0;

	/// b := b - A * x (A = this object)
	virtual bool matmul_minus(IVector&b, const IVector& x) const = 0;

	/// export to ascii - file
	virtual bool printToFile(const char* filename) const = 0;

	/// virtual destructor
	virtual ~IMatrix()
	{}

	private:
	IMatrix(const IMatrix& m);
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
