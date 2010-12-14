/**
 * \file sparsematrix.h
 *
 * \author Martin Rupp
 *
 * \date 26.11.2009
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__CPU_ALGEBRA__SPARSEMATRIX__
#define __H__UG__CPU_ALGEBRA__SPARSEMATRIX__

#include "math.h"


#include "../common/template_expressions.h"
#include "vector.h"


namespace ug{

/// \addtogroup lib_algebra
///	@{

template<typename value_type> class matrixrow;
template<typename vec_type> class Vector;

/** SparseMatrix
 *  \brief sparse matrix for big, variable sparse matrices.
 *
 *  matrix is stored independent row-wise
 *  When doing discretisation, use the add set and get methods
 *  for dealing with submatrices of A.
 *  For other things you can use the row iterators or
 *  operator()-methods.
 *
 * \sa matrixrow, CreateAsMultiplyOf
 * \param T blocktype
 */
template<typename T>
class SparseMatrix
{
public:
	typedef T value_type;
	enum {rows_sorted=true};

	typedef matrixrow<value_type> row_type;
	typedef matrixrow<value_type> matrixrow_type;


public:
	struct connection
	{
		size_t iIndex;		// index to
		value_type dValue; // smallmatrix value;

		void print(){cout << *this;}
		friend ostream &operator<<(ostream &output, const connection &c)
		{
			output << "(" << c.iIndex << "-> ";
			cout << c.dValue;
			cout << ")";
			return output;
		}

		void operator = (const connection &other)
		{
			iIndex = other.iIndex;
			dValue = other.dValue;
		}

		int operator < (const connection &c) const
		{
			return iIndex < c.iIndex;
		}
	};

public:
	// construction etc
	//----------------------

	// constructor for empty SparseMatrix
	SparseMatrix();
	// destructor
	~SparseMatrix ();


	bool create(size_t _rows, size_t _cols);
	bool resize(size_t newRows, size_t newCols);
	bool destroy();

	//! create this as a transpose of SparseMatrix B
	bool create_as_transpose_of(const SparseMatrix &B, double scale=1.0);

	//! create/recreate this as a copy of SparseMatrix B
	bool create_as_copy_of(const SparseMatrix<T> &B, double scale=1.0);
	SparseMatrix<T> &operator = (const SparseMatrix<T> &B)
	{
		create_as_copy_of(B);
		return *this;
	}


public:
	// finalizing functions
	//----------------------
	void defrag();
	void definalize();
	void finalize()
	{
		defrag();
	}
	bool is_finalized() const;

private:
	// safe_set_connections
	/**
	 * "safe" way to set a connection, since when cons[row] is in the big consecutive consmem-array,
	 * you mustnt delete[] mem.
	 */
	void safe_set_connections(size_t row, connection *mem) const;
	//! returns true if the row is stored inside the consecutive mem.
	bool in_consmem(size_t row) const;

private:
	// disallowed operations (not defined):
	//---------------------------------------
	SparseMatrix(SparseMatrix&); ///< disallow copy operator

public:
	// general functions
	//----------------------
	bool set_dirichlet_rows(const size_t *pDirichletRows, size_t iNr);

	//! calculate res = A x
	template<typename Vector_type>
	bool apply(Vector_type &res, const Vector_type &x) const;

	//! calculate res = A.T x
	template<typename Vector_type>
	bool apply_transposed(Vector_type &res, const Vector_type &x) const;

	//! calculate dest = alpha1*v1 + beta1*A*w1 (A = this matrix)
	template<typename vector_t>
	bool axpy(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const vector_t &w1) const;

	//! calculate dest = alpha1*v1 + beta1*A^T*w1 (A = this matrix)
	template<typename vector_t>
	bool axpy_transposed(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const vector_t &w1) const;

	//! calculate res -= A x
	template<typename Vector_type>
	bool matmul_minus(Vector_type &res, const Vector_type &x) const;

	//! get Diagonal A_[i,i] of matrix
	inline const value_type &get_diag(size_t i) const;
	inline value_type &get_diag(size_t i);

	//! isUnconnected: true if only A[i,i] != 0.0.
	inline bool is_isolated(size_t i) const;

	bool scale(double d);
	SparseMatrix<T> &operator *= (double d) { scale(d); return *this; }

public:
	// submatrix set/get functions
	//-------------------------------

	/** Add a local matrix
	 *
	 * The local matrix type must declare the following members:
	 * - num_rows()
	 * - num_cols()
	 * - row_index(size_t i)
	 * - col_index(size_t j)
	 * - operator()(size_t i, size_t j)
	 * so that mat(i,j) will go to SparseMat(mat.row_index(i), mat.col_index(j))
	 * \param mat the whole local matrix type
	 */
	template<typename M>
	void add(const M &mat);
	template<typename M>
	//! set local matrix \sa add
	void set(const M &mat);
	//! get local matrix \sa add
	template<typename M>
	void get(M &mat) const;

	/** Add a local matrix
	 *
	 * The matrix type must declare the following members:
	 * - num_rows()
	 * - num_cols()
	 * - operator()(size_t i, size_t j)
	 * so that mat(i, j) will go to SparseMat(row[i], col[j])
	 * \param mat a small matrix type
	 * \param rows
	 * \param cols
	 */
	template<typename M>
	void add(const M &mat, size_t *rows, size_t *cols);
	//! set local matrix \sa add
	template<typename M>
	void set(const M &mat, size_t *rows, size_t *cols);
	//! get local matrix \sa add
	template<typename M>
	void get(M &mat, size_t *rows, size_t *cols) const;


	//! set matrix to Id*a
	bool set(double a);


	/** operator() (size_t r, size_t c) const
	 * access connection (r, c)
	 * \param r row
	 * \param c column
	 * \note it is assert'ed that connection (r,c) is there
	 * use operator()(r,c,bConnectionFound) to check.
	 * \return SparseMat(r, c)
	 */
	const value_type &operator() (size_t r, size_t c) const;

	/** operator() (size_t r, size_t c) const
	 * access or create connection (r, c)
	 * \param r row
	 * \param c column
	 * \note (r,c) is added to sparsity pattern if not already there
	 * use operator()(r,c,bConnectionFound) to prevent
	 * \return SparseMat(r, c)=0.0 if connection created, otherwise SparseMat(r, c)
	 */
	value_type &operator() (size_t r, size_t c);

	/** operator() (size_t r, size_t c, bool &bConnectionFound) const
	 * access SparseMat(r, c) and check if connection is there
	 * \param r row
	 * \param c column
	 * \param bConnectionFound false if connection couldnt be found
	 * \note the connection (r, c) is not created if not already there.
	 * \return value_type(0.0) if not found, otherwise SparseMat(r, c)
	 */
	const value_type &operator() (size_t r, size_t c, bool &bConnectionFound) const;
	value_type &operator() (size_t r, size_t c, bool &bConnectionFound);

	// for other manipulation/accessor functions see matrixrow functions,
	// that is A[i].matrixrowfunction(params).

public:
	// accessor functions
	//----------------------

	size_t num_rows() const { return rows; }
	size_t num_cols() const { return cols; }

	size_t total_num_connections() const { return iTotalNrOfConnections; }

private:
	enum get_connection_nr_flag
	{
		// find the first (in terms of distance) which is
		LESS =1,
		LESS_EQUAL = 2,
		EQUAL = 3,
		GREATER_EQUAL = 4,
		GREATER = 5
		// than the connection (r,c)
	};

	bool get_connection_nr(size_t r, size_t c, size_t &nr, get_connection_nr_flag flag=EQUAL) const;

	template<get_connection_nr_flag flag>
	bool get_connection_nr_templ(size_t r, size_t c, size_t &nr) const;
public:
	// row functions
	//----------------------

	//! accessor functions for artificial matrixrow-object (= just wrapper with A and row)
	inline const matrixrow_type getrow(size_t i) const;
	inline const matrixrow_type operator [] (size_t i) const;

	//! remove zero entries of SparseMatrix (experimental)
	void remove_zeros(size_t row);

	/** set a row of the matrix. all previous content in this row is destroyed (@sa add_matrix_row).
	 * \param row index of the row to set
	 * \param c pointer to a array of sorted connections of size nr
	 * \param nr number of connections in c
	 * \remark connections have to be sorted
	 */
	void set_matrix_row(size_t row, connection *c, size_t nr);

	/** add_matrix_row
	 *  add a row to a matrix row.
	 * \param row index of the row to set
	 * \param c pointer to a array of sorted connections of size nr
	 * \param nr number of connections in c
	 * \remark if we get new connections, matrix is definalized.
	 * \remark connections have to be sorted
	 */
	void add_matrix_row(size_t row, connection *c, size_t nr);

	//! returns number of connections of row row.
	inline size_t num_connections(size_t row) const;

	template<typename vector_t>
	inline void mat_mult_add_row(size_t row, typename vector_t::value_type &dest, double alpha, const vector_t &v) const;

public:
	// output functions
	//----------------------

	void print(const char * const name = NULL) const;
	void printtype() const;

	void print_to_file(const char *filename) const;
	void printrow(size_t row) const;

	friend ostream &operator<<(ostream &out, const SparseMatrix &m)
	{
		out << "SparseMatrix " //<< m.name
		<< " [ " << m.rows << " x " << m.cols << " ]";
		return out;
	}


	void p() const { print(); } // for use in gdb
	void pr(size_t row) const {printrow(row); } // for use in gdb

public:

	// Iterators
	//---------------------------

	// const_RowIterator

	/** cRowIterator
	 * const iterator over a row
	 * \note due to cLowerLeftIterator, there is no endRow or similar, use cRowIterator::isEnd()
	 */
	class cRowIterator
	{
	public:
		//const SparseMatrix<value_type> &A;
		const connection * pEnd;
		const connection * p;
	public:
		inline cRowIterator(const SparseMatrix<value_type> &A, size_t row)
			: pEnd(A.pRowEnd[row]), p(A.pRowStart[row])
			  { }
		inline cRowIterator(const cRowIterator &other)
			: pEnd(other.pEnd), p(other.p)
			{ }

		inline const connection &operator *() const {return *p;}

		inline size_t index() const { return p->iIndex; }
		inline const value_type &value() const { return p->dValue; }

		inline void operator ++() {	++p; }
		inline void operator += (int nr) { p+=nr;}

		inline bool isEnd() const { return p >= pEnd; } // remove this
	};

	/** rowIterator
	 * iterator over a row
	 */
	class rowIterator
	{
	public:
		connection * pEnd;
		connection * p;
	public:
		inline rowIterator(SparseMatrix<value_type> &A, size_t row_)
		: pEnd(A.pRowEnd[row_]), p(A.pRowStart[row_])
		  { }
		inline rowIterator(const cRowIterator &other)
		: pEnd(other.pEnd), p(other.p)
		  { }

		inline connection &operator *() const {return *p;}

		inline size_t index() const { return p->iIndex; }
		inline value_type &value() const { return p->dValue; }

		inline void operator ++() {	++p; }
		inline void operator += (int nr) { p+=nr;}

		inline bool isEnd() const { return p >= pEnd; } // remove this
	};

	/** cLowerLeftIterator
	 * iterator only over connections with c->iIndex < row
	 */
	class cLowerLeftIterator : public cRowIterator
	{
	private:
		size_t row;
	public:
		cLowerLeftIterator(const SparseMatrix<value_type> &A, size_t row_)
		: cRowIterator(A, row_), row(row_)
		  {	}

		inline bool isEnd() const { return this->p >= this->pEnd || this->p->iIndex >= row; }
	};

	/** cUpperRightIterator
	 * iterator only over connections with c->iIndex > row
	 */
	class cUpperRightIterator : public cRowIterator
	{
	public:
		cUpperRightIterator(const SparseMatrix<value_type> &A, size_t row_) : cRowIterator(A, row_)
		{
			size_t nr=0;
			if(A.get_connection_nr(row_, row_, nr, GREATER))
				this->p+=nr;
			else
				this->p = this->pEnd;
		}
	};

	cRowIterator beginRow(size_t row) const
	{
		return cRowIterator(*this, row);
	}

	rowIterator beginRow(size_t row)
	{
		return rowIterator(*this, row);
	}

	cLowerLeftIterator beginLowerLeftRow(size_t row)  const
	{
		return cLowerLeftIterator(*this, row);
	}

	cUpperRightIterator beginUpperRightRow(size_t row) const
	{
		return cUpperRightIterator(*this, row);
	}

public:
	// connectivity functions
	//-------------------------
	cRowIterator get_connection(size_t r, size_t c, bool &bFound) const;
	rowIterator get_connection(size_t r, size_t c, bool &bFound);

	cRowIterator get_connection(size_t r, size_t c) const;
	rowIterator get_connection(size_t r, size_t c);

public:
	//     data
	//----------------------

private:
	size_t rows;						//!< nr of rows
	size_t cols;						//!< nr of cols
	connection **pRowStart;				//< pointers to array of beginning of connections of each row
	connection **pRowEnd;				//< pointers to array of ends of connections of each row

	size_t iTotalNrOfConnections;		//!< number of connections ("non-zeros")
	size_t bandwidth;					//!< bandwidth (experimental)

	size_t estimatedRowSize;			//!< estimated length of each row
	size_t *iMaxNrOfConnections;		//!< max nr of connections for row [i]. TODO.


	connection *consmem;				//!< consecutive memory for the connections
	size_t consmemsize;					//!< size of the consecutive memory for connections
	size_t iFragmentedMem;				//!< size of connections memory not in consmem

	bool m_bIgnoreZeroes;				//!< if set, add and set wont set connections which are zero

	friend class matrixrow<value_type>;
};


template<typename T>
struct matrix_algebra_type_traits<SparseMatrix<T> >
{
	static const int type = MATRIX_USE_ROW_FUNCTIONS;
};

//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultTransposedAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const SparseMatrix<matrix_t> &A1, const vector_t &w1)
{
	A1.axpy_transposed(dest, alpha1, v1, beta1, w1);
}

///	@}
} // namespace ug

#include "matrixrow.h"
#include "sparsematrix_impl.h"
#include "sparsematrix_util.h"
#include "sparsematrix_print.h"

#endif
