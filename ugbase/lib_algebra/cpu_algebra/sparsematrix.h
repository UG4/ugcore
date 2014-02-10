/**
 * \file sparsematrix.h
 *
 * \author Martin Rupp
 *
 * \date 29.10.2012
 *
 * Goethe-Center for Scientific Computing 2012
 */

#ifndef __H__UG__CPU_ALGEBRA__SPARSEMATRIX__
#define __H__UG__CPU_ALGEBRA__SPARSEMATRIX__

#include "math.h"
#include "common/common.h"
#include "../algebra_common/sparsematrix_util.h"
#include <iostream>
#include <algorithm>
#include "common/util/ostream_util.h"

#include "../algebra_common/connection.h"
#include "../algebra_common/matrixrow.h"
#include "../common/operations_mat/operations_mat.h"

#define PROFILE_SPMATRIX(name) PROFILE_BEGIN_GROUP(name, "SparseMatrix algebra")

#ifndef NDEBUG
#define CHECK_ROW_ITERATORS
#endif

namespace ug{

/// \addtogroup cpu_algebra
///	@{


// example for the variable CRS storage structure:
// say we have:
// rowStart = 0 3 8
// rowEnd = 3 6 11
// rowMax = 3 8 11
// cols ( | marking end of row): 2 5 6 | 2 6 7 x x| 8 9 10

// now insert (0 3): row 0 is full (rowEnd[0]==rowMax[0]), copy it to the end, and insert index
// rowStart = 11 3 8
// rowEnd = 15 6 11
// rowMax = 17 8 11
// cols ( | marking end of row): x x x | 2 6 7 x x| 8 9 10 | 2 3 5 6 x x |

// now insert (1 3): row 1 not full, we can add it
// rowStart 11 3 8
// rowEnd 15 7 11
// rowMax = 17 8 11
// cols : x x x | 2 3 6 7 x | 8 9 10 | 2 3 5 6 x x |

// defragment:
// rowStart 0 4 8
// rowEnd 4 8 11
// rowMax = 4 8 11
// cols : 2 3 5 6 | 2 3 6 7 | 8 9 10


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
 * \param T blocktype
 */
template<typename TValueType> class SparseMatrix
{
public:
	typedef TValueType value_type;
	enum {rows_sorted=true};

	typedef SparseMatrix<value_type> this_type;

public:
	typedef AlgebraicConnection<TValueType> connection;
	typedef MatrixRow<this_type> row_type;
	typedef ConstMatrixRow<this_type> const_row_type;

public:
	// construction etc
	//----------------------

	/// constructor for empty SparseMatrix
	SparseMatrix();
	/// destructor
	virtual ~SparseMatrix () {}


	/**
	 * \brief resizes the SparseMatrix
	 * \param newRows new nr of rows
	 * \param newCols new nr of cols
	 * \return
	 */
	void resize_and_clear(size_t newRows, size_t newCols);
	void resize_and_keep_values(size_t newRows, size_t newCols);

	/**
	 * \brief write in a empty SparseMatrix (this) the transpose SparseMatrix of B.
	 * \param B			the matrix of which to create the transpose of
	 * \param scale		an optional scaling
	 * \return			true on success
	 */
	void set_as_transpose_of(const SparseMatrix<value_type> &B, double scale=1.0);

	/**
	 * \brief create/recreate this as a copy of SparseMatrix B
	 * \param B			the matrix of which to create a copy of
	 * \param scale		an optional scaling
	 * \return			true on success
	 */
	void set_as_copy_of(const SparseMatrix<value_type> &B, double scale=1.0);
	SparseMatrix<value_type> &operator = (const SparseMatrix<value_type> &B)
	{
		set_as_copy_of(B);
		return *this;
	}


public:
	//! calculate dest = alpha1*v1 + beta1*A*w1 (A = this matrix)
	template<typename vector_t>
	void axpy(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const vector_t &w1) const;

	//! calculate dest = alpha1*v1 + beta1*A^T*w1 (A = this matrix)
	template<typename vector_t>
	void axpy_transposed(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const vector_t &w1) const;

	//! calculated dest = beta1*A*w1 . For empty rows, dest will not be changed
	template<typename vector_t>
	void apply_ignore_zero_rows(vector_t &dest,
			const number &beta1, const vector_t &w1) const;

	//! calculated dest = beta1*A*w1 . For empty cols of A (=empty rows of A^T), dest will not be changed
	template<typename vector_t>
	void apply_transposed_ignore_zero_rows(vector_t &dest,
			const number &beta1, const vector_t &w1) const;

	// DEPRECATED!
	//! calculate res = A x
		// apply is deprecated because of axpy(res, 0.0, res, 1.0, beta, w1)
		template<typename Vector_type>
		bool apply(Vector_type &res, const Vector_type &x) const
		{
			axpy(res, 0.0, res, 1.0, x);
			return true;
		}

		//! calculate res = A.T x
		// apply is deprecated because of axpy(res, 0.0, res, 1.0, beta, w1)
		template<typename Vector_type>
		bool apply_transposed(Vector_type &res, const Vector_type &x) const
		{
			axpy_transposed(res, 0.0, res, 1.0, x);
			return true;
		}

		// matmult_minus is deprecated because of axpy(res, 1.0, res, -1.0, x);
		//! calculate res -= A x
		template<typename Vector_type>
		bool matmul_minus(Vector_type &res, const Vector_type &x) const
		{
			axpy(res, 1.0, res, -1.0, x);
			return true;
		}



	/**
	 * \brief check for isolated condition of an index
	 * \param i
	 * \return true if only A[i,i] != 0.0
	 */
	inline bool is_isolated(size_t i) const;

	void scale(double d);
	SparseMatrix<value_type> &operator *= (double d) { scale(d); return *this; }

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

	// finalizing functions
	//----------------------



	inline void check_rc(size_t r, size_t c) const
	{
		UG_ASSERT(r < num_rows() && c < num_cols(), "tried to access element (" << r << ", " << c << ") of " << num_rows() << " x " << num_cols() << " matrix.");
	}

	inline void check_row_modifiable(size_t r) const
	{
#ifdef CHECK_ROW_ITERATORS
		UG_ASSERT(nrOfRowIterators[r] == 0, "row " << r << " is used by an iterator and should not be modified.")
#endif
	}


	//! set matrix to Id*a
	void set(double a);

	/** operator() (size_t r, size_t c) const
	 * access connection (r, c)
	 * \param r row
	 * \param c column
	 * \note if connection (r, c) is not there, returns 0.0
	 * \return SparseMat(r, c)
	 */
	const value_type &operator () (size_t r, size_t c)  const
    {
		check_rc(r, c);
        int j=get_index_const(r, c);
		if(j == -1)
		{
			static value_type v(0.0);
			return v;
		}
        UG_ASSERT(cols[j]==(int)c && j >= rowStart[r] && j < rowEnd[r], "");
        return values[j];
    }

	/** operator() (size_t r, size_t c) const
	 * access or create connection (r, c)
	 * \param r row
	 * \param c column
	 * \note (r,c) is added to sparsity pattern if not already there
	 * use operator()(r,c,bConnectionFound) to prevent
	 * \return SparseMat(r, c)=0.0 if connection created, otherwise SparseMat(r, c)
	 */
	value_type &operator() (size_t r, size_t c)
	{
		check_rc(r, c);
		int j=get_index(r, c);
        UG_ASSERT(j != -1 && cols[j]==(int)c && j >= rowStart[r] && j < rowEnd[r], "");
        return values[j];
    }

public:
	// row functions

	/**
	 * set a row of the matrix. all previous content in this row is destroyed (@sa add_matrix_row).
	 * \param row index of the row to set
	 * \param c pointer to a array of sorted connections of size nr
	 * \param nr number of connections in c
	 * \remark will sort array c
	 */
	void set_matrix_row(size_t row, connection *c, size_t nr);

	/**
	 * adds the connections c to the matrixrow row.
	 * if c has a connection con with con.iIndex=i, and the matrix already has a connection (row, i),
	 * the function will set A(row,i) += con.dValue. otherwise the connection A(row, i) is created
	 * and set to con.dValue.
	 * \param row row to add to
	 * \param c connections ("row") to be added the row.
	 * \param nr number of connections in array c.
	 * \return true on success.
	 * \note you may use double connections in c.
	 */
	void add_matrix_row(size_t row, connection *c, size_t nr);


	//! calculates dest += alpha * A[row, .] v;
	template<typename vector_t>
	inline void mat_mult_add_row(size_t row, typename vector_t::value_type &dest, double alpha, const vector_t &v) const;
public:
	// accessor functions
	//----------------------

	//! returns number of connections of row row.
	inline size_t num_connections(size_t i) const
	{
		if(rowStart[i] == -1) return 0;
		else return rowEnd[i]-rowStart[i];
	}

	//! returns number of rows
	size_t num_rows() const { return rowEnd.size(); }

	//! returns the number of cols
	size_t num_cols() const { return m_numCols; }

	//! returns the total number of connections
	size_t total_num_connections() const { return nnz; }

public:

	// Iterators
	//---------------------------

	// a row_iterator has to suppport
	// operator ++, operator +=, index() const, const value_type &value() const, value_type &value()
	// a const_row_iterator has to suppport
	// operator ++, operator +=, index() const, const value_type &value() const


	/**
	 *  row_iterator
	 *  iterator over a row
	 */
	class row_iterator
    {
        SparseMatrix &A;
        size_t row;
        size_t i;
    public:
        inline void check() const {A.check_row(row, i); }
        row_iterator(SparseMatrix &_A, size_t _row, size_t _i) : A(_A), row(_row), i(_i) { A.add_iterator(row); }
        row_iterator(const row_iterator &other) : A(other.A), row(other.row), i(other.i) { A.add_iterator(row); }
        ~row_iterator() { A.remove_iterator(row); }
        row_iterator *operator ->() { return this; }
        value_type &value() { check(); return A.values[i];   }
        size_t index() const { check(); return A.cols[i]; }
        bool operator != (const row_iterator &o) const { return i != o.i;  }
        void operator ++ () { ++i; }
		void operator += (int nr) { i+=nr; }
		bool operator == (const row_iterator &other) const { return other.i == i; check(); }
    };
    class const_row_iterator
    {
        const SparseMatrix &A;
        size_t row;
        size_t i;
    public:
        inline void check() const {A.check_row(row, i); }
        const_row_iterator(const SparseMatrix &_A, size_t _row, size_t _i) : A(_A), row(_row), i(_i) {A.add_iterator(row);}
        const_row_iterator(const const_row_iterator &other) : A(other.A), row(other.row), i(other.i) { A.add_iterator(row); }
        ~const_row_iterator() { A.remove_iterator(row); }
        const_row_iterator *operator ->() { return this; }
        const value_type &value() const { check(); return A.values[i];   }
        size_t index() const { check(); return A.cols[i];     }
        bool operator != (const const_row_iterator &o) const { return i != o.i; }
        void operator ++ () { ++i; }
        void operator += (int nr) { i+=nr; }
		bool operator == (const const_row_iterator &other) const { return other.i == i; }
    };




	row_iterator         begin_row(size_t r)         { return row_iterator(*this, r, rowStart[r]);  }
    row_iterator         end_row(size_t r)           { return row_iterator(*this, r, rowEnd[r]);  }
    const_row_iterator   begin_row(size_t r) const   { return const_row_iterator(*this, r, rowStart[r]);  }
    const_row_iterator   end_row(size_t r)   const   { return const_row_iterator(*this, r, rowEnd[r]);  }

    row_type 		get_row(size_t r) 		{ return row_type(*this, r); }
    const_row_type 	get_row(size_t r) const { return const_row_type(*this, r); }

public:
	// connectivity functions
	//-------------------------

    bool has_connection(size_t r, size_t c) const
    {
    	check_rc(r, c);
    	bool bFound;
    	get_connection(r, c, bFound);
    	return bFound;
    }

	/**
	 * \param r index of the row
	 * \param c index of the column
	 * \return a const_row_iterator to the connection A(r,c) if existing, otherwise end_row(row)
	 */
	row_iterator get_iterator_or_next(size_t r, size_t c)
	{
		check_rc(r, c);
		if(rowStart[r] == -1 || rowStart[r] == rowEnd[r])
        	return end_row(r);
        else
        {
        	int j=get_index_internal(r, c);
        	if(j > maxValues) return end_row(r);
        	else return row_iterator(*this, r, j);
        }
    }

	/**
	 * \param r index of the row
	 * \param c index of the column
	 * \return a const_row_iterator to the connection A(r,c) if existing, otherwise end_row(row)
	 */
	const_row_iterator get_connection(size_t r, size_t c, bool &bFound) const
	{
		check_rc(r, c);
        int j=get_index_const(r, c);
		if(j != -1)
		{
			bFound = true;
			return const_row_iterator(*this, r, j);
		}
		else
		{
			bFound = false;
			return end_row(r);
		}
    }
	/**
	 * \param r index of the row
	 * \param c index of the column
	 * \return a row_iterator to the connection A(r,c) if existing, otherwise end_row(row)
	 */
	row_iterator get_connection(size_t r, size_t c, bool &bFound)
	{
		check_rc(r, c);
		int j=get_index_const(r, c);
		if(j != -1)
		{
			bFound = true;
			return row_iterator(*this, r, j);
		}
		else
		{
			bFound = false;
			return end_row(r);
		}
	}

	/**
	 * \param r index of the row
	 * \param c index of the column
	 * \return a const_row_iterator to the connection A(r,c) if existing, otherwise end_row(row)
	 */
	const_row_iterator get_connection(size_t r, size_t c) const
	{
		bool b;
		return get_connection(r, c, b);
	}
	/**
	 * \param r index of the row
	 * \param c index of the column
	 * \return a row_iterator to the connection A(r,c)
	 * \remark creates connection if necessary.
	 */
	row_iterator get_connection(size_t r, size_t c)
	{
		check_rc(r, c);
		assert(bNeedsValues);
        int j=get_index(r, c);
		return row_iterator(*this, r, j);
	}


	void defragment()
    {
		if(num_rows() != 0 && num_cols() != 0)
			copyToNewSize(nnz);
    }

	void defragment() const
	{
		(const_cast<this_type*>(this))->defragment();
	}

	/**
	 * copies the matrix to the standard CRS format
	 * @param numRows   	(out) num rows of A
	 * @param numCols		(out) num rows of A
	 * @param argValues		(out) value_type vector with non-zero values
	 * @param argRowStart   (out) row i is from argRowStart[i] to argRowStart[i+1]
	 * @param argColInd		(out) argColInd[i] is colum index of nonzero i
	 */
	void copy_crs(size_t &numRows, size_t &numCols,
			std::vector<value_type> &argValues, std::vector<int> &argRowStart,
			std::vector<int> &argColInd) const
	{
		numRows = num_rows();
		numCols = num_cols();
		defragment();
		argValues = values;
		argRowStart = rowStart;
		argColInd = cols;
	}

	/**
	 * returns pointers to CRS format. note that these are only valid as long
	 * as the matrix is not modified.
	 * @param numRows   	(out) num rows of A
	 * @param numCols		(out) num rows of A
	 * @param pValues		(out) value_type vector with non-zero values
	 * @param pRowStart   (out) row i is from pRowStart[i] to pRowStart[i+1]
	 * @param pColInd		(out) pColInd[i] is colum index of nonzero i
	 */
	void get_crs(size_t &numRows, size_t &numCols,
			value_type *&pValues, size_t *pRowStart, size_t *pColInd, size_t &nnz) const
	{
		defragment();
		pValues = &values[0];
		pRowStart = &rowStart[0];
		pColInd = &cols[0];
		numRows = num_rows();
		numCols = num_cols();
		nnz = total_num_connections();
	}


public:
	// output functions
	//----------------------

	void print(const char * const name = NULL) const;
	void printtype() const;

	void print_to_file(const char *filename) const;
	void printrow(size_t row) const;

	friend std::ostream &operator<<(std::ostream &out, const SparseMatrix &m)
	{
		out << "SparseMatrix " //<< m.name
		<< " [ " << m.num_rows() << " x " << m.num_cols() << " ]";
		return out;
	}


	void p() const { print(); } // for use in gdb
	void pr(size_t row) const {printrow(row); } // for use in gdb





private:
	// private functions

	void add_iterator(size_t row) const
	{
#ifdef CHECK_ROW_ITERATORS
		nrOfRowIterators[row]++;
#endif
		iIterators++;
	}
	void remove_iterator(size_t row) const
	{
#ifdef CHECK_ROW_ITERATORS
		nrOfRowIterators[row]--;
		UG_ASSERT(nrOfRowIterators[row] >= 0, row);
#endif
		iIterators--;
		UG_ASSERT(iIterators >= 0, row);

	}
	inline void check_row(size_t row, int i) const
	{
		UG_ASSERT(i < rowEnd[row] && i >= rowStart[row], "row iterator row " << row << " pos " << i << " out of bounds [" << rowStart[row] << ", " << rowEnd[row] << "]");
	}
    void assureValuesSize(size_t s);
    size_t get_nnz() const { return nnz; }

private:
	// disallowed operations (not defined):
	//---------------------------------------
	SparseMatrix(SparseMatrix&); ///< disallow copy operator


protected:
	int get_index_internal(size_t row, int col) const;
    int get_index_const(int r, int c) const;
    int get_index(int r, int c);
    void copyToNewSize(size_t newSize)
    {
    	copyToNewSize(newSize, num_cols());
    }
    void copyToNewSize(size_t newSize, size_t maxCols);
	void check_fragmentation() const;
	int get_nnz_max_cols(size_t maxCols);


protected:
    std::vector<int> rowStart;
    std::vector<int> rowEnd;
    std::vector<int> rowMax;
    std::vector<int> cols;
    size_t fragmented;
    size_t nnz;
    bool bNeedsValues;

    std::vector<value_type> values;
    int maxValues;
    int m_numCols;
    mutable int iIterators;

#ifdef CHECK_ROW_ITERATORS
    mutable std::vector<int> nrOfRowIterators;
#endif
};


template<typename T>
struct matrix_algebra_type_traits<SparseMatrix<T> >
{
	enum{
		type=MATRIX_USE_ROW_FUNCTIONS
	};
};

//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultTransposedAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const SparseMatrix<matrix_t> &A1, const vector_t &w1)
{
	A1.axpy_transposed(dest, alpha1, v1, beta1, w1);
}








// end group cpu_algebra
/// \}

} // namespace ug

//#include "matrixrow.h"
#include "sparsematrix_impl.h"
#include "sparsematrix_print.h"

#endif
