/**
 * \file sparsematrix.h
 *
 * \author Martin Rupp
 *
 * \date 26.11.2009
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */

#ifndef __H__UG__CRS_ALGEBRA__SPARSEMATRIX__
#define __H__UG__CRS_ALGEBRA__SPARSEMATRIX__

#include "math.h"
#include "common/common.h"

#include "../common/operations_mat/operations_mat.h"

namespace ug{

/// \addtogroup lib_algebra
///	@{


/** SparseMatrix
 *  \brief sparse matrix for big, variable sparse matrices.
 *
 * \sa matrixrow, CreateAsMultiplyOf
 * \param T blocktype
 */
template<typename TValueType> class CRSSparseMatrix
{
public:
	typedef TValueType value_type;
	enum {rows_sorted=true};

	typedef SparseMatrix<value_type> this_type;
	typedef matrixrow<value_type> row_type;
	typedef matrixrow<value_type> matrixrow_type;


	
public:
	struct connection
	{
		size_t iIndex;		// index to
		value_type dValue; // smallmatrix value;
		connection() {}
		connection(size_t i, const value_type &v)
		: iIndex(i), dValue(v) {}
				
		void print(){std::cout << *this;}
		friend std::ostream &operator<<(std::ostream &output, const connection &c)
		{
			output << "(" << c.iIndex << "-> ";
			output << c.dValue;
			output << ")";
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

	/// constructor for empty SparseMatrix
	CRSSparseMatrix();
	/// destructor
	virtual ~CRSSparseMatrix () {}


	/**
	 * \brief resizes the SparseMatrix
	 * \param newRows new nr of rows
	 * \param newCols new nr of cols
	 * \return
	 */
	bool resize(size_t newRows, size_t newCols);


	/**
	 * \brief write in a empty SparseMatrix (this) the transpose SparseMatrix of B.
	 * \param B			the matrix of which to create the transpose of
	 * \param scale		an optional scaling
	 * \return			true on success
	 */
	bool set_as_transpose_of(const SparseMatrix<value_type> &B, double scale=1.0);

	/**
	 * \brief create/recreate this as a copy of SparseMatrix B
	 * \param B			the matrix of which to create a copy of
	 * \param scale		an optional scaling
	 * \return			true on success
	 */
	bool set_as_copy_of(const CRSSparseMatrix<value_type> &B, double scale=1.0);
	SparseMatrix<value_type> &operator = (const SparseMatrix<value_type> &B)
	{
		set_as_copy_of(B);
		return *this;
	}


public:
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


	// DEPRECATED!
	//! calculate res = A x
		// apply is deprecated because of axpy(res, 0.0, res, 1.0, beta, w1)
		template<typename Vector_type>
		bool apply(Vector_type &res, const Vector_type &x) const
		{
			return axpy(res, 0.0, res, 1.0, x);
		}

		//! calculate res = A.T x
		// apply is deprecated because of axpy(res, 0.0, res, 1.0, beta, w1)
		template<typename Vector_type>
		bool apply_transposed(Vector_type &res, const Vector_type &x) const
		{
			return axpy_transposed(res, 0.0, res, 1.0, x);
		}

		// matmult_minus is deprecated because of axpy(res, 1.0, res, -1.0, x);
		//! calculate res -= A x
		template<typename Vector_type>
		bool matmul_minus(Vector_type &res, const Vector_type &x) const
		{
			return axpy(res, 1.0, res, -1.0, x);
		}



	/**
	 * \brief check for isolated condition of an index
	 * \param i
	 * \return true if only A[i,i] != 0.0
	 */
	inline bool is_isolated(size_t i) const;

	bool scale(double d);
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
	const value_type &operator () (size_t r, size_t c)  const
    {
        int j=get_index_const(r, c);
        UG_ASSERT(j != -1, "");
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
		int j=get_index(r, c);
        UG_ASSERT(j != -1, "");
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

	//! returns number of connections of row row.
	inline size_t num_connections(size_t i) const
	{
		if(rowStart[i] == -1) return 0;
		else return rowEnd[i]-rowStart[i];
	}

	//! calculates dest += alpha * A[row, .] v;
	template<typename vector_t>
	inline void mat_mult_add_row(size_t row, typename vector_t::value_type &dest, double alpha, const vector_t &v) const;
public:
	// accessor functions
	//----------------------

	//! returns number of rows
	size_t num_rows() const { return rowEnd.size(); }
	
	//! returns the number of cols
	size_t num_cols() const { return m_numCols; }

	//! returns the total number of connections
	size_t total_num_connections() const { return nnz; }

public:

	// Iterators
	//---------------------------

	// const_row_iterator


	//typedef const connection * const_row_iterator;
	//typedef connection * const_row_iterator;
	/** const_row_iterator
	 * const iterator over a row
	 */

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
        CRSSparseMatrix &A;
        size_t i;
    public:
        row_iterator(CRSSparseMatrix &_A, size_t _i) : A(_A), i(_i) {}
        value_type &value() { return A.values[i];   }
        size_t index() const { return A.cols[i];     }
        bool operator != (const row_iterator &o) const { return i != o.i; }
        void operator ++ () { ++i; }
		void operator += (int nr) { i+=nr; }
		bool operator == (const row_iterator &other) const { return other.i == i; }
    };
    class const_row_iterator
    {
        const CRSSparseMatrix &A;
        size_t i;
    public:
        const_row_iterator(const CRSSparseMatrix &_A, size_t _i) : A(_A), i(_i) {}
        const value_type &value() const { return A.values[i];   }
        size_t index() const { return A.cols[i];     }
        bool operator != (const const_row_iterator &o) const { return i != o.i; }
        void operator ++ () { ++i; }
		bool operator == (const const_row_iterator &other) const { return other.i == i; }
    };
	
	


	row_iterator         begin_row(size_t r)         { return row_iterator(*this, rowStart[r]);  }
    row_iterator         end_row(size_t r)           { return row_iterator(*this, rowEnd[r]);  }
    const_row_iterator   begin_row(size_t r) const   { return const_row_iterator(*this, rowStart[r]);  }
    const_row_iterator   end_row(size_t r)   const   { return const_row_iterator(*this, rowEnd[r]);  }
	
public:
	// connectivity functions
	//-------------------------


	/**
	 * \param r index of the row
	 * \param c index of the column
	 * \return a const_row_iterator to the connection A(r,c) if existing, otherwise end_row(row)
	 */
	const_row_iterator get_connection(size_t r, size_t c, bool &bFound) const
	{
        int j=get_index_const(r, c);
		if(j != -1)
		{
			bFound = true;
			return const_row_iterator(*this, j);
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
		int j=get_index_const(r, c);
		if(j != -1)
		{
			bFound = true;
			return row_iterator(*this, j);
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
		assert(bNeedsValues);
        int j=get_index(r, c);
		return row_iterator(*this, j);
	}
	
	
	void defragment()
    {
        copyToNewSize(nnz);
    }

public:
	// output functions
	//----------------------

	void print(const char * const name = NULL) const;
	void printtype() const;

	void print_to_file(const char *filename) const;
	void printrow(size_t row) const;

	friend std::ostream &operator<<(std::ostream &out, const CRSSparseMatrix &m)
	{
		out << "CRSSparseMatrix " //<< m.name
		<< " [ " << m.rows << " x " << m.cols << " ]";
		return out;
	}


	void p() const { print(); } // for use in gdb
	void pr(size_t row) const {printrow(row); } // for use in gdb

private:
	// disallowed operations (not defined):
	//---------------------------------------
	CRSSparseMatrix(CRSSparseMatrix&); ///< disallow copy operator

	
	
    
    void assureValuesSize(int s)
    {
        if(s < cols.size()) return;
        int newSize = nnz*2;
        if(newSize < s) newSize = s;
        copyToNewSize(newSize);
        
    }
 size_t get_nnz() const { return nnz; }

protected:
	int get_index_internal(size_t row, size_t col) const
    {
        assert(rowStart[row] != -1);
        int l = rowStart[row], r = rowEnd[row], mid=0;
        while(l < r)
        {
            mid = (l+r)/2;
            if(cols[mid] < col)
                l = mid+1;
            else if(cols[mid] > col)
                r = mid-1;
            else
                return mid;
        }
        mid = (l+r)/2;
        if(mid < rowStart[row])
            return rowStart[row];
        if(mid == rowEnd[row] || col <= cols[mid])
            return mid;
        else return mid+1;
    }
    
    
    int get_index_const(size_t r, size_t c) const
    {
        if(rowStart[r] == -1 || rowStart[r] == rowEnd[r]) return -1;
        size_t index=get_index_internal(r, c);
        if(index < maxValues && cols[index] == c)
            return index;
        else
            return -1;
    }
    
    
    int get_index(size_t r, size_t c)
    {
        if(rowStart[r] == -1 || rowStart[r] == rowEnd[r])
        {
            assureValuesSize(maxValues+1);
            rowStart[r] = maxValues;
            rowEnd[r] = maxValues+1;
            rowMax[r] = maxValues+1;
            if(bNeedsValues) values[maxValues] = 0.0;
            cols[maxValues] = c;
            maxValues++;
            nnz++;
            return maxValues-1;
        }
        
        /*    for(int i=rowStart[r]; i<rowEnd[r]; i++)
         if(cols[i] == c)
         return i;*/
        size_t index=get_index_internal(r, c);
        if(index < maxValues && cols[index] == c)
            return index;
        assert(index == rowEnd[r] || cols[index] > c);
        
        if(rowEnd[r] == rowMax[r])
        {
            int newSize = (rowEnd[r]-rowStart[r])*2;
            if(maxValues+newSize > cols.size())
            {
                assureValuesSize(maxValues+newSize);
                index=get_index_internal(r, c);
            }
            fragmented += rowEnd[r]-rowStart[r];
            index = index-rowStart[r]+maxValues;
            size_t j=rowEnd[r]-rowStart[r]+maxValues;
            if(rowEnd[r] != 0)
                for(int i=rowEnd[r]-1; i>=rowStart[r]; i--, j--)
                {
                    if(j==index) j--;
                    if(bNeedsValues) values[j] = values[i];
                    cols[j] = cols[i];
                    if(i==rowStart[r]) break;
                }
            rowEnd[r] = maxValues+rowEnd[r]-rowStart[r]+1;
            rowStart[r] = maxValues;
            rowMax[r] = maxValues+newSize;
            maxValues += newSize;
        }
        else
        {
            if(rowEnd[r] != 0)
                for(int i=rowEnd[r]-1; i>=index; i--)
                {
                    if(bNeedsValues) values[i+1] = values[i];
                    cols[i+1] = cols[i];
                    if(i==index) break;
                }
            rowEnd[r]++;
        }
        if(bNeedsValues) values[index] = 0.0;
        cols[index] = c;
        assert(index >= rowStart[r] && index < rowEnd[r]);
        nnz++;
        return index;
	    
    }
    void copyToNewSize(size_t newSize)
    {
        //std::cout << "increasing from " << cols.size() << " to " << newSize << "\n";
        std::vector<double> v(newSize);
        std::vector<int> c(newSize);
        size_t j=0;
        for(size_t r=0; r<num_rows(); r++)
        {
            if(rowStart[r] == -1)
				rowStart[r] = rowEnd[r] = rowMax[r] = j;
			else
			{
				size_t start=j;
				for(size_t k=rowStart[r]; k<rowEnd[r]; k++, j++)
				{
					if(bNeedsValues) v[j] = values[k];
					c[j] = cols[k];
				}
				rowStart[r] = start;
				rowEnd[r] = rowMax[r] = j;
			}
        }
        rowStart[num_rows()] = rowEnd[num_rows()-1];
        fragmented = 0;
        maxValues = j;
        if(bNeedsValues) values = v;
        cols = c;
    }

	
protected:
    std::vector<int> rowStart;
    std::vector<int> rowEnd;
    std::vector<int> rowMax;
    std::vector<int> cols;
    size_t fragmented;
    int nnz;
    bool bNeedsValues;
    
    std::vector<value_type> values;
    size_t maxValues;
    size_t m_numCols;
};


template<typename T>
struct matrix_algebra_type_traits<CRSSparseMatrix<T> >
{
	static const int type = MATRIX_USE_ROW_FUNCTIONS;
};

//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultTransposedAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const CRSSparseMatrix<matrix_t> &A1, const vector_t &w1)
{
	A1.axpy_transposed(dest, alpha1, v1, beta1, w1);
}

///	@}
} // namespace ug

//#include "matrixrow.h"
#include "crssparsematrix_impl.h"
#include "crssparsematrix_print.h"

#endif
