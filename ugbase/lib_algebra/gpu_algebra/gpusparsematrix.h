/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__CPU_ALGEBRA__GPUSparseMatrix__
#define __H__UG__CPU_ALGEBRA__GPUSparseMatrix__



#include "math.h"
#include "common/common.h"
#include "../algebra_common/sparsematrix_util.h"
#include <iostream>
#include <algorithm>
#include "common/util/ostream_util.h"

#include "../algebra_common/connection.h"
#include "../algebra_common/matrixrow.h"
#include "../common/operations_mat/operations_mat.h"

#include "cuda/cuda_manager.h"
#include "common/debug_print.h"

#define PROFILE_GPUMATRIX(name) PROFILE_BEGIN_GROUP(name, "GPUSparseMatrix algebra")


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


/** GPUSparseMatrix
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
template<typename TValueType>
class GPUSparseMatrix
{
public:
	using value_type = TValueType;
	enum {rows_sorted=true};

	using this_type = GPUSparseMatrix<value_type>;

public:
	using connection = AlgebraicConnection<TValueType>;
	using row_type = MatrixRow<this_type>;
	using const_row_type = ConstMatrixRow<this_type>;

public:
	// construction etc
	//----------------------

	/// constructor for empty GPUSparseMatrix
	GPUSparseMatrix();
	/// destructor
	virtual ~GPUSparseMatrix ()
	{
		freeGPU();
	}


	/**
	 * \brief resizes the GPUSparseMatrix
	 * \param newRows new nr of rows
	 * \param newCols new nr of cols
	 * \return
	 */
	bool resize_and_clear(size_t newRows, size_t newCols);
	bool resize_and_keep_values(size_t newRows, size_t newCols);

	/**
	 * \brief write in a empty GPUSparseMatrix (this) the transpose GPUSparseMatrix of B.
	 * \param B			the matrix of which to create the transpose of
	 * \param scale		an optional scaling
	 * \return			true on success
	 */
	bool set_as_transpose_of(const GPUSparseMatrix<value_type> &B, double scale=1.0);

	/**
	 * \brief create/recreate this as a copy of GPUSparseMatrix B
	 * \param B			the matrix of which to create a copy of
	 * \param scale		an optional scaling
	 * \return			true on success
	 */
	bool set_as_copy_of(const GPUSparseMatrix<value_type> &B, double scale=1.0);
	GPUSparseMatrix<value_type> &operator = (const GPUSparseMatrix<value_type> &B)
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

	//! calculate x = alpha*x + beta*A*y (A = this matrix)
	template<typename vector_t>
	void axpy(double alpha, vector_t &x, double beta, const vector_t &y) const;


	//! calculate dest = alpha1*v1 + beta1*A^T*w1 (A = this matrix)
	template<typename vector_t>
	bool axpy_transposed(vector_t &dest,
			const number &alpha1, const vector_t &v1,
			const number &beta1, const vector_t &w1) const;

	//! calculated dest = beta1*A*w1 . For empty rows, dest will not be changed
	template<typename vector_t>
	void apply_ignore_zero_rows(vector_t &dest,
			const number &beta1, const vector_t &w1) const { assert(0); }

	//! calculated dest = beta1*A*w1 . For empty cols of A (=empty rows of A^T), dest will not be changed
	template<typename vector_t>
	void apply_transposed_ignore_zero_rows(vector_t &dest,
			const number &beta1, const vector_t &w1) const { assert(0); }



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
	GPUSparseMatrix<value_type> &operator *= (double d) { scale(d); return *this; }

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

	//! set matrix to Id*a
	bool set(double a);

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


	//using const_row_iterator = const connection * ;
	//using const_row_iterator = connection * ;
	/** const_row_iterator
	 * const iterator over a row
	 */


	void add_iterator() const
	{
		iIterators++;
	}
	void remove_iterator() const
	{
		iIterators--;
	}
	// a row_iterator has to suppport
	// operator ++, operator +=, index() const, const value_type &value() const, value_type &value()
	// a const_row_iterator has to suppport
	// operator ++, operator +=, index() const, const value_type &value() const

	inline void check_row(size_t row, int i) const
	{
		UG_ASSERT(i < rowEnd[row] && i >= rowStart[row], "row iterator row " << row << " pos " << i << " out of bounds [" << rowStart[row] << ", " << rowEnd[row] << "]");
	}

	/**
	 *  row_iterator
	 *  iterator over a row
	 */
	class row_iterator
    {
        GPUSparseMatrix &A;
        size_t row;
        size_t i;
    public:
        inline void check() const {A.check_row(row, i); }
        row_iterator(GPUSparseMatrix &_A, size_t _row, size_t _i) : A(_A), row(_row), i(_i) { A.add_iterator(); }
        ~row_iterator() { A.remove_iterator(); }
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
        const GPUSparseMatrix &A;
        size_t row;
        size_t i;
    public:
        inline void check() const {A.check_row(row, i); }
        const_row_iterator(const GPUSparseMatrix &_A, size_t _row, size_t _i) : A(_A), row(_row), i(_i) {A.add_iterator();}
        ~const_row_iterator() { A.remove_iterator(); }
        const_row_iterator *operator ->() { return this; }
        const value_type &value() const { check(); return A.values[i];   }
        size_t index() const { check(); return A.cols[i];     }
        bool operator!=(const const_row_iterator &o) const { return i != o.i; }
        void operator++() { ++i; }
        void operator+=(int nr) { i+=nr; }
		bool operator== (const const_row_iterator &other) const { return other.i == i; }
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

public:
	// output functions
	//----------------------

	void print(const char * const name = nullptr) const;
	void printtype() const;

	void print_to_file(const char *filename) const;
	void printrow(size_t row) const;

	friend std::ostream &operator<<(std::ostream &out, const GPUSparseMatrix &m)
	{
		out << "GPUSparseMatrix " //<< m.name
		<< " [ " << m.num_rows() << " x " << m.num_cols() << " ]";
		return out;
	}


	void p() const { print(); } // for use in gdb
	void pr(size_t row) const {printrow(row); } // for use in gdb

private:
	// disallowed operations (not defined):
	//---------------------------------------
	GPUSparseMatrix(GPUSparseMatrix&); ///< disallow copy operator


    void assureValuesSize(size_t s);
    size_t get_nnz() const { return nnz; }

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







    ///////////////////////
public:
   void initGPU()
   {
	   d_cols = d_rowStart = nullptr;
	   d_values = nullptr;
	   descr = 0;
	   bOnDevice = false;
   }
   void freeGPU()
   {
	   cudaFree(d_cols);
	   cudaFree(d_rowStart);
	   cudaFree(d_values);
   }

   const int *get_device_cols() const { check_device(); return d_cols; }
   const int *get_device_rowStart() const { check_device(); return d_rowStart; }
   const double *get_device_value_ptr() const { check_device(); return d_values; }
   cusparseMatDescr_t get_matrix_descr() const { check_device(); return descr; }

   void copy_to_device()
   {
	   CUDAManager::get_instance();
	   descr = 0;
	   cusparseStatus_t cusparseStatus = cusparseCreateMatDescr(&descr);

	   if (checkCudaErrors(cusparseStatus))
	   {
		   exit(EXIT_FAILURE);
	   }

	   cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
	   cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

	   defragment();
	   d_values = &values[0];

	   assert(cols.size() == values.size() && cols.size() == nnz);

	   UG_LOG("cols.size = " << cols.size()*sizeof(int) << " values.size() == " << values.size()*sizeof(value_type) << " rowStart.size = " << rowStart.size()*sizeof(int) << "\n");

	   UG_LOG("gpusparsematrix.h:"<<__LINE__ << "\n")
	   d_cols = CudaCreateAndCopyToDevice(cols);
	   UG_LOG("gpusparsematrix.h:"<<__LINE__ << "\n")
	   d_rowStart = CudaCreateAndCopyToDevice(rowStart);
	   UG_LOG("gpusparsematrix.h:"<<__LINE__ << "\n")
//	   PrintVector(values, "GPUVector::values");
	   d_values = CudaCreateAndCopyToDevice(values);
   }


   void check_device() const
   {
	   if(bOnDevice==true) return;
	   GPUSparseMatrix<value_type>* c = const_cast<GPUSparseMatrix<value_type>*>(this);
	   c->bOnDevice=true;
	   c->copy_to_device();
   }
private:
   //using CRSSparseMatrix::nnz;

   int *d_cols, *d_rowStart;
   double *d_values;
   cusparseMatDescr_t descr;

   bool bOnDevice;
};


template<typename T>
struct matrix_algebra_type_traits<GPUSparseMatrix<T> >
{
	enum{
		type=MATRIX_USE_ROW_FUNCTIONS
	};
};

//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultTransposedAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const GPUSparseMatrix<matrix_t> &A1, const vector_t &w1)
{
	A1.axpy_transposed(dest, alpha1, v1, beta1, w1);
}

// end group cpu_algebra
/// \}

} // namespace ug

//#include "matrixrow.h"
#include "gpusparsematrix_impl.h"
#include "gpusparsematrix_print.h"

#endif
