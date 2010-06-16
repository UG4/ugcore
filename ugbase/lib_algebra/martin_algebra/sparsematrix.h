/*
 *  SparseMatrix.h
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */


#pragma once

#include "math.h"

#include <pmmintrin.h>

#ifdef FLEXAMG
#include "blocks.h"
#include "blockVector.h"
#include "TemplateExpressions.h"

#else
#include "template_expressions.h"
#include "blocks/blocks.h"
#include "blocks/blockVector.h"
#endif

//#define SPECIALIZE_EXPRESSION_TEMPLATES

///////////////////////////////////////////////////////////////////
//							connection
///////////////////////////////////////////////////////////////////
//#include "submatrix.h"

#include "Vector.h"

namespace ug{

template<typename entry_type> class matrixrow;
template<typename vec_type> class Vector;
///////////////////////////////////////////////////////////////////
//							SparseMatrix
///////////////////////////////////////////////////////////////////

//!
//! SparseMatrix
//! template parameter T: blocktype
//! SparseMatrix for big, variable sparse matrices.
//! matrix is stored independent row-wise
//! \sa matrixRow
template<typename T> 
class SparseMatrix : public TE_MAT<SparseMatrix<T> >
{
public:
#ifndef FLEXAMG
	typedef MultiIndex<1> index_type;
	typedef FlexLocalMatrix local_matrix_type;
	typedef std::vector<index_type> local_index_type;
#endif
	// functions
public:
	typedef T entry_type;
	typedef matrixrow<entry_type> row_type;
	typedef matrixrow<entry_type> matrixrow_type;
	//typedef submatrix<entry_type> submatrix_type;

	struct connection
	{
		size_t iIndex;		// index to
		entry_type dValue; // smallmatrix value;
		
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
	
public: // construction etc
	
	// constructor for empty SparseMatrix
	SparseMatrix();
	// destructor
	~SparseMatrix ();	
	

	bool create(size_t _rows, size_t _cols);
	bool destroy();
	
	//! create this as a transpose of SparseMatrix B
	void createAsTransposeOf(const SparseMatrix &B);
	

private: // disallowed operations (not defined):
	SparseMatrix(SparseMatrix&); ///< disallow copy operator
	void operator = (const SparseMatrix &v); ///< disallow assignment
	
public:	// general functions	
	template<typename Vector_type>
	void eliminateDirichletValues(Vector_type &b);
	
	void setDirichletRow(size_t row);
	void setDirichletRows(size_t *pRows, size_t nrows);
#ifndef FLEXAMG
	bool set_dirichlet_rows(const local_index_type &ind);
#endif
	

	//! calculate res = A x
	template<typename Vector_type>
	bool apply(Vector_type &res, const Vector_type &x) const;
	//! calculate res = A.T x
	template<typename Vector_type>
	bool applyTransposed(Vector_type &res, const Vector_type &x) const;
	//! calculate res -= A x
	template<typename Vector_type>
	bool matmul_minus(Vector_type &res, const Vector_type &x) const;


	//! accessor functions for artificial matrixrow-object (= just wrapper with A and row)	
	inline const matrixrow_type getrow(size_t i) const;
	inline const matrixrow_type operator [] (size_t i) const;
	
	
	//! get Diagonal A_[i,i] of matrix
	inline const entry_type &getDiag(size_t i) const;
	inline entry_type &getDiag(size_t i);
	
	//! isUnconnected: true if only A[i,i] != 0.0.
	inline bool isUnconnected(size_t i) const;

	//! adds the submatrix mat to A. 	
	
	template<typename M>
	void add(const M &mat, size_t *rows, size_t *cols);
	template<typename M>
	void set(const M &mat, size_t *rows, size_t *cols);
	template<typename M>
	void get(M &mat, size_t *rows, size_t *cols) const;
	
	template<typename M>
	void add(const M &mat, vector<size_t> &rows, vector<size_t> &cols);
	template<typename M>
	void set(const M &mat, vector<size_t> &rows, vector<size_t> &cols);
	template<typename M>
	void get(M &mat, vector<size_t> &rows, vector<size_t> &cols) const;

#ifndef FLEXAMG
	bool add(const local_matrix_type &mat, const local_index_type &I, const local_index_type &J);
	bool set(const local_matrix_type &mat, const local_index_type &I, const local_index_type &J);
	bool get(local_matrix_type &mat, const local_index_type &I, const local_index_type &J) const;
#endif
	
	
	void add(const entry_type &d, size_t row, size_t col);
	void set(const entry_type &d, size_t row, size_t col);
	void get(entry_type &d, size_t row, size_t col) const;
	
	bool set(double a);

	const entry_type &operator() (int r, int c) const;
	entry_type &operator() (int r, int c);

	// for other manipulation/accessor functions see matrixrow functions,
	// that is A[i].matrixrowfunction(params).

public: // accessor functions
	size_t row_size() const { return rows; } // deprecated
	size_t col_size() const { return cols; } // deprecated
	
	size_t num_rows() const { return rows; }
	size_t num_cols() const { return cols; }

	size_t getTotalNrOfConnections() const { return iTotalNrOfConnections; }
	
	
	
public:	// row functions
	//! remove zero entries of SparseMatrix (experimental)
	void removezeros(size_t row);
	
	void setMatrixRow(size_t row, connection *c, size_t nr);
	void addMatrixRow(size_t row, connection *c, size_t nr);
	inline size_t getNrOfConnections(size_t row) const;

public: // output functions

	void print(const char * const name = NULL) const;
	void printToFile(const char *filename) const;
	void printrow(size_t row) const;
	void p() const; // for use in gdb
	void pr(size_t i) const; // for use in gdb

	friend ostream &operator<<(ostream &output, const SparseMatrix &m)
	{
		output << "SparseMatrix " //<< m.name 
		<< " [ " << m.rows << " x " << m.cols << " ]";
		return output;
	}
	void printtype() const; 
	
	//! writes Matrix into file filename in ConnectionViewer format.
	void writeToFile(const char *filename) const;
	
	
public:
	// finalizing functions
	//----------------------
	void defrag();
	void definalize();
	void finalize()
	{
		defrag();
	}


public:
	
	// Iterators
	//---------------------------
	
	// const_RowIterator

	class cRowIterator 
	{
	public:
		//const SparseMatrix<entry_type> &A;
		const connection * const pStart;
		const connection * const pEnd;
		const connection * p;
	public:
		inline cRowIterator(const SparseMatrix<entry_type> &A, size_t row) : pStart(A.pRowStart[row]), pEnd(A.pRowEnd[row]) { rewind();	}
		inline cRowIterator(const cRowIterator &other) : pStart(other.pStart), pEnd(other.pEnd), p(other.p)	{ }

		inline const connection &operator *() const {return *p;}

		inline void operator ++() {	++p; }

		inline void rewind() { p = pStart;}

		inline bool isEnd() const { return p >= pEnd; }
	};
	
	// unconst row iterator
	class rowIterator 
	{
	public:
		connection * const pStart;
		connection * const pEnd;
		connection * p;
	public:
		inline rowIterator(SparseMatrix<entry_type> &A, size_t row_) : pStart(A.pRowStart[row_]), pEnd(A.pRowEnd[row_]) { rewind(); }
		inline rowIterator(const cRowIterator &other) : pStart(other.pStart), pEnd(other.pEnd), p(other.p) { }

		inline connection &operator *() const {return *p;}

		inline void operator ++() {	++p; }

		inline void rewind() { p = pStart;}

		inline bool isEnd() const { return p >= pEnd; }
	};

	class cLowerLeftIterator : public cRowIterator
	{
	private:
		int row;
	public:
		cLowerLeftIterator(const SparseMatrix<entry_type> &A, size_t row_) : cRowIterator(A, row_), row(row_) { cRowIterator::rewind();	}

		inline bool isEnd() const { return this->p >= this->pEnd || this->p->iIndex >= row; }
	};

	class cUpperRightIterator : public cRowIterator
	{
	private:
		int row;
	public:
		cUpperRightIterator(const SparseMatrix<entry_type> &A, size_t row_) : cRowIterator(A, row_), row(row_) { 	rewind();	}
		inline void rewind()
		{
			int left = 0, right = this->pEnd-this->pStart;
			int mid;
			while(left < right)
			{
				mid = (left+right)/2;
				if(this->pStart[mid].iIndex <= row)
					left = mid+1;
				else
					right = mid-1;
			}
			this->p = this->pStart+mid;
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
	
	cUpperRightIterator beginUpperRightRow(size_t row)  const
	{
		return cUpperRightIterator(*this, row);
	}
	
	
private:	
	//! "safe" way to set a connection, since when cons[row] is in the big consecutive consmem-array,
	void safeSetConnections(size_t row, connection *mem) const;

	bool isFinalized() const;
	bool isInConsMem(size_t row) const;
	size_t getConnection(size_t r, size_t c) const;
public:
	void setEstimatedRowSize(size_t estimatedRowSize_)
	{
		estimatedRowSize = estimatedRowSize_;
	}
	//     data
	//----------------
	
public:
#ifdef FLEXAMG
	int tolevel, fromlevel;
#endif
	const char *name;					//!< name of the SparseMatrix for debuging / printing.
	
private:	

	size_t rows;						//!< nr of rows
	size_t cols;						//!< nr of cols
	connection **pRowStart;				//< pointers to array of connections of each row
	connection **pRowEnd;				//< pointers to array of connections of each row
	
	size_t iTotalNrOfConnections;		//!< number of non-zeros
	size_t bandwidth;					//!< bandwidth (experimental)
	
	size_t estimatedRowSize;			//!< estimated length of each row
	size_t *iMaxNrOfConnections;		//!< max nr of connections for row [i]. TODO.
	

	connection *consmem;				//!< consecutive memory for the connections
	size_t consmemsize;					//!< size of the consecutive memory for connections
	size_t iFragmentedMem;				//!< size of connections memory not in consmem
	
	friend class matrixrow<entry_type>;
	
};

} // namespace ug

#include "matrixRow.h"
#include "sparseMatrix_impl.h"
