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
#include "template_expressions.h"
#include "blocks/blocks.h"
#include "blocks/blockVector.h"


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
	typedef MultiIndex<1> index_type;
	typedef FlexLocalMatrix local_matrix_type;
	typedef std::vector<index_type> local_index_type;
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
	bool set_dirichlet_rows(const local_index_type &ind);

	//! calculate res = A x
	template<typename Vector_type>
	bool apply(Vector_type &res, const Vector_type &x) const;

	//! calculate res = A.T x
	template<typename Vector_type>
	bool applyTransposed(Vector_type &res, const Vector_type &x) const;

	// wrapper for applyTransposed
	template<typename Vector_type>
	bool apply_transposed(Vector_type &res, const Vector_type &x) const;

	//! calculate res -= A x
	template<typename Vector_type>
	bool matmul_minus(Vector_type &res, const Vector_type &x) const;

	//! accessor functions for artificial matrixrow-object (= just wrapper with A and row)
	inline const matrixrow_type getrow(size_t i) const;
	const matrixrow_type operator [] (size_t i) const;


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

	bool add(const local_matrix_type &mat, const local_index_type &I, const local_index_type &J);
	bool set(const local_matrix_type &mat, const local_index_type &I, const local_index_type &J);
	bool get(local_matrix_type &mat, const local_index_type &I, const local_index_type &J) const;


	void add(const entry_type &d, size_t row, size_t col);
	void set(const entry_type &d, size_t row, size_t col);
	void get(entry_type &d, size_t row, size_t col) const;

	bool set(double a);


	// for other manipulation/accessor functions see matrixrow functions,
	// that is A[i].matrixrowfunction(params).

public: // accessor functions

	size_t getLength() const { return rows; }
	size_t getRows() const { return rows; }
	size_t getCols() const { return cols; }
	size_t row_size() const { return rows; }
	size_t col_size() const { return rows; }

	size_t getTotalNrOfConnections() const { return iTotalNrOfConnections; }



public:	// row functions
	//! remove zero entries of SparseMatrix (experimental)
	void removezeros(size_t row);

	void setMatrixRow(size_t row, connection *c, size_t nr);
	void addMatrixRow(size_t row, connection *c, size_t nr);
	size_t getNrOfConnections(size_t row) const { return iNrOfConnections[row]; }

public: // output functions

	void print(const char * const name = NULL) const;
	void printToFile(const char *filename) const;
	void printrow(size_t row) const;
	void p() const; // for use in gdb
	void pr(size_t i) const; // for use in gdb

	friend ostream &operator<<(ostream &output, const SparseMatrix &m)
	{
		output << "SparseMatrix " //<< m.name
		<< " [ " << m.rows << " x " << m.cols << ", level " << m.fromlevel << " to " << m.tolevel << " ]";
		return output;
	}
	void printtype() const;

	//! writes Matrix into file filename in ConnectionViewer format.
	void writeToFile(const char *filename) const;


public:
	// finalizing functions
	//----------------------
	void defrag();
	void finalize()
	{
		defrag();
	}

	// diagcomponent
	// ermöglicht mit Template Expressions Dinge wie x = A.Diag() * y (=> x[i] = A.getDiag(i)*y[i])
	// !!! NICHT x = A.Diag() *x ( == GS) !!!
	class diagcomponent
	{
	public:
		diagcomponent(const SparseMatrix &A_) : A(A_) {}
		entry_type operator [] (size_t i) const { return A.getDiag(i); }
		size_t getLength() const { return A.getLength(); }
		friend ostream &operator<<(ostream &output, const diagcomponent &m)
		{ output << "diagonal of " << m.A; return output; }
		const SparseMatrix &A;
	};
	diagcomponent Diag() const
	{
		diagcomponent a(*this);
		return a;
	}

	//! prefetch (EXPERIMENTAL)
	inline void prefetch(size_t i) const
	{
		prefetchRead(cons[i]);
	}

public:

	// Iterators
	//---------------------------

	// const_RowIterator
#ifdef OLDITERATOR
	class cRowIterator
	{
	public:
		//const SparseMatrix<entry_type> &A;
		const SparseMatrix<entry_type> *A;
		size_t m_position;
		size_t row;
	public:
		inline cRowIterator(const SparseMatrix<entry_type> &A_, size_t row_)
		{
			A = &A_;
			row = row_;
			rewind();

		}
		inline cRowIterator(const cRowIterator &other) : A(other.A) { row = other.row; m_position = other.m_position; }

		inline const connection &operator *() const {return A->cons[row][m_position];}

		inline void operator ++() {	m_position++; }

		inline void rewind() { m_position = 0;}
		inline size_t getPos() const{	return m_position;}

		inline bool isEnd() const { return m_position >= A->getNrOfConnections(row); }
	};
#else

	/*class cRowIterator
	{
	public:
		//const SparseMatrix<entry_type> &A;
		const connection * const pEnd;
		const connection * const pStart;
		const connection *p;
		const size_t row;
	public:
		inline cRowIterator(const SparseMatrix<entry_type> &A, size_t row_) :
			pStart(A.cons[row_]),
			pEnd( A.cons[row_]+A.getNrOfConnections(row_) ),
			row(row_)
		{
			rewind();
		}
		inline cRowIterator(const cRowIterator &other) : pStart(other.pStart), pEnd(other.pEnd), row(other.row)
		{
			p = other.p;
		}

		inline const connection &operator *() const {return *p;}

		inline void operator ++() {	p++; }

		inline void rewind() { p = pStart;}
		inline size_t getPos() const{	return p-pStart;}

		inline bool isEnd() const { return p >= pEnd; }
	};*/

class cRowIterator
{
public:
	//const SparseMatrix<entry_type> &A;
	const connection * const pStart;
	const size_t nrOfConnections;
	const size_t row;
	size_t pos;
public:
	inline cRowIterator(const SparseMatrix<entry_type> &A, size_t row_) :	pStart(A.cons[row_]), nrOfConnections( A.getNrOfConnections(row_) ),
		row(row_)
	{
		rewind();
	}
	inline cRowIterator(const cRowIterator &other) : pStart(other.pStart), nrOfConnections(other.nrOfConnections), row(other.row), pos(other.pos)
	{
	}

	inline const connection &operator *() const {return pStart[pos];}

	inline void operator ++() {	pos++; }

	inline void rewind() { pos = 0;}
	inline size_t getPos() const{	return pos;}

	inline bool isEnd() const { return pos >= nrOfConnections; }
};
#endif

	cRowIterator beginRow(size_t row) const
	{
		return cRowIterator(*this, row);
	}


	// unconst row iterator
	class rowIterator
	{
	public:
		SparseMatrix<entry_type> &A;
		size_t m_position;
		size_t row;
	public:
		inline rowIterator(SparseMatrix<entry_type> &A_, size_t row_) : A(A_) { row = row_; rewind(); }
		inline rowIterator(const cRowIterator &other) : A(other.A) { row = other.row; m_position = other.m_position; }

		inline connection &operator *() const {return A.cons[row][m_position];}

		inline void operator ++() {	m_position++; }

		inline void rewind() { m_position = 0;}
		inline size_t getPos() const {	return m_position;}

		inline bool isEnd() const { return m_position >= A.getNrOfConnections(row); }
	};

	rowIterator beginRow(size_t row)
	{
		return rowIterator(*this, row);
	}


	// make this with templates?
/*	class rowIterator
	{
	protected:
		matrixrow<entry_type> &row;
		size_t m_position;
	public:
		inline iterator(matrixrow<entry_type> &r) : row(r) { m_position = 0;}
		inline connection &operator *(){	return row[m_position];	}
		inline void operator ++(){		m_position++;			}
		inline void rewind(){m_position = 0;}
		inline size_t getPos() const{	return m_position;	}
		inline const matrixrow &getMatrixRow() const{return row;}
		inline bool isEnd() const { return m_position >= row.getNrOfConnections(); }
	}
 ;*/

#ifdef OLDITERATOR
	class cLowerLeftIterator : public cRowIterator
	{
	public:
		cLowerLeftIterator(const SparseMatrix<entry_type> &A_, size_t row_) : cRowIterator(A_, row_) { 	rewind();	}
		inline void operator ++()
		{
			cRowIterator::m_position ++;
			while(cRowIterator::m_position < cRowIterator::A->getNrOfConnections(cRowIterator::row)
				  && (  cRowIterator::operator *().iIndex >= cRowIterator::row)) // || cRowIterator::operator *().dValue == 0.0))
				cRowIterator::m_position++;
		}
		inline void rewind()
		{
			cRowIterator::m_position = -1;
			operator ++ ();
		}
	};

	cLowerLeftIterator beginLowerLeftRow(size_t row)  const
	{
		return cLowerLeftIterator(*this, row);
	}

	class cUpperRightIterator : public cRowIterator
	{
	public:
		cUpperRightIterator(const SparseMatrix<entry_type> &A_, size_t row_) : cRowIterator(A_, row_) { 	rewind();	}
		inline void operator ++()
		{
			cRowIterator::m_position ++;
			while(cRowIterator::m_position < cRowIterator::A->getNrOfConnections(cRowIterator::row)
				  && (  cRowIterator::operator *().iIndex <= cRowIterator::row))// || cRowIterator::operator *().dValue == 0.0))
				cRowIterator::m_position++;
		}
		inline void rewind()
		{
			cRowIterator::m_position = -1;
			operator ++ ();
		}
	};


	cUpperRightIterator beginUpperRightRow(size_t row)  const
	{
		return cUpperRightIterator(*this, row);
	}
#else
	class cLowerLeftIterator : public cRowIterator
	{
	public:
		cLowerLeftIterator(const SparseMatrix<entry_type> &A, size_t row_) : cRowIterator(A, row_) { 	rewind();	}
		inline void operator ++()
		{
			cRowIterator::operator ++();
			while(cRowIterator::operator *().iIndex >= cRowIterator::row &&
				  !cRowIterator::isEnd()) // || cRowIterator::operator *().dValue == 0.0))
				cRowIterator::operator ++();
		}
		inline void rewind()
		{
			cRowIterator::rewind();
			while(cRowIterator::operator *().iIndex >= cRowIterator::row &&
				  !cRowIterator::isEnd()) // || cRowIterator::operator *().dValue == 0.0))
				cRowIterator::operator ++();
		}
	};

	cLowerLeftIterator beginLowerLeftRow(size_t row)  const
	{
		return cLowerLeftIterator(*this, row);
	}

	class cUpperRightIterator : public cRowIterator
		{
		public:
			cUpperRightIterator(const SparseMatrix<entry_type> &A, size_t row_) : cRowIterator(A, row_) { 	rewind();	}
			inline void operator ++()
			{
				cRowIterator::operator ++();
				while(cRowIterator::operator *().iIndex <= cRowIterator::row && !cRowIterator::isEnd()) // || cRowIterator::operator *().dValue == 0.0))
					cRowIterator::operator ++();
			}
			inline void rewind()
			{
				cRowIterator::rewind();
				while(cRowIterator::operator *().iIndex <= cRowIterator::row && !cRowIterator::isEnd()) // || cRowIterator::operator *().dValue == 0.0))
					cRowIterator::operator ++();
					cRowIterator::operator ++();
			}
		};


	cUpperRightIterator beginUpperRightRow(size_t row)  const
	{
		return cUpperRightIterator(*this, row);
	}
#endif


private:
	//! "safe" way to set a connection, since when cons[row] is in the big consecutive consmem-array,
	void safeSetConnections(size_t row, connection *mem) const;
	void recreateWithMaxNrOfConnections(size_t newMax) const;

public:
	void setEstimatedRowSize(size_t estimatedRowSize_)
	{
		estimatedRowSize = estimatedRowSize_;
	}
//     data
//----------------

public:
	size_t fromlevel;						//!< fromlevel
	size_t tolevel;						//!< tolevel
	const char *name;					//!< name of the SparseMatrix for debuging / printing.

private:

	size_t rows;							//!< nr of rows
	size_t cols;							//!< nr of cols
	size_t *iNrOfConnections;				//< array nr of Connections for each row
	connection **cons;			//< pointers to array of connections of each row

	size_t iTotalNrOfConnections;	//!< number of non-zeros
	size_t bandwidth;				//!< bandwidth (experimental)

	size_t estimatedRowSize;				//!< estimated length of each row
	size_t *iMaxNrOfConnections;			//!< max nr of connections for row [i]. TODO.


	connection *consmem;		//!< consecutive memory for the connections
	size_t consmemsize;					//!< size of the consecutive memory for connections
	size_t iFragmentedMem;			//!< size of connections memory not in consmem

	friend class matrixrow<entry_type>;

};

} // namespace ug

#include "matrixrow.h"
#include "sparsematrix_impl.h"
