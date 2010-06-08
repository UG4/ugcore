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
//#include "blockMatrix.h"


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
		int iIndex;		// index to
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
	

	bool create(int _rows, int _cols);
	bool destroy();
	
	//! create this as a transpose of SparseMatrix B
	void createAsTransposeOf(const SparseMatrix &B);

	
private: // disallowed operations (not defined):
	SparseMatrix(SparseMatrix&); ///< disallow copy operator
	void operator = (const SparseMatrix &v); ///< disallow assignment
	
public:	// general functions	
	template<typename Vector_type>
	void eliminateDirichletValues(Vector_type &b);
	
	void setDirichletRow(int row);
	void setDirichletRows(int *pRows, int nrows);
	bool set_dirichlet_rows(const local_index_type &ind);
	
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
	inline const matrixrow_type getrow(int i) const;
	const matrixrow_type operator [] (int i) const;
	
	
	//! get Diagonal A_[i,i] of matrix
	inline const entry_type &getDiag(int i) const;
	inline entry_type &getDiag(int i);
	
	//! isUnconnected: true if only A[i,i] != 0.0.
	inline bool isUnconnected(int i) const;

	//! adds the submatrix mat to A. 	
	
	template<typename M>
	void add(const M &mat, int *rows, int *cols);
	template<typename M>
	void set(const M &mat, int *rows, int *cols);
	template<typename M>
	void get(M &mat, int *rows, int *cols) const;
	
	template<typename M>
	void add(const M &mat, vector<int> &rows, vector<int> &cols);
	template<typename M>
	void set(const M &mat, vector<int> &rows, vector<int> &cols);
	template<typename M>
	void get(M &mat, vector<int> &rows, vector<int> &cols) const;
	
	bool add(const local_matrix_type &mat, const local_index_type &I, const local_index_type &J);
	bool set(const local_matrix_type &mat, const local_index_type &I, const local_index_type &J);
	bool get(local_matrix_type &mat, const local_index_type &I, const local_index_type &J) const;
	
	
	void add(const entry_type &d, int row, int col);
	void set(const entry_type &d, int row, int col);
	void get(entry_type &d, int row, int col) const;
	
	bool set(double a);


	// for other manipulation/accessor functions see matrixrow functions,
	// that is A[i].matrixrowfunction(params).
		
public: // accessor functions
	
	int getLength() const { return rows; }
	int getRows() const { return rows; }
	int getCols() const { return cols; }
	int row_size() const { return rows; }
	int col_size() const { return rows; }
	
	int getTotalNrOfConnections() const { return iTotalNrOfConnections; }
	
	
	
public:	// row functions
	//! remove zero entries of SparseMatrix (experimental)
	void removezeros(int row);
	
	void setMatrixRow(int row, connection *c, int nr);
	void addMatrixRow(int row, connection *c, int nr);
	int getNrOfConnections(int row) const { return iNrOfConnections[row]; }

public: // output functions

	void print(const char * const name = NULL) const;
	void printToFile(const char *filename) const;
	void printrow(int row) const;
	void p() const; // for use in gdb
	void pr(int i) const; // for use in gdb

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
		entry_type operator [] (int i) const { return A.getDiag(i); }
		int getLength() const { return A.getLength(); }
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
	inline void prefetch(int i) const
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
		int m_position;
		int row;
	public:
		inline cRowIterator(const SparseMatrix<entry_type> &A_, int row_)
		{
			A = &A_;
			row = row_; 
			rewind(); 
			
		}
		inline cRowIterator(const cRowIterator &other) : A(other.A) { row = other.row; m_position = other.m_position; }
		
		inline const connection &operator *() const {return A->cons[row][m_position];}
		
		inline void operator ++() {	m_position++; }
		
		inline void rewind() { m_position = 0;}
		inline int getPos() const{	return m_position;}
		
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
		const int row;
	public:
		inline cRowIterator(const SparseMatrix<entry_type> &A, int row_) :
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
		inline int getPos() const{	return p-pStart;}
		
		inline bool isEnd() const { return p >= pEnd; }
	};*/
	
class cRowIterator 
{
public:
	//const SparseMatrix<entry_type> &A;
	const connection * const pStart;
	const int nrOfConnections;
	const int row;
	int pos;
public:
	inline cRowIterator(const SparseMatrix<entry_type> &A, int row_) :	pStart(A.cons[row_]), nrOfConnections( A.getNrOfConnections(row_) ),
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
	inline int getPos() const{	return pos;}
	
	inline bool isEnd() const { return pos >= nrOfConnections; }
};
#endif
	
	cRowIterator beginRow(int row) const
	{
		return cRowIterator(*this, row);
	}
	
	
	// unconst row iterator
	class rowIterator 
	{
	public:
		SparseMatrix<entry_type> &A;
		int m_position;
		int row;
	public:
		inline rowIterator(SparseMatrix<entry_type> &A_, int row_) : A(A_) { row = row_; rewind(); }
		inline rowIterator(const cRowIterator &other) : A(other.A) { row = other.row; m_position = other.m_position; }
		
		inline connection &operator *() const {return A.cons[row][m_position];}
		
		inline void operator ++() {	m_position++; }
		
		inline void rewind() { m_position = 0;}
		inline int getPos() const {	return m_position;}
		
		inline bool isEnd() const { return m_position >= A.getNrOfConnections(row); }
	};
	
	rowIterator beginRow(int row)
	{
		return rowIterator(*this, row);
	}
	
	
	// make this with templates?
/*	class rowIterator 
	{
	protected:
		matrixrow<entry_type> &row;
		int m_position;
	public:
		inline iterator(matrixrow<entry_type> &r) : row(r) { m_position = 0;}
		inline connection &operator *(){	return row[m_position];	}
		inline void operator ++(){		m_position++;			}
		inline void rewind(){m_position = 0;}
		inline int getPos() const{	return m_position;	}
		inline const matrixrow &getMatrixRow() const{return row;}
		inline bool isEnd() const { return m_position >= row.getNrOfConnections(); }
	}
 ;*/	
		
#ifdef OLDITERATOR
	class cLowerLeftIterator : public cRowIterator
	{
	public:
		cLowerLeftIterator(const SparseMatrix<entry_type> &A_, int row_) : cRowIterator(A_, row_) { 	rewind();	}
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
				 
	cLowerLeftIterator beginLowerLeftRow(int row)  const
	{
		return cLowerLeftIterator(*this, row);
	}
	
	class cUpperRightIterator : public cRowIterator
	{
	public:
		cUpperRightIterator(const SparseMatrix<entry_type> &A_, int row_) : cRowIterator(A_, row_) { 	rewind();	}
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
				
				
	cUpperRightIterator beginUpperRightRow(int row)  const
	{
		return cUpperRightIterator(*this, row);
	}
#else
	class cLowerLeftIterator : public cRowIterator
	{
	public:
		cLowerLeftIterator(const SparseMatrix<entry_type> &A, int row_) : cRowIterator(A, row_) { 	rewind();	}
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
	
	cLowerLeftIterator beginLowerLeftRow(int row)  const
	{
		return cLowerLeftIterator(*this, row);
	}
	
	class cUpperRightIterator : public cRowIterator
		{
		public:
			cUpperRightIterator(const SparseMatrix<entry_type> &A, int row_) : cRowIterator(A, row_) { 	rewind();	}
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
	
	
	cUpperRightIterator beginUpperRightRow(int row)  const
	{
		return cUpperRightIterator(*this, row);
	}
#endif
	
	
private:	
	//! "safe" way to set a connection, since when cons[row] is in the big consecutive consmem-array,
	void safeSetConnections(int row, connection *mem) const;
	void recreateWithMaxNrOfConnections(int newMax) const;
	
public:
	void setEstimatedRowSize(int estimatedRowSize_)
	{
		estimatedRowSize = estimatedRowSize_;
	}
//     data
//----------------
	
public:
	int fromlevel;						//!< fromlevel
	int tolevel;						//!< tolevel
	const char *name;					//!< name of the SparseMatrix for debuging / printing.
	
private:	

	int rows;							//!< nr of rows
	int cols;							//!< nr of cols
	int *iNrOfConnections;				//< array nr of Connections for each row
	connection **cons;			//< pointers to array of connections of each row
	
	int iTotalNrOfConnections;	//!< number of non-zeros
	int bandwidth;				//!< bandwidth (experimental)
	
	int estimatedRowSize;				//!< estimated length of each row
	int *iMaxNrOfConnections;			//!< max nr of connections for row [i]. TODO.

	
	connection *consmem;		//!< consecutive memory for the connections
	int consmemsize;					//!< size of the consecutive memory for connections
	int iFragmentedMem;			//!< size of connections memory not in consmem
	
	friend class matrixrow<entry_type>;
	
};

} // namespace ug

#include "matrixrow.h"
#include "sparsematrix_impl.h"
