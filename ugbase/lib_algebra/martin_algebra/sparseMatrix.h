#pragma once

#include "math.h"
#include "TemplateExpressions.h"
#include <pmmintrin.h>
#include "blockMatrix.h"

///////////////////////////////////////////////////////////////////

// ForwardExpression (fuer i=0; i<n; i++)
// BackwardExpression
// RandomExpression



/*
template <int n, int nnz>
struct matrix_trait<sparseMatrix<n, nnz> >
{
	typedef blockVector<n> vec_type;
};*/


template<typename vec_type> class Vector;

//#define SPECIALIZE_EXPRESSION_TEMPLATES

///////////////////////////////////////////////////////////////////
//							connection
///////////////////////////////////////////////////////////////////
#include "submatrix.h"

#include "Vector.h"

template<typename mat_type> class matrixrow;
template<typename vec_type> class Vector;
///////////////////////////////////////////////////////////////////
//							SparseMatrix
///////////////////////////////////////////////////////////////////

//!
//! class SparseMatrix
//! templated class, parameter is a blockmatrix type like double or blockDenseMatrix
//! trough matrix_trait<mat_type>::vec_type corresponding blockvector is determined
//! 
template<typename templ_mat_type> 
class SparseMatrix
{
	// functions
public:
	typedef templ_mat_type mat_type;
	typedef typename matrix_trait<mat_type>::vec_type vec_type;
	typedef matrixrow<mat_type> row_type;
	typedef Vector<vec_type> Vector_type;
	typedef matrixrow<mat_type> matrixrow_type;

	struct connection
	{
		int iIndex;		// index to
		mat_type dValue; // smallmatrix value;
		
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
	};
	
public:
	// constructor for empty SparseMatrix
	SparseMatrix();
	// destructor
	~SparseMatrix ();	
	
private:
	// disallowed operations (not defined):
	SparseMatrix(SparseMatrix&); // disallow copy operator
	void operator = (const SparseMatrix &v); // disallow assignment
	
public:	
	//! create
	//! creates an empty SparseMatrix with lenght length_.
	void create(int _rows, int _cols);
	
	//!
	//! createAsTransposeOf
	//! create this as a transpose of SparseMatrix B
	void createAsTransposeOf(const SparseMatrix &B);
	
	//!
	//! eliminateDirichletValues
	//! Dirichlet rows are rows i with only A[i,i] = 1.0 and A[i,j] = 0.0 for all i != j
	//! these rows are added to other rows j with A[j,i] != 0.0, until A[j,i] = 0.0,
	//! and corresponding entries in rhs vector b are changed.
	void eliminateDirichletValues(Vector<vec_type> &b);
	
	void setDirichletRow(int row);
	void setDirichletRows(int *pRows, int nrows);
	
	
	void apply(Vector_type &res, const Vector_type &x) const;
	void applyTransposed(Vector_type &res, const Vector_type &x) const;
	void matmul_minus(Vector_type &res, const Vector_type &x) const;

	// accessor functions for artificial matrixrow-object (= just wrapper with A and row)	
	inline const matrixrow_type getrow(int i) const;
	const matrixrow_type operator [] (int i) const;
	
	
	inline const mat_type getDiag(int i) const;
	inline mat_type &getDiag(int i);	
	
	//!
	//! isUnconnected: true if only A[i,i] != 0.0.
	inline bool isUnconnected(int i) const;	

	//!
	//! add (submatrix)
	//! adds the submatrix mat to A. 
	//! @see submatrix
	void add(const submatrix<mat_type> &mat);
	void set(const submatrix<mat_type> &mat);
	void get(submatrix<mat_type> &mat) const;

	// for other manipulation/accessor functions see matrixrow functions,
	// that is A[i].matrixrowfunction(params).
	
	// accessor functions
public:
	int getLength() const { return rows; }
	int getRows() const { return rows; }
	int getCols() const { return cols; }
	int getTotalNrOfConnections() const { return iTotalNrOfConnections; }
	
	
	// row functions
public:	
	//!
	//! removezeros
	//! remove zero entries of SparseMatrix (experimental)
	void removezeros(int row);
	
	void setMatrixRow(int row, connection *c, int nr);
	void addMatrixRow(int row, connection *c, int nr);	
	int getNrOfConnections(int row) const { return iNrOfConnections[row]; }	

public:
	// print functions
	//------------------
	void print(const char * const name = NULL) const;
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
	
	
	//!
	//! writeToFile
	//! writes Matrix into file filename in ConnectionViewer format.
	void writeToFile(const char *filename) const;
	
	
public:
	// finalizing functions
	//----------------------
	void defrag();
	void finish()
	{
		//defrag();
	}
	
	// diagcomponent
	// ermöglicht mit Template Expressions Dinge wie x = A.Diag() * y (=> x[i] = A.getDiag(i)*y[i])
	// !!! NICHT x = A.Diag() *x ( == GS) !!!
	class diagcomponent
	{	
	public:
		diagcomponent(const SparseMatrix &A_) : A(A_) {}
		mat_type operator [] (int i) const { return A.getDiag(i); }
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
	
	//!
	//! prefetch (EXPERIMENTAL)
	inline void prefetch(int i) const
	{
		prefetchRead(cons[i]);
	}
	
public:
	
	// Iterators
	//---------------------------
	
	// const_RowIterator
	
	class cRowIterator 
	{
	public:
		const SparseMatrix<templ_mat_type> &A;
		int m_position;
		int row;
	public:
		inline cRowIterator(const SparseMatrix<templ_mat_type> &A_, int row_) : A(A_) { row = row_; rewind(); }
		inline cRowIterator(const cRowIterator &other) : A(other.A) { row = other.row; m_position = other.m_position; }
		
		inline const connection &operator *() const {return A.cons[row][m_position];}
		
		inline void operator ++() {	m_position++; }
		
		inline void rewind() { m_position = 0;}
		inline int getPos() const{	return m_position;}
		
		inline bool isEnd() const { return m_position >= A.getNrOfConnections(row); }
	};
	
	cRowIterator beginRow(int row) const
	{
		return cRowIterator(*this, row);
	}
	
	
	// unconst row iterator
	class rowIterator 
	{
	public:
		SparseMatrix<templ_mat_type> &A;
		int m_position;
		int row;
	public:
		inline rowIterator(SparseMatrix<templ_mat_type> &A_, int row_) : A(A_) { row = row_; rewind(); }
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
		matrixrow<mat_type> &row;
		int m_position;
	public:
		inline iterator(matrixrow<mat_type> &r) : row(r) { m_position = 0;}
		inline connection &operator *(){	return row[m_position];	}
		inline void operator ++(){		m_position++;			}
		inline void rewind(){m_position = 0;}
		inline int getPos() const{	return m_position;	}		
		inline const matrixrow &getMatrixRow() const{return row;}
		inline bool isEnd() const { return m_position >= row.getNrOfConnections(); }
	}
 ;*/	
		
	
	class cLowerLeftIterator : public cRowIterator
	{
	public:
		cLowerLeftIterator(const SparseMatrix<templ_mat_type> &A_, int row_) : cRowIterator(A_, row_) { 	rewind();	}
		inline void operator ++()
		{
			cRowIterator::m_position ++;
			while(cRowIterator::m_position < cRowIterator::A.getNrOfConnections(cRowIterator::row) 
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
		cUpperRightIterator(const SparseMatrix<templ_mat_type> &A_, int row_) : cRowIterator(A_, row_) { 	rewind();	}
		inline void operator ++()
		{
			cRowIterator::m_position ++;
			while(cRowIterator::m_position < cRowIterator::A.getNrOfConnections(cRowIterator::row) 
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
	
	
private:	
	void safeSetConnections(int row, connection *mem) const;
	void recreateWithMaxNrOfConnections(int newMax) const;
	
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
	int *iMaxNrOfConnections;			//< max nr of connections for row [i]. TODO.
	connection **cons;			//< pointers to array of connections of each row
	
	int iTotalNrOfConnections;	//!< number of non-zeros
	int bandwidth;				//!< bandwidth (experimental)
	
	connection *consmem;		//!< consecutive memory for the connections
	int consmemsize;					//!< size of the consecutive memory for connections
	int iFragmentedMem;			//!< size of connections memory not in consmem 	
	
	friend class matrixrow<mat_type>;
	
};

#include "matrixRow.h"

#include "Vector.hpp"
#include "SparseMatrix.hpp"
