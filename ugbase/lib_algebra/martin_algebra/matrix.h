#pragma once

#include "math.h"
#include "TemplateExpressions.h"
#include <pmmintrin.h>
#include "smallMatrix.h"

///////////////////////////////////////////////////////////////////

// ForwardExpression (fuer i=0; i<n; i++)
// BackwardExpression
// RandomExpression



/*
template <int n, int nnz>
struct matrix_trait<sparseMatrix<n, nnz> >
{
	typedef fixedVector<n> vec_type;
};*/


template<typename vec_type> class Vector;

//#define SPECIALIZE_EXPRESSION_TEMPLATES

///////////////////////////////////////////////////////////////////
//							connection
///////////////////////////////////////////////////////////////////

template <typename mat_type>
class submatrix
{
public:
	submatrix(int *indices_, int nr_)
	{
		nr = nr_;
		values = new mat_type[nr*nr];
		
		for(int a=0; a<nr; a++) 
			for(int b=0; b<nr; b++)
				values[a*nr + b] = 0.0;
		
		indices = new int[nr];
		memcpy(indices, indices_, sizeof(int)*nr);
	}
	
	submatrix(int *indices_, int *unknowns, int nr_)
	{
		nr = nr_;
		values = new mat_type[nr*nr];
		
		for(int a=0; a<nr; a++) 
		{
			for(int b=0; b<nr; b++)
			{
				// could be flipped
				setSize(values[a*nr + b], unknowns[a], unknowns[b]);
				values[a*nr + b] = 0.0;
			}			
		}
		
		indices = new int[nr];
		memcpy(indices, indices_, sizeof(int)*nr);
	}
	~submatrix()
	{
		delete[] values;
		delete[] indices;
	}
	
	mat_type &operator () (int from, int to)
	{
		ASSERT(from < nr && to < nr && to >= 0 && from >= 0);
		return values[from*nr + to];
	}
	
	int getNr() { return nr; }
	int getIndex(int localIndex) { return indices[localIndex]; }
	int getLocalFromIndex(int index)
	{
		for(int i=0; i<nr; i++)
			if(indices[i] == index) 
				return i;
		
		return -1;
	}
	
	
private:	
	int *indices;
	int nr;
	mat_type *values;
};


#include "Vector.h"

template<typename mat_type> class matrixrow;
template<typename vec_type> class Vector;
///////////////////////////////////////////////////////////////////
//							matrix
///////////////////////////////////////////////////////////////////

//!
//! class matrix
//! templated class, parameter is a blockmatrix type like double or fixedMatrix
//! trough matrix_trait<mat_type>::vec_type corresponding blockvector is determined
//! 
//! this is a little bit stupid with mutable and const bcs of the matrixrow class...
template<typename templ_mat_type> 
class matrix
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
	
	// constructor for empty matrix
	matrix(const char *_name = "");
	// constructor for matrix with length and name
	matrix (int _length, const char *name = "");
	// destructor
	~matrix ();	
	
	// disallowed operations (not defined):
	matrix(matrix&); // disallow copy operator
	void operator = (const matrix &v); // disallow assignment
	
	//!
	//! createAsTransposeOf
	//! create this as a transpose of matrix B (needs length of this, since B has only nr of rows)
	void createAsTransposeOf(const matrix &B, int length);
	
	//! create
	//! creates an empty matrix with lenght length_.
	void create(int _length);
	
	//!
	//! eliminateDirichletValues
	//! Dirichlet rows are rows i with only A[i,i] = 1.0 and A[i,j] = 0.0 for all i != j
	//! these rows are added to other rows j with A[j,i] != 0.0, until A[j,i] = 0.0,
	//! and corresponding entries in rhs vector b are changed.
	void eliminateDirichletValues(Vector<vec_type> &b);

	// accessor functions for rows and elements	
	inline matrixrow_type getrow(int i);
	inline const matrixrow_type getrow(int i) const;
	
	matrixrow_type operator [] (int i) ;	
	const matrixrow_type operator [] (int i) const;
	
	inline const mat_type getDiag(int i) const;
	inline mat_type &getDiag(int i);	
	
	//!
	//! isUnconnected: true if only A[i,i] != 0.0.
	inline bool isUnconnected(int i) const;
	

	//!
	//! addSubmatrix
	//! adds the submatrix mat to A. 
	//! @see submatrix
	void addSubmatrix(submatrix<mat_type> &mat);

	// for other manipulation/accessor functions see matrixrow functions,
	// that is A[i].matrixrowfunction(params).

	// print functions
	//------------------
	void print(const char * const name = NULL) const;
	void p(); // gdb
	friend ostream &operator<<(ostream &output, const matrix &m)
	{
		output << "matrix " << m.name << "[" << m.length << "]";
		return output;
	}
	void printtype() const; 

	// finalizing functions
	void defrag();
	void finish()
	{
		//defrag();
	}
	
	//!
	//! writeToFile
	//! writes Matrix into file filename in ConnectionViewer format.
	void writeToFile(const char *filename) const;
	
	// diagcomponent
	// ermöglicht mit Template Expressions Dinge wie x = A.Diag() * y (=> x[i] = A.getDiag(i)*y[i])
	// !!! NICHT x = A.Diag() *x ( == GS) !!!
	class diagcomponent
	{	
	public:
		diagcomponent(const matrix &A_) : A(A_) {}
		mat_type operator [] (int i) const { return A.getDiag(i); }
		int getLength() const { return A.getLength(); }
		friend ostream &operator<<(ostream &output, const diagcomponent &m)
		{ output << "diagonal of " << m.A; return output; }
		const matrix &A;
	};
	diagcomponent Diag() const
	{
		diagcomponent a(*this);
		return a;
	}
	
	//!
	//! prefetch
	inline void prefetch(int i) const
	{
		prefetchRead(cons[i]);
	}
	
	// accessor functions
	int getLength() const { return length; }
	int getTotalNrOfConnections() const { return iTotalNrOfConnections; }
	

private:	
	void saveSetConnections(int row, connection *mem) const;
	void recreateWithMaxNrOfConnections(int newMax) const;
	
//     data
//----------------

public:
	int fromlevel;						//!< fromlevel
	int tolevel;						//!< tolevel
	const char *name;					//!< name of the matrix for debuging / printing.

private:	

	int length;							//!< nr of rows
	int *iNrOfConnections;				//< array nr of Connections for each row
	mutable connection **cons;			//< pointers to array of connections of each row

	mutable int iTotalNrOfConnections;	//!< number of non-zeros
	mutable int bandwidth;				//!< bandwidth (experimental)
	
	mutable connection *consmem;		//!< consecutive memory for the connections
	int consmemsize;					//!< size of the consecutive memory for connections
	mutable int iFragmentedMem;			//!< size of connections memory not in consmem 
	
	
	friend class matrixrow<mat_type>;
};

///////////////////////////////////////////////////////////////////

//!
//! class matrixrow
//! templated class, parameter is a blockmatrix type like double or fixedMatrix
//! you get a matrixrow variable when A is a matrix and you do A[i] or A.getrow(i)
//! you can do A[i].getDiag() or A[i].addMatrixRow(c, nr);
//! iterators: matrixrow<mat_type>::citerator it(A[i]), ++it, (*it).iIndex, dValue.
//! also possible: double a = A[i]*x. (thats why i like this concept)
//! its not 100% atm, because that const matrix_type *A would rather be a const matrix_type &A,
//! but this gives all problems with const and non-const objects (like in iterator and citerator). Ideas?
template<typename mat_type>
class matrixrow
{
	typedef matrix<mat_type> matrix_type;
	//typedef typename matrix<mat_type>::vec_type vec_type;
	typedef typename matrix<mat_type>::connection connection;
	//typedef Vector<vec_type> Vector_type;
	// functions
public:
	matrixrow(const matrix_type *_A, const int _row)
	: row(_row), A(_A)
	{
		ASSERT2(row < A->length, *this);
	}
	
	matrixrow(matrix_type *_A, const int _row)
	: row(_row), A(_A)
	{
		ASSERT2(row < A->length, *this);
	}
	
	// get the sum of the entries of the matrix row
	//mat_type sum() const;
	
	//!
	//! removezeros
	//! remove zero entries of matrix (experimental)
	void removezeros();
	
	void setMatrixRow(connection *c, int nr);
	void addMatrixRow(connection *c, int nr);
	
	//!
	//! operator []:
	//! get connection nr i.
	inline const connection &operator [] (int i) const;	
	inline connection &operator [] (int i);
	
	inline int getConNr(int index) const;
	
	inline mat_type getDiag() const
	{
		ASSERT2(A->cons[row][0].iIndex == row, *this << " first entry has to be diagonal");
		return A->cons[row][0].dValue;
	}
	
	/*double& getCon(int index)
	 {
		 int i=getConNr(index);
		 if(i == -1)
		 {
		 addConnection(index, 0.0);
		 i = getNrOfConnections()-1;
		 }
		 
		 return A->cons[row][i].dValue;
	 }*/		

		
	//!
	//! operator * with Vector. This Vector is templated, since then one can do
	//! fixedMatrix<3> = matrixrow<double> * Vector<fixedMatrix<3> >, or
	//! Vector<fixedMatrix<3> > = matrix<double> * Vector<fixedMatrix<3> >
	//! like the prolongation/restriction in MG.
	template<typename vec_type>
	inline vec_type operator * (const Vector<vec_type> &x) const;
	
	template<typename vec_type>
	inline void copyToMult(vec_type &d, const Vector<vec_type> &x) const;
	template<typename vec_type>
	inline void addToMult(vec_type &d, const Vector<vec_type> &x) const;
	template<typename vec_type>
	inline void substractFromMult(vec_type &d, const Vector<vec_type> &x) const;
	
	// accessor functions
	inline bool indexWithinBounds(int i) const	{	return i < getNrOfConnections(); }
	inline int getRow() const { return row; }	
	inline int getNrOfConnections() const {	return A->iNrOfConnections[row];	}
	inline bool isUnconnected() const {	return A->isUnconnected(row); }
	
	// printing
	void printtype() const;
	friend ostream &operator<<(ostream &output, const matrixrow &r)
	{
		output << "matrixrow[row = " << r.row << "] of " << *(r.A);
		return output;
	}
	void print() const;
	void p(); // gdb

	// Iterators
	//---------------------------
	class citerator 
	{
	public:
		const matrixrow<mat_type> &row;
		int m_position;
	public:
		inline citerator(const matrixrow<mat_type> &r) : row(r) { rewind(); }
		//iterator(const iterator &other);
		
		inline const connection &operator *() const {return row[m_position];}
				
		inline void operator ++() {	m_position++; }
		
		inline void rewind() { m_position = 0;}
		inline int getPos() const{	return m_position;}
		
		inline const matrixrow<mat_type> &getMatrixRow() const { return row;}
		inline bool isEnd() const { return m_position >= row.getNrOfConnections(); }
	};
	
	class cLowerLeftIterator : public citerator
	{
	public:
		cLowerLeftIterator(const matrixrow<mat_type> &r) : citerator(r) { 	rewind();	}
		inline void operator ++()
		{
			citerator::m_position ++;
			while(citerator::m_position < citerator::row.getNrOfConnections() && (citerator::row[citerator::m_position].iIndex >= citerator::row.getRow() || citerator::row[citerator::m_position].dValue == 0.0)) 
				citerator::m_position++;
		}
		inline void rewind()
		{
			citerator::m_position = -1;
			operator ++ ();
		}
	};
	
	class cUpperRightIterator : public citerator
	{
	public:
		cUpperRightIterator(const matrixrow<mat_type> &r) : citerator(r) { 	rewind();	}
		inline void operator ++()
		{
			citerator::m_position ++;
			while(citerator::m_position < citerator::row.getNrOfConnections() && (citerator::row[citerator::m_position].iIndex <= citerator::row.getRow() || citerator::row[citerator::m_position].dValue == 0.0)) 
				citerator::m_position++;
		}
		inline void rewind()
		{
			citerator::m_position = -1;
			operator ++ ();
		}
	};
	
	// make this with templates?
	class iterator 
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
	};
	
		
//     data
//----------------
	
private:
	const matrix_type *A;
	const int row;
};


#include "Vector.hpp"
#include "matrix.hpp"