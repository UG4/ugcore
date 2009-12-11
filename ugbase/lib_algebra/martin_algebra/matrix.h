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
class untermatrix
{
public:
	untermatrix(int *indices_, int nr_)
	{
		nr = nr_;
		values = new mat_type[nr*nr];
		for(int i=0; i<nr*nr; i++) values[i] = 0.0;
		
		indices = new int[nr];
		memcpy(indices, indices_, sizeof(int)*nr);
	}
	~untermatrix()
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
	};
	
	// constructor for empty matrix
	matrix(const char *_name = "");
	// constructor for matrix with length
	matrix (int _length, const char *name = "");	
	~matrix ();	
	
	matrix(matrix&); // disallow copy operator
	void operator = (const matrix &v); // disallow assignment
	
	void eliminateDirichletValues(Vector<vec_type> &b); 
	void createAsTransposeOf(const matrix &B, int length);
	
	void create(int _length);
	void saveSetConnections(int row, connection *mem) const;
	void defrag();
	
	inline matrixrow_type getrow(int i);
	inline const matrixrow_type getrow(int i) const;
	
	matrixrow_type operator [] (int i) ;	
	const matrixrow_type operator [] (int i) const;

	void print(const char * const name = NULL) const;
	void p(); // gdb
	
	void addUnterMatrix(untermatrix<mat_type> &mat);
	
	friend ostream &operator<<(ostream &output, const matrix &m)
	{
		output << "matrix " << m.name << "[" << m.length << "]";
		return output;
	}
	void printtype() const; 

	void recreateWithMaxNrOfConnections(int newMax) const;
	void finish()
	{
	//	defrag();
	}
	
	inline const mat_type getDiag(int i) const;
	inline mat_type &getDiag(int i);	
	inline bool isUnconnected(int i) const;
	
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
	
	inline void prefetch(int i) const
	{
		prefetchRead(cons[i]);
	}
	
	int getLength() const { return length; }
	int getTotalNrOfConnections() const { return iTotalNrOfConnections; }
	
	// data
	// this is a little bit stupid with mutable bcs of the matrixrow class
	
public:
	int fromlevel;
	int tolevel;
	const char *name;

private:	
	int length;	
	mutable int iTotalNrOfConnections;
	mutable int iFragmentedMem;
	mutable connection **cons;

	mutable connection *consmem;
	int consmemsize;
	
	int *iNrOfConnections;
	mutable int bandwidth;
	
	friend class matrixrow<mat_type>;
};

///////////////////////////////////////////////////////////////////
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
	
	inline const connection &operator [] (int i) const;	
	inline connection &operator [] (int i);
	

	class citerator 
	{
	public:
		const matrixrow<mat_type> &row;
		int m_position;
	public:
		inline citerator(const matrixrow<mat_type> &r) : row(r) { rewind(); }
		//iterator(const iterator &other);
		
		inline const connection &operator *() const
		{
			return row[m_position];
		}
				
		inline void operator ++()
		{
			m_position++;			
		}
		
		inline void rewind()
		{
			m_position = 0;
		}
		inline int getPos() const
		{
			return m_position;
		}
		
		inline const matrixrow<mat_type> &getMatrixRow() const
		{
			return row;
		}
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
	
	// make this with templates:
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
	
		
	inline bool indexWithinBounds(int i) const	{	return i < getNrOfConnections(); }
	
	inline mat_type getDiag() const
	{
		ASSERT2(A->cons[row][0].iIndex == row, *this << " first entry has to be diagonal");
		return A->cons[row][0].dValue;
	}
	
	inline int getNrOfConnections() const {	return A->iNrOfConnections[row];	}
	inline bool isUnconnected() const {	return A->isUnconnected(row); }
	
	// get the sum of the entries of the matrix row
	//mat_type sum() const;
	
	void printtype() const;
	
	friend ostream &operator<<(ostream &output, const matrixrow &r)
	{
		output << "matrixrow[row = " << r.row << "] of " << *(r.A);
		return output;
	}
	void removezeros();
	
	void setMatrixRow(connection *c, int nr);
	void addMatrixRow(connection *c, int nr);
	
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
	
	// other
	
	inline int getConNr(int index) const;		
	template<typename vec_type>
	inline vec_type operator * (const Vector<vec_type> &x) const;
	
	void print() const;
	void p(); // gdb
	
	inline int getRow() const { return row; }	
	
	// data
	
private:
	const matrix_type *A;
	const int row;
};


#include "Vector.hpp"
#include "matrix.hpp"