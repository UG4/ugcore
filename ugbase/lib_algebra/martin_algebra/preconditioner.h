/*
 *  preconditioner.h
 *  flexamg
 *
 *  Created by Martin Rupp on 29.09.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#pragma once
//#if UNKNOWN_NR > 1
#define USE_PREPED
//#endif

#include "sparseMatrix.h"

//!
//! preconditioner class
//! has following functions
//! init: initialize with Matrix A
//! precond(px, b). Apply x = B^{-1} b
//! iterate(px, b). Apply x = (I - B^{-1}A)x - B^{-1} b, that is x += B^{-1} (b - Ax)
//! question: is it faster to calculate r = b-Ax, precond(c, b), x+=c, or directly x += B^{-1} (b - Ax)?
//!
//! pre-calculation of diag(A)^{-1} speeds up algorithm considerably (2x + for 6x6 matrices) 
//!
template<typename Matrix_type, typename Vector_type>
class preconditioner
{
public:
	preconditioner() {}
	virtual ~preconditioner()  { destroy(); }
	virtual bool init(const Matrix_type &A)
	{
		pA = &A;
		return true;
	}

	//! calc x = (I - B^{-1}A)x - B^{-1} b, that is x += B^{-1} (b - Ax)
	virtual double iterate(Vector_type &x, const Vector_type &b)
	{
		ASSERT2(0, "preconditioner::iterate: dont use this function");
		return 0; // dont use
	}
	
	//! precond is like iterate with zero starting vector, calcs x = B^{-1} b	
	virtual void precond(Vector_type &x, const Vector_type &b)
	{
		x = b;		
	}
	
	virtual bool destroy() { return true; }
	
protected:
	const Matrix_type *pA;
};

//! class diagonalInversePreconditioner
//! base class for preconditioners like gs and jacobi, which need the inverse of diagonal elements of A
template<typename Matrix_type, typename Vector_type>
class diagonalInversePreconditioner : public preconditioner<Matrix_type, Vector_type>
{
public:
	typedef typename Matrix_type::entry_type entry_type;
	typedef typename matrix_trait<entry_type>::inverse_type inverse_type;
	typedef preconditioner<Matrix_type, Vector_type> super;
	
#ifdef USE_PREPED
	virtual bool init(const Matrix_type &A)
	{
		preconditioner<Matrix_type, Vector_type>::init(A);
		diagonal = new inverse_type[A.getLength()];
		for(int i=0; i<A.getLength(); i++)
			diagonal[i].setAsInverseOf(A.getDiag(i));
		return true;
	}
	
	virtual bool destroy() 
	{ 
		delete[] diagonal;
		return preconditioner<Matrix_type, Vector_type>::destroy(); 
	}
	
	inverse_type &getDiagInverse(int i)
	{
		return diagonal[i];
	}	
public:
	inverse_type *diagonal;
	
#else
	inverse_type &getDiagInverse(int i)
	{
		const Matrix_type &A = *preconditioner<Matrix_type, Vector_type>::pA;
		inv.setAsInverseOf(A.getDiag(i));
		return inv;
	}
	
	inverse_type inv;
#endif
	
};

//! specialization of diagonalInversePreconditioner for doubles: not pre-calculation needed.
template<typename Vector_type>
class diagonalInversePreconditioner<SparseMatrix<double>, Vector_type> : public preconditioner<SparseMatrix<double>, Vector_type>
{
public:
	double getDiagInverse(int i)
	{
		return 1.0 / preconditioner<SparseMatrix<double>, Vector_type>::pA->getDiag(i);
	}
};

//!
//! class jacobi
template<typename Matrix_type, typename Vector_type>
class jacobi : public diagonalInversePreconditioner<Matrix_type, Vector_type>
{
public:	
	typedef typename Matrix_type::entry_type entry_type;
	typedef typename matrix_trait<entry_type>::inverse_type inverse_type;
	
	virtual void precond(Vector_type &x, const Vector_type &b)
	{
		const Matrix_type &A = *preconditioner<Matrix_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		//	x = b / A.Diag();
		for(int i=0; i<x.getLength(); i++)
			x[i] = getDiagInverse(i) * b[i];
	}
	
	//! thats actually gauss-seidel o_O !!!
	double iterate(Vector_type &x, const Vector_type &b)
	{
		const Matrix_type &A = *preconditioner<Matrix_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		typename Vector_type::entry_type d;
		for(int j=0; j < A.getLength(); j++)
		{
			d = b[j];
			d -= A[j] * x;
			x[j] += getDiagInverse(j) * d;
		}
		return 1.0;
	}
	
	inverse_type &getDiagInverse(int i)
	{
		return diagonalInversePreconditioner<Matrix_type, Vector_type>::getDiagInverse(i);
	}
};

//!
//! class gs
template<typename Matrix_type, typename Vector_type>
class gs : public diagonalInversePreconditioner<Matrix_type, Vector_type>
{
public:
	typedef typename Matrix_type::entry_type entry_type;;
	typedef typename matrix_trait<entry_type>::inverse_type inverse_type;
	
	virtual void precond(Vector_type &x, const Vector_type &b)
	{			
		const Matrix_type &A = *preconditioner<Matrix_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		entry_type d;
		for(int i=0; i < x.getLength(); i++)
		{
			x[i] = b[i];
			for(typename matrixrow<entry_type>::cLowerLeftIterator it(A[i]); !it.isEnd(); ++it)
				x[i] -= ((*it).dValue * x[(*it).iIndex]);
			x[i] = getDiagInverse(i) * x[i];
		}
	}
	
	double iterate(Vector_type &x, const Vector_type &b)
	{
		const Matrix_type &A = *preconditioner<Matrix_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		typename Vector_type::entry_type d;
		for(int j=0; j < A.getLength(); j++)
		{
			d = b[j];
			d -= A[j] * x;
			x[j] += getDiagInverse(j) * d;
		}
		return 1.0;
	}
	 
	inverse_type &getDiagInverse(int i)
	{
		return diagonalInversePreconditioner<Matrix_type, Vector_type>::getDiagInverse(i);
	}
};

template <typename T>
inline void copyTo(T &dest, const T &src)
{
	dest = src;
}

//!
//! class sgs
template<typename Matrix_type, typename Vector_type>
class sgs : public diagonalInversePreconditioner<Matrix_type, Vector_type>
{
public:
	typedef typename Matrix_type::entry_type entry_type;;
	typedef typename matrix_trait<entry_type>::inverse_type inverse_type;
	
	// B = (D-R)^{-1} D (D-L)^{-1}
	// x = x + B^{-1}(b-Ax)

	virtual void precond(Vector_type &x, const Vector_type &b)
	{			
		const Matrix_type &A = *preconditioner<Matrix_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");

#if 1	// without temporaries
		typename Vector_type::entry_type s;
		for(int i=0; i < x.getLength(); i++)
		{
			//entry_type d= A.getDiag(i);
			s = b[i];
			for(typename SparseMatrix<entry_type>::cLowerLeftIterator it(A, i); !it.isEnd(); ++it)
				sub_mult(s, (*it).dValue, x[(*it).iIndex]);
			assign_mult(x[i], getDiagInverse(i), s);
		}
		// x = A.Diag() * x;
		
		for(int i = x.getLength()-1; i >= 0; i--)
		{
			//entry_type d= A.getDiag(i);
			//x[i] = x[i] / d;
			s = 0.0;
			for(typename SparseMatrix<entry_type>::cUpperRightIterator it(A, i); !it.isEnd(); ++it)
				sub_mult(s, (*it).dValue, x[(*it).iIndex]);
			add_mult(x[i], getDiagInverse(i), s);
		}		
#else	// with temporaries
		typename Vector_type::entry_type s;
		for(int i=0; i < x.getLength(); i++)
		{
			//entry_type d= A.getDiag(i);
			s = b[i];
			for(typename SparseMatrix<entry_type>::cLowerLeftIterator it(A, i); !it.isEnd(); ++it)
				s -= (*it).dValue * x[(*it).iIndex];
			x[i] = getDiagInverse(i) * s;
		}
		// x = A.Diag() * x;
		
		for(int i = x.getLength()-1; i >= 0; i--)
		{
			//entry_type d= A.getDiag(i);
			//x[i] = x[i] / d;
			s = 0.0;
			for(typename SparseMatrix<entry_type>::cUpperRightIterator it(A, i); !it.isEnd(); ++it)
				s -= (*it).dValue * x[(*it).iIndex];
			x[i] += getDiagInverse(i) * s;
		}
#endif

	}
	

	double iterate(Vector_type &x, const Vector_type &b)
	{
		const Matrix_type &A = *preconditioner<Matrix_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
#if 1   // without temporaries (10 secs compared to 16 secs)
		// todo: replace add_mult etc. with template expressions
		typename Vector_type::entry_type d;
		for(int j=0; j < A.getLength(); j++)
		{
			d = b[j];
			sub_mult(d, A[j], x);
			//x[j] += getDiagInverse(j)*d;
			add_mult(x[j], getDiagInverse(j), d);
		}
		
		//double res = 0;
		for(int j=A.getLength()-1; j >= 0; j--)
		{
			d = b[j];
			sub_mult(d, A[j], x);
			//x[j] += getDiagInverse(j) *d;
			add_mult(x[j], getDiagInverse(j), d);
		}
#else
		typename Vector_type::entry_type d;
		for(int j=0; j < A.getLength(); j++)
		{
			d = b[j];
			d -= A[j] * x;
			x[j] += getDiagInverse(j) * d;
		}
		
		//double res = 0;
		for(int j=A.getLength()-1; j >= 0; j--)
		{
			d = b[j];
			d -= A[j] * x;
			x[j] += getDiagInverse(j) * d;
		}
#endif		
		
		return 1.0; //sqrt(res);
	}
	
	inline inverse_type getDiagInverse(int i)
	{
		return diagonalInversePreconditioner<Matrix_type, Vector_type>::getDiagInverse(i);
	}
};