/*
 *  preconditioner.h
 *  flexamg
 *
 *  Created by Martin Rupp on 29.09.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#pragma once
#if UNKNOWN_NR > 1
#define USE_PREPED
#endif

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
template<typename entry_type, typename Vector_type>
class preconditioner
{
	typedef SparseMatrix<entry_type> matrix_type;
	
public:
	preconditioner() {}
	virtual ~preconditioner()  { destroy(); }
	virtual bool init(const matrix_type &A)
	{
		pA = &A;
		return true;
	}
	
	virtual double iterate(Vector_type &x, const Vector_type &b)
	{
		ASSERT2(0, "preconditioner::iterate: dont use this function");
		return 0; // dont use
	}
	
	// precond is like iterate with zero starting Vector
	virtual void precond(Vector_type &x, const Vector_type &b)
	{
		x = b;		
	}
	
	virtual bool destroy() { return true; }
	
protected:
	const matrix_type *pA;
};

//! class diagonalInversePreconditioner
//! base class for preconditioners like gs and jacobi, which need the inverse of diagonal elements of A
template<typename entry_type, typename Vector_type>
class diagonalInversePreconditioner : public preconditioner<entry_type, Vector_type>
{
public:
	typedef SparseMatrix<entry_type> matrix_type;
	typedef typename matrix_trait<entry_type>::inverse_type inverse_type;
	
#ifdef USE_PREPED
	virtual bool init(const matrix_type &A)
	{
		preconditioner<entry_type>::init(A);
		diagonal = new inverse_type[A.getLength()];
		for(int i=0; i<A.getLength(); i++)
			diagonal[i].setAsInverseOf(A.getDiag(i));
		return true;
	}
	
	virtual bool destroy() 
	{ 
		delete[] diagonal;
		return preconditioner<entry_type>::destroy(); 
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
		const matrix_type &A = *preconditioner<entry_type, Vector_type>::pA;
		inv.setAsInverseOf(A.getDiag(i));
		return inv;
	}
	
	inverse_type inv;
#endif
	
};

template<typename Vector_type>
class diagonalInversePreconditioner<double, Vector_type> : public preconditioner<double, Vector_type>
{
public:
	double inv;
	double &getDiagInverse(int i)
	{
		inv = 1/ preconditioner<double, Vector_type>::pA->getDiag(i);
		return inv;
	}
};

//!
//! class jacobi
template<typename entry_type, typename Vector_type>
class jacobi : public diagonalInversePreconditioner<entry_type, Vector_type>
{
public:
	typedef SparseMatrix<entry_type> matrix_type;
	typedef typename matrix_trait<entry_type>::inverse_type inverse_type;
	
	virtual void precond(Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<entry_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		//	x = b / A.Diag();
		for(int i=0; i<x.getLength(); i++)
			x[i] = getDiagInverse(i) * b[i];
	}
	
	// thats actually gauss-seidel o_O
	double iterate(Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<entry_type, Vector_type>::pA;
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
		return diagonalInversePreconditioner<entry_type, Vector_type>::getDiagInverse(i);
	}
};

//!
//! class gs
template<typename entry_type, typename Vector_type>
class gs : public diagonalInversePreconditioner<entry_type, Vector_type>
{
public:
	typedef SparseMatrix<entry_type> matrix_type;
	typedef typename matrix_trait<entry_type>::inverse_type inverse_type;
	
	virtual void precond(Vector_type &x, const Vector_type &b)
	{			
		const matrix_type &A = *preconditioner<entry_type, Vector_type>::pA;
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
		const matrix_type &A = *preconditioner<entry_type, Vector_type>::pA;
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
		return diagonalInversePreconditioner<entry_type, Vector_type>::getDiagInverse(i);
	}
};

template <typename T>
inline void copyTo(T &dest, const T &src)
{
	dest = src;
}

//!
//! class sgs
template<typename entry_type, typename Vector_type>
class sgs : public diagonalInversePreconditioner<entry_type, Vector_type>
{
public:
	typedef SparseMatrix<entry_type> matrix_type;
	typedef typename matrix_trait<entry_type>::inverse_type inverse_type;
	
	// B = (D-R)^{-1} D (D-L)^{-1}
	// x = x + B^{-1}(b-Ax)

	virtual void precond(Vector_type &x, const Vector_type &b)
	{			
		const matrix_type &A = *preconditioner<entry_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");

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
	}
	

	double iterate(Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<entry_type, Vector_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
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
		
		return 1.0; //sqrt(res);
	}
	
	inline inverse_type &getDiagInverse(int i)
	{
		return diagonalInversePreconditioner<entry_type, Vector_type>::getDiagInverse(i);
	}
};