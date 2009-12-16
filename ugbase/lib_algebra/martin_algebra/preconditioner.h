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
template<typename mat_type>
class preconditioner
	{
		typedef matrix<mat_type> matrix_type;
		typedef typename matrix_type::vec_type vec_type;
		typedef Vector< typename matrix_type::vec_type> Vector_type;
		
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
template<typename mat_type>
class diagonalInversePreconditioner : public preconditioner<mat_type>
{
public:
	typedef matrix<mat_type> matrix_type;
	typedef typename matrix_type::vec_type vec_type;
	typedef Vector< typename matrix_type::vec_type> Vector_type;
	typedef typename matrix_trait<mat_type>::inverse_type inverse_type;
	
#ifdef USE_PREPED
	virtual bool init(const matrix_type &A)
	{
		preconditioner<mat_type>::init(A);
		diagonal = new inverse_type[A.getLength()];
		for(int i=0; i<A.getLength(); i++)
			diagonal[i].setAsInverseOf(A.getDiag(i));
		return true;
	}
	
	virtual bool destroy() 
	{ 
		delete[] diagonal;
		return preconditioner<mat_type>::destroy(); 
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
		const matrix_type &A = *preconditioner<mat_type>::pA;
		inv.setAsInverseOf(A.getDiag(i));
		return inv;
	}
	
	inverse_type inv;
#endif
	
};

template<>
class diagonalInversePreconditioner<double> : public preconditioner<double>
{
public:
	double inv;
	double &getDiagInverse(int i)
	{
		inv = 1/ preconditioner<double>::pA->getDiag(i);
		return inv;
	}
};

//!
//! class jacobi
template<typename mat_type>
class jacobi : public diagonalInversePreconditioner<mat_type>
{
public:
	typedef matrix<mat_type> matrix_type;
	typedef typename matrix_type::vec_type vec_type;
	typedef Vector< typename matrix_type::vec_type> Vector_type;
	typedef typename matrix_trait<mat_type>::inverse_type inverse_type;
	
	virtual void precond(Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		//	x = b / A.Diag();
		for(int i=0; i<x.getLength(); i++)
			x[i] = getDiagInverse(i) * b[i];
	}
	
	// thats actually gauss-seidel o_O
	double iterate(Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		double res = 0;
		vec_type d;
		for(int j=0; j < A.getLength(); j++)
		{
			d = (b[j] - A[j] * x);
			x[j] +=  getDiagInverse(j) * d;
			res += d*d;
		}
		return sqrt(res);
	}
	
	inverse_type &getDiagInverse(int i)
	{
		return diagonalInversePreconditioner<mat_type>::getDiagInverse(i);
	}
};

//!
//! class gs
template<typename mat_type>
class gs : public diagonalInversePreconditioner<mat_type>
{
public:
	typedef matrix<mat_type> matrix_type;
	typedef typename matrix_type::vec_type vec_type;
	typedef Vector< typename matrix_type::vec_type> Vector_type;
	typedef typename matrix_trait<mat_type>::inverse_type inverse_type;
	
	double dampediter(double damp, Vector_type *pc, Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA; Vector_type &c = *pc;
		double res=0;
		vec_type d;
		int j;
		
		c = x;
		for(j=0; j < A.getLength(); j++)
		{			
			d = (b[j] - A[j] * x);
			x[j] += getDiagInverse(j) *  d;
			res += d*d;
		}		
		x -= (1-damp)*(x-c);	
		
		return sqrt(res);
	}
	
	
	virtual void precond(Vector_type &x, const Vector_type &b)
	{			
		const matrix_type &A = *preconditioner<mat_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		mat_type d;
		for(int i=0; i < x.getLength(); i++)
		{
			x[i] = b[i];
			for(typename matrixrow<mat_type>::cLowerLeftIterator it(A[i]); !it.isEnd(); ++it)
				x[i] -= ((*it).dValue * x[(*it).iIndex]);
			x[i] = getDiagInverse(i) * x[i];
		}
	}
	double iterate(Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		for(int j=0; j < A.getLength(); j++)
			x[j] += getDiagInverse(j) * (b[j] - A[j] * x);			
		
		return 0.0;
	}
	 
	inverse_type &getDiagInverse(int i)
	{
		return diagonalInversePreconditioner<mat_type>::getDiagInverse(i);
	}
};

//!
//! class sgs
template<typename mat_type>
class sgs : public diagonalInversePreconditioner<mat_type>
{
public:
	typedef matrix<mat_type> matrix_type;
	typedef typename matrix_type::vec_type vec_type;
	typedef Vector< typename matrix_type::vec_type> Vector_type;
	typedef typename matrix_trait<mat_type>::inverse_type inverse_type;
	
	// B = (D-R)^{-1} D (D-L)^{-1}
	// x = x + B^{-1}(b-Ax)
	virtual void precond(Vector_type &x, const Vector_type &b)
	{			
		const matrix_type &A = *preconditioner<mat_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		vec_type s;
		for(int i=0; i < x.getLength(); i++)
		{
			//mat_type d= A.getDiag(i);
			s = b[i];
			for(typename matrixrow<mat_type>::cLowerLeftIterator it(A[i]); !it.isEnd(); ++it)
				s -= ((*it).dValue * x[(*it).iIndex]);
			x[i] = getDiagInverse(i) * s;
		}
		// x = A.Diag() * x;
		
		for(int i = x.getLength()-1; i >= 0; i--)
		{
			//mat_type d= A.getDiag(i);
			//x[i] = x[i] / d;
			s = 0.0;
			for(typename matrixrow<mat_type>::cUpperRightIterator it(A[i]); !it.isEnd(); ++it)
				s -= ((*it).dValue *x[(*it).iIndex]);
			x[i] += getDiagInverse(i) * s;
		}		
	}
	
	double iterate(Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		vec_type d;
		for(int j=0; j < A.getLength(); j++)
			x[j] += getDiagInverse(j) * (b[j] - A[j] * x);
		
		//double res = 0;
		for(int j=A.getLength()-1; j >= 0; j--)
		{
			d = (b[j] - A[j] * x);
			x[j] += getDiagInverse(j) *d;
			//x[j] += d / A.getDiag(j);
			//res += mnorm2(d);
		}
		
		return 1.0; //sqrt(res);
	}
	
	inline inverse_type &getDiagInverse(int i)
	{
		return diagonalInversePreconditioner<mat_type>::getDiagInverse(i);
	}
};