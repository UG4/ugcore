/*
 *  preconditioner.h
 *  flexamg
 *
 *  Created by Martin Rupp on 29.09.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#pragma once

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
	

	virtual double iterate(Vector_type *px, const Vector_type &b)
	{
		ASSERT2(0, "preconditioner::iterate: dont use this function");
		return 0; // dont use
	}

	// precond is like iterate with zero starting Vector
	virtual void precond2(Vector_type *px, const Vector_type &b)
	{
		(*px) = b;		
	}
	
	class preconditioner_parameter : public FunctionExpression<vec_type>
	{
	public:
		preconditioner_parameter(preconditioner *pre_, Vector_type &v) : vec(v) { pre = pre_;}  
		preconditioner *pre;
		Vector_type &vec;
		virtual void applyto(Vector_type &x)
		{
			pre->precond2(&x, vec);
		}
	};	
	
	virtual FunctionExpression<vec_type> &precond(Vector_type &v)
	{
		static preconditioner_parameter par(this, v);
		return par;
	}
	
	virtual bool destroy() { return true; }
	
protected:
	const matrix_type *pA;
};

template<typename mat_type>
class smoother : public preconditioner<mat_type>
	{
		typedef matrix<mat_type> matrix_type;
		typedef typename matrix_type::vec_type vec_type;
		typedef Vector< typename matrix_type::vec_type> Vector_type;
		
		smoother() {}
		virtual ~smoother()  {}
		
		virtual void smoothRes(Vector_type *x, const Vector_type &res) = 0;
	};

template<typename mat_type>
class jacobi : public preconditioner<mat_type>
{
public:
	typedef matrix<mat_type> matrix_type;
	typedef typename matrix_type::vec_type vec_type;
	typedef Vector< typename matrix_type::vec_type> Vector_type;
	
	virtual void precond2(Vector_type *px, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA; Vector_type &x = *px;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		//	x = b / A.Diag();
		for(int i=0; i<x.getLength(); i++)
			x[i] = b[i] / A.getDiag(i);
	}
	
		
	double iterate(Vector_type *px, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA; Vector_type &x = *px;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		double res = 0;
		vec_type d;
		for(int j=0; j < A.getLength(); j++)
		{
			d = (b[j] - A[j] * x);
			x[j] +=  d / A.getDiag(j);			
			res += d*d;
		}
		return sqrt(res);
	}
};

template<typename mat_type>
class gs : public preconditioner<mat_type>
{
public:
	typedef matrix<mat_type> matrix_type;
	typedef typename matrix_type::vec_type vec_type;
	typedef Vector< typename matrix_type::vec_type> Vector_type;
	
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
			x[j] += d/ A.getDiag(j);			
			res += d*d;
		}		
		x -= (1-damp)*(x-c);	
		
		return sqrt(res);
	}

	double dampediter2(double damp, Vector_type *pc, Vector_type &x, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA; Vector_type &c = *pc;
		mat_type d;
		c = b - A*x;
		//double res = norm(c);
		for(int i=0; i < x.getLength(); i++)
		{
			d = A.getDiag(i);
			c[i] = c[i] / d;
			for(typename matrixrow<mat_type>::cLowerLeftIterator it(A[i]); !it.isEnd(); ++it)
				c[i] -= (c[(*it).iIndex] * (*it).dValue) / d;
			//x[i] += damp*c[i];
		}		
		x += damp * c;
		return 1.0; //res;
	}
	
	virtual void precond2(Vector_type *px, const Vector_type &b)
	{			
		const matrix_type &A = *preconditioner<mat_type>::pA; Vector_type &x = *px;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		mat_type d;
		for(int i=0; i < x.getLength(); i++)
		{
			d= A.getDiag(i);
			x[i] = b[i] / d;
			for(typename matrixrow<mat_type>::cLowerLeftIterator it(A[i]); !it.isEnd(); ++it)
				x[i] -= x[(*it).iIndex] * (*it).dValue / d;
		}
	}
	double iterate(Vector_type *px, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA; Vector_type &x = *px;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");
		
		for(int j=0; j < A.getLength(); j++)
			x[j] += (b[j] - A[j] * x) / A.getDiag(j);
		
		return 0.0;
	}
};

template<typename mat_type>
class sgs : public preconditioner<mat_type>
{
public:
	typedef matrix<mat_type> matrix_type;
	typedef typename matrix_type::vec_type vec_type;
	typedef Vector< typename matrix_type::vec_type> Vector_type;
	typedef typename matrix_trait<mat_type>::inverse_type inverse_type;
	//typedef mat_type inverse_type;
	
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
	
	// B = (D-R)^{-1} D (D-L)^{-1}
	// x = x + B^{-1}(b-Ax)
	virtual void precond2(Vector_type *px, const Vector_type &b)
	{			
		const matrix_type &A = *preconditioner<mat_type>::pA; Vector_type &x = *px;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");

		for(int i=0; i < x.getLength(); i++)
		{
			mat_type d= A.getDiag(i);
			x[i] = b[i] / d;
			for(typename matrixrow<mat_type>::cLowerLeftIterator it(A[i]); !it.isEnd(); ++it)
				x[i] -= ((*it).dValue * x[(*it).iIndex]) / d;
		}
		// x = A.Diag() * x;

		for(int i = x.getLength()-1; i >= 0; i--)
		{
			mat_type d= A.getDiag(i);
			//x[i] = x[i] / d;
			for(typename matrixrow<mat_type>::cUpperRightIterator it(A[i]); !it.isEnd(); ++it)
				x[i] -= ((*it).dValue *x[(*it).iIndex]) / d;
		}		
	}
	
	double iterate(Vector_type *px, const Vector_type &b)
	{
		const matrix_type &A = *preconditioner<mat_type>::pA; Vector_type &x = *px;
		ASSERT2(x.getLength() == b.getLength() && x.getLength() == A.getLength(), x << ", " << b << " and " << A << " need to have same size.");

		vec_type d;
		for(int j=0; j < A.getLength(); j++)
			x[j] += diagonal[j] * (b[j] - A[j] * x);
		
		//double res = 0;
		for(int j=A.getLength()-1; j >= 0; j--)
		{
			d = (b[j] - A[j] * x);
			x[j] += diagonal[j] * d;
			//res += mnorm2(d);
		}
		
		return 1.0; //sqrt(res);
	}
	
private:
	inverse_type *diagonal;
};
