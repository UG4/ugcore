/*
 *  TemplateExpressions.h
 *  flexamg
 *
 *  Created by Martin Rupp on 29.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */

#pragma once
#include "smallMatrix.h"

//typedef double MAT_TYPE;
//typedef double VEC_TYPE;


template<class A> class XD
{
public:
	const A& cast() const {return static_cast<const A&>(*this); }
};

/*
class myClass : public XD<myClass>
{
	// implementation
}

template<typename tempType>
void doit(XD<tempType> bla)
{
	tempType c = bla.cast();
}
*/

// TemplateExpressions Assignment
//---------------------------------------------
// zB. v = x + y, v = x + y + z
// v = A * x, v = x + y - A*x, v -= A*b
// !NICHT! x = A*x

/*
 möglich ist zB. auch 
 (x - A*x + b).printtype() (ausgabe der internen Expression)
 norm2(x - A*b)  (entspricht for(...) d += (x[i] - A[i]*b)^2 )
 norm(ex), maxabs(ex)
 (x - A*b)[5]
 x = A.Diag() * b; (sh. diagcomponent)
 x = 5.0*y; (mit spezialisierung Expression<double, Multiply_Exp, Vector>
 
 */


/* Function Expression
 ermöglicht Dinge wie zB. x = interpolate(xH);
 zB.
 class assignall_par 
 {
 public:
	assignall_par(int i) {n=i;} 
	void applyto(Vector &v) { for(int j=0; j<v.length; j++) v[j] = n;}
	int n; 
 }
 assignall_par assignall(int i) 
 {
	return assignall_par(i); 
 }
 
 x = assignall(5);  
*/

/*
 intern Funktioniert das so, zB.:
 x = y + z;
 x.operator = ( Expression(Vector, Add_Exp, Vector)(y, z);)
 das dröselt sich dann auf zu
 x.operator = (expression)
 {
	for(int i=0; i<x.length; i++)
	{
		x[i] = expression[i];
	}
 }
 
 wobei expression = Expression(Vector, Add_Exp, Vector)(y, z), also
 double Expression::operator [] (i)
 {
	return Operator::apply(l[i], r[i]) 
	=> Add_Exp::apply(x[i], y[i]) ( da Operator = Add_Exp, l = x und r = y)
		={ return a + b; }	
 }
 
 bei komplizierteren genauso, zB. b - A*x
 -> Expression(Vector, Minus_Exp, Expression(matrix, Mult_Exp, Vector))
 und exp[i] wird dann zu b[i] - A[i]*x 
 (sh. spezialisierung template<> struct Expression<matrix, Multiply_Exp, Vector> in matrix.h) 
 */


template<typename L, typename Operator, typename R> 
class Expression : public XD<Expression<L, Operator, R> >
{ 
public:
	const L& l; 
	const R& r; 
	inline Expression(const L& l_, const R& r_) : l(l_),r(r_) 
	{ ASSERT2(l.getLength() == r.getLength(), l << " has different length as " <<  r); } 
		
	inline typename Operator::ReturnType operator [] (int i) const
	{
		return Operator::apply( l[i], r[i] );
	}
	inline int getLength() const	{	return l.getLength();	}

	
	// print functions
	friend ostream &operator<<(ostream &output, const Expression<L, Operator, R>  &ex)
	{
		output << "(" << ex.l << Operator::cTyp() << ex.r << ")"; 
		return output;
	}
	inline void printtype() const	{	cout << *this; }
}; 


///////////////////////////////////////////////////////////////////
template<typename vec_type> class Vector;
template<typename vec_type>
class FunctionExpression
{
public:
	virtual void applyto(Vector<vec_type> &x) {}
};


template<typename Operator, typename R> 
class Expression<double, Operator, R> : public XD<Expression<double, Operator, R> >
{ 
public:
	const double ld; 
	const R& r; 
	inline Expression(const double& ld_, const R& r_) : ld(ld_),r(r_) {} 
	
	inline typename Operator::ReturnType operator [] (int i) const
	{
		return Operator::apply( ld, r[i] );
	} 
	
	inline int getLength() const	{	return r.getLength();	}
	
	
	// print functions
	friend ostream &operator<<(ostream &output, const Expression<double, Operator, R>  &ex)
	{
		output << "( double " << ex.ld << Operator::cTyp() << ex.r << ")"; 
		return output;
	}
	inline void printtype() const	{	cout << *this; }
}; 


template<typename vec_type>
struct Add_Exp
{ 
	typedef vec_type ReturnType;
	static inline vec_type apply(vec_type a, vec_type b) {return a + b;} 
	static inline const char *cTyp() { return " + "; }
}; 

template<typename vec_type>
struct Minus_Exp
{ 
	typedef vec_type ReturnType;
	static inline vec_type apply(vec_type a, vec_type b) {return a - b;} 
	static inline const char *cTyp() { return " - "; }
}; 

template<typename mat_type, typename vec_type>
struct Multiply_Exp
{ 
	typedef mat_type MultType;
	typedef typename Mult_Traits<mat_type, vec_type>::ReturnType ReturnType;
	static inline ReturnType apply(mat_type a, vec_type b) {return a * b;} 
	static inline const char *cTyp() { return " * "; }
};
	
/*struct Divide_Exp
{ 
	static inline VEC_TYPE apply(VEC_TYPE a, MAT_TYPE b) {return a / b;} 
	static inline const char *cTyp() { return " / "; }
};*/
	
// allow + for all
template<typename L, typename R> Expression< L, Add_Exp< typename L::vec_type >, R> operator+(const XD<L> &l,const XD<R> &r)
{ 
	return Expression<L, Add_Exp< typename L::vec_type >, R> (l.cast(), r.cast()); 
}


// allow - for all
template<typename L, typename R> Expression< L, Minus_Exp< typename L::vec_type >, R> operator-(const XD<L> &l,const XD<R> &r)
{ 
	return Expression<L, Minus_Exp< typename L::vec_type >, R> (l.cast(), r.cast()); 
}

// allow * for doubles and all.
template<typename R> Expression<double, Multiply_Exp<double, typename R::vec_type >, R> operator*(double l,const XD<R> &r)
{ 
	return Expression<double, Multiply_Exp<double, typename R::vec_type>, R> (l, r.cast()); 
}
// * and / only for special

/*inline double norm2(double a)
{
	return a*a;
}
inline double norm(double a)
{
	return (a);
}*/

template<typename Type> 
inline double norm2(const XD<Type> &t_)
{
	const Type &t = t_.cast();
	double sum=0;
	for(int i=0; i < t.getLength(); i++)	
	{		
		sum += mnorm2(t[i]);
	}
	return sum;
}

template<typename Type> 
inline double norm(const XD<Type> &t)
{
	return sqrt(norm2<Type>(t));
}


