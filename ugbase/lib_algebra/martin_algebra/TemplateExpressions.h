/*
 *  TemplateExpressions.h
 *  flexamg
 *
 *  Created by Martin Rupp on 29.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */

#pragma once
#include "blockMatrix.h"

/////////////////////////////////////////////////////////////////

//!
//! class XD: class for template Expressions.
//! only classes which inherit from XD (via class myClass : public XD<myClass> )
//! can use templateExpressions used in this file
template<class A> class XD
{
public:
	const A& cast() const {return static_cast<const A&>(*this); }
};

/////////////////////////////////////////////////////////////////
//!
//! general template Expression: L Operator R.
//! all Expressions support
//! - operator []
//! ostream operator <<, printtype
//! getLength
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
		// dont use
		return Operator::apply( l[i], r[i] );
	}
	
	
	inline void copyTo(typename Operator::ReturnType &d, int i) const
	{
		Operator::copyTo(d, l[i], r[i]);
	}

	inline void addTo(typename Operator::ReturnType &d, int i) const
	{
		Operator::addTo(d, l[i], r[i]);
	}
	
	inline void substractFrom(typename Operator::ReturnType &d, int i) const
	{
		Operator::substractFrom(d, l[i], r[i]);
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


//!
//! template Expression: double Operator R.
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
	
	inline void copyTo(typename Operator::ReturnType &d, int i) const
	{
		Operator::copyTo(d, ld, r[i]);
	}
	
	inline void addTo(typename Operator::ReturnType &d, int i) const
	{
		Operator::addTo(d, ld, r[i]);
	}
	
	inline void substractFrom(typename Operator::ReturnType &d, int i) const
	{
		Operator::substractFrom(d, ld, r[i]);
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

// Operator Structs, for defining Operator behaviour
/////////////////////////////////////////////////////////////////

template<typename vec_type>
struct Add_Operator
{ 
	typedef vec_type ReturnType;
	static inline vec_type apply(const vec_type &a, const vec_type &b) {return a + b;} 
	
	static inline void copyTo(vec_type &r, const vec_type &a, const vec_type &b) {r = a; r += b;} 
	static inline void addTo(vec_type &r, const vec_type &a, const vec_type &b) {r += a; r += b;} 
	static inline void substractFrom(vec_type &r, const vec_type a, const vec_type &b) {r = a; r += b; r *= -1.0;} 
	
	static inline const char *cTyp() { return " + "; }
}; 

template<typename vec_type>
struct Minus_Operator
{ 
	typedef vec_type ReturnType;
	static inline vec_type apply(const vec_type &a, const vec_type &b) {return a - b;} 
	
	static inline void copyTo(vec_type &r, const vec_type &a, const vec_type &b) {r = a; r -= b;} 
	static inline void addTo(vec_type &r, const vec_type &a, const vec_type &b) {r += a; r -= b;} 
	static inline void substractFrom(vec_type &r, const vec_type a, const vec_type &b) {r = a; r -= b; r *= -1.0;} 
	
	static inline const char *cTyp() { return " - "; }
}; 

template<typename mat_type, typename vec_type>
struct Multiply_Operator
{ 
	typedef mat_type MultType;
	typedef typename Mult_Traits<mat_type, vec_type>::ReturnType ReturnType;
	static inline ReturnType apply(const vec_type &a, const vec_type &b) {return a * b;} 
	
	static inline void copyTo(ReturnType &r, const vec_type &a, const vec_type &b) {r = a * b;} 
	static inline void addTo(ReturnType &r, const vec_type &a, const vec_type &b) {r += a * b;} 
	static inline void substractFrom(ReturnType &r, const vec_type &a, const vec_type &b) {r -= a * b;} 
	
	static inline const char *cTyp() { return " * "; }
};
	
/*struct Divide_Operator
{ 
	static inline VEC_TYPE apply(VEC_TYPE a, MAT_TYPE b) {return a / b;} 
	static inline const char *cTyp() { return " / "; }
};*/

// operator x -> X_Operator
/////////////////////////////////////////////////////////////////

//
// allow + for all XDs
template<typename L, typename R> Expression< L, Add_Operator< typename L::vec_type >, R> operator+(const XD<L> &l,const XD<R> &r)
{ 
	return Expression<L, Add_Operator< typename L::vec_type >, R> (l.cast(), r.cast()); 
}

//
// allow - for all XDs
template<typename L, typename R> Expression< L, Minus_Operator< typename L::vec_type >, R> operator-(const XD<L> &l,const XD<R> &r)
{ 
	return Expression<L, Minus_Operator< typename L::vec_type >, R> (l.cast(), r.cast()); 
}

//
// allow * for doubles and all XDs
template<typename R> Expression<double, Multiply_Operator<double, typename R::vec_type >, R> operator*(double l,const XD<R> &r)
{ 
	return Expression<double, Multiply_Operator<double, typename R::vec_type>, R> (l, r.cast()); 
}
// * and / only for special



/////////////////////////////////////////////////////////////////

//!
//! template expression norm2
template<typename Type> 
inline double norm2(const XD<Type> &t_)
{
	const Type &t = t_.cast();
	double sum=0;
	for(int i=0; i < t.getLength(); i++)	
		sum += mnorm2(t[i]);
	return sum;
}

//!
//! template expression norm
template<typename Type> 
inline double norm(const XD<Type> &t)
{
	return sqrt(norm2<Type>(t));
}




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
 x = 5.0*y; (mit spezialisierung Expression<double, Multiply_Operator, Vector>
 
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
 x.operator = ( Expression(Vector, Add_Operator, Vector)(y, z);)
 das dröselt sich dann auf zu
 x.operator = (expression)
 {
 for(int i=0; i<x.length; i++)
 {
 x[i] = expression[i];
 }
 }
 
 wobei expression = Expression(Vector, Add_Operator, Vector)(y, z), also
 double Expression::operator [] (i)
 {
 return Operator::apply(l[i], r[i]) 
 => Add_Operator::apply(x[i], y[i]) ( da Operator = Add_Operator, l = x und r = y)
 ={ return a + b; }	
 }
 
 bei komplizierteren genauso, zB. b - A*x
 -> Expression(Vector, Minus_Operator, Expression(SparseMatrix, Mult_Operator, Vector))
 und exp[i] wird dann zu b[i] - A[i]*x 
 (sh. spezialisierung template<> struct Expression<SparseMatrix, Multiply_Operator, Vector> in SparseMatrix.h) 
 */