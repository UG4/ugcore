/*
 *  template_expressions.h
 *
 *  Created by Martin Rupp on 08.02.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__MARTIN_ALGEBRA__TEMPLATE_EXPRESSIONS__
#define __H__UG__MARTIN_ALGEBRA__TEMPLATE_EXPRESSIONS__

//#include "blockMatrix.h"



using namespace std;

namespace ug{

// the problem:
// template<typename T> void bla(T &t) { }
// can be used for all kinds of T. if we do not want this, we can use CRTP, the Curiously recurring template pattern
// (http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern)
// example:
// class myVectorClass : public TE_VEC<myVectorClass> { ... }
// now one can write functions
// template<typename T> void bla(TE_VEC<T> &ttvec) { T &t = ttvec.cast(); // see cast() }
// which only work for functions with the "attribute" TE_VEC.

// if you want to use the template expression mechanism for your vector class myVectorClass, then write
// class myVectorClass : public TE_VEC<myVectorClass> { ... }
// and add functions
// template<typename T> myVectorClass& operator = (const TE_AMV_X<T> &t) { VectorAssign(&this, t); return *this; }
// template<typename T> myVectorClass& operator += (const TE_AMV_X<T> &t) { VectorAdd(&this, t); return *this; }
// template<typename T> myVectorClass& operator -= (const TE_AMV_X<T> &t) { VectorSub(&this, t); return *this; }
// if your class supports operator [] and int size(), then you can use the standart VecScaleAdd functions of operations_vec.h
// otherwise, declare your own
// void VecScale(myVectorClass &dest, double alpha1, myVectorClass &v1)
// { // dest = alpha1*v1; // .. }
// void VecScaleAdd(myVectorClass &dest, double alpha1, myVectorClass &v1, double alpha2, myVectorClass &v2)
// { // dest = alpha1*v1 + alpha2*v2; // .. }
// void VecScaleAdd(myVectorClass &dest, double alpha1, myVectorClass &v1, double alpha2, myVectorClass &v2, double alpha3, myVectorClass &v3)
// { // dest = alpha1*v1 + alpha2*v2 + alpha3*v3; // .. }
// and you are ready to go. an expression like x += 0.6*y - z; is then transformed to VecScaleAdd(x, 0.6, y, -1.0, z, 1.0, x);


//! only classes which inherit from TE_AMV_X, TE_MAT or TE_VEC (via class myClass : public TE_VEC<myClass> )
//! can use templateExpressions used in this file

////////////////////////////////////////////////////////////////////////////////

//! this helper class is a transposed of class A
template<class A> class TRANSPOSED
{
public:
	TRANSPOSED(const A &a_) : a(a_) {}
	//! get untransposed original class A.
	const A& T() const {return a; }
private:
	const A &a;
};

////////////////////////////////////////////////////////////////////////////////
//!
//! class TE_AMV_X: class for template Expressions.
//! this is the CRTP base class for all expressions with
//! alpha Mat Vec + alpha Mat Vec used in this file.
template<class A> class TE_AMV_X
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
	const TRANSPOSED<TE_AMV_X<A> > T() const { return TRANSPOSED<TE_AMV_X<A> > (*this); }
};


//! CRTP base class representing scaled and non-scaled vectors
template<class A> class TE_SCALED_VEC : public TE_AMV_X<A>
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
};

//! CRTP base class representing scaled and non-scaled matrices
template<class A> class TE_SCALED_MAT
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
};

////////////////////////////////////////////////////////////////////////////////

//! CRTP base class representing a normal Matrix class
template<class A> class TE_MAT : public TE_SCALED_MAT<A>
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
};

//! CRTP base class representing a matrix class which supports mat_mult_add_row
template<class A> class TE_MAT_RowApplicable : public TE_MAT<A>
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
};

//! CRTP base class representing a normal vector
template<class A> class TE_VEC : public TE_SCALED_VEC<A>
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
};


////////////////////////////////////////////////////////////////////////////////
//! AlphaVec_Expression
//! class for Template Expressions of the form
//! double * Vector.
//! \attention to case x = alpha1*y + alpha1*x
template<typename R>
class AlphaVec_Expression : public TE_SCALED_VEC<AlphaVec_Expression<R> >
{ 
public:
	double alpha; 
	const R& r; 
	inline AlphaVec_Expression(double alpha_, const R & r_) : alpha(alpha_), r(r_) 
	{ }
};  

// a normal vector is a vector which is scaled by 1.0.

template<typename T>
double getScaling(const T &t)
{
	return 1.0;
}

template<typename T>
const T &getVector(const TE_SCALED_VEC<T> &t)
{
	return t.cast();
}

template<typename T>
double getScaling(const TE_SCALED_VEC<AlphaVec_Expression<T> > &t)
{
	return t.cast().alpha;
}

template<typename T>
const T &getVector(const TE_SCALED_VEC< AlphaVec_Expression<T> > &t)
{
	return t.cast().r;
}

////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//! MatVec_Expression
//! class for Template Expressions of the form
//! Matrix * Vector.
//! todo \attention x = A*x, this would be GaussSeidel
template<typename L, typename R>
class MatVec_Expression : public TE_AMV_X<MatVec_Expression<L, R> >
{ 
public:
	double alpha;
	const L& l; 
	const R& r; 
	inline MatVec_Expression(double alpha_, const L & l_, const R & r_) : alpha(alpha_), l(l_),r(r_)
	{  }
}; 

////////////////////////////////////////////////////////////////////////////////
struct operation_add
{
	static inline bool is_add() { return true; }
};

struct operation_sub
{
	static bool inline is_add() { return false; }
};

//! AlphaMatVec_X_Expression
//! class for nested Template Expression
//! use this class to do TE_AMV + TE_AMV or TE_AMV - TE_AMV,
//! like A*x - b, or (x-0.3*b)-dAlpha*c.
template<typename L, typename operation, typename R>
class AlphaMatVec_X_Expression : public TE_AMV_X<AlphaMatVec_X_Expression<L, operation, R> >
{ 
public:
	const L& l; 
	const R& r; 
	inline AlphaMatVec_X_Expression(const L& l_, const R& r_) : l(l_),r(r_) {}
	//{ UG_ASSERT(l.size() == r.size(), l << " has different length as " <<  r); }

	static inline bool is_add() { return operation::is_add(); }
}; 

////////////////////////////////////////////////////////////////////////////////
//! AlphaMat_Expression
//! use this class to do double * TE_MAT
template<typename R>
class AlphaMat_Expression
{
public:
	double alpha;
	const R& r;
	inline AlphaMat_Expression(double d, const R & r_) : alpha(d), r(r_)
	{  }
};

//#pragma mark operator x -> X_Operator
////////////////////////////////////////////////////////////////////////////////
// operators


//! create AlphaMatVec_X_Expression<L, operation_add, R> by conjunction of TE_AMV_X<L> + TE_AMV_X<R>
template<typename L, typename R>
AlphaMatVec_X_Expression<L, operation_add, R> operator + (const TE_AMV_X<L> &l, const TE_AMV_X<R> &r)
{ 
	return AlphaMatVec_X_Expression<L, operation_add, R> (l.cast(), r.cast());
}

//! create AlphaMatVec_X_Expression<L, operation_minus, R> by conjunction of TE_AMV_X<L> + TE_AMV_X<R>
template<typename L, typename R>
AlphaMatVec_X_Expression<L, operation_sub, R> operator - (const TE_AMV_X<L> &l, const TE_AMV_X<R> &r)
{ 
	return AlphaMatVec_X_Expression<L, operation_sub, R> (l.cast(), r.cast());
}

//! create a MatVec_Expression by TE_MAT * TE_VEC
template<typename L, typename R>
MatVec_Expression <L, R>  operator * (const TE_MAT<L> &l, const TE_SCALED_VEC<R> &r)
{ 
	return MatVec_Expression<L, R> (getScaling(r.cast()), l.cast(), getVector(r.cast()));
}

//! create a MatVec_Expression by (alpha*TE_MAT) * TE_VEC
template<typename L, typename R>
MatVec_Expression <L, R>  operator * (const AlphaMat_Expression<L> &l, const TE_SCALED_VEC<R> &r)
{
	return MatVec_Expression<L, R> (l.alpha*getScaling(r.cast()), l.r, getVector(r.cast()));
}


//! create a AlphaVec_Expression by double * TE_VEC
template<typename R>
AlphaVec_Expression <R>  operator * (double d, const TE_VEC<R> &r)
{ 
	return AlphaVec_Expression<R> (d, r.cast()); 
}

//! create a AlphaMat_Expression by double * TE_MAT
template<typename R>
AlphaMat_Expression <R > operator * (double d, const TE_MAT<R> &r)
{ 
	return AlphaMat_Expression <R> (d, r.cast());
}


} // namespace ug

#include "operations.h"

#endif
