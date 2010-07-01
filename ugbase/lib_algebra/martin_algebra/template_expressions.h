/*
 *  NewTemplateExpressions.h
 *  flexamg
 *
 *  Created by Martin Rupp on 08.02.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__MARTIN_ALGEBRA__TEMPLATE_EXPRESSIONS__
#define __H__UG__MARTIN_ALGEBRA__TEMPLATE_EXPRESSIONS__

//#include "blockMatrix.h"

// misc stuff remove this!!
//!
//! template struct for sorting some keys after values
//! for example, sorting a vector of ints and know original pos
template<typename T>
struct sortStruct
{
	size_t index; // for example "original" position.
	T sortValue;

	bool operator < (const sortStruct<T> &other) const
	{
		return sortValue < other.sortValue;
	}
};

using namespace std;

namespace ug{

//#pragma mark template expression CRTP base classes
//! only classes which inherit from TE_AMV_X, TE_MAT or TE_VEC (via class myClass : public XD<myClass> )
//! can use templateExpressions used in this file
//! CRTP = Curiously recurring template pattern

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
//! need functions assign, addTo, substractFrom, prevent, size and <<.
template<class A> class TE_AMV_X
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
	const TRANSPOSED<TE_AMV_X<A> > T() const { return TRANSPOSED<TE_AMV_X<A> > (*this); }
};

////////////////////////////////////////////////////////////////////////////////
//! class TE_MAT: class for template Expressions.
//! this is the CRTP base class for all Sparse Matrices used by
//! Template Expressions in this file.
//! need functions assign, addTo, substractFrom, prevent, size and <<.
template<class A> class TE_MAT
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
};

////////////////////////////////////////////////////////////////////////////////
//! class TE_MAT: class for template Expressions.
//! this is the CRTP base class for all Vectors used by
//! Template Expressions in this file.
//! need functions assign, addTo, substractFrom, prevent, size and <<.
template<class A> class TE_VEC : public TE_AMV_X<A>
{
public:
	//! cast this class down to original class A.
	const A& cast() const {return static_cast<const A&>(*this); }
};






////////////////////////////////////////////////////////////////////////////////
//! AlphaVec_Expression
//! class for Template Expressions of the form
//! double * Vector.
//! \attention you cannot use x = d*y + e*x. instead use x = e*x + d*a.
//! this is because of the mechanism to prevent local variables.
template<typename R>
class AlphaVec_Expression : public TE_AMV_X<AlphaVec_Expression<R> >
{ 
public:
	typedef typename R::entry_type entry_type;

	double alpha; 
	const R& r; 
	inline AlphaVec_Expression(double alpha_, const R & r_) : alpha(alpha_), r(r_) 
	{ } 
	
	//! calcs d = expression[i]
	inline void assign(entry_type &d, size_t i) const
	{
		AssignMult(d, alpha, r[i]);
	}
	
	//! calcs d += expression[i]
	inline void addTo(entry_type &d, size_t i) const
	{
		AddMult(d, alpha, r[i]);
	}
	
	//! calcs d -= expression[i]
	inline void substractFrom(entry_type &d, size_t i) const
	{
		SubMult(d, alpha, r[i]);
	}
	
	//! as only argument on right side, this is always ok.
	void preventForbiddenDestination(void *p) const
	{
	}
	
	//! x = a+x is forbidden, x = x + etc  ok.
	void preventForbiddenDestination(void *p, bool &bFirst) const
	{
		UG_ASSERT(bFirst == true || p != &r, r << " is on left and right side of template expression. only possible as first summand on right side.");
		bFirst = false;
	}	
	
	inline size_t size() const	{	return r.size();	}
	
// print functions
	friend ostream &operator<<(ostream &out, const AlphaVec_Expression<R>  &ex)
	{
		out << "( double " << ex.alpha << " * " << ex.r << " )";
		return out;
	}
	inline void printtype() const	{	cout << *this; }
};  


// L = TE_MAT<L>, R = TE_VEC<R>

////////////////////////////////////////////////////////////////////////////////
//! MatVec_Expression
//! class for Template Expressions of the form
//! Matrix * Vector.
//! \attention x = A*x, this would be GaussSeidel (use prevent to check).
template<typename L, typename R>
class MatVec_Expression : public  TE_AMV_X<MatVec_Expression<L, R> >
{ 
#define LEFT_RIGHT_CHECK(a) UG_ASSERT(a, r << " is on left and right side of template expression, not possible for matrix multiplication (otherwise GS).");
public:
	typedef typename R::entry_type entry_type;
	const L& l; 
	const R& r; 
	inline MatVec_Expression(const L & l_, const R & r_) : l(l_),r(r_) 
	{ UG_ASSERT(l.num_cols() == r.size(), l << " has different length as " <<  r); }
	
	//! calcs d = expression[i]
	inline void assign(entry_type &d, size_t i) const
	{
		//LEFT_RIGHT_CHECK(&d != &r[i]);
		l[i].assign_mult(d, r);
	}
	
	//! calcs d += expression[i]
	inline void addTo(entry_type &d, size_t i) const
	{
		//LEFT_RIGHT_CHECK(&d != &r[i]);
		l[i].add_mult(d, r);
	}
	
	//! calcs d -= expression[i]
	inline void substractFrom(entry_type &d, size_t i) const
	{
		//LEFT_RIGHT_CHECK(&d != &r[i]);
		l[i].sub_mult(d, r);
	}	
	
	//! ASSERTs that r is not on left and right side (as r = A*r)
	void preventForbiddenDestination(void *p) const
	{
		
		bool bFirst = true;
		preventForbiddenDestination(p, bFirst);
	}

	//! ASSERTs that r is not on left and right side (as r = A*r)
	void preventForbiddenDestination(void *p, bool &bFirst) const
	{
		LEFT_RIGHT_CHECK(p != &r);
		bFirst = false;		
	}	
	
	inline size_t size() const	{	return l.num_rows();	}
	
// print functions
	friend ostream &operator<<(ostream &out, const MatVec_Expression<L, R>  &ex)
	{
		out << "( " << ex.l << " * " << ex.r << " )";
		return out;
	}
	inline void printtype() const	{	cout << *this; }	
}; 

/*
template<typename L, typename R>
class AlphaMatVec_Expression
{ 
public:
	double ld;
	const L& l; 
	const R& r; 
	inline Expression(double ld_, const L& l_, const R& r_) : ld(ld_), l(l_), r(r_) 
	{ UG_ASSERT(l.size() == r.size(), l << " has different length as " <<  r); }
	
	
	inline void assign(vector_type &d, size_t i) const
	{			
		l[i].assign_mult_scale(d, r);
	}
	
	inline void addTo(vector_type &d, size_t i) const
	{
		l[i].add_mult_scale(d, r, ld);
	}
	
	inline void substractFrom(vector_type &d, size_t i) const
	{
		l[i].sub_mult_scale(d, r, ld);
	}	
	
	inline size_t size() const	{	return l.size();	}
}; */

////////////////////////////////////////////////////////////////////////////////
//! AlphaMatVec_Add_Expression
//! class for nested Template Expressions of 
//! MatVec_Expression, AlphaVec_Expression and Vectors. ADD Version
template<typename L, typename R> 
class AlphaMatVec_Add_Expression : public TE_AMV_X<AlphaMatVec_Add_Expression<L, R> >
{ 
public:
	typedef typename L::entry_type entry_type;
	const L& l; 
	const R& r; 
	inline AlphaMatVec_Add_Expression(const L& l_, const R& r_) : l(l_),r(r_) 
	{ UG_ASSERT(l.size() == r.size(), l << " has different length as " <<  r); }
	
	//! calcs d = expression[i]
	inline void assign(entry_type &d, size_t i) const
	{
		l.assign(d, i);
		r.addTo(d, i);
	}
	
	//! calcs d += expression[i]
	inline void addTo(entry_type &d, size_t i) const
	{
		l.addTo(d, i);
		r.addTo(d, i);
	}
	
	//! calcs d -= expression[i]
	inline void substractFrom(entry_type &d, size_t i) const
	{
		l.substractFrom(d, i);
		r.substractFrom(d, i);
	}
	
	//! throws an assert if one argument of expression is forbidden
	void preventForbiddenDestination(void *p) const
	{
		// perhaps this could be done by static asserts (boost)
		bool bFirst = true;
		preventForbiddenDestination(p, bFirst);
	}
	
	//! throws an assert if one argument of expression is forbidden
	void preventForbiddenDestination(void *p, bool &bFirst) const
	{
		// delegate check to arguments
		l.preventForbiddenDestination(p, bFirst);
		r.preventForbiddenDestination(p, bFirst);
	}
	
	inline size_t size() const	{	return l.size();	}
	
// print functions
	friend ostream &operator<<(ostream &out, const AlphaMatVec_Add_Expression<L, R>   &ex)
	{
		out << "( " << ex.l << " + " << ex.r << " )";
		return out;
	}
	inline void printtype() const	{	cout << *this; }	
}; 

////////////////////////////////////////////////////////////////////////////////
//! AlphaMatVec_Sub_Expression
//! class for nested Template Expressions of 
//! MatVec_Expression, AlphaVec_Expression and Vectors. SUBSTRACT Version
template<typename L, typename R> 
class AlphaMatVec_Sub_Expression : public TE_AMV_X<AlphaMatVec_Sub_Expression<L, R> >
{ 
public:
	typedef typename L::entry_type entry_type;
	const L& l; 
	const R& r; 
	inline AlphaMatVec_Sub_Expression(const L& l_, const R& r_) : l(l_),r(r_) 
	{ UG_ASSERT(l.size() == r.size(), l << " has different length as " <<  r); }
	
	//! calcs d = expression[i]
	inline void assign(entry_type &d, size_t i) const
	{
		l.assign(d, i);
		r.substractFrom(d, i);
	}
	
	//! calcs d += expression[i]
	inline void addTo(entry_type &d, size_t i) const
	{
		l.substractFrom(d, i);
		r.substractFrom(d, i);
	}
	
	//! calcs d -= expression[i]
	inline void substractFrom(entry_type &d, size_t i) const
	{
		l.substractFrom(d, i);
		r.addTo(d, i);
	}	
	
	//! throws an assert if one argument of expression is forbidden
	void preventForbiddenDestination(void *p) const
	{
		bool bFirst = true;
		preventForbiddenDestination(p, bFirst);
	}
	
	//! throws an assert if one argument of expression is forbidden
	void preventForbiddenDestination(void *p, bool &bFirst) const
	{
		// delegate check to arguments
		l.preventForbiddenDestination(p, bFirst);
		r.preventForbiddenDestination(p, bFirst);
	}

	inline size_t size() const	{	return l.size();	}

// print functions
	friend ostream &operator<<(ostream &out, const AlphaMatVec_Sub_Expression<L, R>   &ex)
	{
		out << "( " << ex.l << " - " << ex.r << " )";
		return out;
	}
	inline void printtype() const	{	cout << *this; }	
}; 


//#pragma mark operator x -> X_Operator
////////////////////////////////////////////////////////////////////////////////

//! create AlphaMatVec_Add_Expression by conjunction of two TE_AMV_X
template<typename L, typename R>
AlphaMatVec_Add_Expression<L, R> operator + (const TE_AMV_X<L> &l, const TE_AMV_X<R> &r)
{ 
	return AlphaMatVec_Add_Expression<L, R> (l.cast(), r.cast()); 
}

//! create AlphaMatVec_Sub_Expression by conjunction of two TE_AMV_X
template<typename L, typename R>
AlphaMatVec_Sub_Expression<L, R> operator - (const TE_AMV_X<L> &l, const TE_AMV_X<R> &r)
{ 
	return AlphaMatVec_Sub_Expression<L, R> (l.cast(), r.cast()); 
}

//! create a MatVec_Expression by TE_MAT * TE_VEC
template<typename L, typename R>
MatVec_Expression <L, R>  operator * (const TE_MAT<L> &l, const TE_VEC<R> &r)
{ 
	return MatVec_Expression<L, R> (l.cast(), r.cast()); 
}

//! create a AlphaVec_Expression by double * TE_VEC
template<typename R>
AlphaVec_Expression <R>  operator * (double d, const TE_VEC<R> &r)
{ 
	return AlphaVec_Expression<R> (d, r.cast()); 
}

/*
//! create a AlphaMatVec_Expression by double * MatVec_Expression
//! problem here is add_mult_scale and so on could be done by ONE temporal variable. 
template<typename L, typename R>
AlphaMatVec_Expression <L, R>  operator * (double d, const MatVec_Expression <L, R> &v)
{ 
	return AlphaMatVec_Expression<L, R> (d, l, r); 
}*/

// * and / only for special


//#pragma mark Templated Functions
////////////////////////////////////////////////////////////////////////////////

double mnorm2(double &a)
{
	return a*a;
}

//! template expression norm2
template<typename X> 
inline double norm2(const TE_AMV_X<X> &ex_)
{
	const X &ex = ex_.cast();
	double sum=0;
	typename X::entry_type t;
	for(size_t i=0; i < ex.size(); i++)
	{
		ex.assign(t, i);
		sum += mnorm2(t);
	}
	return sum;
}

//! template expression norm
template<typename X> 
inline double norm(const TE_AMV_X<X> &t)
{
	return sqrt(norm2<X>(t));
}

template<typename X> 
inline double absmax2(const TE_AMV_X<X> &ex_)
{
	const X &ex = ex_.cast();
	double absmax=-1;
	typename X::entry_type t;
	for(size_t i=0; i < ex.size(); i++)
	{
		ex.assign(t, i);
		double d= mnorm2(t);
		if(d > absmax)
			absmax = d;
	}
	return absmax;
}

template<typename X> 
inline double absmax(const TE_AMV_X<X> &ex_)
{
	return sqrt( absmax2<X> (ex_) );
}



//! template expression dotProd between two TE_AMV_Xs.
template<typename L, typename R> 
inline double dotProd(const TE_AMV_X<L> &l_, const TE_AMV_X<R> &r_)
{
	const L &l = l_.cast();
	const R &r = r_.cast();
	UG_ASSERT(l.size() == r.size(), "expression L=" << l << " and R=" << r << " differ in length.");
	
	double sum=0;
	typename L::entry_type block_l;
	typename R::entry_type block_r;
	
	size_t N = l.size();
	for(size_t i=0; i < N; i++)
	{
		l.assign(block_l, i);
		r.assign(block_r, i);
		sum += block_l * block_r;
	}
	return sum;
}


template<typename L, typename R> 
inline double operator *(const TRANSPOSED<TE_AMV_X<L> > &l, const TE_AMV_X<R> &r)
{
	return dotProd<L, R> (l.T() ,r);
}

} // namespace ug

#endif
