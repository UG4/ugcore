/*
 *  AATemplateExpressions.h
 *  flexamg
 *
 *  Created by Martin Rupp on 02.03.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#pragma once

template<typename T> class VectorAccessorBase;
template<typename T> class MatrixAccessorBase;

#pragma mark template expression Barton-Nackman base classes

////////////////////////////////////////////////////////////////////////////////
//! AlphaVec_Expression
//! class for Template Expressions of the form
//! double * Vector.
//! \attention you cannot use x = d*y + e*x. instead use x = e*x + d*a.
//! this is because of the mechanism to prevent local variables.
template<typename T>
class AA_AlphaVec_Expression
{
public:
	double l;
	const VectorAccessorBase<T> *r;
	inline AA_AlphaVec_Expression(double alpha_, const VectorAccessorBase<T> *r_) : l(alpha_), r(r_) 
	{ } 
};  

template<typename T>
class AA_MatVec_Expression
{ 
public:
	const MatrixAccessorBase<T> *l; 
	const VectorAccessorBase<T> *r; 
	inline AA_MatVec_Expression(const MatrixAccessorBase<T> *l_, const VectorAccessorBase<T> *r_) : l(l_),r(r_) {}
}; 


////////////////////////////////////////////////////////////////////////////////
//! AlphaMatVec_Add_Expression
//! class for nested Template Expressions of 
//! MatVec_Expression, AlphaVec_Expression and Vectors. ADD Version
template<typename L, typename R> 
class AA_AlphaMatVec_Add_Expression
{ 
public:
	const L& l; 
	const R& r; 
	inline AA_AlphaMatVec_Add_Expression(const L& l_, const R& r_) : l(l_),r(r_) 
	{ } 	
	
}; 

////////////////////////////////////////////////////////////////////////////////
//! AlphaMatVec_Sub_Expression
//! class for nested Template Expressions of 
//! MatVec_Expression, AlphaVec_Expression and Vectors. SUBSTRACT Version
template<typename L, typename R> 
class AA_AlphaMatVec_Sub_Expression
{ 
public:
	typedef typename L::entry_type entry_type;
	const L& l; 
	const R& r; 
	inline AA_AlphaMatVec_Sub_Expression(const L& l_, const R& r_) : l(l_),r(r_) 
	{  } 		
}; 


#pragma mark operator x -> X_Operator
////////////////////////////////////////////////////////////////////////////////


template<typename T>
AA_AlphaMatVec_Add_Expression<AA_MatVec_Expression<T>, AA_MatVec_Expression<T> > operator + (const AA_MatVec_Expression<T> &l, const AA_MatVec_Expression<T> &r)
{ 
	return AA_AlphaMatVec_Add_Expression<AA_MatVec_Expression<T>, AA_MatVec_Expression<T> > (l, r); 
}

template<typename T>
AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > operator + (const AA_MatVec_Expression<T> &l, const AA_AlphaVec_Expression<T> &r)
{ 
	return AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > (r, l); 
}

template<typename T>
AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > operator + (const AA_AlphaVec_Expression<T> &l, const AA_MatVec_Expression<T> &r)
{ 
	return AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > (l, r); 
}


// avec avec

template<typename T>
AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > operator + (const AA_AlphaVec_Expression<T> &l, const AA_AlphaVec_Expression<T> &r)
{ 
	return AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > (l, r); 
}


template<typename T>
AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > operator - (const AA_AlphaVec_Expression<T> &l, const AA_AlphaVec_Expression<T> &r)
{ 
	return AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > (l, AA_AlphaVec_Expression<T>(-r.l, r.r)); 
}

/////

template<typename T>
AA_AlphaVec_Expression<T> v(const VectorAccessorBase<T> *r)
{
	return AA_AlphaVec_Expression<T> (1.0, r);
}



template<typename T>
AA_MatVec_Expression<T> mult (const MatrixAccessorBase<T> *l, const VectorAccessorBase<T> *r)
{
	return AA_MatVec_Expression<T> (l, r);
}

template<typename T>
AA_AlphaVec_Expression<T> operator * (double d, const AA_AlphaVec_Expression<T> &r)
{ 
	return AA_AlphaVec_Expression<T> (d*r.l, r.r); 
}

