/*
 *  AlgebraAccessor.h
 *  flexamg
 *
 *  Created by Martin Rupp on 23.02.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#pragma once
#include "interface.h"



typedef blockDenseMatrix<double, variableStorage> DenseMatrix;
typedef blockVector<double, variableStorage> DenseVector;

#include "AATemplateExpressions.h"

template<typename T>
class MatrixAccessorBase
{
public:
	typedef T multi_index_type;
	
	MatrixAccessorBase() {}
	virtual ~MatrixAccessorBase() {}
	virtual bool Create(const IndexInfo &i) = 0;
	
	virtual bool add(DenseMatrix &M, vector<T> &I, vector<T> &J) = 0;
	//	virtual bool set(DenseMatrix &M, vector<T> &I, vector<T> &J) = 0;
	//	virtual bool get(DenseMatrix &M, vector<T> &I, vector<T> &J) = 0;
	virtual void print() {return;};
};


template<typename T>
class VectorAccessorBase
{
public:
	typedef T multi_index_type;
	
	VectorAccessorBase() {}
	virtual ~VectorAccessorBase() {}
	virtual bool Create(const IndexInfo &i) = 0;
	virtual VectorAccessorBase<T> *newClone() const = 0;
	
	virtual bool add(DenseVector &v, vector<T> &I) = 0;
	/*	virtual bool set(DenseVector &v, vector<T> &I) = 0;
	 virtual bool get(DenseVector &v, vector<T> &I) = 0;*/
	 

	// = vec
	virtual bool apply(Operation_type operation, const AA_AlphaVec_Expression<T> &ex ) = 0;	
	bool operator = (const AA_AlphaVec_Expression<T> &ex ) { return apply(OPERATION_SET, ex); }
	bool operator += (const AA_AlphaVec_Expression<T> &ex ) { return apply(OPERATION_ADD, ex); }
	bool operator -= (const AA_AlphaVec_Expression<T> &ex ) { return apply(OPERATION_SUB, ex); }

	// = vec + vec
	virtual bool apply(Operation_type operation, const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > &ex) = 0;
	bool operator = (const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > &ex) { return apply(OPERATION_SET, ex); }
	bool operator += (const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > &ex) { return apply(OPERATION_ADD, ex); }
	bool operator -= (const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_AlphaVec_Expression<T> > &ex) { return apply(OPERATION_SUB, ex); }
	
	// = vec + matvec
	virtual bool apply(Operation_type operation, const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > &ex) = 0;
	bool operator = (const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > &ex ) { return apply(OPERATION_SET, ex); }
	bool operator += (const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > &ex ) { return apply(OPERATION_ADD, ex); }
	bool operator -= (const AA_AlphaMatVec_Add_Expression<AA_AlphaVec_Expression<T>, AA_MatVec_Expression<T> > &ex ) { return apply(OPERATION_SUB, ex); }

	// = matvec
	virtual bool apply(Operation_type operation, const AA_MatVec_Expression<T> &ex) = 0;
	bool operator = (const AA_MatVec_Expression<T> &ex) { return apply(OPERATION_SET, ex); }
	bool operator += (const AA_MatVec_Expression<T> &ex ) { return apply(OPERATION_ADD, ex); }
	bool operator -= (const AA_MatVec_Expression<T> &ex ) { return apply(OPERATION_SUB, ex); }
	
	/*bool operator = ( ) { return apply(OPERATION_SET, ex); }
	bool operator += ( ) { return apply(OPERATION_ADD, ex); }
	bool operator -= ( ) { return apply(OPERATION_SUB, ex); }*/
	
	virtual double scal_prod(const VectorAccessorBase<T> *v) const = 0;
	virtual double norm2() const = 0;
	double norm() const { return sqrt(norm2()); }
	
	virtual void print() {return;};
};

template<typename T>
double scal_prod(const VectorAccessorBase<T> *x, const VectorAccessorBase<T> *y)
{
	return x->scal_prod(y);
}
	
template<typename T>
class AlgebraAccessorBase
{
public:
	typedef T multi_index_type;
	
	AlgebraAccessorBase() {}
	virtual ~AlgebraAccessorBase() {}
	
	virtual MatrixAccessorBase<T> *newMatrix() = 0;
	virtual VectorAccessorBase<T> *newVector() = 0;	
};


/////////////////////////////
template<typename T, typename M, typename V>
class templateAlgebraAccessorBase : public AlgebraAccessorBase<T>
{
public:
	typedef T multi_index_type;
	
	templateAlgebraAccessorBase() {}
	virtual ~templateAlgebraAccessorBase() {}
	
	virtual MatrixAccessorBase<T> *newMatrix()
	{
		return new M;
	}
	virtual VectorAccessorBase<T> *newVector()
	{
		return new V;
	}
};

///////////////////////////////////////////////////////////////////////////////

#include "UnknownWise.h"
#include "PointBlock.h"