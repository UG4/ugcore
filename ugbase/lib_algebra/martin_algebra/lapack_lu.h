/*
 *  LapackLU.h
 *  flexamg
 *
 *  Created by Martin Rupp on 26.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#ifndef __H__UG__MARTIN_ALGEBRA__LAPACK_LU__
#define __H__UG__MARTIN_ALGEBRA__LAPACK_LU__


#include "sparsematrix.h"
#include "vector.h"

#include <cblas.h>
#include <clapack.h>


// TODO: for smallmatrix tasks like this, better use sth. like FLENS
namespace ug{

//template<typename T>
class LapackLU
{
	//typedef T entry_type;
public:
	LapackLU()
	{
		densemat = NULL;
		interchange = NULL;
	}

	~LapackLU()
	{
		if(densemat) delete[] densemat;
		if(interchange) delete[] interchange;
	}

	// todo: for non-double
	void init(const SparseMatrix<double> &A);

	void prepare(const Vector<double> &b, Vector<double> &x) {}
	void apply(const Vector<double> &b, Vector<double> &x);

private:
	size_t size;
	double *densemat;
	int *interchange;
};

} // namespace ug
#include "lapack_lu_impl.h"

#endif /* __H__UG__MARTIN_ALGEBRA__LAPACK_LU__ */
