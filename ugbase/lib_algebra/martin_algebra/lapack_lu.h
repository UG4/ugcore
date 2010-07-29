/*
 *  Lapack_lu.h
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

class LapackLU
{
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
	template<typename matrix_type>
	void init(const matrix_type &A);

	template<typename vec_type>
	void prepare(const vec_type &b, vec_type &x) {}
	template<typename vec_type>
	void apply(const vec_type &b, vec_type &x);

private:
	size_t size;
	double *densemat;
	int *interchange;
};

} // namespace ug
#include "lapack_lu_impl.h"

#endif /* __H__UG__MARTIN_ALGEBRA__LAPACK_LU__ */
