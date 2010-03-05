/*
 *  CoarseSolver.h
 *  flexamg
 *
 *  Created by Martin Rupp on 26.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#pragma once

#include "sparseMatrix.h"
//typedef long int __CLPK_integer;

// TODO: for smallmatrix tasks like this, better use sth. like FLENS
namespace ug{
	
class CoarseSolver
{
public:
	CoarseSolver()
	{
		densemat = NULL;
		interchange = NULL;
	}
	
	~CoarseSolver()
	{	
		if(densemat) delete[] densemat;
		if(interchange) delete[] interchange;
		if(vec) delete[] vec;		
	}
	template <typename entry_type>
	void create(const SparseMatrix<entry_type> &A);
	
	template<typename vec_type>
	void solve(const Vector<vec_type> &b, Vector<vec_type> &x);
	
	int size;
	double *densemat;
	double *vec;
	__CLPK_integer *interchange;
};
	
} // namespace ug
#include "CoarseSolver.hpp"