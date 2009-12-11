/*
 *  CoarseSolver.h
 *  flexamg
 *
 *  Created by Martin Rupp on 26.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */
#pragma once

typedef long int __CLPK_integer;

class CoarseSolver
{
public:
	CoarseSolver()
	{
		densemat = NULL;
		interchange = NULL;
	}
	
	~CoarseSolver();	
	template <typename mat_type>
	void create(const matrix<mat_type> &A);
	
	template<typename vec_type>
	void solve(const Vector<vec_type> &b, Vector<vec_type> &x);
	
	int size;
	double *densemat;
	double *vec;
	__CLPK_integer *interchange;
};
#include "CoarseSolver.hpp"