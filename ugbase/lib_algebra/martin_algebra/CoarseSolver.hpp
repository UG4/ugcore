/*
 *  CoarseSolver.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 26.11.09.
 *  Copyright 2009 . All rights reserved.
 *
 */
/*#include <iostream>

#include <veclib/cblas.h>
#include <veclib/clapack.h>

#include "misc.h"
#include "matrix.h"

#include "CoarseSolver.h" */

#include <veclib/cblas.h>
#include <veclib/clapack.h>

using namespace std;

CoarseSolver::~CoarseSolver()
{
	if(densemat) delete[] densemat;
	if(interchange) delete[] interchange;
	if(vec) delete[] vec;
}

template<typename vec_type>
void CoarseSolver::solve(const Vector<vec_type> &b, Vector<vec_type> &x)
{
	// TODO: Variing nr of unknowns
	int nrOfUnknowns = vec_traits<vec_type>::nrOfUnknowns;
	ASSERT2(size == b.getLength() * nrOfUnknowns && size == x.getLength() * nrOfUnknowns, " wrong size! has to be " << size << ", but is " << b << " and " << x);
	
	for(int i=0; i < b.getLength(); i++)
		for(int j=0; j<nrOfUnknowns; j++)
			vec[i*nrOfUnknowns+j] = getAt(b[i], j);
	
	// solve system
	char trans ='N';
	__CLPK_integer dim = size;
	__CLPK_integer nrhs = 1;
	__CLPK_integer info;
	
	dgetrs_(&trans, &dim, &nrhs, densemat, &dim, interchange, vec, &dim, &info);
	ASSERT2(info == 0, "info is " << info);
	for(int i=0; i<x.getLength(); i++)
		for(int j=0; j<nrOfUnknowns; j++)
			setAt(x[i], j, vec[i*nrOfUnknowns+j]);

}


template<typename mat_type>
void CoarseSolver::create(const matrix<mat_type> &A)
{
	int nrOfUnknowns = matrix_trait<mat_type>::nrOfUnknowns;
	size = (__CLPK_integer) A.getLength() * nrOfUnknowns;
	
	densemat = new double[size*size];
	interchange = new __CLPK_integer[size];
	
	memset(densemat, 0, sizeof(double)*size*size);
	
	for(int r=0; r<A.getLength(); r++)
		for(typename matrixrow<mat_type>::citerator it(A[r]); !it.isEnd(); ++it)
		{
			int rr = r*nrOfUnknowns;
			int cc = (*it).iIndex*nrOfUnknowns;
			for(int r2=0; r2<nrOfUnknowns; r2++)
					for(int c2=0; c2<nrOfUnknowns; c2++)
						densemat[(rr+r2) + (cc+c2)*size] = getAt((*it).dValue, r2, c2);
		}


	vec = new double[size];	
	
	
/*
 cout << endl << endl;
	 for(int j=0; j<A.getLength(); j++)
	 {
	 for(int i=0; i<A.getLength(); i++)		
	 cout << cut(densemat[j + i*size], 1e-7) << "\t";
	 cout << endl;
	 }
	//*/
	
	__CLPK_integer info = 0;
	__CLPK_integer dim = size;
	dgetrf_(&dim, &dim, densemat, &dim, interchange, &info);
	ASSERT2(info == 0, "info is " << info << ( info > 0 ? ": matrix singular in U(i,i)" : ": i-th argument had had illegal value"));
}