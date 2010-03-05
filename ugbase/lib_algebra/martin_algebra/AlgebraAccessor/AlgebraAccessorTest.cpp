#if 0
/*
 *  AlgebraAccessorTest.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 24.02.10.
 *  Copyright 2010 G-CSC, University of Frankfurt. All rights reserved.
 *
 */


#include "algebra.h"

#include "AlgebraAccessor.h"
int flexamg_dimensions;

extern std::vector<postype> positions;

typedef MultiIndex<4> index4;

void disc2d(MatrixAccessorBase<index4> *A, VectorAccessorBase<index4> *x, VectorAccessorBase<index4> *b, int NX, int NY, int ALPHA);
void disc1d(MatrixAccessorBase<index4> *A, VectorAccessorBase<index4> *x, VectorAccessorBase<index4> *b, int N, int ALPHA);

//!
//! CG-Solver
//! CG Solver on system Ax = b with Preconditioner P.
template<typename index_type>
void CG(VectorAccessorBase<index_type> *x, MatrixAccessorBase<index_type> *A, const VectorAccessorBase<index_type> *b, int maxit)
{
	VectorAccessorBase<index_type> *r = x->newClone();
	*r = v(b);
	*r -= mult(A, x);

	
	VectorAccessorBase<index_type> *d = x->newClone();

	*d = v(r); //P.precond(d, r);
	
	VectorAccessorBase<index_type> *z = x->newClone();
	*z = v(d);
	
	VectorAccessorBase<index_type> *t = x->newClone();
	
	int i=0;
	double alpha, rz_new, rz_old;
	
	rz_old = scal_prod(r, z);
	
	double res, oldres;
	
	res = r->norm();
	cout << "[0] res: " << res << endl;
	
	for(i=1; i <= maxit && (res > 1e-10); i++)
	{
		
		*t = mult(A,d);
		alpha = rz_old/scal_prod(d, t);
		*x += alpha * v(d);
		*r -= alpha * v(t);
		
		*z = v(r); // P.precond(z, r);
		
		rz_new = scal_prod(r, z);
		
		*d = (rz_new/rz_old) *v(d) + v(z);
		
		rz_old = rz_new;
		
		oldres = res;
		res = r->norm();
		//if(i %10 == 0)
		cout << "[" << i << "] res: " << res << " conv.: " << res/oldres << endl;		
		cout.flush();
	}
	//if(--i %10 != 0) 
	cout << "[" << i << "] res: " << res << " conv.: " << res/oldres << endl;
	cout << i << " Iterations." << endl;
}


int main()
{
	const int ALPHA = 2;


	AlgebraAccessorBase<index4> *alg;
	if(0)
		alg = new AlgebraAccessor_UnknownWise<index4>;
	else
		alg = new AlgebraAccessor_PointBlock<index4, ALPHA>;
	
	const int dimensions =2;
	int NX=10;
	int NY=10;
	int N;
	if(dimensions==1)
		N = NX;
	else
		N = NX*NY;
	positions.resize(N);
	
	IndexInfo tree;
	
	
	IndexInfo indexBranch;
	indexBranch.set_num_index(N);
	indexBranch.set_num_comp(ALPHA);
	
	IndexInfo dofGroupBranch;
	dofGroupBranch.set_num_index(1);
	dofGroupBranch.get_index_info(0) = indexBranch;
	
	IndexInfo &subdomainBranch = tree;	
	subdomainBranch.set_num_index(1); // subdomain
	subdomainBranch.get_index_info(0) = dofGroupBranch;
	
	MatrixAccessorBase<index4> *A = alg->newMatrix();
	A->Create(tree);
	
	VectorAccessorBase<index4> *b = alg->newVector();	
	b->Create(tree);
	VectorAccessorBase<index4> *x = b->newClone();

	if(dimensions==1)
		disc1d(A, x, b, N, ALPHA);
	else
		disc2d(A, x, b, NX, NY, ALPHA);

	/*A->print();
	x->print();
	b->print();*/
	
	
	CG<index4> (x, A, b, 1000);
	//x->print();

	
	return 3;
	
}
#endif