/*
 *  linearsolver.h
 *  flexamg
 *
 *  Created by Martin Rupp on 11.12.09.
 *  Copyright 2009 . All rights reserved.
 *
 */

#include "sparseMatrix.h"

//!
//! Linear Solver
//! Performs maxit steps of iteration x -> x + P (b-Ax), where P is a preconditioner
// this one uses preconditioner::iterate
template<typename entry_type, typename Vector_type>
void LinearSolver(Vector_type &x, const SparseMatrix<entry_type> &A, 
				  const Vector_type &b, preconditioner<entry_type, Vector_type> &P, int maxit)
{
	P.init(A);
	stopwatch SW;
	SW.start();
	
	
	double res=1, oldres = 0;
	int i;
	
	//	Vector_type // 0.00001
	for(i=0; i< maxit && (res > 1e-10); i++)
	{
		P.iterate(x, b);
		res = norm(b-A*x);		
		cout << "[" << i << "] res: " << res << " conv.: " << res/oldres << endl;
		cout.flush();
		oldres = res;
	}
	cout << i << " Iterations." << endl;
	cout << "res: " << norm(b-A* x) << endl;
	
	SW.printTimeDiff();
}


//!
//! Linear Solver
//! Performs maxit steps of iteration x -> x + P (b-Ax), where P is a preconditioner
// this one uses preconditioner::precond
template<typename entry_type, typename Vector_type>
void LinearSolver2(typename SparseMatrix<entry_type>::Vector_type &x, const SparseMatrix<entry_type> &A, 
				  const Vector_type &b, preconditioner<entry_type, Vector_type> &P, int maxit)
{	
	
	P.init(A);
	stopwatch SW;
	SW.start();
	
	Vector_type r("r"), c("c");
	r.create(A.getLength());
	c.create(A.getLength());
	r = b-A*x;
	double res=1, oldres = norm(r);	
	int i;
	//	Vector_type
	
	for(i=0; i< maxit && (res > 1e-10); i++)
	{
		P.precond(c, r);
		x += c;
		r = b-A*x;
		res = norm(r);
		cout << "[" << i << "] res: " << res << " conv.: " << res/oldres << endl;
		cout.flush();
		oldres = res;
		
	}
	cout << i << " Iterations." << endl;
	cout << "res: " << norm(b-A*x) << endl;
	
	SW.printTimeDiff();
}

//!
//! CG-Solver
//! CG Solver on system Ax = b with Preconditioner P.
template<typename entry_type, typename Vector_type>
void CG(Vector_type &x, const SparseMatrix<entry_type> &A, const Vector_type &b, preconditioner<entry_type, Vector_type> &P, int maxit)
{
	P.init(A);
	
	
	Vector_type r(x.getLength(), "CG:r");
	r = b - A*x;
	
	Vector_type d(x.getLength(), "CG:d");
	
	P.precond(d, r);
	
	Vector_type z(x.getLength(), "CG:z");
	z = d;
	
	Vector_type t(x.getLength(), "CG:t");
	int i=0;
	double alpha, rz_new, rz_old;
	rz_old = r*z;
	double res, oldres;
	
	res = norm(r);
	cout << "[0] res: " << res << endl;
	
	for(i=1; i <= maxit && (res > 0.00001); i++)
	{
		
		t = A*d;
		alpha = rz_old/(d*t);
		x += alpha * d;
		r -= alpha * t;
		
		P.precond(z, r);
		
		rz_new = r*z;
		
		d = z + (rz_new/rz_old) *d;
		
		rz_old = rz_new;
		
		oldres = res;
		res = norm(r);
		//if(i %10 == 0)
		cout << "[" << i << "] res: " << res << " conv.: " << res/oldres << endl;		
		cout.flush();
	}
	//if(--i %10 != 0) 
	cout << "[" << i << "] res: " << res << " conv.: " << res/oldres << endl;
	cout << i << " Iterations." << endl;
}
