/*
 *  howtouse.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 07.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */
#if 0

#include "algebra.h"

#if 1
double solution(double x, double y)
{
	return 1;
}

double phi(double x, double y)
{	
	return 1;
}
double f(double x, double y)
{
	return 0;
}
#else
double solution(double x, double y)
{
	return x*x+y*y;
}

double phi(double x, double y)
{	
	return x*x+y*y;
}
double f(double x, double y)
{
	return 2+2;	
}
#endif

// since this is 4-point
#define AGGRESSIVE_COARSENING

extern std::vector<postype> positions;

////////// main
////////////////////////////////////////////////////////////////
int main(int argc, char **argv) 
{	
	int NX = 513;
	int NY = 513; // 9
	int n = NX*NY;
	positions.resize(n);
	double hX = 1/((double)(NX-1));
	double hY = 1/((double)(NY-1));
#ifndef EASY_MATRIX
	double hx2 = hX*hX, hy2 = hY*hY;
#else
	double hx2 = 1.0; double hy2 = hY*hY/(hX*hX);
#endif
	double h2 = hX*hY;
	double res;	
	
	typedef SparseMatrix<MAT_TYPE> myMatrix;
	typedef Vector<VEC_TYPE> myVector;	
	
	
	cout << "grid " << NX << " x " << NY << ", " << UNKNOWN_NR << " unknowns." << endl << endl;
	
	// for debug output: save positions.
	for(int y=0; y<NY; y++)
		for(int x=0; x<NX; x++)
		{
			int i = x + y*NX;
			positions[i].x = x*hX;
			positions[i].y = y*hY;
		}

	
	// create matrix
	//----------------
	
	myMatrix A; 
	A.name = "A";
	A.fromlevel = 0;
	A.tolevel = 0;
	A.create(n, n);	
	// first do finite differences in x-direction
	for(int y=0; y<NY; y++)
	{
		for(int x=0; x<NX; x++)
		{	
			int i = x + y*NX;
			if(x == 0 || x == NX-1 || y == 0 || y == NY-1)
			{
				myMatrix::connection con; 
				con.iIndex = i; 
				setSize(con.dValue, UNKNOWN_NR, UNKNOWN_NR);  // use set size for variable blockvector length
				con.dValue = 1.0;
				A.setMatrixRow(i, &con, 1);
			}
			else
			{
				int indices[3] = {i-1, i, i+1};
				int unknowns[3] = {UNKNOWN_NR, UNKNOWN_NR, UNKNOWN_NR}; 
				
				// get submatrix of indices, specify nr of unknowns of node (needed since if matrix empty informations is not known).
				submatrix<MAT_TYPE> UM(indices, unknowns, 3);
				UM(1, 0) = -1.0 / hx2;				
				UM(1, 1) = 2.0 / hx2;
				UM(1, 2) = -1.0 / hx2;
				A.add(UM);
			}
		}
	}
	if(n < 50)	A.print();
	// then do finite differences in y-direction
	for(int y=0; y<NY; y++)
	{
		for(int x=0; x<NX; x++)
		{	
			if(x == 0 || x == NX-1 || y == 0 || y == NY-1)
				continue;
			else
			{
				int indices[3] = {x+(y-1)*NX, x+y*NX, x+(y+1)*NX};
				int unknowns[3] = {UNKNOWN_NR, UNKNOWN_NR, UNKNOWN_NR};
				submatrix<MAT_TYPE> UM(indices, unknowns, 3);
				UM(1, 0) = -1.0 / hy2;
				UM(1, 1) = 2.0 / hy2;
				UM(1, 2) = -1.0 / hy2;
				A.add(UM);
			}
		}
	} 	
	
	if(n < 50) A.print();
	
	// create rhs
	//-------------
	myVector b(n, "b");
	b = 0.0;
	for(int i=0; i<n; i++)
	{
		postype pos = GetPosForIndex(i);
		setSize(b[i], UNKNOWN_NR); // use set size for variable blockvector length
		if(pos.x == 0 || pos.x == 1.0 || pos.y == 0 || pos.y == 1)
			b[i] = phi(pos.x, pos.y);
		else
			b[i] = -f(pos.x, pos.y);
	}
	
	// eliminate dirichlet values
	A.eliminateDirichletValues(b);
	A.finish();	
	
	// create solution vector
	//-------------------------
	myVector vx(n, "x");
	vx = 0.0;	
	for(int i=0; i<n; i++)
	{
		setSize(vx[i], UNKNOWN_NR); // use set size for variable blockvector length
		if(A.isUnconnected(i))
			vx[i] = b[i]; // dirichlet values
		else
			vx[i] = 0.0;
	}
	
	/*A.print();
	 b.print();
	 vx.print();*/
	
	// solve
	//-----------
	
	stopwatch SW;
	// ANG TESTS
	
	amg<MAT_TYPE> AMG;
	//AMG.max_levels = 20;
#ifdef AGGRESSIVE_COARSENING
	AMG.setAggressiveCoarsening_A_2();
#endif

	cout << " sol: "; 
	LinearSolver(vx, A, b, AMG, 100);
	
	if(n < 10003)
	{
		writeToPosFile("/tmp/grid");
		AMG.writeMatrices("/Users/mrupp/Documents/Xcode/flexamg/matrices/matrix");
		vx.printtofile("/Users/mrupp/Documents/Xcode/flexamg/matrices/matrix/vx.dat"); 
	}

	
		
	cout << endl << endl << "Check: " << endl;
	res = 0;
	// calculate discretization error
	for(int i=0; i<n; i++)
	{			
		double d = solution(GetPosForIndex(i).x, GetPosForIndex(i).y)- getAt(vx[i], 0);
		res += d*d *h2;
	}
	cout << "res:       " << norm(b - A*vx) << endl;
	cout << "disc. err: " << sqrt(res) << endl;	
	
    return 0;
}

#endif