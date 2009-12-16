/*
 *  howtouse.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 07.09.09.
 *  Copyright 2009 G-CSC. All rights reserved.
 *
 */
#include <iostream>
using namespace std;
#include <math.h>


#include "misc.h"

const double sigma = 0.3;
const double theta = 0.3;

#include "amg.h"
#include "ls.h"



#define UNKNOWN_NR 3

//#define EASY_MATRIX // set hx2 = 1.0


#if UNKNOWN_NR > 1
typedef fixedMatrix<UNKNOWN_NR> MAT_TYPE;
typedef fixedVector<UNKNOWN_NR> VEC_TYPE;
#else
typedef double VEC_TYPE;
typedef double MAT_TYPE;
#endif



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

////////// main
////////////////////////////////////////////////////////////////
int main(int argc, char **argv) 
{	
	int NX = 257;
	int NY = 513; // 9
	int n = NX*NY;
	iNrOfPositions=n;
	double hX = 1/((double)(NX-1));
	double hY = 1/((double)(NY-1));
#ifndef EASY_MATRIX
	double hx2 = hX*hX, hy2 = hY*hY;
#else
	double hx2 = 1.0; double hy2 = hY*hY/(hX*hX);
#endif
	double h2 = hX*hY;
	double res;	
	
	typedef matrix<MAT_TYPE> myMatrix;
	typedef Vector<VEC_TYPE> myVector;	
	
	
	cout << "grid " << NX << " x " << NY << ", " << UNKNOWN_NR << " unknowns." << endl << endl;
	
	// for debug output: save positions.
	positions = new pos2d[n];
	for(int y=0; y<NY; y++)
		for(int x=0; x<NX; x++)
		{
			int i = x + y*NX;
			positions[i].x = x*hX;
			positions[i].y = y*hY;
		}

	
	// create matrix
	//----------------
	
	myMatrix A("A");
	A.create(n);	
	// first do finite differences in x-direction
	for(int y=0; y<NY; y++)
	{
		for(int x=0; x<NX; x++)
		{	
			int i = x + y*NX;
			if(x == 0 || x == NX-1 || y == 0 || y == NY-1)
			{
				myMatrix::connection con; con.iIndex = i; con.dValue = 1.0;
				A[i].setMatrixRow(&con, 1);
			}
			else
			{
				int indices[3] = {i-1, i, i+1};
				submatrix<MAT_TYPE> UM(indices, 3);
				UM(1, 0) = -1.0 / hx2;
				UM(1, 1) = 2.0 / hx2;
				UM(1, 2) = -1.0 / hx2;
				A.addSubmatrix(UM);
			}
		}
	}
	if(n < 50)	A.print();
	// then in do finite differences in y-direction
	for(int y=0; y<NY; y++)
	{
		for(int x=0; x<NX; x++)
		{	
			if(x == 0 || x == NX-1 || y == 0 || y == NY-1)
				continue;
			else
			{
				int indices[3] = {x+(y-1)*NX, x+y*NX, x+(y+1)*NX};
				submatrix<MAT_TYPE> UM(indices, 3);
				UM(1, 0) = -1.0 / hy2;
				UM(1, 1) = 2.0 / hy2;
				UM(1, 2) = -1.0 / hy2;
				A.addSubmatrix(UM);
			}
		}
	} 	
	
	// create rhs
	//-------------
	myVector b(n, "b");
	b = 0.0;
	for(int i=0; i<n; i++)
	{
		pos2d pos = GetPosForIndex(i);
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

