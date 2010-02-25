/*
 *  disc.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 26.01.10.
 *  Copyright 2010 . All rights reserved.
 *
 */
#include "algebra.h"

extern std::vector<postype> positions;


class cProblem
{
public:
	virtual double solution(double x, double y, double z) = 0;
	virtual double phi(double x, double y, double z) = 0;
	virtual double f(double x, double y, double z) = 0;
};

class Laplace : public cProblem
{
public:
	virtual double solution(double x, double y, double z) 
	{
		return 1;
	}
	
	virtual double phi(double x, double y, double z)
	{
		return 1;
	}
	
	virtual double f(double x, double y, double z)
	{
		return 0;
	}	
};

class Problemf6_2d : public cProblem
{
public:
	virtual double solution(double x, double y, double z) 
	{
		return x*x+y*y;
	}
	
	virtual double phi(double x, double y, double z)
	{
		return x*x+y*y;
	}
	
	virtual double f(double x, double y, double z)
	{
		return 2+2;
	}	
};

class Problemf6_3d : public cProblem
{
public:
	virtual double solution(double x, double y, double z) 
	{
		return x*x+y*y+z*z;
	}
	
	virtual double phi(double x, double y, double z)
	{
		return x*x+y*y+z*z;
	}
	
	virtual double f(double x, double y, double z)
	{
		return 2+2+2;
	}	
};

#define GETINDEX2D(x, y) (x + NX*(y))

//!
//! do finite-difference discretisation in structured grid on [0..1]x[0..1]
//! @param A SparseMatrix A for discretisation of the matrix
//! @param vx vector which is set to b[i] for dirichlet nodes, otherwise 0.0
//! @param b vector of the right hand side in Ax=b.
//! @param NX number of nodes in X-direction
//! @param NY number of nodes in Y-direction
void discretization2d(myMatrix &A, myVector &vx, myVector &b, int NX, int NY, int NZ, bool full_star, int problemtype)
{
	cProblem *problem;
	if(problemtype == 0)
		problem = new Laplace;
	else
		problem = new Problemf6_2d;
	
	int n = NX*NY;
	positions.resize(n);
	double hX = 1/((double)(NX-1));
	double hY = 1/((double)(NY-1));
#ifndef EASY_MATRIX
	double hx2 = hX*hX, hy2 = hY*hY;
#else
	double hx2 = 1.0; double hy2 = hY*hY/(hX*hX);
#endif

	stopwatch SW;
	SW.start();
	
	A.create(n, n);
	A.fromlevel = A.tolevel = 0;
	for(int y=0; y<NY; y++)
		for(int x=0; x<NX; x++)
		{	
			int i = GETINDEX2D(x, y);
			positions[i].x = x*hX;
			positions[i].y = y*hY;
			positions[i].z = 0;
			if(x == 0 || x == NX-1 || y == 0 || y == NY-1)
			{
				myMatrix::connection con; con.iIndex = i; 
				setSize(con.dValue, UNKNOWN_NR, UNKNOWN_NR);
				con.dValue = 1.0;
				A.setMatrixRow(i, &con, 1);
			}
			else
			{
				if(full_star)
				{
					// 9-point
					myMatrix::connection con[9];
					for(int k=0; k<9; k++) setSize(con[k].dValue, UNKNOWN_NR, UNKNOWN_NR);					
					double diag=0;
					
					int k=1;
					//x 1 2
					con[k].iIndex = GETINDEX2D(x-1, y);	con[k++].dValue = -1.0 / hx2;
					con[k].iIndex = GETINDEX2D(x+1, y);	con[k++].dValue = -1.0 / hx2;
					diag += 2.0/hx2;
					//y 3 4
					con[k].iIndex = GETINDEX2D(x, y-1);	con[k++].dValue = -1.0 / hy2;
					con[k].iIndex = GETINDEX2D(x, y+1);	con[k++].dValue = -1.0 / hy2;
					diag += 2.0/hy2;
					
					//xy 
					con[k].iIndex = GETINDEX2D(x-1, y+1);	con[k++].dValue = -1.0 / (hx2+hy2);
					con[k].iIndex = GETINDEX2D(x+1, y-1);	con[k++].dValue = -1.0 / (hx2+hy2);
					diag += 2.0/(hx2+hy2);
					con[k].iIndex = GETINDEX2D(x-1, y-1);	con[k++].dValue = -1.0 / (hx2+hy2);
					con[k].iIndex = GETINDEX2D(x+1, y+1);	con[k++].dValue = -1.0 / (hx2+hy2);
					diag += 2.0/(hx2+hy2);
					
					con[0].iIndex = i;
					con[0].dValue = diag;
					A.setMatrixRow(i, con, 9);				
				}
				else
				{
					// 5-point (needs agressive coarsening)
					myMatrix::connection con[5];
					for(int k=0; k<5; k++) setSize(con[k].dValue, UNKNOWN_NR, UNKNOWN_NR);
					con[0].iIndex = i;				con[0].dValue = 2.0 / hx2 + 2.0 / hy2;
					con[1].iIndex = GETINDEX2D(x-1, y);	con[1].dValue = -1.0/hx2;
					con[2].iIndex = GETINDEX2D(x+1, y);	con[2].dValue = -1.0/hx2;
					con[3].iIndex = GETINDEX2D(x, y-1);	con[3].dValue = -1.0/hy2;
					con[4].iIndex = GETINDEX2D(x, y+1);	con[4].dValue = -1.0/hy2;
					A.setMatrixRow(i, con, 5);
				}
			}
		}
	
	cout << "assemble took ";
	SW.printTimeDiff();
	
	
	if(n<10003)
		A.writeToFile("/Users/mrupp/matrices/A.mat");
	
	// create rhs
	b.create(n);
	b = 0.0;
	for(int i=0; i<n; i++)
	{
		postype pos = GetPosForIndex(i);
		setSize(b[i], UNKNOWN_NR);
		if(pos.x == 0 || pos.x == 1.0 || pos.y == 0 || pos.y == 1 || ((flexamg_dimensions == 3) && (pos.z == 0 || pos.z == 1)))
			b[i] = problem->phi(pos.x, pos.y, pos.z);
		else
			b[i] = - problem->f(pos.x, pos.y, pos.z);
	}
	
	/*A.print();
	b.print();*/
	// eliminate dirichlet values
	A.eliminateDirichletValues(b);
	A.finish();	
	
	vx.create(n);
	vx = 0.0;
	// init randwerte
	for(int i=0; i<n; i++)
	{
		postype pos = GetPosForIndex(i);
		setSize(vx[i], UNKNOWN_NR);
		if(pos.x == 0 || pos.x == 1.0 || pos.y == 0 || pos.y == 1 || ((flexamg_dimensions == 3) && (pos.z == 0 || pos.z == 1)))
			vx[i] = problem->phi(pos.x, pos.y, pos.z);
		else
			vx[i] = 0.0;
	}	
	
	delete problem;
}


#define GETINDEX3D(x, y, z) (x + NX*(y) + NX*NY*(z))

//!
//! do finite-difference discretisation in structured grid on [0..1]x[0..1]x[0..1]
//! @param A SparseMatrix A for discretisation of the matrix
//! @param vx vector which is set to b[i] for dirichlet nodes, otherwise 0.0
//! @param b vector of the right hand side in Ax=b.
//! @param NX number of nodes in X-direction
//! @param NY number of nodes in Y-direction
//! @param NZ number of nodes in Z-direction
void discretization3d(myMatrix &A, myVector &vx, myVector &b, int NX, int NY, int NZ, bool full_star, int problemtype)
{
	cProblem *problem;
	if(problemtype == 0)
		problem = new Laplace;
	else
		problem = new Problemf6_3d;
	
	int n = NX*NY*NZ;

	positions.resize(n);
	double hX = 1/((double)(NX-1));
	double hY = 1/((double)(NY-1));
	double hZ = 1/((double)(NZ-1));
#ifndef EASY_MATRIX
	double hx2 = hX*hX, hy2 = hY*hY, hz2 = hZ*hZ;
#else
	double hx2 = 1.0; double hy2 = hY*hY/(hX*hX), hz2 = hZ*hZ/(hX*hX);
#endif
	
	positions.resize(n);
	stopwatch SW;
	SW.start();
	
	A.create(n, n);
	A.fromlevel = A.tolevel = 0;
	for(int y=0; y<NY; y++)
		for(int x=0; x<NX; x++)
			for(int z=0; z<NZ; z++)				
			{	
				int i = GETINDEX3D(x,y,z);
				positions[i].x = x*hX;
				positions[i].y = y*hY;
				positions[i].z = z*hZ;
				if(x == 0 || x == NX-1 || y == 0 || y == NY-1 || z == 0 || z == NZ-1)
				{
					myMatrix::connection con; con.iIndex = i; 
					setSize(con.dValue, UNKNOWN_NR, UNKNOWN_NR);
					con.dValue = 1.0;
					A.setMatrixRow(i, &con, 1);
				}
				else
				{
					if(full_star)
					{
						myMatrix::connection con[27];
						for(int k=0; k<27; k++) setSize(con[k].dValue, UNKNOWN_NR, UNKNOWN_NR);
						double diag = 0;
						
						int k=1;
						//x 1 2
						con[k].iIndex = GETINDEX3D(x-1, y, z);	con[k++].dValue = -1.0 / hx2;
						con[k].iIndex = GETINDEX3D(x+1, y, z);	con[k++].dValue = -1.0 / hx2;
						diag += 2.0/hx2;
						//y 3 4
						con[k].iIndex = GETINDEX3D(x, y-1, z);	con[k++].dValue = -1.0 / hy2;
						con[k].iIndex = GETINDEX3D(x, y+1, z);	con[k++].dValue = -1.0 / hy2;
						diag += 2.0/hy2;
						//z 5 6
						con[k].iIndex = GETINDEX3D(x, y, z-1);	con[k++].dValue = -1.0 / hz2;
						con[k].iIndex = GETINDEX3D(x, y, z+1);	con[k++].dValue = -1.0 / hz2;
						diag += 2.0/hz2;
						
						//xy 7 8 9 10
						con[k].iIndex = GETINDEX3D(x-1, y+1, z);	con[k++].dValue = -1.0 / (hx2+hy2);
						con[k].iIndex = GETINDEX3D(x+1, y-1, z);	con[k++].dValue = -1.0 / (hx2+hy2);
						diag += 2.0/(hx2+hy2);
						con[k].iIndex = GETINDEX3D(x-1, y-1, z);	con[k++].dValue = -1.0 / (hx2+hy2);
						con[k].iIndex = GETINDEX3D(x+1, y+1, z);	con[k++].dValue = -1.0 / (hx2+hy2);
						diag += 2.0/(hx2+hy2);
						
						//xz 11 12 13 14			
						con[k].iIndex = GETINDEX3D(x-1, y, z+1);	con[k++].dValue = -1.0 / (hx2+hz2);
						con[k].iIndex = GETINDEX3D(x+1, y, z-1);	con[k++].dValue = -1.0 / (hx2+hz2);
						diag += 2.0/(hx2+hz2);					
						con[k].iIndex = GETINDEX3D(x-1, y, z-1);	con[k++].dValue = -1.0 / (hx2+hz2);
						con[k].iIndex = GETINDEX3D(x+1, y, z+1);	con[k++].dValue = -1.0 / (hx2+hz2);
						diag += 2.0/(hx2+hz2);
						
						//yz 15 16 17 18
						con[k].iIndex = GETINDEX3D(x, y-1, z+1);	con[k++].dValue = -1.0 / (hy2+hz2);
						con[k].iIndex = GETINDEX3D(x, y+1, z-1);	con[k++].dValue = -1.0 / (hy2+hz2);
						diag += 2.0/(hy2+hz2);
						con[k].iIndex = GETINDEX3D(x, y-1, z-1);	con[k++].dValue = -1.0 / (hy2+hz2);
						con[k].iIndex = GETINDEX3D(x, y+1, z+1);	con[k++].dValue = -1.0 / (hy2+hz2);
						diag += 2.0/(hy2+hz2);
						
						//xyz 19 20 21 22, 23 24 25 26
						con[k].iIndex = GETINDEX3D(x-1, y-1, z-1);	con[k++].dValue = -1.0 / (hx2+hy2+hz2);
						con[k].iIndex = GETINDEX3D(x+1, y+1, z+1);	con[k++].dValue = -1.0 / (hx2+hy2+hz2);
						diag += 2.0/(hx2+hy2+hz2);
						con[k].iIndex = GETINDEX3D(x-1, y+1, z-1);	con[k++].dValue = -1.0 / (hx2+hy2+hz2);
						con[k].iIndex = GETINDEX3D(x+1, y-1, z+1);	con[k++].dValue = -1.0 / (hx2+hy2+hz2);
						diag += 2.0/(hx2+hy2+hz2);
						con[k].iIndex = GETINDEX3D(x-1, y+1, z+1);	con[k++].dValue = -1.0 / (hx2+hy2+hz2);
						con[k].iIndex = GETINDEX3D(x+1, y-1, z-1);	con[k++].dValue = -1.0 / (hx2+hy2+hz2);
						diag += 2.0/(hx2+hy2+hz2);
						con[k].iIndex = GETINDEX3D(x-1, y-1, z+1);	con[k++].dValue = -1.0 / (hx2+hy2+hz2);
						con[k].iIndex = GETINDEX3D(x+1, y+1, z-1);	con[k++].dValue = -1.0 / (hx2+hy2+hz2);
						diag += 2.0/(hx2+hy2+hz2);					
						
						con[0].iIndex = i;
						con[0].dValue = diag;
						//assert(k==27);
						A.setMatrixRow(i, con, k);
					}
					else
					{
						// 7-point (needs agressive coarsening)
						myMatrix::connection con[7];
						for(int k=0; k<7; k++) setSize(con[k].dValue, UNKNOWN_NR, UNKNOWN_NR);
						con[0].iIndex = i;						con[0].dValue = 2.0 / hx2 + 2.0 / hy2 + 2.0/hz2;
						con[1].iIndex = GETINDEX3D(x-1,y,z);	con[1].dValue = -1.0/hx2;
						con[2].iIndex = GETINDEX3D(x+1,y,z);	con[2].dValue = -1.0/hx2;
						con[3].iIndex = GETINDEX3D(x,y-1,z);	con[3].dValue = -1.0/hy2;
						con[4].iIndex = GETINDEX3D(x,y+1,z);	con[4].dValue = -1.0/hy2;
						con[5].iIndex = GETINDEX3D(x,y,z-1);	con[5].dValue = -1.0/hz2;
						con[6].iIndex = GETINDEX3D(x,y,z+1);	con[6].dValue = -1.0/hz2;
						A.setMatrixRow(i, con, 7);
					}				
				}
			}
	
	cout << "assemble took ";
	SW.printTimeDiff();
	
	
	if(n<10003)
		A.writeToFile("/Users/mrupp/matrices/A.mat");
	
	// create rhs
	b.create(n);
	b = 0.0;
	for(int i=0; i<n; i++)
	{
		postype pos = GetPosForIndex(i);
		setSize(b[i], UNKNOWN_NR);
		if(pos.x == 0 || pos.x == 1.0 || pos.y == 0 || pos.y == 1 || pos.z == 0 || pos.z == 1)
			b[i] = problem->phi(pos.x, pos.y, pos.z);
		else
			b[i] = - problem->f(pos.x, pos.y, pos.z);
	}
	
	// eliminate dirichlet values
	A.eliminateDirichletValues(b);
	A.finish();	
	
	vx.create(n);
	vx = 0.0;
	// init randwerte
	for(int i=0; i<n; i++)
	{
		postype pos = GetPosForIndex(i);
		setSize(vx[i], UNKNOWN_NR);
		if(pos.x == 0 || pos.x == 1.0 || pos.y == 0 || pos.y == 1 || pos.z == 0 || pos.z == 1)
			vx[i] = problem->phi(pos.x, pos.y, pos.z);
		else
			vx[i] = 0.0;
	}	
	
	delete problem;
}

//! returns the error of vx to the known solution
double getDiscErr(myVector &vx, int NX, int NY, int NZ, int problemtype)
{
	cProblem *problem;
	if(problemtype == 0)
		problem = new Laplace;
	else if(flexamg_dimensions == 2)
		problem = new Problemf6_2d;
	else
		problem = new Problemf6_3d;
	
	double hX = 1/((double)(NX-1));
	double hY = 1/((double)(NY-1));
	double hZ = 1/((double)(NZ-1));
	double h2;
	if(flexamg_dimensions == 2)
		h2= hX*hY;
	else
		h2 = hX*hY*hZ;
	
	double res = 0;
	// check discretization error
	for(int i=0; i < vx.getLength(); i++)
	{			
		double d = problem->solution(GetPosForIndex(i).x, GetPosForIndex(i).y, GetPosForIndex(i).z) - getAt(vx[i], 0);
		res += d*d *h2;
	}
	delete problem;
	return sqrt(res);
}

void discretization(myMatrix &A, myVector &vx, myVector &b, int NX, int NY, int NZ, bool full_star, int problemtype)
{
	if(flexamg_dimensions == 2)
		discretization2d(A, vx, b, NX, NY, NZ, full_star, problemtype);
	else
		discretization3d(A, vx, b, NX, NY, NZ, full_star, problemtype);
}