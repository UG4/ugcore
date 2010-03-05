/*
 *  AADisc.cpp
 *  flexamg
 *
 *  Created by Martin Rupp on 05.03.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#if 0
#include "algebra.h"
#include "AlgebraAccessor.h"


#define EASY_MATRIX
#define GETINDEX2D(x, y) (x + NX*(y))

typedef MultiIndex<4> index4;

extern std::vector<postype> positions;

void disc2d(MatrixAccessorBase<index4> *A, VectorAccessorBase<index4> *x, VectorAccessorBase<index4> *b, int NX, int NY, int ALPHA)
{
	DenseMatrix m(1,1);
	DenseMatrix mat(1,5);
	
	double hX = 1/((double)(NX-1));
	double hY = 1/((double)(NY-1));
	
#ifndef EASY_MATRIX
	double hx2 = hX*hX, hy2 = hY*hY;
#else
	double hx2 = 1.0; double hy2 = hY*hY/(hX*hX);
#endif
	
	
	vector<index4> I(1);
	vector<index4> J(5);
	for(int j=0; j<ALPHA; j++)
		for(int y=0; y<NY; y++)
			for(int x=0; x<NX; x++)
			{	
				int i = GETINDEX2D(x, y);	
				I[0][2] = i;
				I[0][3] = j;
				positions[i].x = x*hX;
				positions[i].y = y*hY;
				positions[i].z = 0;
				if(x == 0 || x == NX-1 || y == 0 || y == NY-1)
				{
					m(0,0) = 1.0;
					A->add(m, I, I);
				}
				else
				{
					// 5-point (needs agressive coarsening)
					J[0][3] = j; J[0][2] = i;					mat(0,0) = 2.0/hx2 + 2.0/hy2;
					J[1][3] = j; J[1][2] = GETINDEX2D(x-1, y);	mat(0,1) = -1.0/hx2;
					J[2][3] = j; J[2][2] = GETINDEX2D(x+1, y);	mat(0,2) = -1.0/hx2;
					J[3][3] = j; J[3][2] = GETINDEX2D(x, y+1);	mat(0,3) = -1.0/hy2;
					J[4][3] = j; J[4][2] = GETINDEX2D(x, y-1);	mat(0,4) = -1.0/hy2;
					A->add(mat, I, J);
				}
			}
	
	DenseVector d(1);
	
	// create rhs
	for(int j=0; j<ALPHA; j++)
	{
		I[0][3] = j;
		for(int i=0; i<NX*NY; i++)
		{
			postype pos = GetPosForIndex(i);
			if(pos.x == 0 || pos.x == 1.0 || pos.y == 0 || pos.y == 1 || ((flexamg_dimensions == 3) && (pos.z == 0 || pos.z == 1)))
				d(0) = 1.0;
			else
				d(0) = 0.0;
			I[0][2] = i;
			b->add(d, I);
		}
	}
	
	// init randwerte
	for(int j=0; j<ALPHA; j++)
	{
		I[0][3] = j;
		for(int i=0; i<NX*NY; i++)
		{
			postype pos = GetPosForIndex(i);
			if(pos.x == 0 || pos.x == 1.0 || pos.y == 0 || pos.y == 1 || ((flexamg_dimensions == 3) && (pos.z == 0 || pos.z == 1)))
				d(0) = 1.0;
			else
				d(0) = 0.0;
			I[0][2] = i;
			x->add(d, I);
		}	
	}
	
}


void disc1d(MatrixAccessorBase<index4> *A, VectorAccessorBase<index4> *x, VectorAccessorBase<index4> *b, int N, int ALPHA)
{
	DenseMatrix m(1,1);
	DenseMatrix mat(1,3);
	DenseVector vec(1);
	
	vector<index4> I(1), I2(1);
	vector<index4> J(3);
	for(int j=0; j<ALPHA; j++)
	{		
		index4 ind;
		ind[0] = 0;
		ind[1] = 0;
		ind[3] = j;
		
		// linkes ende
		m(0,0) = 1.0;
		
		ind[2] = 0;
		I[0] = ind;
		A->add(m, I,I);
		
		vec(0) = 3.0;	
		b->add(vec, I);
		
		// rechtes ende
		ind[2] = N-1;
		I[0] = ind;
		A->add(m, I, I);
		
		vec(0) = 1.0;	
		b->add(vec, I);
		
		
		for(int i=1; i<N-1; i++)
		{		
			ind[2] = i;
			I[0] = J[0] = ind;
			
			ind[2] = i-1;
			J[1] = ind;
			ind[2] = i+1;
			J[2] = ind;
			
			mat(0,0) = 2;
			mat(0,1) = -1;
			mat(0,2) = -1;
			
			A->add(mat, I, J);
		}	
	}	
	//A->print();
	//b->print();
}
#endif