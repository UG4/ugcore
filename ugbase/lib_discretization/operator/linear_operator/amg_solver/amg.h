/*
 *  amg.h
 *  flexamg
 *
 *  Created by Martin Rupp on 03.12.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 *  class declaration for amg
 */
#pragma once
#include <vector>
#include <iostream>

using namespace std;


#ifdef FLEXAMG
#include "sparseMatrix.h"
#include "CoarseSolver.h"
#include "preconditioner.h"
#endif

#include "maxheap.h"



#ifndef FLEXAMG
struct position2d
{
	double x, y, z;
	friend std::ostream &operator << (std::ostream &out, const position2d &p)
	{
		out << p.x << " " << p.y;
		//out << p.x << " " << p.y << " " << p.z;
		return out;
	}
};

struct cAMG_helper
{
	int GetOriginalIndex(int level, int i) const
	{
		while(level > 0)
			i = parentIndex[level--][i];
		return i;		
	}
	
	
	position2d GetPosForIndexAtLevel(int level, int i) const
	{
		return positions[GetOriginalIndex(level, i)];
	}
	
	
	void writePosToStream(std::ostream &out) const
	{
		out << size << endl;
		for(int i=0; i< size ; i++)
			out << positions[i] << endl;
	}
	
	const position2d *positions;
	int size;
	int **parentIndex;
};

#endif
#include "graph.h"



namespace ug{
#define AMG_MAX_LEVELS 32

struct amg_nodeinfo;
	
// AMG
//---------------------------------
//! algebraic multigrid class.
//!
template<typename Matrix_type, typename Vector_type>
#ifdef FLEXAMG
class amg : public preconditioner<Matrix_type, Vector_type>
#else
class amg
#endif
{
public:
	typedef typename Matrix_type::entry_type entry_type;
	
//  functions
	void writeMatrices(const char *pathAndName);
	amg() ;
	virtual ~amg();
	virtual bool init(const Matrix_type& A);

#ifdef FLEXAMG
	virtual void precond2(Vector_type &x, const Vector_type &b)
	{
		for(int i=0; i<gamma; i++)
			MGCycle(x, b, 0);
	}
	virtual double iterate(Vector_type &x, const Vector_type &b)
	{
		for(int i=0; i<gamma; i++)
			return MGCycle(x, b, 0);
		return 0.0;
	}

	
	double MGCycle(Vector_type &x, const Vector_type &b, int i=0);	
	void printCoarsening(int level, amg_nodeinfo *grid);
	
	bool onlyOneLevel(const Matrix_type& A_);
	void amgTest(const Matrix_type& A_, Vector_type &x, const Vector_type &b);
	void amgTestLevel(Vector_type &x, const Vector_type &b, int level);
#else
	bool get_correction_and_update_defect(Vector_type &d, Vector_type &c, int level=0);
#endif
	
	int getNrOfCoarse(int level)
	{
		assert(level+1 < used_levels);
		return A[level+1]->length;
	}
	
	int getNrOfUsedLevels() { return used_levels; }
//  data
	int max_levels;
	int aggressiveCoarsening;
	int aggressiveCoarseningNrOfPaths;
	
	void setAggressiveCoarsening_A_2() { aggressiveCoarsening = true; aggressiveCoarseningNrOfPaths = 2;}
	void setAggressiveCoarsening_A_1() { aggressiveCoarsening = true; aggressiveCoarseningNrOfPaths = 1;}
	
private:
	int	nu1;		///< nu_1 : nr. of pre-smoothing steps
	int nu2;		///< nu_2: nr. of post-smoothing steps
	int gamma;		///< gamma: cycle type (1 = V-Cycle, 2 = W-Cycle)
	double eps_truncation_of_interpolation;
	double theta; ///< measure for strong connectivity
	double sigma; ///< dunno
	
private:
//  functions
	void createAMGLevel(Matrix_type &AH, SparseMatrix<double> &R, const Matrix_type &A, SparseMatrix<double> &P, int level);
	//void createGalerkinMatrix(Matrix_type &AH, const SparseMatrix<double> &R, const Matrix_type &A, const SparseMatrix<double> &P, int *posInConnections);
	
	void CreateProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int iNrOfCoarse, int &unassigned, amg_nodeinfo *grid);
	void CreateIndirectProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int *posInConnections, int unassigned, amg_nodeinfo *grid);
	
	//! creates the graph
	void CreateGraph(const Matrix_type &A, cgraph &graph, maxheap<amg_nodeinfo> &PQ, int &unassigned, amg_nodeinfo *grid);
	void CreateGraph2(cgraph &graph, cgraph &graph2, maxheap<amg_nodeinfo> &PQ, int &unassigned, int &iNrOfCoarse, int *posInConnections, int *newIndex, amg_nodeinfo *grid);
	int Coarsen(cgraph &graph, maxheap<amg_nodeinfo> &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, const Matrix_type &A, amg_nodeinfo *grid);
	

//	data
#ifdef FLEXAMG
	CoarseSolver coarseSolver;	///< the coarse(st) grid solver
#else
	LapackLU coarseSolver;
#endif

	int used_levels;			///< nr of AMG levels used

	Vector_type *vec1[AMG_MAX_LEVELS]; ///< temporary Vector for storing r = Ax -b
	Vector_type *vec2[AMG_MAX_LEVELS]; ///< temporary Vector for storing rH
	Vector_type *vec3[AMG_MAX_LEVELS]; ///< temporary Vector for storing eH
	
	SparseMatrix<double> R[AMG_MAX_LEVELS]; ///< R Restriction Matrices
	SparseMatrix<double> P[AMG_MAX_LEVELS]; ///< P Prolongation Matrices
	Matrix_type *A[AMG_MAX_LEVELS+1];		///< A Matrices
	

#ifdef FLEXAMG
	sgs<Matrix_type, Vector_type> smoother[AMG_MAX_LEVELS];  ///< smoother for each level
#endif

public:
	int *parentIndex[AMG_MAX_LEVELS];
	cAMG_helper amghelper;

#ifndef FLEXAMG
	const position2d *positions;
#endif
};
	
	


} // namespace ug
#include "amg_impl.h"
