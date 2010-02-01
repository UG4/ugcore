#pragma once
#include <vector>
#include <iostream>
//using namespace std;

#include "sparseMatrix.h"

#include "preconditioner.h"
#include "maxheap.h"
#include "graph.h"


#include "CoarseSolver.h"

#define AMG_MAX_LEVELS 32

// AMG
//---------------------------------
//! algebraic multigrid class.
//!
template<typename entry_type, typename Vector_type>
class amg : public preconditioner<entry_type, Vector_type>
{
public:
	typedef SparseMatrix<entry_type> matrix_type;
	
//  functions
	void writeMatrices(const char *pathAndName);
	amg() ;
	~amg();
	virtual bool init(const matrix_type& A);

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
	void printCoarsening(int level);
	
	void interpolate(Vector_type *pto, const Vector_type &from, int tolevel); // 2h -> h prolongate
	void restriction(Vector_type *pto, const Vector_type &from, int tolevel); // h -> 2h
	
	
	bool onlyOneLevel(const matrix_type& A_);
	
	//void restrictioninjection(Vector_type &to, const Vector_type &from, int tolevel);
	
	int getNrOfCoarse(int level)
	{
		ASSERT1(level+1 < used_levels);
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
	//  structs
#define ASSIGNED_RATING			(-2000000000)
#define COARSE_RATING			(-2000000001)
#define FINE_RATING_INDIRECT	(-2000000002)
#define FINE_RATING				(-2000000003)
	struct nodeinfo
	{
		int rating;
		//int newIndex;		
		inline void setAssigned(){rating = ASSIGNED_RATING;}		
		inline void setCoarse()	{rating = COARSE_RATING;}
		inline void setFine(){rating = FINE_RATING;	}
		inline void setFineIndirect(){rating = FINE_RATING_INDIRECT;}
		inline bool isCoarse(){	return rating == COARSE_RATING;}
		inline bool isFine(){return (rating == FINE_RATING) || (rating == FINE_RATING_INDIRECT);}
		inline bool isFineNotIndirect()	{return (rating == FINE_RATING);}
		inline bool isFineIndirect(){ return rating == FINE_RATING_INDIRECT;}		
		inline bool isAssigned(){return (rating <= ASSIGNED_RATING);}
	
		friend ostream &operator << (ostream &out, nodeinfo &n)
		{
			out << "Rating: " << n.rating;
			if(n.isFineNotIndirect()) out << " (fine not ind)";
			else if(n.isFineIndirect()) out << " (fine ind)";
			else if(n.isCoarse()) out << " (coarse)";
			out << " ";
			return out;
		}
		void print()
		{
			cout << *this << endl;
		} // << " newindex: " << newIndex << endl;
		
		inline bool operator > (const nodeinfo &other)
		{
			if(rating == other.rating)
				return this < &other; // we somehow want a STABLE sort, for that coarsening is in the direction of the numbering of the elements
			else
				return rating > other.rating;
		}
	};
	
	int	nu1;		///< nu_1 : nr. of pre-smoothing steps
	int nu2;		///< nu_2: nr. of post-smoothing steps
	int gamma;		///< gamma: cycle type (1 = V-Cycle, 2 = W-Cycle)
private:
//  functions
	int getNodeWithBestRating(int n);
		
	void createAMGLevel(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A, SparseMatrix<double> &P, int level);
	void createGalerkinMatrix(matrix_type &AH, const SparseMatrix<double> &R, const matrix_type &A, const SparseMatrix<double> &P, int *posInConnections);
	
	void CreateProlongation(SparseMatrix<double> &P, const matrix_type &A, int *newIndex, int iNrOfCoarse);
	void CreateIndirectProlongation(SparseMatrix<double> &P, const matrix_type &A, int *newIndex, int *posInConnections);
	
	//! creates the graph
	void CreateGraph(const matrix_type &A, cgraph &graph, maxheap<nodeinfo> &PQ, int &unassigned);
	void CreateGraph2(cgraph &graph, cgraph &graph2, maxheap<nodeinfo> &PQ, int &unassigned, int *posInConnections);
	int Coarsen(cgraph &graph, maxheap<nodeinfo> &PQ, int *newIndex, int unassigned, bool bIndirect, const matrix_type &A);
	

//	data
	CoarseSolver coarseSolver;	///< the coarse(st) grid solver
	nodeinfo *grid;				///< needed for construction
	int used_levels;			///< nr of AMG levels used
	
	Vector_type *vec1[AMG_MAX_LEVELS]; ///< temporary Vector for storing r = Ax -b
	Vector_type *vec2[AMG_MAX_LEVELS]; ///< temporary Vector for storing rH
	Vector_type *vec3[AMG_MAX_LEVELS]; ///< temporary Vector for storing eH
	
	SparseMatrix<double> R[AMG_MAX_LEVELS]; ///< R Restriction Matrices
	SparseMatrix<double> P[AMG_MAX_LEVELS]; ///< P Restriction Matrices
	matrix_type *A[AMG_MAX_LEVELS+1];		///< A Matrices
	
	sgs<entry_type, Vector_type> smoother[AMG_MAX_LEVELS];  ///< smoother for each level
};


#include "amg.hpp"