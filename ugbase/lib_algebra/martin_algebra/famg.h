#pragma once
#include <vector>
#include <iostream>
//using namespace std;

#include "sparseMatrix.h"

#include "preconditioner.h"
#include "maxheap.h"
#include "graph.h"


#include "CoarseSolver.h"

#define FAMG_MAX_LEVELS 32

// AMG
//---------------------------------
//! algebraic multigrid class.
//!
template<typename Matrix_type, typename Vector_type>
class famg : public preconditioner<Matrix_type, Vector_type>
{
public:
	typedef typename Matrix_type::entry_type entry_type;
	
	//  functions
	void writeMatrices(const char *pathAndName);
	famg() ;
	~famg();
	virtual bool init(const Matrix_type& A);
	
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
	
	bool onlyOneLevel(const Matrix_type& A_);
	void amgTest(const Matrix_type& A_, Vector_type &x, const Vector_type &b);
	void amgTestLevel(Vector_type &x, const Vector_type &b, int level);
	
	int getNrOfCoarse(int level)
	{
		ASSERT1(level+1 < used_levels);
		return A[level+1]->length;
	}
	
	int getNrOfUsedLevels() { return used_levels; }
	//  data
	int max_levels;

	
private:
	//  structs
#define ASSIGNED_RATING					(-1000000000)
#define COARSE_RATING					(-1000000001)
#define FINE_RATING						(-2000000000)
#define FINE_RATING_INDIRECT_UNASSIGNED	(-1000000002)
	struct nodeinfo
	{
		int rating;
		//int newIndex;		
		inline void setAssigned(){rating = ASSIGNED_RATING;}		
		inline void setCoarse()	{rating = COARSE_RATING;}
		inline void setFineDirect(){rating = FINE_RATING;	}
		
		inline void setFineIndirect(){rating = FINE_RATING_INDIRECT_UNASSIGNED;}
		
		inline bool isCoarse(){	return rating == COARSE_RATING;}
		inline bool isFineDirect()	{return (rating == FINE_RATING);}
		
		inline bool isAssigned(){return (rating <= ASSIGNED_RATING);}
		
		friend ostream &operator << (ostream &out, nodeinfo &n)
		{
			out << "Rating: " << n.rating;
			if(n.isCoarse()) out << " (coarse)";
			else if(n.isFineDirect()) out << " (fine)";
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
	double eps_truncation_of_interpolation;
	double theta; ///< measure for strong connectivity
	double sigma; ///< dunno
private:
	//  functions
	void createAMGLevel(Matrix_type &AH, SparseMatrix<double> &R, const Matrix_type &A, SparseMatrix<double> &P, int level);
	void createGalerkinMatrix(Matrix_type &AH, const SparseMatrix<double> &R, const Matrix_type &A, const SparseMatrix<double> &P, int *posInConnections);
	
	void CreateProlongation(SparseMatrix<double> &P, const Matrix_type &A, int *newIndex, int iNrOfCoarse, int &unassigned);
	
	//! creates the graph
	void CreateGraph(const Matrix_type &A, cgraph &graph, maxheap<nodeinfo> &PQ, int &unassigned);
	int Coarsen(cgraph &graph, maxheap<nodeinfo> &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, const Matrix_type &A);
	
	
	//	data
	CoarseSolver coarseSolver;	///< the coarse(st) grid solver
	nodeinfo *grid;				///< needed for construction
	int used_levels;			///< nr of AMG levels used
	
	Vector_type *vec1[FAMG_MAX_LEVELS]; ///< temporary Vector for storing r = Ax -b
	Vector_type *vec2[FAMG_MAX_LEVELS]; ///< temporary Vector for storing rH
	Vector_type *vec3[FAMG_MAX_LEVELS]; ///< temporary Vector for storing eH
	
	SparseMatrix<double> R[FAMG_MAX_LEVELS]; ///< R Restriction Matrices
	SparseMatrix<double> P[FAMG_MAX_LEVELS]; ///< P Restriction Matrices
	Matrix_type *A[FAMG_MAX_LEVELS+1];		///< A Matrices
	
	sgs<Matrix_type, Vector_type> smoother[FAMG_MAX_LEVELS];  ///< smoother for each level
};


#include "amg.hpp"