/**
 * \file amg_rs_prolongation.h
 *
 * \author Martin Rupp
 *
 * \date 06.08.2010
 *
 * class declaration for amg
 *
 * Goethe-Center for Scientific Computing 2009-2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__

#include <vector>
#include <iostream>

#include "maxheap.h"

using namespace std;

namespace ug{

template<typename T>
std::string ToString(const T &t)
{
	std::stringstream out;
	out << t;
	return out.str();
}


struct cAMG_helper
{
	int GetOriginalIndex(int level, int i) const
	{
		while(level > 0)
			i = parentIndex[level--][i];
		return i;		
	}
	
	
	MathVector<3> GetPosForIndexAtLevel(int level, int i) const
	{
		return positions[GetOriginalIndex(level, i)];
	}
	
	
	void writePosToStream(std::ostream &out) const
	{
		out << dimension << endl;
		out << size << endl;
		for(int i=0; i< size ; i++)
		{
			out << positions[i][0] << " " << positions[i][1];
			if(dimension == 3)
				out << " " << positions[i][2];
			out << endl;
		}
	}
	
	const MathVector<3> *positions;
	int size;
	int dimension;
	int **parentIndex;
};

}
#include "graph.h"

#include "amg_rs_prolongation.h"
#include "amg_debug.h"

namespace ug{

#define AMG_MAX_LEVELS 32

struct amg_nodeinfo;


template<typename Matrix_type>
void
CreateStrongConnectionGraph(const Matrix_type &A, cgraph &graph, double theta=0.25);

void CreateMeasureOfImportancePQ(cgraph &strong, cgraph &strongT, maxheap<amg_nodeinfo> &PQ, int &unassigned, amg_nodeinfo *nodes);

void CreateAggressiveCoarseningGraph(cgraph &graph, cgraph &graph2, amg_nodeinfo *nodes,
		int nrOfPaths, int *posInConnections);


void CreateMeasureOfImportanceAggressiveCoarseningPQ(cgraph &graphAC, maxheap<amg_nodeinfo> &PQ, int &unassigned, int &iNrOfCoarse, int *newIndex, amg_nodeinfo *nodes);

int Coarsen(cgraph &graph, maxheap<amg_nodeinfo> &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, amg_nodeinfo *nodes);

void PreventFFConnections(cgraph &graphS, cgraph &graphST, amg_nodeinfo *nodes, int *newIndex, int &nrOfCoarse);
	

// AMG
//---------------------------------
//! algebraic multigrid class.
//!
template<typename Matrix_type, typename Vector_type>
class amg
{
public:
	typedef typename Matrix_type::entry_type entry_type;
	
//  functions
	void writeMatrices(const char *pathAndName);
	amg() ;
	virtual ~amg();
	virtual bool init(const Matrix_type& A);

	bool get_correction_and_update_defect(Vector_type &d, Vector_type &c, int level=0);
	
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

//	data
    LapackLU coarseSolver;

	int used_levels;			///< nr of AMG levels used

	Vector_type *vec1[AMG_MAX_LEVELS]; ///< temporary Vector for storing r = Ax -b
	Vector_type *vec2[AMG_MAX_LEVELS]; ///< temporary Vector for storing rH
	Vector_type *vec3[AMG_MAX_LEVELS]; ///< temporary Vector for storing eH
	
	SparseMatrix<double> R[AMG_MAX_LEVELS]; ///< R Restriction Matrices
	SparseMatrix<double> P[AMG_MAX_LEVELS]; ///< P Prolongation Matrices
	Matrix_type *A[AMG_MAX_LEVELS+1];		///< A Matrices
	
public:
	int *parentIndex[AMG_MAX_LEVELS];
	cAMG_helper amghelper;

	void set_debug_positions(const MathVector<2> *pos, size_t size)
	{
		dbg_positions.resize(size);
		for(size_t i=0; i<size; ++i)
		{
			dbg_positions[i].x = pos[i].x;
			dbg_positions[i].y = pos[i].y;
			dbg_positions[i].z = 0.0;
		}
		dbg_dimension = 2;
	}
	void set_debug_positions(const MathVector<3> *pos, size_t size)
	{
		dbg_positions.resize(size);
		for(size_t i=0; i<size; ++i)
			dbg_positions[i] = pos[i];
		dbg_dimension = 3;
	}

private:
	vector<MathVector<3> > dbg_positions;
	int dbg_dimension;
};
	
	


} // namespace ug
#include "amg_impl.h"



#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__
