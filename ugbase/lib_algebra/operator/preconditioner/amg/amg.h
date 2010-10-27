/**
 * \file amg.h
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

using namespace std;

template<typename T>
std::string ToString(const T &t)
{
	std::stringstream out;
	out << t;
	return out.str();
}

#include "amg_debug_helper.h"
#include "graph.h"

#include "amg_rs_prolongation.h"
#include "amg_debug.h"
#include "amg_nodeinfo.h"

/// \defgroup AMG

/**
 * \brief Algebraic Multigrid Functions.
 *
 *
 * \defgroup lib_algebra_AMG AMG
 * \ingroup lib_algebra
 */



namespace ug{

/// \addtogroup lib_algebra_AMG
///	@{

#define AMG_MAX_LEVELS 32


template<typename Matrix_type>
void
CreateStrongConnectionGraph(const Matrix_type &A, cgraph &graph, double theta=0.25);

void CreateMeasureOfImportancePQ(cgraph &strong, cgraph &strongT, nodeinfo_pq_type &PQ, int &unassigned, amg_nodeinfo *nodes);

void CreateAggressiveCoarseningGraph(cgraph &graph, cgraph &graph2, amg_nodeinfo *nodes,
		int nrOfPaths, int *posInConnections);


void CreateMeasureOfImportanceAggressiveCoarseningPQ(cgraph &graphAC, nodeinfo_pq_type &PQ, int &unassigned, int &iNrOfCoarse, int *newIndex, amg_nodeinfo *nodes);

int Coarsen(cgraph &graph, nodeinfo_pq_type &PQ, int *newIndex, int unassigned, int &iNrOfCoarse, amg_nodeinfo *nodes);

void PreventFFConnections(cgraph &graphS, cgraph &graphST, amg_nodeinfo *nodes, int *newIndex, int &nrOfCoarse);
	

template <typename matrix_type, typename vector_type>
class SparseMatrixOperator : public virtual IMatrixOperator<vector_type, vector_type, matrix_type >
{
	public:
		SparseMatrixOperator() {}
		SparseMatrixOperator(matrix_type *pMatrix) : m_pMatrix(pMatrix) { }
		void setmatrix(matrix_type *pMatrix) { m_pMatrix = pMatrix; }
		bool init(const vector_type &) { return true; }
		bool init() { return true; }
		bool apply(vector_type &c, const vector_type &d)
		{
			MatMult(c, 1.0, *m_pMatrix, d);
			return true;
		}
		bool apply_sub(vector_type &c, const vector_type &d)
		{
			MatMultAdd(c, -1.0, *m_pMatrix, d, 1.0, c);
			return true;
		}

		virtual ~SparseMatrixOperator(){}

	public:
	// 	Access to matrix
		virtual matrix_type& get_matrix() { return *m_pMatrix; };

	private:
		matrix_type *m_pMatrix;
};

// AMG
//---------------------------------
//! algebraic multigrid class.
//!

template <typename TAlgebra>
class amg:
	public IPreconditioner<	TAlgebra >
{
public:
//	Algebra type
	typedef TAlgebra algebra_type;

//	Vector type
	typedef typename TAlgebra::vector_type vector_type;

//	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;

	typedef typename matrix_type::entry_type entry_type;
	
//  functions
	amg() ;
	virtual ILinearIterator<vector_type,vector_type>* clone()
	{
		amg<algebra_type>* clone = new amg<algebra_type>();
		return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
	}
	//	Name of preconditioner
	virtual ~amg();
	void cleanup();

protected:
	virtual const char* name() const {return "AMGPreconditioner";}

//	Preprocess routine
	virtual bool preprocess(matrix_type& mat);

//	Postprocess routine
	virtual bool postprocess() {return true;}

//	Stepping routine
	virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
	{
		if(m_bInited == false)
		{
#ifdef UG_PARALLEL
			// set level 0 communicator
			com[0] = &c.get_communicator();
#endif
			init_amg();
			m_bInited = true;
		}
		return get_correction(c, d);
	}

public:
	bool get_correction_and_update_defect(vector_type &c, vector_type &d, int level=0);
	bool get_correction(vector_type &c, const vector_type &d);
	
	int get_nr_of_coarse(int level)
	{
		assert(level+1 < used_levels);
		return A[level+1]->length;
	}
	
	int get_nr_of_used_levels() { return used_levels; }

//  data
	

	void set_nu1(int new_nu1) { nu1 = new_nu1; }
	void set_nu2(int new_nu2) { nu2 = new_nu2; }
	void set_gamma(int new_gamma) { gamma = new_gamma; }
	void set_theta(double new_theta) { theta = new_theta; }
	void set_sigma(double new_sigma) { sigma = new_sigma; }
	void set_max_levels(int new_max_levels)
	{
		max_levels = new_max_levels;
		UG_ASSERT(max_levels <= AMG_MAX_LEVELS, "Currently only " << AMG_MAX_LEVELS << " level supported.\n");
	}
	void set_aggressive_coarsening_A_2() { aggressiveCoarsening = true; aggressiveCoarseningNrOfPaths = 2;}
	void set_aggressive_coarsening_A_1() { aggressiveCoarsening = true; aggressiveCoarseningNrOfPaths = 1;}
	void set_presmoother(ILinearIterator<vector_type, vector_type> *presmoother) {	m_presmoother = presmoother; }
	void set_postsmoother(ILinearIterator<vector_type, vector_type> *postsmoother) { m_postsmoother = postsmoother; }
	void set_base_solver(ILinearOperatorInverse<vector_type, vector_type> *basesolver) { m_basesolver = basesolver; }

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

	template <typename TGridFunction>
	bool set_debug(	TGridFunction& u)
	{
		static const int dim = TGridFunction::domain_type::dim;

		vector<MathVector<dim> > positions;
		ExtractPositions(u, positions);
		set_debug_positions(&positions[0], positions.size());
		UG_LOG("successfully set " << positions.size() << " positions.\n");
		return true;
	}

	void tostring() const;
private:
//  functions
	void create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
							SparseMatrix<double> &P, int level);

	bool init_amg();

private:
// data
	
	int	nu1;								///< nu_1 : nr. of pre-smoothing steps
	int nu2;								///< nu_2: nr. of post-smoothing steps
	int gamma;								///< gamma: cycle type (1 = V-Cycle, 2 = W-Cycle)

	double eps_truncation_of_interpolation;	///< parameter used for truncation of interpolation
	double theta; 							///< measure for strong connectivity
	double sigma;

	int max_levels;							///< max. nr of levels used for AMG
	int used_levels;						///< nr of AMG levels used
	bool aggressiveCoarsening;				///< true if aggressive coarsening is used on first level
	int aggressiveCoarseningNrOfPaths;  	///<

	vector_type *vec1[AMG_MAX_LEVELS]; 		///< temporary Vector for storing r = Ax -b
	vector_type *vec2[AMG_MAX_LEVELS]; 		///< temporary Vector for storing rH
	vector_type *vec3[AMG_MAX_LEVELS]; 		///< temporary Vector for storing eH
	vector_type *vec4;						///< temporary Vector for defect (in get_correction)
	
	SparseMatrix<double> R[AMG_MAX_LEVELS]; ///< R Restriction Matrices
	SparseMatrix<double> P[AMG_MAX_LEVELS]; ///< P Prolongation Matrices
	matrix_type *A[AMG_MAX_LEVELS+1];		///< A Matrices
	SparseMatrixOperator<matrix_type, vector_type> SMO[AMG_MAX_LEVELS];

	pcl::ParallelCommunicator<IndexLayout>
		*com[AMG_MAX_LEVELS]; 				///< the communicator objects on the levels
	IndexLayout pseudoLayout;				///< Pseudo-IndexLayout for the created ParallelVectors.

	int *parentIndex[AMG_MAX_LEVELS];		///< parentIndex[L][i] is the index of i on level L-1
	cAMG_helper amghelper;					///< helper struct for viewing matrices (optional)
	vector<MathVector<3> > dbg_positions;	///< positions of geometric grid (optional)
	int dbg_dimension;						///< dimension of geometric grid (optional)


	ILinearIterator<vector_type, vector_type> *m_presmoother;	///< presmoother template
	ILinearIterator<vector_type, vector_type> *m_postsmoother;	///< postsmoother template \note: may be pre=post, is optimized away.

	ILinearIterator<vector_type, vector_type> *m_presmoothers[AMG_MAX_LEVELS];	///< presmoothers for each level
	ILinearIterator<vector_type, vector_type> *m_postsmoothers[AMG_MAX_LEVELS];	///< postsmoothers for each level

	ILinearOperatorInverse<vector_type, vector_type> *m_basesolver; ///< the base solver


	bool m_bInited;				///< true if inited. needed since preprocess doesnt give us a ParallelCommunicator atm.
};
	
	
///	@}

} // namespace ug




#include "amg_impl.h"



#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_H__
