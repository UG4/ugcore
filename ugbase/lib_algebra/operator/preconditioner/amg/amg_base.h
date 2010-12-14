/**
 * \file amg_base.h
 *
 * \author Martin Rupp
 *
 * \date 01.12.2010
 *
 * class declaration for amg base
 *
 * Goethe-Center for Scientific Computing 2010.
 */


#ifndef __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_BASE_H__
#define __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_BASE_H__

#include <vector>
#include <iostream>

#include "amg_debug_helper.h"
#include "sparsematrix_operator.h"


template<typename T>
std::string ToString(const T &t)
{
	std::stringstream out;
	out << t;
	return out.str();
}

namespace ug{

#define AMG_MAX_LEVELS 32


template <typename TAlgebra>
class amg_base:
	public IPreconditioner<	TAlgebra >
{
public:
//	Algebra type
	typedef TAlgebra algebra_type;

//	Vector type
	typedef typename TAlgebra::vector_type vector_type;

//	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;

	typedef typename matrix_type::value_type value_type;

//  functions
	amg_base();
	virtual ILinearIterator<vector_type,vector_type>* clone() = 0;

	//	Name of preconditioner
	virtual ~amg_base();
	void cleanup();

protected:
	virtual const char* name() const = 0;

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
			init();
			m_bInited = true;
		}
		return get_correction(c, d);
	}

public:
	bool get_correction_and_update_defect(vector_type &c, vector_type &d, int level=0);
	bool get_correction(vector_type &c, const vector_type &d);
/*
	int get_nr_of_coarse(int level)
	{
		assert(level+1 < used_levels);
		return A[level+1]->length;
	}
*/
	int get_nr_of_used_levels() { return used_levels; }

//  data

	void set_num_presmooth(int new_presmooth) { m_numPreSmooth = new_presmooth; } // nu1
	void set_num_postsmooth(int new_postsmooth) { m_numPostSmooth = new_postsmooth; } // nu2
	void set_cycle_type(int new_cycletype) { m_cycleType = new_cycletype; UG_ASSERT(m_cycleType == 1, "only cycle type 1 supported"); } // gamma

	void set_max_levels(int new_max_levels)
	{
		max_levels = new_max_levels;
		UG_ASSERT(max_levels <= AMG_MAX_LEVELS, "Currently only " << AMG_MAX_LEVELS << " level supported.\n");
	}
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

protected:
	bool create_level_vectors(int level);
	virtual void create_AMG_level(matrix_type &AH, SparseMatrix<double> &R, const matrix_type &A,
								SparseMatrix<double> &P, int level) = 0;
	virtual bool init();

// data
	int	m_numPreSmooth;						///< nu_1 : nr. of pre-smoothing steps
	int m_numPostSmooth;					///< nu_2: nr. of post-smoothing steps
	int m_cycleType;						///< gamma: cycle type (1 = V-Cycle, 2 = W-Cycle)

	int max_levels;							///< max. nr of levels used for FAMG
	int used_levels;						///< nr of FAMG levels used

	vector_type *vec1[AMG_MAX_LEVELS]; 	///< temporary Vector for storing r = Ax -b
	vector_type *vec2[AMG_MAX_LEVELS]; 	///< temporary Vector for storing rH
	vector_type *vec3[AMG_MAX_LEVELS]; 	///< temporary Vector for storing eH
	vector_type *vec4;						///< temporary Vector for defect (in get_correction)

	SparseMatrix<double> R[AMG_MAX_LEVELS]; ///< R Restriction Matrices
	SparseMatrix<double> P[AMG_MAX_LEVELS]; ///< P Prolongation Matrices
	matrix_type *A[AMG_MAX_LEVELS+1];		///< A Matrices
	SparseMatrixOperator<matrix_type, vector_type> SMO[AMG_MAX_LEVELS];

#ifdef UG_PARALLEL
	pcl::ParallelCommunicator<IndexLayout>
		*com[AMG_MAX_LEVELS]; 				///< the communicator objects on the levels
	IndexLayout pseudoLayout;				///< Pseudo-IndexLayout for the created ParallelVectors.
#endif

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



#include "amg_base_impl.h"


#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_BASE_H__
