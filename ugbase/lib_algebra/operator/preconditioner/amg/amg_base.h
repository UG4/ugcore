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
	typedef typename TAlgebra::matrix_type prolongation_matrix_type;

	typedef typename matrix_type::value_type value_type;

	class LevelInformation
	{
	public:
		LevelInformation(double creationTimeMS, size_t iNrOfNodes) : m_dCreationTimeMS(creationTimeMS), m_iNrOfNodes(iNrOfNodes) { }
		double get_creation_time_ms() { return m_dCreationTimeMS; }
		double get_nr_of_nodes() { return m_iNrOfNodes; }
		bool is_valid() { return this != NULL; }
	private:
		double m_dCreationTimeMS;
		size_t m_iNrOfNodes;
	};


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
			com = &c.get_communicator();
#endif
			init();
			m_bInited = true;
		}
		return get_correction(c, d);
	}

public:
	bool get_correction_and_update_defect(vector_type &c, vector_type &d, size_t level=0);
	bool get_correction(vector_type &c, const vector_type &d);
/*
	size_t get_nr_of_coarse(size_t level)
	{
		assert(level+1 < m_usedLevels);
		return A[level+1]->length;
	}
*/
	size_t get_nr_of_used_levels() { return m_usedLevels; }

	bool check_level(vector_type &c, vector_type &d, size_t level);
	bool check(const vector_type &const_c, const vector_type &const_d);
//  data

	void 	set_num_presmooth(size_t new_presmooth) 				{ m_numPreSmooth = new_presmooth; }
	size_t	get_num_presmooth()	const				 				{ return m_numPreSmooth; }

	void 	set_num_postsmooth(size_t new_postsmooth) 				{ m_numPostSmooth = new_postsmooth; }
	size_t	get_num_postsmooth() const					 			{ return m_numPostSmooth; }

	void 	set_cycle_type(size_t new_cycletype)
	{
		UG_ASSERT(new_cycletype > 0, "cannot set cycle type " << new_cycletype << ", has to be > 0");
		m_cycleType = new_cycletype;
	}
	size_t	get_cycle_type() const									{ return m_cycleType; }

	void 	set_max_levels(size_t new_max_levels)
	{
		UG_ASSERT(new_max_levels > 0, "cannot set max levels " << new_max_levels << ", has to be > 0");
		m_maxLevels = new_max_levels;
	}
	size_t	get_max_levels() const									{ return m_maxLevels; }

	void 	set_fsmoothing(double fdamping) 						{ m_fDamp = fdamping; }
	double	get_fsmoothing() const									{ return m_fDamp; }

	void 	set_max_nodes_for_base(size_t newMaxNodesForBase) 	{ m_maxNodesForBase = newMaxNodesForBase; }
	size_t 	get_max_nodes_for_base() const							{ return m_maxNodesForBase; }

	void 	set_max_fill_before_base(double newMaxFillBeforeBase) { m_dMaxFillBeforeBase = newMaxFillBeforeBase;	}
	double	get_max_fill_before_base() const						{ return m_dMaxFillBeforeBase;	}


	void 	set_presmoother(ILinearIterator<vector_type, vector_type> *presmoother) {	m_presmoother = presmoother; }
	void 	set_postsmoother(ILinearIterator<vector_type, vector_type> *postsmoother) { m_postsmoother = postsmoother; }
	void 	set_base_solver(ILinearOperatorInverse<vector_type, vector_type> *basesolver) { m_basesolver = basesolver; }



	void set_matrix_write_path(const char *path)
	{
		m_writeMatrixPath = path;
		m_writeMatrices = true;
	}

	void set_debug_positions(const MathVector<2> *pos, size_t size)
	{
		m_dbgPositions.resize(size);
		for(size_t i=0; i<size; ++i)
		{
			m_dbgPositions[i].x = pos[i].x;
			m_dbgPositions[i].y = pos[i].y;
			m_dbgPositions[i].z = 0.0;
		}
		m_dbgDimension = 2;
	}
	void set_debug_positions(const MathVector<3> *pos, size_t size)
	{
		m_dbgPositions.resize(size);
		for(size_t i=0; i<size; ++i)
			m_dbgPositions[i] = pos[i];
		m_dbgDimension = 3;
	}

	template <typename TGridFunction>
	bool set_debug(	TGridFunction& u)
	{
		static const int dim = TGridFunction::domain_type::dim;

		stdvector<MathVector<dim> > positions;
		ExtractPositions(u, positions);
		set_debug_positions(&positions[0], positions.size());
		UG_LOG("successfully set " << positions.size() << " positions.\n");
		return true;
	}

	void tostring() const;

protected:
	bool create_level_vectors(size_t level);
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
			prolongation_matrix_type &P, size_t level) = 0;
	virtual bool init();
	bool do_f_smoothing(vector_type &corr, vector_type &d, size_t level);

// data
	size_t 	m_numPreSmooth;						///< nu_1 : nr. of pre-smoothing steps
	size_t 	m_numPostSmooth;					///< nu_2: nr. of post-smoothing steps
	int 	m_cycleType;						///< gamma: cycle type (1 = V-Cycle, 2 = W-Cycle)

	size_t 	m_maxLevels;						///< max. nr of levels used for FAMG
	size_t	m_usedLevels;						///< nr of FAMG levels used

	size_t 	m_maxNodesForBase;					///< max nr of coarse nodes before Base solver is used
	double 	m_dMaxFillBeforeBase;				///< max fill rate before Base solver is used

	std::string m_writeMatrixPath;
	bool 	m_writeMatrices;
	double 	m_fDamp;

	stdvector<vector_type*> m_vec1; 			///< temporary Vector for storing r = Ax -b
	stdvector<vector_type*> m_vec2; 			///< temporary Vector for storing rH
	stdvector<vector_type*> m_vec3; 			///< temporary Vector for storing eH
	vector_type *m_vec4;						///< temporary Vector for defect (in get_correction)

	stdvector<stdvector<bool> > is_fine;

	stdvector<prolongation_matrix_type *> m_R; 	///< R Restriction Matrices
	stdvector<prolongation_matrix_type *> m_P; 	///< P Prolongation Matrices
	stdvector<matrix_type *> m_A;				///< A Matrices
	stdvector< SparseMatrixOperator<matrix_type, vector_type> > m_SMO;

#ifdef UG_PARALLEL
	pcl::ParallelCommunicator<IndexLayout> * com;  ///< the communicator object on the levels
	stdvector<IndexLayout> slaveLayouts, masterLayouts;				///< Pseudo-IndexLayout for the created ParallelVectors.

#endif

	stdvector< stdvector<int> > m_parentIndex;	///< parentIndex[L][i] is the index of i on level L-1
	cAMG_helper m_amghelper;					///< helper struct for viewing matrices (optional)
	stdvector<MathVector<3> > m_dbgPositions;	///< positions of geometric grid (optional)
	int m_dbgDimension;							///< dimension of geometric grid (optional)


	ILinearIterator<vector_type, vector_type> *m_presmoother;	///< presmoother template
	ILinearIterator<vector_type, vector_type> *m_postsmoother;	///< postsmoother template \note: may be pre=post, is optimized away.

	stdvector<ILinearIterator<vector_type, vector_type> *> m_presmoothers;	///< presmoothers for each level
	stdvector<ILinearIterator<vector_type, vector_type> *> m_postsmoothers;	///< postsmoothers for each level

	ILinearOperatorInverse<vector_type, vector_type> *m_basesolver; ///< the base solver


	bool m_bInited;				///< true if inited. needed since preprocess doesnt give us a ParallelCommunicator atm.
	double m_dOperatorComplexity;
	double m_dNodesComplexity;
	double m_dTimingWholeSetupMS;
	double m_dTimingCoarseSolverMS;

	stdvector<LevelInformation> m_levelInformation;

public:
	double get_operator_complexity() { return m_dOperatorComplexity; }
	double get_nodes_complexity() { return m_dNodesComplexity; }
	double get_timing_whole_setup_ms() { return m_dTimingWholeSetupMS; }
	double get_timing_coarse_solver_ms() { return m_dTimingCoarseSolverMS; }
	LevelInformation *get_level_information(size_t i)
	{
		if(i < m_levelInformation.size())
			return &m_levelInformation[i];
		else return NULL;
	}
};


///	@}

} // namespace ug



#include "amg_base_impl.h"


#endif // __H__LIB_DISCRETIZATION__AMG_SOLVER__AMG_BASE_H__
