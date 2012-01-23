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


#ifndef __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_H__
#define __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_H__

#include "lib_algebra/lib_algebra.h"

#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/operator/operator_iterator_interface.h"
#include "lib_algebra/operator/vector_writer.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "amg_debug_helper.h"
#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallel_nodes.h"
#endif

template<typename T>
std::string ToString(const T &t)
{
	std::stringstream out;
	out << t;
	return out.str();
}
#define PRINTLAYOUT(com, Layout1, Layout2) MyPrintLayout(com, (Layout1), (Layout2), #Layout1, #Layout2)

namespace ug{
#ifdef UG_PARALLEL
void MyPrintLayout(pcl::ParallelCommunicator<IndexLayout> &communicator, IndexLayout &layout1, IndexLayout &layout2, const char *name1, const char *name2);
#endif


template <typename TAlgebra>
class AMGBase:
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

///	Matrix Operator type
	typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	typedef typename matrix_type::value_type value_type;

	class LevelInformation
	{
	public:
		LevelInformation(int _level)
		{
			level = _level;
		}
		double get_creation_time_ms() const { return m_dCreationTimeMS; }

		size_t get_nr_of_nodes_min() const { return m_iNrOfNodesMin; }
		size_t get_nr_of_nodes_max() const { return m_iNrOfNodesMax; }
		size_t get_nr_of_nodes() const { return m_iNrOfNodesSum; }
		size_t get_nnz() const { return m_iNNZsSum; }
		size_t get_nnz_min() const { return m_iNNZsMin; }
		size_t get_nnz_max() const { return m_iNNZsMax; }

		// nr of elements in master/slave interfaces (including multiplicities)
		size_t get_nr_of_interface_elements() const { return m_iInterfaceElements; }

		double get_fill_in() const { return ((double)m_iNNZsSum)/(((double)m_iNrOfNodesSum)*((double)m_iNrOfNodesSum)); }
		double get_avg_nnz_per_row() const { return m_iNNZsSum/(double)m_iNrOfNodesSum; }
		size_t get_max_connections() const { return m_connectionsMax; }
		bool is_valid() const { return this != NULL; }

	public:
		void set_nr_of_nodes(size_t nrOfNodesMin, size_t nrOfNodesMax, size_t nrOfNodesSum)
		{
			m_iNrOfNodesMin = nrOfNodesMin;	m_iNrOfNodesMax = nrOfNodesMax;	m_iNrOfNodesSum = nrOfNodesSum;
		}
		void set_nnz(size_t nnzMin, size_t nnzMax, size_t nnzSum)
		{
			m_iNNZsMin = nnzMin; m_iNNZsMax = nnzMax; m_iNNZsSum= nnzSum;
		}
		void set_max_connections(size_t connectionsMax)
		{
			m_connectionsMax = connectionsMax;
		}

		void set_coarsening_rate(double rate)
		{
			m_dCoarseningRate = rate;
		}

		double get_coarsening_rate() const
		{
			return m_dCoarseningRate;
		}

		std::string tostring() const
		{
			std::stringstream ss;
			ss << "Level " << level << ": creation time: " << get_creation_time_ms() << " ms. number of nodes: " << get_nr_of_nodes() << ". fill in " <<
					get_fill_in()*100 << "%. ";
			if(level != 0) ss << "coarsening rate: " << get_coarsening_rate()*100 << "%. ";
			ss << "nr of interface elements: " << get_nr_of_interface_elements() << " (" << (double)get_nr_of_interface_elements()/(double)get_nr_of_nodes()*100.0 << "%) "
					<< "nnzs: " << get_nnz() << " avgNNZs/row: " << get_avg_nnz_per_row() << " maxCon: " << get_max_connections();
			return ss.str();
		}

	public:
		int level;
		double m_dCoarseningRate;
		double m_dCreationTimeMS;
		size_t m_iNrOfNodesMin;
		size_t m_iNrOfNodesMax;
		size_t m_iNrOfNodesSum;
		size_t m_iNNZsMin;
		size_t m_iNNZsMax;
		size_t m_iNNZsSum;
		size_t m_iInterfaceElements;
		size_t m_connectionsMax;
	};


//  functions
	AMGBase();
	virtual ILinearIterator<vector_type,vector_type>* clone() = 0;

	//	Name of preconditioner
	virtual ~AMGBase();
	void cleanup();


protected:
	virtual const char* name() const = 0;

//	Preprocess routine
	virtual bool preprocess(matrix_operator_type& mat);

//	Postprocess routine
	virtual bool postprocess() {return true;}

//	Stepping routine
	virtual bool step(matrix_operator_type& mat, vector_type& c, const vector_type& d)
	{
		UG_ASSERT(m_bInited, "not inited?");
		return get_correction(c, d);
	}

	void create_direct_solver(size_t level);
	bool solve_on_base(vector_type &c, vector_type &d, size_t level);

	template<typename TAMGNodes>
	void create_fine_marks(int level, TAMGNodes &amgnodes, size_t N);

	void create_parent_index(int level, stdvector<int> newIndex, size_t nrOfCoarse);

	template<typename TAMGNodes>
	void create_new_indices(prolongation_matrix_type &PoldIndices, prolongation_matrix_type &PnewIndices,
				size_t N, TAMGNodes &amgnodes, stdvector<int> &newIndex, double dEpsilonTruncation);
#ifndef UG_PARALLEL
	template<typename TAMGNodes>
	void serial_process_prolongation(prolongation_matrix_type &PoldIndices, prolongation_matrix_type &PnewIndices, double dEpsilonTruncation, int level, TAMGNodes &amgnodes);
#else
	template<typename TAMGNodes>
	void parallel_process_prolongation(prolongation_matrix_type &PoldIndices, prolongation_matrix_type &PnewIndices, double dEpsilonTruncation, int level, TAMGNodes &amgnodes,
			ParallelNodes &PN, bool bCreateNewNodes, IndexLayout &nextLevelMasterLayout, IndexLayout &nextLevelSlaveLayout);

	void communicate_prolongation(ParallelNodes &PN, prolongation_matrix_type &PoldIndices, bool bCreateNewNodes);

	void create_minimal_layout_for_prolongation(ParallelNodes &PN, prolongation_matrix_type &P, IndexLayout &newMasterLayout, IndexLayout &newSlaveLayout);

	template<typename TAMGNodes>
	void postset_coarse(ParallelNodes &PN, prolongation_matrix_type &PoldIndices, TAMGNodes &nodes);

#endif

public:
	void write_interfaces();
	bool add_correction_and_update_defect(vector_type &c, vector_type &d, size_t level, size_t exactLevel);
	bool check2(const vector_type &const_c, const vector_type &const_d);
	bool check_fsmoothing();

	bool add_correction_and_update_defect(vector_type &c, vector_type &d, size_t level=0);
	bool add_correction_and_update_defect2(vector_type &c, vector_type &d, size_t level=0);
	bool get_correction(vector_type &c, const vector_type &d);

	bool injection(vector_type &vH, const vector_type &v, size_t level);
/*
	size_t get_nr_of_coarse(size_t level)
	{
		assert(level+1 < m_usedLevels);
		return A[level+1]->length;
	}
*/
	size_t get_used_levels() const { return m_usedLevels; }

	bool check_level(vector_type &c, vector_type &d, size_t level);
//	bool check(IMatrixOperator const vector_type &const_c, const vector_type &const_d);
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

	void 	set_fsmoothing(bool enable) 							{ m_bFSmoothing = enable; }
	bool 	get_fsmoothing() const									{ return m_bFSmoothing; }

	void 	set_max_nodes_for_base(size_t newMaxNodesForBase) 		{ m_maxNodesForBase = newMaxNodesForBase; }
	size_t 	get_max_nodes_for_base() const							{ return m_maxNodesForBase; }

	void 	set_max_fill_before_base(double newMaxFillBeforeBase)	{ m_dMaxFillBeforeBase = newMaxFillBeforeBase;	}
	double	get_max_fill_before_base() const						{ return m_dMaxFillBeforeBase;	}


	void 	set_min_nodes_on_one_processor(size_t newMinNodes)		{ m_minNodesOnOneProcessor = newMinNodes; }
	size_t	get_min_nodes_on_one_processor() const					{ return m_minNodesOnOneProcessor; }

	void 	set_preferred_nodes_on_one_processor(size_t i)			{ m_preferredNodesOnOneProcessor = i; }
	size_t	get_preferred_nodes_on_one_processor() const			{ return m_preferredNodesOnOneProcessor; }


	void 	set_presmoother(ILinearIterator<vector_type, vector_type> *presmoother) {	m_presmoother = presmoother; }
	void 	set_postsmoother(ILinearIterator<vector_type, vector_type> *postsmoother) { m_postsmoother = postsmoother; }
	void 	set_base_solver(ILinearOperatorInverse<vector_type, vector_type> *basesolver) { m_basesolver = basesolver; }


	void set_position_provider(IPositionProvider<2> *ppp2d)
	{
		m_pPositionProvider2d = ppp2d;
	}

	void set_position_provider(IPositionProvider<3> *ppp3d)
	{
		m_pPositionProvider3d = ppp3d;
	}

	void set_matrix_write_path(const char *path)
	{
		m_writeMatrixPath = path;
		m_writeMatrices = true;
	}

	void set_positions(const MathVector<2> *pos, size_t size)
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
	void set_positions(const MathVector<3> *pos, size_t size)
	{
		m_dbgPositions.resize(size);
		for(size_t i=0; i<size; ++i)
			m_dbgPositions[i] = pos[i];
		m_dbgDimension = 3;
	}


	virtual void tostring() const;

public:
	void write_debug_matrices(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
			prolongation_matrix_type &P, size_t level);

	template<typename TMatrix>
	void write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name);

	template<typename TNodeType>
	void write_debug_matrix_markers(size_t level, const TNodeType &nodes);

	void set_one_init(bool b)
	{
		m_bOneInit = b;
	}

protected:
	void init_fsmoothing();
	bool writevec(std::string filename, const vector_type &d, size_t level);
	void update_positions();

	bool create_level_vectors(size_t level);
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
			prolongation_matrix_type &P, size_t level) = 0;
	bool f_smoothing(vector_type &corr, vector_type &d, size_t level);



#ifdef UG_PARALLEL
	bool agglomerate(size_t level);
#endif

// data
	size_t 	m_numPreSmooth;						///< nu_1 : nr. of pre-smoothing steps
	size_t 	m_numPostSmooth;					///< nu_2: nr. of post-smoothing steps
	int 	m_cycleType;						///< gamma: cycle type (1 = V-Cycle, 2 = W-Cycle)

	size_t 	m_maxLevels;						///< max. nr of levels used for FAMG
	size_t	m_usedLevels;						///< nr of FAMG levels used

	size_t 	m_maxNodesForBase;					///< max nr of coarse nodes before Base solver is used
	double 	m_dMaxFillBeforeBase;				///< max fill rate before Base solver is used
	size_t	m_minNodesOnOneProcessor;			///< min nr of nodes on one processor (for agglomeration)
	size_t	m_preferredNodesOnOneProcessor;		///< preferred nr of nodes on one processor (for agglomeration)
	bool 	m_bUseCollectedSolver;

	bool	m_bFSmoothing;
	bool	m_bOneInit;


	vector_type *m_vec4;						///< temporary Vector for defect (in get_correction)

	int iteration_glboal;


	bool 	m_writeMatrices;

	std::string m_writeMatrixPath;
	cAMG_helper m_amghelper;					///< helper struct for viewing matrices (optional)
	stdvector<MathVector<3> > m_dbgPositions;	///< positions of geometric grid (optional)
	int m_dbgDimension;							///< dimension of geometric grid (optional)


	ILinearIterator<vector_type, vector_type> *m_presmoother;	///< presmoother template
	ILinearIterator<vector_type, vector_type> *m_postsmoother;	///< postsmoother template \note: may be pre=post, is optimized away.

	ILinearOperatorInverse<vector_type, vector_type> *m_basesolver; ///< the base solver


	bool m_bInited;					///< true if inited. needed since preprocess doesnt give us a ParallelCommunicator atm.
	double m_dOperatorComplexity;
	double m_dGridComplexity;
	double m_dTimingWholeSetupMS;
	double m_dTimingCoarseSolverSetupMS;


	IPositionProvider<2> *m_pPositionProvider2d;
	IPositionProvider<3> *m_pPositionProvider3d;

	stdvector<stdvector<int> > m_parentIndex;		///< parentIndex[i] is the index of i on the finer level

	struct AMGLevel
	{
		AMGLevel(int level) : m_levelInformation(level)
		{
			pA = NULL;
			presmoother = NULL;
			postsmoother = NULL;
#ifdef UG_PARALLEL
			slaveLayout.clear();
			masterLayout.clear();
			agglomerateMasterLayout.clear();
			bHasBeenMerged = false;
#endif
		}
		vector_type corr;					///< temporary Vector for storing the correction made on this level
		vector_type cH;						///< temporary Vector for storing rH
		vector_type dH; 					///< temporary Vector for storing eH

		ILinearIterator<vector_type, vector_type> *presmoother;	///< presmoothers for each level
		ILinearIterator<vector_type, vector_type> *postsmoother;	///< postsmoothers for each level

		stdvector<bool> is_fine;

		prolongation_matrix_type R; 	///< R Restriction Matrices
		prolongation_matrix_type P; 	///< P Prolongation Matrices
		matrix_operator_type *pA;				///< A Matrices

#ifdef UG_PARALLEL
		pcl::ParallelCommunicator<IndexLayout> com; ///< the communicator object on this level
		IndexLayout slaveLayout, masterLayout;

		stdvector< typename block_traits<typename matrix_type::value_type>::inverse_type > m_diagInv;

		// agglomeration
		bool bHasBeenMerged;
		// level 0 - m_agglomerateLevel
		pcl::ProcessCommunicator processCommunicator;

		// level 0 - m_agglomerateLevel-1
		IndexLayout agglomerateMasterLayout;
		vector_type collC, collD;
		matrix_operator_type uncollectedA;
#endif


		LevelInformation m_levelInformation;
	};

	stdvector<AMGLevel*> levels;

#ifdef UG_PARALLEL
	// on level m_agglomerateLevel
	IndexLayout agglomerateSlaveLayout;
	size_t m_agglomerateLevel;

	// stuff from old level 0 agglomeration
	pcl::ProcessCommunicator m_emptyPC;
	IndexLayout m_emptyLayout;
#endif

	void calculate_level_information(size_t level, double createAMGlevelTiming);

public:
	//! \return c_A = total nnz of all matrices divided by nnz of matrix A
	double get_operator_complexity() const { return m_dOperatorComplexity; }

	//! \return c_G = total number of nodes of all levels divided by number of nodes on level 0
	double get_grid_complexity() const { return m_dGridComplexity; }

	//! \return the time spent on the whole setup in ms
	double get_timing_whole_setup_ms() const { return m_dTimingWholeSetupMS; }

	//! \return the time spent in the coarse solver setup in ms
	double get_timing_coarse_solver_setup_ms() const { return m_dTimingCoarseSolverSetupMS; }

	void print_level_information() const
	{
		std::cout.setf(std::ios::fixed);
		UG_LOG("Operator Complexity: " << get_operator_complexity() << ", grid complexity: " << get_grid_complexity() << ".\n");
		UG_LOG("Whole setup took " << get_timing_whole_setup_ms() << " ms, coarse solver setup took " << get_timing_coarse_solver_setup_ms() << " ms.\n");
		/*for(int i = 0; i<get_used_levels(); i++)
			UG_LOG(get_level_information(i)->tostring() << "\n");*/

	}

	const LevelInformation *get_level_information(size_t i) const
	{
		if(i < levels.size())
			return &levels[i]->m_levelInformation;
		else return NULL;
	}
};


///	@}

} // namespace ug



#include "amg_base_impl.h"


#endif // __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_H__
