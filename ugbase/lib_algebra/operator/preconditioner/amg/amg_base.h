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

//#include "lib_algebra/lib_algebra.h"

#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/operator_iterator.h"
#include "lib_algebra/operator/vector_writer.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "amg_debug_helper.h"

#include "lib_algebra/common/stl_debug.h" // stdvector

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
#define PRINTLAYOUT(pc, com, Layout1, Layout2) MyPrintLayout(pc, com, (Layout1), (Layout2), #Layout1, #Layout2)

namespace ug{
#ifdef UG_PARALLEL
void MyPrintLayout(pcl::InterfaceCommunicator<IndexLayout> &communicator, IndexLayout &layout1, IndexLayout &layout2, const char *name1, const char *name2);
#endif


template <typename TAlgebra>
class AMGBase:
	public IPreconditioner<	TAlgebra >
{
public:

///	Algebra type
	typedef TAlgebra algebra_type;

///	Vector type
	typedef typename TAlgebra::vector_type vector_type;

///	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;
///	Prolongation Matrix type
	typedef typename TAlgebra::matrix_type prolongation_matrix_type;

///	Matrix Operator type
	typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	typedef typename matrix_type::value_type value_type;


	class LevelInformation
	{
	public:
		LevelInformation(int _level=0)
		{
			level = _level;
			m_dCoarseningRate = 0.0;
			m_dCreationTimeMS = 0.0;
			m_iNrOfNodesMin=0;
			m_iNrOfNodesMax=0;
			m_iNrOfNodesSum=0;
			m_iNNZsMin=0;
			m_iNNZsMax=0;
			m_iNNZsSum=0;
			m_iInterfaceElements=0;
			m_connectionsMax=0;
			m_iParticipating=0;
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

		int get_participating()
		{
			return m_iParticipating;
		}

		std::string tostring() const
		{
			std::stringstream ss;
			ss << "Level " << level << ": creation time: " << get_creation_time_ms() << " ms. number of nodes: " << get_nr_of_nodes() << ". fill in " <<
					get_fill_in()*100 << "%. ";
			if(level != 0) ss << "coarsening rate: " << get_coarsening_rate()*100 << "%. ";
			ss << "nr of interface elements: " << get_nr_of_interface_elements() << " (" << (double)get_nr_of_interface_elements()/(double)get_nr_of_nodes()*100.0 << "%) "
					<< "nnzs: " << get_nnz() << " avgNNZs/row: " << get_avg_nnz_per_row() << " maxCon: " << get_max_connections();
#ifdef UG_PARALLEL
			ss << " participating: " << m_iParticipating;
#endif
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
		size_t m_iParticipating;
	};

	struct checkResult
	{
		double preSmoothing;
		double preFSmoothing;
		size_t iInnerIterations;
		double lastCoarseReduction;
		double coarseDefect;
		double coarseCorrection;
		double reduction;
		double postFSmoothing;
		double postSmoothing;
	};

//  functions
	AMGBase();
	virtual SmartPtr<ILinearIterator<vector_type> > clone() = 0;

	//	Name of preconditioner
	virtual ~AMGBase();
	void cleanup();


protected:
	virtual const char* name() const = 0;

	/**
	 * creates MG Hierachy for with matrix_operator_type A and temporary vectors for higher levels
	 * @param A matrix A.
	 * @return true
	 */
	virtual bool preprocess(matrix_operator_type& A);

	///	nothing to do here
	virtual bool postprocess() {return true;}

//	Stepping routine
	virtual bool step(matrix_operator_type& mat, vector_type& c, const vector_type& d)
	{
		UG_ASSERT(m_bInited, "not inited?");
		return get_correction(c, d);
	}

	/**
	 * create direct solver on a specified level
	 * @param level
	 */
	void create_direct_solver(size_t level);

	/**
	 * solves A[level]c = d with the direct solver,
	 * then calculates d -= A[level]c
	 * @param c output correction
	 * @param d defect. will be updated
	 * @param level
	 * @return true
	 */
	bool solve_on_base(vector_type &c, vector_type &d, size_t level);

	/**
	 * create fine marks for F-Smoothing
	 * @param level		current AMG level
	 * @param amgnodes 	for getting amgnodes[i].is_fine()
	 * @param N			number of elements in amgnodes with is_fine() == true.
	 */
	template<typename TAMGNodes>
	void create_fine_marks(int level, TAMGNodes &amgnodes, size_t N);

	/**
	 * calculate so that parentIndex[level+1][newIndex[i]] = i
	 * @param level		 current AMG level
	 * @param newIndex	 new indices of the nodes. -1 means node is fine.
	 * @param nrOfCoarse nr of Coarse nodes.
	 */
	void create_parent_index(int level, stdvector<int> newIndex, size_t nrOfCoarse);

	/**
	 * every coarse node c gets a unique new index newIndex[c]
	 * calculate PnewIndices(r, newIndex[c]) = PoldIndex(r, c)
	 *
	 * @param PoldIndices			prolongation matrix with old indices
	 * @param PnewIndices			prolongation matrix with new indices
	 * @param N						overlap 0 size (size of master+slave+inner nodes)
	 * @param amgnodes				amgnodes used for coarse/fine
	 * @param newIndex				new indices are stored here
	 * @param dEpsilonTruncation	connections in P with P(i,j) < dEpsilonTruction * max_k P(i,k) are dropped.
	 */
	template<typename TAMGNodes>
	void create_new_indices(prolongation_matrix_type &PoldIndices, prolongation_matrix_type &PnewIndices,
				size_t N, TAMGNodes &amgnodes, stdvector<int> &newIndex, double dEpsilonTruncation);
#ifndef UG_PARALLEL
	/**
	 * this function processes the prolongation in old indices and returns a prolongation with new indices,
	 * computing the new indices and and everything associated with the new indices (fine marks, parent indices etc.)
	 * in detail, this is
	 * - create the new indices by counting used coarse nodes and creation of PnewIndices @sa create_new_indices
	 * - create parent index (for debugging) @sa create_parent_index
	 * - create the coarser level of positions (for debugging) @sa make_coarse_level
	 * - create fine marks for f-smoothing @sa create_fine_marks
	 *
	 * @note this function is the serial version of @sa parallel_process_prolongation
	 *
	 * @param 	PoldIndices				used Prolongation matrix with old (fine grid) indices. in/out.
	 * @param	PnewIndices				Prolongation matrix with new indices.
	 * @param	dEpsilonTruncation		used to trucation the prolongation /sa create_new_indices
	 * @param	level					the amg level
	 * @param	amgnodes				coarse/fine ratings of the nodes
	 * @param	nextLevelMasterLayout	the master layout on the next coarser level
	 * @param	nextLevelSLaveLayout	the slave layout on the next coarser level
	 */
	template<typename TAMGNodes>
	void serial_process_prolongation(prolongation_matrix_type &PoldIndices, prolongation_matrix_type &PnewIndices, double dEpsilonTruncation, int level, TAMGNodes &amgnodes);
#else

	/**
	 * this function processes the prolongation in old indices and returns a prolongation with new indices,
	 * computing the new indices and and everything associated with the new indices (fine marks, parent indices etc.)
	 * in detail, this is
	 * - since we calculate the prolongation only in the master/inner nodes, we need to send this communication to the slave nodes (so we get additive AH=RAP again)
	 * 	 @sa communicate_prolongation
	 * - post set nodes which are interpolated from because of newly received prolongation rows @sa postset_coarse
	 * - create a minimal layout for the next level: we start with the total overlap layout we have in ParallelNodes PN, which also includes some node which might have
	 *   been added because of the sending of the prolongation rows, and then send a 'delete' to the master proc @sa create_minimal_layout_for_prolongation
	 * - create the new indices by counting used coarse nodes and creation of PnewIndices @sa create_new_indices
	 * - Replace Indices in the layout so we have our next level layouts with new indices @sa ReplaceIndicesInLayout
	 * - create parent index (for debugging) @sa create_parent_index
	 * - create the coarser level of positions (for debugging) @sa make_coarse_level
	 * - create fine marks for f-smoothing @sa create_fine_marks
	 *
	 * \note this function is the parallel version of @sa serial_process_prolongation
	 *
	 * @param 	PoldIndices				used Prolongation matrix with old (fine grid) indices. in/out.
	 * @param	PnewIndices				Prolongation matrix with new indices.
	 * @param	dEpsilonTruncation		used to trucation the prolongation /sa create_new_indices
	 * @param	level					the amg level
	 * @param	amgnodes				coarse/fine ratings of the nodes
	 * @param	nextLevelMasterLayout	the master layout on the next coarser level
	 * @param	nextLevelSLaveLayout	the slave layout on the next coarser level
	 */
	template<typename TAMGNodes>
	void parallel_process_prolongation(prolongation_matrix_type &PoldIndices, prolongation_matrix_type &PnewIndices, double dEpsilonTruncation, int level, TAMGNodes &amgnodes,
			ParallelNodes &PN, bool bCreateNewNodes, IndexLayout &nextLevelMasterLayout, IndexLayout &nextLevelSlaveLayout);

	/** since we calculate the prolongation only in the master/inner nodes, we need to send this communication to the slave nodes (so we get additive AH=RAP again)
	 *
	 * @param	PN				ParallelNodes. here we use master/slaveLayouts
	 * @param 	PoldIndices		used Prolongation matrix with old (fine grid) indices. in/out.
	 * @note this function can create new nodes in ParallelNodes. This is because for example a Slave node may be interpolated by a node which is in Overlap1.
	 */
	void communicate_prolongation(ParallelNodes &PN, prolongation_matrix_type &PoldIndices, bool bCreateNewNodes);

	/** post set nodes which are interpolated from because of newly received prolongation rows
	 *
	 * @param	PN				ParallelNodes. here we use slaveLayout
	 * @param 	PoldIndices		used Prolongation matrix with old (fine grid) indices. in.
	 * @param	amgnodes		struct to set fine/coarse.
	 */
	template<typename TAMGNodes>
	void postset_coarse(ParallelNodes &PN, prolongation_matrix_type &PoldIndices, TAMGNodes &nodes);

	/**
	 * P is a matrix which has nodes from ParallelNodes PN. We want to get a minimal layout for P.
	 * For this we first mark all slaves j for which there exist i with P(i, j) != 0.
	 * Since the P is P : Coarse -> Fine, this means we mark all slave elements we use on this processor.
	 * Then we create a slave layout out of it an communicate it to the master. This is necessary
	 * since the master cannot know which nodes we use.
	 *
	 * @param PN				Parallel Nodes structure
	 * @param P					used prolongation matrix
	 * @param newMasterLayout	new master layout for P, condensed from PN.get_total_master_layout();
	 * @param newSlaveLayout	new slave layout for P, condensed from PN.get_total_slave_layout();
	 */
	void create_minimal_layout_for_prolongation(ParallelNodes &PN, prolongation_matrix_type &P, IndexLayout &newMasterLayout, IndexLayout &newSlaveLayout);


#endif

#ifdef UG_PARALLEL
public:

	/**
	 * since we can turn off cores, we need to collect the data from them
	 * here we collect the data from vec into collectedVec.
	 * @param vec			uncollected vec
	 * @param collectedVec	collectedVec (result)
	 * @param level			level we operate on
	 * @param type			the type we want the result to be in (PST_ADDITIVE or PST_CONSISTENT)
	 */
	bool gather_vertical(vector_type &vec, vector_type &collectedVec, size_t level, ParallelStorageType type);

	/**
	 * in this function, we broadcast the data to cores we have shut off
	 * i.e. we write the data to vec
	 * @param vec			uncollected vec (result)
	 * @param collectedVec	collectedVec
	 * @param level			level we operate on
	 * @param type			the type we want the result to be in (PST_ADDITIVE or PST_CONSISTENT)
	 */
	bool broadcast_vertical(vector_type &vec, vector_type &collectedVec, size_t level, ParallelStorageType type);
#endif


public:
	/// debug output of interfaces written in ConnectionViewer-format
	void write_interfaces();

	bool check2(const vector_type &const_c, const vector_type &const_d);
	bool check_fsmoothing();

	/**
	 * handles agglomeration, and then calls
	 * add_correction_and_update_defect2 with collected vectors
	 * @param c			result: computed correction is added
	 * @param d			defect, result: updated defect.
	 * @param level		AMG level we operate on
	 * @return			true
	 */
	bool add_correction_and_update_defect(vector_type &c, vector_type &d, size_t level=0);
	/**
	 * adds the correction dc to c and updates the defect d with d -= A dc.
	 * uses the normal recursive AMG cycle for the coarse problem
	 * actual AMG cycle.
	 * @param c			result: computed correction is added
	 * @param d			defect, result: updated defect.
	 * @param A			matrix A to use
	 * @return			true
	 */
	bool add_correction_and_update_defect2(vector_type &c, vector_type &d,  matrix_operator_type &A, size_t level=0);

	/**
	 * for the purposes: like a normal AMG cycle, except that it will try to solve
	 * the coarse problem at exactLevel exactly by linear iteration
	 * @param c				result: computed correction is added
	 * @param d				defect, result: updated defect.
	 * @param level			current level
	 * @param exactLevel	level to solve problem exactly
	 * @return				true
	 */
	bool add_correction_and_update_defect(vector_type &c, vector_type &d, size_t level, size_t exactLevel);

	bool get_correction(vector_type &c, const vector_type &d);

	/// calculates vH[j] = v[m_parentIndex[level+1][j]];
	bool injection(vector_type &vH, const vector_type &v, size_t level);
/*
	size_t get_nr_of_coarse(size_t level)
	{
		assert(level+1 < m_usedLevels);
		return A[level+1]->length;
	}
*/
	/// returns the number of AMG levels
	size_t get_used_levels() const { return m_totalUsedLevels; }

	bool check_level(vector_type &c, vector_type &d, matrix_type &A, size_t level,
			checkResult &res, const vector_type *solution=NULL);
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

	void 	set_checkLevel_post_iterations(size_t i)				{ m_checkLevelPostIterations = i; }
	size_t	get_checkLevel_post_iterations() const					{ return m_checkLevelPostIterations; }


	void 	set_presmoother(SmartPtr<ILinearIterator<vector_type> > presmoother) {	m_presmoother = presmoother; }
	void 	set_postsmoother(SmartPtr<ILinearIterator<vector_type> > postsmoother) { m_postsmoother = postsmoother; }
	void 	set_base_solver(SmartPtr<ILinearOperatorInverse<vector_type> > basesolver) { m_basesolver = basesolver; }


	void set_position_provider(IPositionProvider<2> *ppp2d)
	{
		m_pPositionProvider2d = ppp2d;
	}

	void set_position_provider(IPositionProvider<3> *ppp3d)
	{
		m_pPositionProvider3d = ppp3d;
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

	void set_matrix_write_path(const char *path)
	{
		m_writeMatrixPath = path;
		m_writeMatrices = true;
	}


	virtual void tostring() const;

	/// if true, do not re-create AMG we preprocess is called
	void set_one_init(bool b)
	{
		m_bOneInit = b;
	}

	void set_nr_of_preiterations_at_check(size_t i)
	{
		m_iNrOfPreiterationsCheck = i;
	}

	void set_preiterations_mimum_defect_at_check(double d)
	{
		m_dPreiterationsMimumDefect = d;
	}

	void set_Y_cycle(int maxIterations, double dYreduce, double dYabs)
	{
		m_dYreduce = dYreduce;
		m_dYabs = dYabs;
		m_iYCycle = maxIterations;
	}

	/// @return c_A = total nnz of all matrices divided by nnz of matrix A
	double get_operator_complexity() const { return m_dOperatorComplexity; }

	/// @return c_G = total number of nodes of all levels divided by number of nodes on level 0
	double get_grid_complexity() const { return m_dGridComplexity; }

	/// @return the time spent on the whole setup in ms
	double get_timing_whole_setup_ms() const { return m_dTimingWholeSetupMS; }

	/// @return the time spent in the coarse solver setup in ms
	double get_timing_coarse_solver_setup_ms() const { return m_dTimingCoarseSolverSetupMS; }

	/// print level informations
	void print_level_information() const
	{
		std::cout.setf(std::ios::fixed);
		UG_LOG("Operator Complexity: " << get_operator_complexity() << ", grid complexity: " << get_grid_complexity() << ".\n");
		UG_LOG("Whole setup took " << get_timing_whole_setup_ms() << " ms, coarse solver setup took " << get_timing_coarse_solver_setup_ms() << " ms.\n");
		/*for(int i = 0; i<get_used_levels(); i++)
			UG_LOG(get_level_information(i)->tostring() << "\n");*/
	}

	const LevelInformation *get_level_information(size_t level) const
	{
		if(level < m_levelInformation.size())
			return &m_levelInformation[level];
		else return NULL;
	}

protected:
	void write_debug_matrices(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
			prolongation_matrix_type &P, size_t level);

	template<typename TMatrix>
	void write_debug_matrix(TMatrix &mat, size_t fromlevel, size_t tolevel, const char *name);

	template<typename TNodeType>
	void write_debug_matrix_markers(size_t level, const TNodeType &nodes);

	/// init f-smoothing (get consistent matrix)
	void init_fsmoothing();
	/// f-smoothing
	bool f_smoothing(vector_type &corr, vector_type &d, size_t level);


	/**
	 * writes the vector d in a connection-viewer-vec format
	 *
	 * @param 	filename
	 * @param	d
	 * @param	level
	 */
	bool writevec(std::string filename, const vector_type &d, size_t level, const vector_type *solution=NULL);

	/**
	 * reads in the positions from the position provider
	 * @sa set_position_provider, set_positions
	 */
	void update_positions();



// pure virtual functions
	virtual void create_AMG_level(matrix_type &AH, prolongation_matrix_type &R, const matrix_type &A,
				prolongation_matrix_type &P, size_t level) = 0;
	virtual void precalc_level(size_t level) = 0;


private:
#ifdef UG_PARALLEL
	/// starts the agglomeration process
	bool agglomerate(size_t level);

	/**
	 * agglomerate the vectors and matrix from all processors in mergeWith
	 * this processor becomes father for this level
	 * @param merge		this map maps neighbors which got agglomerated to their new cores
	 * @param mergeWith the processors to agglomerate with
	 * @param level
	 */
	bool agglomerate_from(const std::map<int, int> &merge, const std::vector<int> &mergeWith, size_t level, ParallelNodes &PN);

	/**
	 * agglomerates this processor on the level to mergeTo
	 * this processor becomes child for this level
	 * @param merge		this map maps neighbors which got agglomerated to their new cores
	 * @param mergeTo 	the processor to merge to
	 * @param level		the AMG level0
	 */
	bool agglomerate_to(const std::map<int, int> &merge, int mergeTo, size_t level, ParallelNodes &PN);
	/**
	 * even though this processor does not get agglomerated and does not agglomerate,
	 * he might need to change his interfaces because neighbors could have gotten
	 * agglomerated
	 * @param merge		this map maps neighbors which got agglomerated to their new cores
	 * @param level		the AMG level0
	 */
	bool agglomerate_adjust_interfaces(const std::map<int, int> &merge, size_t level, ParallelNodes &PN);
#endif

protected:
// data
	size_t 	m_numPreSmooth;						///< nu_1 : nr. of pre-smoothing steps
	size_t 	m_numPostSmooth;					///< nu_2: nr. of post-smoothing steps
	int 	m_cycleType;						///< gamma: cycle type (1 = V-Cycle, 2 = W-Cycle)

	size_t 	m_maxLevels;						///< max. nr of levels used for FAMG
	size_t	m_usedLevels;						///< nr of FAMG levels used
	size_t  m_totalUsedLevels;

	size_t 	m_maxNodesForBase;					///< max nr of coarse nodes before Base solver is used
	double 	m_dMaxFillBeforeBase;				///< max fill rate before Base solver is used
	size_t	m_minNodesOnOneProcessor;			///< min nr of nodes on one processor (for agglomeration)
	size_t	m_preferredNodesOnOneProcessor;		///< preferred nr of nodes on one processor (for agglomeration)
	bool 	m_bUseCollectedSolver;

	bool	m_bFSmoothing;
	bool	m_bOneInit;							///< if true, do not re-create AMG we preprocess is called

	size_t 	m_iNrOfPreiterationsCheck;			///< nr of mg cycles performed before checking
	double m_dPreiterationsMimumDefect;		///< minimum defect for preiterations before checking

	size_t m_checkLevelPostIterations;

	vector_type *m_vec4;						///< temporary Vector for defect (in get_correction)

	int iteration_glboal;


	bool 	m_writeMatrices;

	std::string m_writeMatrixPath;
	cAMG_helper m_amghelper;					///< helper struct for viewing matrices (optional)
	stdvector<MathVector<3> > m_dbgPositions;	///< positions of geometric grid (optional)
	int m_dbgDimension;							///< dimension of geometric grid (optional)


	SmartPtr<ILinearIterator<vector_type> > m_presmoother;	///< presmoother template
	SmartPtr<ILinearIterator<vector_type> > m_postsmoother;	///< postsmoother template \note: may be pre=post, is optimized away.

	SmartPtr<ILinearOperatorInverse<vector_type> > m_basesolver; ///< the base solver


	bool m_bInited;					///< true if inited. needed since preprocess doesnt give us a InterfaceCommunicator atm.
	double m_dOperatorComplexity;
	double m_dGridComplexity;
	double m_dTimingWholeSetupMS;
	double m_dTimingCoarseSolverSetupMS;

	size_t m_iYCycle;
	double m_dYreduce;
	double m_dYabs;

	IPositionProvider<2> *m_pPositionProvider2d;
	IPositionProvider<3> *m_pPositionProvider3d;

	stdvector<stdvector<int> > m_parentIndex;		///< parentIndex[i] is the index of i on the finer level
	stdvector<LevelInformation> m_levelInformation;
	struct AMGLevel
	{
		AMGLevel(int level)
			: pA(NULL),
			  presmoother(NULL),
			  postsmoother(NULL)
#ifdef UG_PARALLEL
			  ,collectedA(new matrix_operator_type)
#endif
		{
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

		stdvector<bool> is_fine;			///< fine marks for f-smoothing

		prolongation_matrix_type R; 		///< R Restriction Matrices
		prolongation_matrix_type P; 		///< P Prolongation Matrices
		SmartPtr<matrix_operator_type> pA;	///< A Matrices

		SmartPtr<ILinearIterator<vector_type> > presmoother;	///< presmoothers for each level
		SmartPtr<ILinearIterator<vector_type> > postsmoother;	///< postsmoothers for each level

#ifdef UG_PARALLEL
		pcl::InterfaceCommunicator<IndexLayout> com; 										///< the communicator object on this level
		IndexLayout slaveLayout, masterLayout;												///< layouts
		IndexLayout slaveLayoutAfterAgglomeration, masterLayoutAfterAgglomeration;			///< since layouts can change due to neighbors merging, these are the new layouts

		stdvector< typename block_traits<typename matrix_type::value_type>::inverse_type > m_diagInv;	///< inverse for diagonal, for f-smoothing

		// agglomeration
		bool bHasBeenMerged;								///< if true, this core is either father or child of a collected group
		// level 0 - m_agglomerateLevel
		pcl::ProcessCommunicator *pProcessCommunicator;		///< process Communicator on this level


		// level 0 - m_agglomerateLevel-1
		pcl::ProcessCommunicator agglomeratedPC;			///< process Communicator if level contains mergers
		IndexLayout agglomerateMasterLayout;

		vector_type collC, collD;
		SmartPtr<matrix_operator_type> collectedA;
		bool bLevelHasMergers;
#endif
		// in contrast to collC, collD, and collectedA, this can also point to other vectors.
		SmartPtr<matrix_operator_type> pAgglomeratedA;
	};

	stdvector<AMGLevel*> levels;

#ifdef UG_PARALLEL
	// on level m_agglomerateLevel
	IndexLayout agglomerateSlaveLayout;
	size_t m_agglomerateLevel;

	// stuff from old level 0 agglomeration
	pcl::ProcessCommunicator m_emptyPC;
	IndexLayout m_emptyLayout;

	/// return true if we have been merged to another core on this level
	bool isMergingSlave(size_t level);
	/// return true if we are the master for agglomeration on this level
	bool isMergingMaster(size_t level);
	/// return true if we are neither horizontal master nor slave on this level
	bool isNotMerging(size_t level);

#endif

	void calculate_level_information(size_t level, double createAMGlevelTiming);
};


///	@}

} // namespace ug



#include "amg_base_impl.h"


#endif // __H__UG__LIB_ALGEBRA__AMG_SOLVER__AMG_BASE_H__
