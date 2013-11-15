/*
 * mg_solver.h
 *
 *  Created on: 07.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER__
#define __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER__

// extern includes
#include <vector>
#include <iostream>

// other ug4 modules
#include "common/common.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallel_storage_type.h"
	#include "lib_grid/parallelization/distributed_grid.h"
#endif
#include "lib_algebra/operator/interface/linear_iterator.h"
#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/preconditioner/jacobi.h"
#include "lib_algebra/operator/linear_solver/lu.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/operator/linear_operator/transfer_interface.h"
//only for debugging!!!
#include "lib_grid/algorithms/debug_util.h"

// library intern headers
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// AssembledMultiGridCycle
////////////////////////////////////////////////////////////////////////////////

/// geometric multi grid preconditioner
/**
 * This class implements one step of the geometric multi grid as a
 * preconditioner for linear iteration schemes such as linear iteration, cg
 * or bicgstab.
 *
 * The coarse grid spaces are build up according to the Approximation Space
 * that is set from outside. In addition an Assembling routine must be
 * specified that is used to assemble the coarse grid matrices.
 *
 * \tparam		TApproximationSpace		Type of Approximation Space
 * \tparam		TAlgebra				Type of Algebra
 */
template <typename TDomain, typename TAlgebra>
class AssembledMultiGridCycle :
 public ILinearIterator<	typename TAlgebra::vector_type>
{
	public:
	///	Domain
		typedef TDomain domain_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	///	Grid Function type
		typedef GridFunction<TDomain, TAlgebra> GF;

	///////////////////////////////////////////////////////////////////////////
	//	Setup
	///////////////////////////////////////////////////////////////////////////

	public:
	/// constructor setting approximation space
		AssembledMultiGridCycle(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	/// sets the assembling procedure that is used to compute coarse grid matrices
		void set_discretization(SmartPtr<IAssemble<TAlgebra> > spAss)
			{m_spAss = spAss;}

	///	sets the level where exact solving is performed in the mg cycle
		void set_base_level(int baseLevel) {m_baseLev = baseLevel;}

	///	sets the surface level (default is top-surface)
		void set_surface_level(int topLevel) {m_surfaceLev = topLevel;}

	///	sets the base solver that is used
		void set_base_solver(SmartPtr<ILinearOperatorInverse<vector_type> > baseSolver)
			{m_spBaseSolver = baseSolver;}

	///	sets if the base solver is applied in parallel
		void set_parallel_base_solver(bool bParallel) {m_bParallelBaseSolverIfAmbiguous = bParallel;}

	///	sets the cycle type (1 = V-cycle, 2 = W-cycle, ...)
		void set_cycle_type(int type) {m_cycleType = type;}

	///	sets the number of pre-smoothing steps to be performed
		void set_num_presmooth(int num) {m_numPreSmooth = num;}

	///	sets the number of post-smoothing steps to be performed
		void set_num_postsmooth(int num) {m_numPostSmooth = num;}

	///	sets the smoother that is used
		void set_smoother(SmartPtr<ILinearIterator<vector_type> > smoother)
			{set_presmoother(smoother); set_postsmoother(smoother);}

	///	sets the pre-smoother that is used
		void set_presmoother(SmartPtr<ILinearIterator<vector_type> > smoother)
			{m_spPreSmootherPrototype = smoother;}

	///	sets the post-smoother that is used
		void set_postsmoother(SmartPtr<ILinearIterator<vector_type> > smoother)
			{m_spPostSmootherPrototype = smoother;}

	///	sets the transfer operator
		void set_transfer(SmartPtr<ITransferOperator<TAlgebra> > P)
			{set_prolongation(P); set_restriction(P);}

	///	sets the prolongation operator
		void set_prolongation(SmartPtr<ITransferOperator<TAlgebra> > P)
			{m_spProlongationPrototype = P;}

	///	sets the restriction operator
		void set_restriction(SmartPtr<ITransferOperator<TAlgebra> > P)
			{m_spRestrictionPrototype = P;}

	///	sets the projection operator
		void set_projection(SmartPtr<ITransferOperator<TAlgebra> > P)
			{m_spProjectionPrototype = P;}

	///	clears all transfer post process
		void clear_transfer_post_process()
			{m_vspProlongationPostProcess.clear(); m_vspRestrictionPostProcess.clear();}

	///	add prolongation post process
		void add_prolongation_post_process(SmartPtr<ITransferPostProcess<TAlgebra> > PP)
			{m_vspProlongationPostProcess.push_back(PP);}

	///	add restriction post process
		void add_restriction_post_process(SmartPtr<ITransferPostProcess<TAlgebra> > PP)
			{m_vspRestrictionPostProcess.push_back(PP);}

	///////////////////////////////////////////////////////////////////////////
	//	Linear Solver interface methods
	///////////////////////////////////////////////////////////////////////////

	///	name
		virtual const char* name() const {return "Geometric MultiGrid";}

	///	returns information about configuration parameters
		virtual std::string config_string() const;

	/// Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

	///	Prepare for Linear Operartor L
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L);

	///	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d);

	///	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d);

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone();

	///	Destructor
		~AssembledMultiGridCycle();

 	protected:
 	/// compute correction on level and update defect
		void lmgc(size_t lev);

	////////////////////////////////////////////////////////////////
	//	The methods in this section rely on each other and should be called in sequence
	///	performs presmoothing on the given level
		void presmooth(size_t lev);

	///	performs restriction on to the level below
		void restriction(size_t lev);

	///	performs prolongation to the level above
		void prolongation(size_t lev);

	///	performs postsmoothin
		void postsmooth(size_t lev);
	//	end of section
	////////////////////////////////////////////////////////////////

	///	compute base solver
		void base_solve(size_t lev);

	/// performs smoothing on level l, nu times
		void smooth(SmartPtr<GF> sc, SmartPtr<GF> sd, SmartPtr<GF> st,
		            MatrixOperator<matrix_type, vector_type>& A,
		            ILinearIterator<vector_type>& S, size_t lev, int nu);

	///	returns the number of allocated levels
		size_t num_levels() const {return m_vLevData.size();}

	///	allocates the memory
		void init_level_memory(int baseLev, int topLev);

	///	initializes common part
		void init();

	///	initializes the smoother and base solver
		void init_smoother();

	///	initializes the coarse grid matrices
		void init_level_operator();

	///	initializes the smoother and base solver
		void init_base_solver();

	///	initializes the prolongation
		void init_transfer();

	///	initializes the prolongation
		void init_projection();

	///	assembles the missing matrix part on the coarse level, that must be
	///	added if the correction has been computed to ensure a correctly updated
	///	defect. (i.e. assembles A^c, with d^f -= A^c * c^c)
		void init_missing_coarse_grid_coupling(const vector_type* u);

	protected:
	/// operator to invert (surface grid)
		ConstSmartPtr<matrix_type> m_spSurfaceMat;

	///	Solution on surface grid
		const vector_type* m_pSurfaceSol;

	///	assembling routine for coarse grid matrices
		SmartPtr<IAssemble<TAlgebra> > m_spAss;

	///	approximation space for level and surface grid
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	top level (i.e. highest level in hierarchy. This is the surface level
	///	in case of non-adaptive refinement)
		int m_topLev;
		int m_surfaceLev;

	///	base level (where exact inverse is computed)
		int m_baseLev;

	///	cylce type (1 = V-cycle, 2 = W-cylcle, ...)
		int m_cycleType;

	///	number of Presmooth steps
		int m_numPreSmooth;

	///	number of Postsmooth steps
		int m_numPostSmooth;

	///	flag indicating if grid is full refined
		bool m_bAdaptive;

	///	Structure used to realize Surface to Level mapping
	/// \{
		struct LevelIndex{
			LevelIndex() : index(-1), level(-1) {}
			LevelIndex(size_t index_, int level_) : index(index_), level(level_) {}
			size_t index;
			int level;
		};
		std::vector<LevelIndex> m_vSurfToLevelMap;
		template <typename TElem>
		void init_surface_to_level_mapping();
		void init_surface_to_level_mapping();
	/// \}

	///	init mapping from noghost -> w/ ghost
	/// \{
		template <typename TElem>
		void init_noghost_to_ghost_mapping(	std::vector<size_t>& vNoGhostToGhostMap,
		                                   	ConstSmartPtr<DoFDistribution> spNoGhostDD,
		                                   	ConstSmartPtr<DoFDistribution> spGhostDD);
		void init_noghost_to_ghost_mapping(int lev);
	/// \}

	///	collects all shadowing indices into the std::vector
	/// \{
		template <typename TElem>
		void collect_shadowing_indices(std::vector<size_t>& vShadowing,
		                               ConstSmartPtr<DoFDistribution> spDD,
		                               std::vector<size_t>& vSurfShadowing,
		                               ConstSmartPtr<DoFDistribution> spSurfDD);
		void collect_shadowing_indices(int lev);
	/// \}

	///	prototype for pre-smoother
		SmartPtr<ILinearIterator<vector_type> > m_spPreSmootherPrototype;

	///	prototype for post-smoother
		SmartPtr<ILinearIterator<vector_type> > m_spPostSmootherPrototype;

	///	prototype for projection operator
		SmartPtr<ITransferOperator<TAlgebra> > m_spProjectionPrototype;

	///	prototype for prolongation operator
		SmartPtr<ITransferOperator<TAlgebra> > m_spProlongationPrototype;

	///	prototype for prolongation operator
		SmartPtr<ITransferOperator<TAlgebra> > m_spRestrictionPrototype;

	///	prototpe for transfer post process
		std::vector<SmartPtr<ITransferPostProcess<TAlgebra> > > m_vspProlongationPostProcess;
		std::vector<SmartPtr<ITransferPostProcess<TAlgebra> > > m_vspRestrictionPostProcess;

	///	base solver for the coarse problem
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spBaseSolver;

		////////////////////////////////////
		// Storage for each grid level
		////////////////////////////////////

		struct LevData
		{
		///	Level matrix operator
			SmartPtr<MatrixOperator<matrix_type, vector_type> > A;

		///	Smoother
			SmartPtr<ILinearIterator<vector_type> > PreSmoother;
			SmartPtr<ILinearIterator<vector_type> > PostSmoother;

		///	Projection operator
			SmartPtr<ITransferOperator<TAlgebra> > Projection;

		///	Transfer operator
			SmartPtr<ITransferOperator<TAlgebra> > Prolongation;
			SmartPtr<ITransferOperator<TAlgebra> > Restriction;

		///	transfer post process
            std::vector<SmartPtr<ITransferPostProcess<TAlgebra> > > vProlongationPP;
            std::vector<SmartPtr<ITransferPostProcess<TAlgebra> > > vRestrictionPP;

		///	vectors needed (sx = no-ghosts [for smoothing], t = for transfer)
			SmartPtr<GF> sc, sd, st, t;

		///	maps global indices (including ghosts) to patch indices (no ghosts included).
			std::vector<size_t> vMapPatchToGlobal;

		/// list of shadowing indices
			std::vector<size_t> vShadowing;

		/// list of corresponding surface index to shadowing indices
			std::vector<size_t> vSurfShadowing;

		///	missing coarse grid correction
			matrix_type CoarseGridContribution;
		};

	///	flag, if to solve base problem in parallel when gathered and (!) parallel possible
		bool m_bParallelBaseSolverIfAmbiguous;

	///	flag if using parallel base solver
		bool m_bGatheredBaseUsed;

	///	Matrix for gathered base solver
		SmartPtr<MatrixOperator<matrix_type, vector_type> > spBaseSolverMat;

	///	vector for gathered base solver
		SmartPtr<GF> spGatheredBaseCorr;

	///	storage for all level
		std::vector<SmartPtr<LevData> > m_vLevData;

	///	copies vector to smoothing patch using cached mapping
		void copy_ghost_to_noghost(SmartPtr<GF> spVecTo,
		                           ConstSmartPtr<GF> spVecFrom,
		                           const std::vector<size_t>& vMapPatchToGlobal);

	///	copies vector from smoothing patch using cached mapping
		void copy_noghost_to_ghost(SmartPtr<GF> spVecTo,
								   ConstSmartPtr<GF> spVecFrom,
								   const std::vector<size_t>& vMapPatchToGlobal);

	/// gathers the vector using vertical interfaces. Entries are summed at vmasters.
		void gather_vertical(vector_type& d);

	/// broadcasts the vector using vertical interfaces.
		void copy_to_vertical_slaves(vector_type& t);

	///	copies values from h-masters to h-slaves
		void copy_to_vertical_masters(vector_type& c);

#ifdef UG_PARALLEL
	/// communicator
		pcl::InterfaceCommunicator<IndexLayout> m_Com;
#endif

	public:
	///	set debug output
	/**
	 * If a DebugWriter is passed by this method, the multi grid cycle writes
	 * the level/surface vectors and matrices for debug purposes.
	 *
	 * \param[in]	debugWriter		Debug Writer to use
	 */
		void set_debug(SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra> > spDebugWriter)
		{
			m_spDebugWriter = spDebugWriter;
		}

	protected:
	///	writes debug output for a level vector only on smooth path
	/**
	 * This method writes the level vector to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		spGF		Level Vector to write for debug purpose
	 * \param[in]		name		Filename
	 */
		void write_debug(ConstSmartPtr<GF> spGF, std::string name);
		void write_debug(const GF& rGF, std::string name);

	///	writes debug output for a level matrix only on smooth path
	/**
	 * This method writes the level matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		name		Filename
	 * \param[in]		level		grid level corresponding to the matrix
	 */
		void write_smooth_level_debug(const matrix_type& mat, std::string name, size_t lev);

	///	writes debug output for a level matrix
	/**
	 * This method writes the level matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		name		Filename
	 * \param[in]		level		grid level corresponding to the matrix
	 */
		void write_level_debug(const matrix_type& mat, std::string name, size_t lev);

	///	writes debug output for a surface matrix
	/**
	 * This method writes the surface matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		name		Filename
	 */
		void write_surface_debug(const matrix_type& mat, std::string name);

	///	logs a level-data-struct to the terminal
		void log_debug_data(int lvl, std::string name);

	///	Debug Writer
		SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra> > m_spDebugWriter;

	///	counter for debug, to distinguish the iterations
		int m_dbgIterCnt;
};

////////////////////////////////////////////////////////////////////////////////
// Selections
////////////////////////////////////////////////////////////////////////////////

/// selects all non-shadows, that are adjacent to a shadow in the multigrid
void SelectNonShadowsAdjacentToShadows(BoolMarker& sel, const SurfaceView& surfView);

/// selects all non-shadows, that are adjacent to a shadow on a grid levels
void SelectNonShadowsAdjacentToShadowsOnLevel(BoolMarker& sel,
										   const SurfaceView& surfView,
										   int level);

} // end namespace ug

// include implementation
#include "mg_solver_impl.hpp"

#endif /* __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER__ */
