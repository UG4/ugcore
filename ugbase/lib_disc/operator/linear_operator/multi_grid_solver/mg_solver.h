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
#include "lib_disc/dof_manager/mg_dof_distribution.h"
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
		void set_parallel_base_solver(bool bParallel) {m_bBaseParallel = bParallel;}

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
		void smooth(vector_type& c, vector_type& d, vector_type& t,
		            MatrixOperator<matrix_type, vector_type>& A,
		            ILinearIterator<vector_type>& S, size_t lev, int nu);

	///	returns the number of allocated levels
		size_t num_levels() const {return m_vLevData.size();}

	///	allocates the memory
		void top_level_required(size_t topLevel);

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
		SmartPtr<matrix_type> m_spSurfaceMat;

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

	///	flag, if to solve base problem in parallel
		bool m_bBaseParallel;

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
		//	constructor
			LevData();

		//	prepares the grid level data for appication
			void update(size_t lev,
			            SmartPtr<DoFDistribution> levelDD,
			            SmartPtr<ApproximationSpace<TDomain> > approxSpace,
			            SmartPtr<IAssemble<TAlgebra> > spAss,
			            ILinearIterator<vector_type>& presmoother,
			            ILinearIterator<vector_type>& postsmoother,
			            ITransferOperator<TAlgebra>& projection,
			            ITransferOperator<TAlgebra>& prolongation,
			            ITransferOperator<TAlgebra>& restriction,
			            std::vector<SmartPtr<ITransferPostProcess<TAlgebra> > >& vprolongationPP,
			            std::vector<SmartPtr<ITransferPostProcess<TAlgebra> > >& vrestrictionPP,
			            BoolMarker& nonGhostMarker);

		//	returns if ghosts are present on the level
			bool has_ghosts() const {return num_smooth_indices() != num_indices();}

		//	number of smoothing indices
			size_t num_smooth_indices() const {return m_numSmoothIndices;}

		//	number of indices on whole level
			size_t num_indices() const
				{UG_ASSERT(spLevDD.valid(), "Missing LevDD"); return spLevDD->num_indices();}

		//	returns the smoothing matrix (depends if smooth patch needed or not)
			SmartPtr<MatrixOperator<matrix_type, vector_type> >
			get_smooth_mat()
			{
				if(has_ghosts()) return spSmoothMat;
				else return spLevMat;
			}

		//	returns the vectors used for smoothing (patch only vectors)
			vector_type& get_smooth_solution() {if(has_ghosts()) return *su; else return *u;}
			vector_type& get_smooth_defect() {if(has_ghosts()) return *sd; else return *d;}
			vector_type& get_smooth_correction(){if(has_ghosts()) return *sc; else return *c;}
			vector_type& get_smooth_tmp(){if(has_ghosts()) return *st; else return *t;}

		//	copies values of defect to smoothing patch
			void copy_defect_to_smooth_patch()
			{
				#ifdef UG_PARALLEL
				if(has_ghosts()) {
					for(size_t i = 0; i < vMapPatchToGlobal.size(); ++i) (*sd)[i] = (*d)[ vMapPatchToGlobal[i] ];
					sd->set_storage_type(d->get_storage_mask());
				}
				#endif
			}

		//	copies values of tmp to smoothing patch
			void copy_tmp_to_smooth_patch()
			{
				#ifdef UG_PARALLEL
				if(has_ghosts()) {
					for(size_t i = 0; i < vMapPatchToGlobal.size(); ++i) (*st)[i] = (*t)[ vMapPatchToGlobal[i] ];
					st->set_storage_type(t->get_storage_mask());
				}
				#endif
			}

		//	copies values of defect from smoothing patch
			void copy_defect_from_smooth_patch(bool clearGhosts = false)
			{
				#ifdef UG_PARALLEL
				if(has_ghosts()) {
					if(clearGhosts) d->set(0.0);
					for(size_t i = 0; i < vMapPatchToGlobal.size(); ++i) (*d)[ vMapPatchToGlobal[i] ] = (*sd)[i];
					d->set_storage_type(sd->get_storage_mask());
				}
				#endif
			}

		//	copies values of defect to smoothing patch
			void copy_correction_from_smooth_patch(bool clearGhosts = false)
			{
				#ifdef UG_PARALLEL
				if(has_ghosts()) {
					if(clearGhosts) c->set(0.0);
					for(size_t i = 0; i < vMapPatchToGlobal.size(); ++i) (*c)[ vMapPatchToGlobal[i] ] = (*sc)[i];
					c->set_storage_type(sc->get_storage_mask());
				}
				#endif
			}

		//	destructor
			~LevData();

			public:
		//	level DoF Distribution
			SmartPtr<DoFDistribution> spLevDD;

		//	Approximation Space
			SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

		//	matrix operator for whole grid level
			SmartPtr<MatrixOperator<matrix_type, vector_type> > spLevMat;

		//	matrix for smoothing on smoothing patch of grid level
			SmartPtr<MatrixOperator<matrix_type, vector_type> > spSmoothMat;

		//	smoother
			SmartPtr<ILinearIterator<vector_type> > PreSmoother;
			SmartPtr<ILinearIterator<vector_type> > PostSmoother;

		//	projection operator
			SmartPtr<ITransferOperator<TAlgebra> > Projection;

		//	transfer operator
			SmartPtr<ITransferOperator<TAlgebra> > Prolongation;
			SmartPtr<ITransferOperator<TAlgebra> > Restriction;

		//	transfer post process
            std::vector<SmartPtr<ITransferPostProcess<TAlgebra> > > vProlongationPP;
            std::vector<SmartPtr<ITransferPostProcess<TAlgebra> > > vRestrictionPP;

		//	vectors needed
			SmartPtr<GridFunction<TDomain, TAlgebra> > u, c, d, t;

		//	vectors needed for smoothing
			SmartPtr<GridFunction<TDomain, TAlgebra> > su, sc, sd, st;

		//	missing coarse grid correction
			matrix_type CoarseGridContribution;

		///	maps global indices (including ghosts) to patch indices (no ghosts included).
		/**	maps are only filled if ghosts are present on the given level.
		 * \{ */
			std::vector<size_t> vMapPatchToGlobal;
			std::vector<int> vMapGlobalToPatch;
		/** \} */

		//	number of smooth indices
			size_t m_numSmoothIndices;

#ifdef UG_PARALLEL
		//	interfaces needed for smoothing
			SmartPtr<AlgebraLayouts> spSmoothLayouts;
#endif
		};

	///	storage for all level
		std::vector<SmartPtr<LevData> > m_vLevData;

	///	bool marker of non-ghosts
		BoolMarker m_NonGhostMarker;

		std::vector<vector_type*> level_defects()
		{
			std::vector<vector_type*> vVec;
			for(size_t i = 0; i < m_vLevData.size(); ++i)
			{
				if(m_vLevData[i]->num_smooth_indices() > 0)
					vVec.push_back(m_vLevData[i]->d.get());
				else vVec.push_back(NULL);
			}
			return vVec;
		}

#ifdef UG_PARALLEL
	/**
	 *	gathers the vector using vertical interfaces.
	 *	Entries are summed at vmasters.
	 */
		void gather_vertical(vector_type& d);

	/**
	 *	gathers the vector using vertical interfaces.
	 *	Entries are copied from vslaves to vmasters
	 */
		void gather_vertical_copy(vector_type& d);

	/**
	 *	gathers the vector using vertical interfaces.
	 *	Only ghost values are adjusted. The specified tmp vector is for internal
	 *	calculations only and has to be of the same size as d!
	 *	mapGlobalToPatch either has to be empty or of the same size as d.
	 */
		void gather_on_ghosts(vector_type& d, vector_type& tmp,
					std::vector<int>& mapGlobalToPatch);

	/**
	 *	broadcasts the vector using vertical interfaces.
	 */
		void broadcast_vertical(vector_type& t);

	/**
	 *	broadcasts and adds the vector using vertical interfaces.
	 */
		void broadcast_vertical_add(vector_type& t);

	///	copies values from h-masters to h-slaves
		void copy_to_horizontal_slaves(vector_type& c);

	///	copies values from h-masters to h-slaves
		void copy_to_vertical_masters(vector_type& c);

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
		void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter)
		{
			m_spDebugWriter = spDebugWriter;
		}

	protected:
	///	writes debug output for a level vector only on smooth path
	/**
	 * This method writes the level vector to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		vec			Level Vector to write for debug purpose
	 * \param[in]		filename	Filename
	 * \param[in]		level		grid level corresponding to the vector
	 */
		void write_smooth_level_debug(const vector_type& vec, const char* filename, size_t lev);

	///	writes debug output for a level matrix only on smooth path
	/**
	 * This method writes the level matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		filename	Filename
	 * \param[in]		level		grid level corresponding to the matrix
	 */
		void write_smooth_level_debug(const matrix_type& mat, const char* filename, size_t lev);

	///	writes debug output for a level vector
	/**
	 * This method writes the level vector to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		vec			Level Vector to write for debug purpose
	 * \param[in]		filename	Filename
	 * \param[in]		level		grid level corresponding to the vector
	 */
		void write_level_debug(const vector_type& vec, const char* filename, size_t lev);

	///	writes debug output for a level matrix
	/**
	 * This method writes the level matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		filename	Filename
	 * \param[in]		level		grid level corresponding to the matrix
	 */
		void write_level_debug(const matrix_type& mat, const char* filename, size_t lev);

	///	writes debug output for a surface vector
	/**
	 * This method writes the surface vector to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		vec			Level Vector to write for debug purpose
	 * \param[in]		filename	Filename
	 */
		void write_surface_debug(const vector_type& vec, const char* filename);

	///	writes debug output for a surface matrix
	/**
	 * This method writes the surface matrix to a debug file, if a debug writer
	 * has been set.
	 *
	 * \param[in]		mat			Level Matrix to write for debug purpose
	 * \param[in]		filename	Filename
	 */
		void write_surface_debug(const matrix_type& mat, const char* filename);

	///	logs a level-data-struct to the terminal
		void log_level_data(size_t lvl);

	///	Debug Writer
		SmartPtr<IDebugWriter<algebra_type> > m_spDebugWriter;

	///	counter for debug, to distinguish the iterations
		int m_dbgIterCnt;
};

////////////////////////////////////////////////////////////////////////////////
// Operation on Shadows/Shadowing
////////////////////////////////////////////////////////////////////////////////

/**
* This functions adds the shadow values from a coarser grid to the shadowing
* DoFs on the finer grid.
*
* \param[out]	fineVec			fine grid vector
* \param[out]	coarseVec		coarse grid vector
* \param[in]	scale			scaling, when adding
* \param[in] 	ddFine			dof distribution on fine space
* \param[in] 	ddCoarse		dof distribution on coarse space
* \param[in]	surfView		surface view
*/
template <typename TBaseElem, typename TVector>
void AddProjectionOfShadows(const std::vector<TVector*>& vFineVector,
						 std::vector<ConstSmartPtr<DoFDistribution> > vDDFine,
						 const TVector& coarseVec,
						 ConstSmartPtr<DoFDistribution> ddCoarse,
						 const int level,
						 const number scale,
						 const SurfaceView& surfView)
{
	PROFILE_FUNC_GROUP("gmg");
	std::vector<size_t> fineInd, coarseInd;

	// 	iterators
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator const_iterator;
	const_iterator iter, iterEnd;

	// 	loop subsets of fine level
	for(int si = 0; si < ddCoarse->num_subsets(); ++si)
	{
		iter = ddCoarse->begin<TBaseElem>(si);
		iterEnd = ddCoarse->end<TBaseElem>(si);

	// 	loop elements of coarse subset
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			TBaseElem* pElem = *iter;

		// 	get child (i.e. shadow)
			TBaseElem* pShadowing = surfView.child_if_copy(pElem);

		//	offset to count which child currently handling
			int offset = 0;

		// 	get global indices
			ddCoarse->inner_algebra_indices(pElem, coarseInd);

		//	skip if not a copy
			while(pShadowing){
			//	increase offset
				++offset;

			// 	get global indices
				vDDFine[level+offset]->inner_algebra_indices(pShadowing, fineInd);

			//	add coarse vector entries to fine vector entries
				for(size_t i = 0; i < coarseInd.size(); ++i)
				{
					VecScaleAdd((*vFineVector[level+offset])[fineInd[i]],
								1.0, (*vFineVector[level+offset])[fineInd[i]],
								scale, coarseVec[coarseInd[i]]);
				}

			//	next child
				pElem = pShadowing;
				pShadowing = surfView.child_if_copy(pElem);
			}
		}
	}
}

template <typename TVector>
void AddProjectionOfShadows(const std::vector<TVector*>& vFineVector,
						 std::vector<ConstSmartPtr<DoFDistribution> > vDDFine,
						 const TVector& coarseVec,
						 ConstSmartPtr<DoFDistribution> ddCoarse,
						 const int level,
						 const number scale,
						 const SurfaceView& surfView)
{
	PROFILE_FUNC_GROUP("gmg");
	//	forward for all BaseObject types
	if(ddCoarse->max_dofs(VERTEX))
		AddProjectionOfShadows<VertexBase, TVector>
					(vFineVector, vDDFine, coarseVec, ddCoarse, level, scale, surfView);
	if(ddCoarse->max_dofs(EDGE))
		AddProjectionOfShadows<EdgeBase, TVector>
					(vFineVector, vDDFine, coarseVec, ddCoarse, level, scale, surfView);
	if(ddCoarse->max_dofs(FACE))
		AddProjectionOfShadows<Face, TVector>
					(vFineVector, vDDFine, coarseVec, ddCoarse, level, scale, surfView);
	if(ddCoarse->max_dofs(VOLUME))
		AddProjectionOfShadows<Volume, TVector>
					(vFineVector, vDDFine, coarseVec, ddCoarse, level, scale, surfView);
}


/**
* This functions sets the values of a vector to zero, where the index
* corresponds to a refine-patch boundary (i.e. the geomeric object is a
* shadowing object) for an element type
*
* \param[out]	vec					grid vector
* \param[in] 	dd					DoFDistribution
* \param[in]	surfView			SurfaceView
* \param[in]	pmapGlobalToPatch	(optional) mapping of global indices to patch indices
*/
template <typename TBaseElem, typename TVector>
void SetZeroOnShadowing(TVector& vec,
					 ConstSmartPtr<DoFDistribution> dd,
					 const SurfaceView& surfView,
					 const std::vector<int>* pmapGlobalToPatch = NULL)
{
	PROFILE_FUNC_GROUP("gmg");
	//	indices
	std::vector<size_t> ind;

	// 	Vertex iterators
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator const_iterator;
	const_iterator iter, iterEnd;

	// 	loop subsets of fine level
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
		iter = dd->begin<TBaseElem>(si);
		iterEnd = dd->end<TBaseElem>(si);

	// 	loop nodes of fine subset
		if(pmapGlobalToPatch){
			const std::vector<int>& mapGlobalToPatch = *pmapGlobalToPatch;

			for(; iter != iterEnd; ++iter)
			{
			//	get vertex
				TBaseElem* vrt = *iter;

				if(!surfView.is_shadowing(vrt))
					continue;

			// 	get global indices
				dd->inner_algebra_indices(vrt, ind);

			//	set vector entries to zero
				for(size_t i = 0; i < ind.size(); ++i)
				{
					UG_ASSERT(ind[i] < mapGlobalToPatch.size(),
							 "mapGlobalToPatch is too small on level " << dd->grid_level() << "."
							 << "size: " << mapGlobalToPatch.size() << ", "
							 << "index: " << ind[i]
							 << ", at " << GetGeometricObjectCenter(
									 *const_cast<MultiGrid*>(surfView.subset_handler()->multi_grid()), vrt)
							 << ", on level " << surfView.subset_handler()->multi_grid()->get_level(vrt));
					UG_ASSERT((mapGlobalToPatch[ind[i]] >= 0) && (mapGlobalToPatch[ind[i]] < (int)vec.size()),
							  "Some problem with mapGlobalToPatch... probably trying to set a ghost to zero? "
							  << "mapGlobalToPatch[ind[i]] = " << mapGlobalToPatch[ind[i]]
							  << ", num patch indices " << vec.size()
							  << ", at " << GetGeometricObjectCenter(
									  *const_cast<MultiGrid*>(surfView.subset_handler()->multi_grid()), vrt)
							  << ", on level " << surfView.subset_handler()->multi_grid()->get_level(vrt));
					vec[mapGlobalToPatch[ind[i]]] = 0.0;
				}
			}
		}
		else{
			for(; iter != iterEnd; ++iter)
			{
			//	get vertex
				TBaseElem* vrt = *iter;

				if(!surfView.is_shadowing(vrt))
					continue;

			// 	get global indices
				dd->inner_algebra_indices(vrt, ind);

			//	set vector entries to zero
				for(size_t i = 0; i < ind.size(); ++i)				{
					vec[ind[i]] = 0.0;
				}
			}
		}
	}
}

/**
* This functions sets the values of a vector to zero, where the index
* corresponds to a refine-patch boundary (i.e. the geomeric object is a
* shadowing object)
*
* \param[out]	vec				grid vector
* \param[in] 	dd				DoFDistribution
* \param[in]	surfView		SurfaceView
*/
template <typename TVector>
void SetZeroOnShadowing(TVector& vec,
					 ConstSmartPtr<DoFDistribution> dd,
					 const SurfaceView& surfView,
					 const std::vector<int>* pmapGlobalToPatch = NULL)
{
	//	forward for all BaseObject types
	if(dd->max_dofs(VERTEX))
		SetZeroOnShadowing<VertexBase, TVector>(vec, dd, surfView, pmapGlobalToPatch);
	if(dd->max_dofs(EDGE))
		SetZeroOnShadowing<EdgeBase, TVector>(vec, dd, surfView, pmapGlobalToPatch);
	if(dd->max_dofs(FACE))
		SetZeroOnShadowing<Face, TVector>(vec, dd, surfView, pmapGlobalToPatch);
	if(dd->max_dofs(VOLUME))
		SetZeroOnShadowing<Volume, TVector>(vec, dd, surfView, pmapGlobalToPatch);
}

////////////////////////////////////////////////////////////////////////////////
// Selections
////////////////////////////////////////////////////////////////////////////////

/// selects all non-shadows, that are adjacent to a shadow in the multigrid
void SelectNonShadowsAdjacentToShadows(BoolMarker& sel, const SurfaceView& surfView);

/// selects all non-shadows, that are adjacent to a shadow on a grid levels
void SelectNonShadowsAdjacentToShadowsOnLevel(BoolMarker& sel,
										   const SurfaceView& surfView,
										   int level);

#ifdef UG_PARALLEL
template <typename TElemBase, typename TIter>
void SelectNonGhosts(BoolMarker& sel,
				  DistributedGridManager& dstGrMgr,
				  TIter iter,
				  TIter iterEnd)
{
	//	loop all base elems
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TElemBase* elem = *iter;

	//	select ghosts
		if(!dstGrMgr.is_ghost(elem)) sel.mark(elem);
	}
}
#endif

////////////////////////////////////////////////////////////////////////////////
// Matrix Copy operations
////////////////////////////////////////////////////////////////////////////////

/// copies a matrix from a larger one into a smaller one
/**
* This function copies a matrix of a larger index set into a matrix with a
* smaller (or equally sized) index set. The copying is performed using a
* mapping between the index set, that returns smallIndex = vMap[largeIndex],
* and a -1 if the largeIndex is dropped.
*/
template <typename TMatrix>
void CopyMatrixByMapping(TMatrix& smallMat,
					  const std::vector<int>& vMap,
					  const TMatrix& origMat)
{
	PROFILE_FUNC_GROUP("gmg");
	//	check size
	UG_ASSERT(vMap.size() == origMat.num_rows(), "Size must match.");
	UG_ASSERT(vMap.size() == origMat.num_cols(), "Size must match.");

	//	type of matrix row iterator
	typedef typename TMatrix::const_row_iterator const_row_iterator;

	//	loop all mapped indices
	for(size_t origInd = 0; origInd < vMap.size(); ++origInd)
	{
	//	get mapped level index
		const int smallInd = vMap[origInd];

	//	skipped non-mapped indices (indicated by -1)
		if(smallInd < 0) continue;

	//	loop all connections of the surface dof to other surface dofs and copy
	//	the matrix coupling into the level matrix

		for(const_row_iterator conn = origMat.begin_row(origInd);
								conn != origMat.end_row(origInd); ++conn)
		{
		//	get corresponding level connection index
			const int smallConnIndex = vMap[conn.index()];

		//	check that index is in small matrix, too
			if(smallConnIndex < 0) continue;

		//	copy connection to smaller matrix
			smallMat(smallInd, smallConnIndex) = conn.value();
		}
	}

	#ifdef UG_PARALLEL
	smallMat.set_storage_type(origMat.get_storage_mask());
	#endif
}

/// copies a matrix into a equally sized second one using a mapping
/**
* This function copies a matrix of a index set into a matrix with a
* equally sized index set. The copying is performed using a
* mapping between the index set, that returns newIndex = vMap[origIndex].
*/
template <typename TMatrix>
void CopyMatrixByMapping(TMatrix& newMat,
					  const std::vector<size_t>& vMap,
					  const TMatrix& origMat)
{
	PROFILE_FUNC_GROUP("gmg");
	//	check size
	UG_ASSERT(vMap.size() <= newMat.num_rows(), "Size must match. Map:"<<vMap.size()<<", mat:"<<newMat.num_rows());
	UG_ASSERT(vMap.size() <= newMat.num_cols(), "Size must match. Map:"<<vMap.size()<<", mat:"<<newMat.num_cols());
	UG_ASSERT(vMap.size() <= origMat.num_rows(), "Size must match. Map:"<<vMap.size()<<", mat:"<<origMat.num_rows());
	UG_ASSERT(vMap.size() <= origMat.num_cols(), "Size must match. Map:"<<vMap.size()<<", mat:"<<origMat.num_cols());

	newMat.resize_and_clear(origMat.num_rows(), origMat.num_cols());

	//	type of matrix row iterator
	typedef typename TMatrix::const_row_iterator const_row_iterator;

	//	loop all mapped indices
	for(size_t origInd = 0; origInd < vMap.size(); ++origInd)
	{
	//	get mapped level index
		const size_t newInd = vMap[origInd];

	//	loop all connections of the surface dof to other surface dofs and copy
	//	the matrix coupling into the level matrix

		for(const_row_iterator conn = origMat.begin_row(origInd);
									conn != origMat.end_row(origInd); ++conn)
		{
			size_t newConnIndex = conn.index();

		//	get corresponding level connection index
			if(conn.index() < vMap.size())
				newConnIndex = vMap[conn.index()];

		//	copy connection to level matrix
			newMat(newInd, newConnIndex) = conn.value();
		}
	}

	#ifdef UG_PARALLEL
	newMat.set_storage_type(origMat.get_storage_mask());
	#endif
}


/// projects surface function to level functions
template <typename TElem>
void CreateSurfaceToLevelMapping(std::vector<std::vector<int> >& vSurfLevelMapping,
                                 const std::vector<ConstSmartPtr<DoFDistribution> >& vLevelDD,
                                 ConstSmartPtr<DoFDistribution> surfaceDD,
                                 const SurfaceView& surfaceView)
{
//	type of element iterator
	typedef typename DoFDistribution::traits<TElem>::const_iterator iter_type;

//	iterators for subset
	iter_type iter = surfaceDD->begin<TElem>(SurfaceView::ALL);
	iter_type iterEnd = surfaceDD->end<TElem>(SurfaceView::ALL);

//	vector of indices
	std::vector<size_t> surfaceInd, levelInd;

//	loop all elements of type
	for( ; iter != iterEnd; ++iter)
	{
	//	get elem
		TElem* elem = *iter;

	//	extract all algebra indices for the element on surface
		surfaceDD->inner_algebra_indices(elem, surfaceInd);

	//	get level of element in hierarchy
		int level = surfaceView.get_level(elem);

	//	get corresponding level matrix for element
		UG_ASSERT(level < (int)vSurfLevelMapping.size(), "Level missing");
		std::vector<int>& levelMapping = vSurfLevelMapping[level];

	//	check that level is correct
		UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, "
				"									that is not passed.");

	//	extract all algebra indices for the element on level
		UG_ASSERT(vLevelDD[level].valid(), "DoF Distribution missing");
		vLevelDD[level]->inner_algebra_indices(elem, levelInd);

	//	check that index sets have same cardinality
		UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

	//	copy all elements of the matrix
		for(size_t i = 0; i < surfaceInd.size(); ++i)
		{
			UG_ASSERT(surfaceInd[i] < levelMapping.size(), "Index to large.");
			levelMapping[surfaceInd[i]] = levelInd[i];
		}
	}
}

/// creates a mapping of indices from the surface dof distribution to the level dof distribution
inline void CreateSurfaceToLevelMapping(std::vector<std::vector<int> >& vSurfLevelMapping,
                                        const std::vector<ConstSmartPtr<DoFDistribution> >& vLevelDD,
                                        ConstSmartPtr<DoFDistribution> surfDD,
                                        const SurfaceView& surfView)
{
//	resize the mapping
	vSurfLevelMapping.clear();
	vSurfLevelMapping.resize(vLevelDD.size());
	for(size_t lev = 0; lev < vSurfLevelMapping.size(); ++lev)
		vSurfLevelMapping[lev].resize(surfDD->num_indices(), -1);

	if(surfDD->max_dofs(VERTEX))
		CreateSurfaceToLevelMapping<VertexBase>(vSurfLevelMapping, vLevelDD, surfDD, surfView);
	if(surfDD->max_dofs(EDGE))
		CreateSurfaceToLevelMapping<EdgeBase>(vSurfLevelMapping, vLevelDD, surfDD, surfView);
	if(surfDD->max_dofs(FACE))
		CreateSurfaceToLevelMapping<Face>(vSurfLevelMapping, vLevelDD, surfDD, surfView);
	if(surfDD->max_dofs(VOLUME))
		CreateSurfaceToLevelMapping<Volume>(vSurfLevelMapping, vLevelDD, surfDD, surfView);
}

/// projects surface function to level functions
template <typename TMatrix>
void CopyMatrixSurfaceToLevel(TMatrix& levelMatrix,
                              const std::vector<int>& surfLevelMapping,
                              const TMatrix& surfMatrix)
{
//	type of matrix row iterator
	typedef typename TMatrix::const_row_iterator const_row_iterator;

//	loop all mapped indices
	for(size_t surfInd = 0; surfInd < surfLevelMapping.size(); ++surfInd)
	{
	//	get mapped level index
		const int levelInd = surfLevelMapping[surfInd];

	//	skipped non-mapped indices (indicated by -1)
		if(levelInd < 0) continue;

	//	loop all connections of the surface dof to other surface dofs and copy
	//	the matrix coupling into the level matrix

		for(const_row_iterator conn = surfMatrix.begin_row(surfInd);
						conn != surfMatrix.end_row(surfInd); ++conn)
		{
		//	get corresponding level connection index
			const int levelConnIndex = surfLevelMapping[conn.index()];

		//	check that index is from same level, too
			if(levelConnIndex < 0) continue;

		//	copy connection to level matrix
			levelMatrix(levelInd, levelConnIndex) = conn.value();
		}
	}

#ifdef UG_PARALLEL
	levelMatrix.set_storage_type(surfMatrix.get_storage_mask());
#endif
}

} // end namespace ug

// include implementation
#include "mg_solver_impl.hpp"

#endif /* __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER__ */
