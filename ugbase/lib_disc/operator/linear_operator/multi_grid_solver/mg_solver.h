/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
	/// constructor without arguments
		AssembledMultiGridCycle();

	/// constructor setting approximation space
		AssembledMultiGridCycle(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	/// sets the approximation space
		void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

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
		void set_gathered_base_solver_if_ambiguous(bool bGathered) {m_bGatheredBaseIfAmbiguous = bGathered;}

	///	sets if copies should be used to emulate a full-refined grid
		void set_emulate_full_refined_grid(bool bEmulate){
			if(bEmulate) m_GridLevelType = GridLevel::SURFACE;
			else m_GridLevelType = GridLevel::LEVEL;
		}

	///	sets if RAP - Product used to build coarse grid matrices
		void set_rap(bool bRAP) {m_bUseRAP = bRAP;}

	///	sets if smoothing is performed on surface rim
		void set_smooth_on_surface_rim(bool bSmooth) {m_bSmoothOnSurfaceRim = bSmooth;}

	///	sets the cycle type (1 = V-cycle, 2 = W-cycle, ...)
		void set_cycle_type(int type) {m_cycleType = type;}

	///	sets the cycle type (1 = V-cycle, 2 = W-cycle, ...)
		void set_cycle_type(const std::string& type) {
			if(TrimString(type) == "V") {m_cycleType = _V_;}
			else if(TrimString(type) == "W") {m_cycleType = _W_;}
			else if(TrimString(type) == "F") {m_cycleType = _F_;}
			else {UG_THROW("GMG::set_cycle_type: option '"<<type<<"' not supported.");}
		}

	///	sets if communication and computation should be overlaped
		void set_comm_comp_overlap(bool bOverlap) {m_bCommCompOverlap = bOverlap;}

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
		void set_transfer(SmartPtr<ITransferOperator<TDomain, TAlgebra> > P)
			{set_prolongation(P); set_restriction(P);}

	///	sets the prolongation operator
		void set_prolongation(SmartPtr<ITransferOperator<TDomain, TAlgebra> > P)
			{m_spProlongationPrototype = P;}

	///	sets the restriction operator
		void set_restriction(SmartPtr<ITransferOperator<TDomain, TAlgebra> > P)
			{m_spRestrictionPrototype = P;}

	///	sets the projection operator
		void set_projection(SmartPtr<ITransferOperator<TDomain, TAlgebra> > P)
			{m_spProjectionPrototype = P;}

	///	clears all transfer post process
		void clear_transfer_post_process()
			{m_vspProlongationPostProcess.clear(); m_vspRestrictionPostProcess.clear();}

	///	add prolongation post process
		void add_prolongation_post_process(SmartPtr<ITransferPostProcess<TDomain, TAlgebra> > PP)
			{m_vspProlongationPostProcess.push_back(PP);}

	///	add restriction post process
		void add_restriction_post_process(SmartPtr<ITransferPostProcess<TDomain, TAlgebra> > PP)
			{m_vspRestrictionPostProcess.push_back(PP);}

	///////////////////////////////////////////////////////////////////////////
	//	Linear Solver interface methods
	///////////////////////////////////////////////////////////////////////////

	///	name
		virtual const char* name() const {return "Geometric MultiGrid";}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
			if(!m_spPreSmootherPrototype->supports_parallel())
				return false;
			if(!m_spPostSmootherPrototype->supports_parallel())
				return false;
			return true;
		}

	///	returns information about configuration parameters
		virtual std::string config_string() const;

	/// Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

	///	Prepare for Linear Operator L
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L);

	///	does not call init on base-solver during initialization
	/**	Use this method with care. It can be useful e.g. during repeated
	 * adaption of a linear problem, where the base-level does not change.
	 * \warning	Currently only works in serial environments
	 * \note	This only affects calls to init(...) and has no effect on
	 *			direct calls to init_base_solver.
	 * \{ */
		virtual void ignore_init_for_base_solver(bool ignore);
		virtual bool ignore_init_for_base_solver() const;
	/** \\ */

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
		void lmgc(int lev, int cycleType);

	////////////////////////////////////////////////////////////////
	//	The methods in this section rely on each other and should be called in sequence
	///	performs presmoothing on the given level
		void presmooth_and_restriction(int lev);

	///	performs prolongation to the level above
		void prolongation_and_postsmooth(int lev);

	///	compute base solver
		void base_solve(int lev);
	//	end of section
	////////////////////////////////////////////////////////////////

	///	allocates the memory
		void init_level_memory(int baseLev, int topLev);

	///	initializes common part
		void init();

	///	initializes the smoother and base solver
		void init_smoother();

	///	initializes the coarse grid matrices
		void assemble_level_operator();
		void init_rap_operator();

	///	initializes the smoother and base solver
		void init_base_solver();

	///	initializes the transfers
		void init_transfer();

	///	initializes the prolongation
		void init_projection();

	///	assembles the missing matrix part on the coarse level, that must be
	///	added if the correction has been computed to ensure a correctly updated
	///	defect. (i.e. assembles A^c, with d^f -= A^c * c^c)
		void assemble_rim_cpl(const vector_type* u);
		void init_rap_rim_cpl();

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
		static const int _V_ = 1;
		static const int _W_ = 2;
		static const int _F_ = -1;

	///	number of Presmooth steps
		int m_numPreSmooth;

	///	number of Postsmooth steps
		int m_numPostSmooth;

	///	lowest level containing surface geom obj (proc-locally)
		int m_LocalFullRefLevel;

	///	grid-view for level vectors
		GridLevel::ViewType m_GridLevelType;

	///	using RAP-Product (assemble coarse-grid matrices otherwise)
		bool m_bUseRAP;

	///	flag if smoothing on surface rim
		bool m_bSmoothOnSurfaceRim;

	///	flag if overlapping communication and computation
		bool m_bCommCompOverlap;

	///	approximation space revision of cached values
		RevisionCounter m_ApproxSpaceRevision;

	///	prototype for pre-smoother
		SmartPtr<ILinearIterator<vector_type> > m_spPreSmootherPrototype;

	///	prototype for post-smoother
		SmartPtr<ILinearIterator<vector_type> > m_spPostSmootherPrototype;

	///	prototype for projection operator
		SmartPtr<ITransferOperator<TDomain, TAlgebra> > m_spProjectionPrototype;

	///	prototype for prolongation operator
		SmartPtr<ITransferOperator<TDomain, TAlgebra> > m_spProlongationPrototype;

	///	prototype for restriction operator
		SmartPtr<ITransferOperator<TDomain, TAlgebra> > m_spRestrictionPrototype;

	///	prototpe for transfer post process
		std::vector<SmartPtr<ITransferPostProcess<TDomain, TAlgebra> > > m_vspProlongationPostProcess;
		std::vector<SmartPtr<ITransferPostProcess<TDomain, TAlgebra> > > m_vspRestrictionPostProcess;

	///	base solver for the coarse problem
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spBaseSolver;

		////////////////////////////////////
		// Storage for each grid level
		////////////////////////////////////

	///	Structure used to realize Surface to Level mapping
	/// \{
		struct LevelIndex{
			LevelIndex() : index(-1), indexLower(-1), level(-1), levelLower(-1) {}
			size_t index, indexLower;
			int level, levelLower;
		};
		std::vector<LevelIndex> m_vSurfToLevelMap;

		template <typename TElem>
		void init_index_mappings();
		void init_index_mappings();
	/// \}

	///	Structure used to realize Surface to Level mapping
	/// \{
		struct SurfLevelMap{
			SurfLevelMap() : levIndex(-1), surfIndex(-1) {}
			SurfLevelMap(size_t levIndex_, size_t surfIndex_)
				: levIndex(levIndex_), surfIndex(surfIndex_) {}
			size_t levIndex, surfIndex;
		};
	/// \}

		struct LevData
		{
		///	Level matrix operator
			SmartPtr<MatrixOperator<matrix_type, vector_type> > A;

		///	Smoother
			SmartPtr<ILinearIterator<vector_type> > PreSmoother;
			SmartPtr<ILinearIterator<vector_type> > PostSmoother;

		///	Projection operator
			SmartPtr<ITransferOperator<TDomain, TAlgebra> > Projection;

		///	Transfer operator
			SmartPtr<ITransferOperator<TDomain, TAlgebra> > Prolongation;
			SmartPtr<ITransferOperator<TDomain, TAlgebra> > Restriction;

		///	vectors needed (sx = no-ghosts [for smoothing], t = for transfer)
			SmartPtr<GF> sc, sd, st, t;

		///	maps global indices (including ghosts) to patch indices (no ghosts included).
			std::vector<size_t> vMapPatchToGlobal;

		/// list of shadowing indices
			std::vector<size_t> vShadowing;

		/// list of corresponding surface index to shadowing indices
			std::vector<size_t> vSurfShadowing;

		///	map surface to level
			std::vector<SurfLevelMap> vSurfLevelMap;

		///	missing coarse grid correction
			matrix_type RimCpl_Fine_Coarse;

		///	missing coarse grid correction
			matrix_type RimCpl_Coarse_Fine;
		};

	///	storage for all level
		std::vector<SmartPtr<LevData> > m_vLevData;

	///	flag, if to solve base problem in parallel when gathered and (!) parallel possible
		bool m_bGatheredBaseIfAmbiguous;

	///	flag if using parallel base solver
		bool m_bGatheredBaseUsed;

	///	flag indicating whether the base-solver shall be initialized during a call to init()
		bool m_ignoreInitForBaseSolver;

	///	Matrix for gathered base solver
		SmartPtr<MatrixOperator<matrix_type, vector_type> > spGatheredBaseMat;

	///	vector for gathered base solver
		SmartPtr<GF> spGatheredBaseCorr;

	///	returns if gathered base master
		bool gathered_base_master() const;

	///	current surface correction
		GF* m_pC;

	///	init mapping from noghost -> w/ ghost
	/// \{
		template <typename TElem>
		void init_noghost_to_ghost_mapping(	std::vector<size_t>& vNoGhostToGhostMap,
											ConstSmartPtr<DoFDistribution> spNoGhostDD,
											ConstSmartPtr<DoFDistribution> spGhostDD);
		void init_noghost_to_ghost_mapping(int lev);
	/// \}

	///	copies vector to smoothing patch using cached mapping
		void copy_ghost_to_noghost(SmartPtr<GF> spVecTo,
		                           ConstSmartPtr<GF> spVecFrom,
		                           const std::vector<size_t>& vMapPatchToGlobal);

	///	copies vector from smoothing patch using cached mapping
		void copy_noghost_to_ghost(SmartPtr<GF> spVecTo,
								   ConstSmartPtr<GF> spVecFrom,
								   const std::vector<size_t>& vMapPatchToGlobal);

	///	copies matrix from smoothing patch using cached mapping
		void copy_noghost_to_ghost(SmartPtr<matrix_type> spMatTo,
		                           ConstSmartPtr<matrix_type> spMatFrom,
		                           const std::vector<size_t>& vMapPatchToGlobal);

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
	 */
	/// \{
		void write_debug(const matrix_type& mat, std::string name,
		                 const GridLevel& glTo, const GridLevel& glFrom);
		void write_debug(const matrix_type& mat, std::string name,
		                 const GF& rTo, const GF& rFrom);
	/// \}

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

/// selects all non-shadows, that are adjacent to a shadow on a grid levels
void SelectNonShadowsAdjacentToShadowsOnLevel(BoolMarker& sel,
										   const SurfaceView& surfView,
										   int level);

} // end namespace ug

// include implementation
#include "mg_solver_impl.hpp"

#endif /* __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER__ */
