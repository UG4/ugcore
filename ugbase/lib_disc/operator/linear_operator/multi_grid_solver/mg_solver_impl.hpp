/*
 * mg_solver_impl.h
 *
 *  Created on: 04.01.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__
#define __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__

#include <iostream>
#include <sstream>
#include <string>
#include "common/profiler/profiler.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/operator/linear_operator/prolongation_operator.h"
#include "lib_disc/operator/linear_operator/projection_operator.h"

#include "mg_solver.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
	#include "pcl/pcl_util.h"
//	the debug barrier is used to eliminate synchronization overhead from
//	profiling stats. Only used for parallel builds.
//	PCL_DEBUG_BARRIER only has an effect if PCL_DEBUG_BARRIER_ENABLED is defined.
	#define GMG_PARALLEL_DEBUG_BARRIER(comm) PCL_DEBUG_BARRIER(comm)

#else
	#define GMG_PARALLEL_DEBUG_BARRIER(comm)
#endif

#define PROFILE_GMG
#ifdef PROFILE_GMG
	#define GMG_PROFILE_FUNC()		PROFILE_FUNC()
	#define GMG_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "gmg")
	#define GMG_PROFILE_END()		PROFILE_END()
#else
	#define GMG_PROFILE_FUNC()
	#define GMG_PROFILE_BEGIN(name)
	#define GMG_PROFILE_END()
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
AssembledMultiGridCycle<TDomain, TAlgebra>::
AssembledMultiGridCycle(SmartPtr<ApproximationSpace<TDomain> > approxSpace) :
	m_spSurfaceMat(NULL), m_spAss(NULL), m_spApproxSpace(approxSpace),
	m_topLev(GridLevel::TOP), m_surfaceLev(GridLevel::TOP),
	m_baseLev(0), m_cycleType(1),
	m_numPreSmooth(2), m_numPostSmooth(2),
	m_LocalFullRefLevel(0), m_GridLevelType(GridLevel::LEVEL), m_bUseRAP(false),
	m_spPreSmootherPrototype(new Jacobi<TAlgebra>()),
	m_spPostSmootherPrototype(m_spPreSmootherPrototype),
	m_spProjectionPrototype(new InjectionTransfer<TDomain,TAlgebra>(m_spApproxSpace)),
	m_spProlongationPrototype(new StdTransfer<TDomain,TAlgebra>(m_spApproxSpace)),
	m_spRestrictionPrototype(m_spProlongationPrototype),
	m_spBaseSolver(new LU<TAlgebra>()),
	m_bGatheredBaseIfAmbiguous(true),
	m_spDebugWriter(NULL), m_dbgIterCnt(0)
{};

template <typename TDomain, typename TAlgebra>
SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >
AssembledMultiGridCycle<TDomain, TAlgebra>::
clone()
{
	SmartPtr<AssembledMultiGridCycle<TDomain, TAlgebra> > clone(
		new AssembledMultiGridCycle<TDomain, TAlgebra>(m_spApproxSpace));

	clone->set_base_level(m_baseLev);
	clone->set_base_solver(m_spBaseSolver);
	clone->set_cycle_type(m_cycleType);
	clone->set_debug(m_spDebugWriter);
	clone->set_discretization(m_spAss);
	clone->set_num_postsmooth(m_numPostSmooth);
	clone->set_num_presmooth(m_numPreSmooth);
	clone->set_projection(m_spProjectionPrototype);
	clone->set_prolongation(m_spProlongationPrototype);
	clone->set_restriction(m_spRestrictionPrototype);
	clone->set_presmoother(m_spPreSmootherPrototype);
	clone->set_postsmoother(m_spPostSmootherPrototype);
	clone->set_surface_level(m_surfaceLev);

	for(size_t i = 0; i < m_vspProlongationPostProcess.size(); ++i)
		clone->add_prolongation_post_process(m_vspProlongationPostProcess[i]);

	for(size_t i = 0; i < m_vspRestrictionPostProcess.size(); ++i)
		clone->add_restriction_post_process(m_vspRestrictionPostProcess[i]);

	return clone;
}

template <typename TDomain, typename TAlgebra>
AssembledMultiGridCycle<TDomain, TAlgebra>::
~AssembledMultiGridCycle()
{};

////////////////////////////////////////////////////////////////////////////////
// apply and init
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
bool AssembledMultiGridCycle<TDomain, TAlgebra>::
apply(vector_type& rC, const vector_type& rD)
{
	PROFILE_FUNC_GROUP("gmg");
	GF* pC = dynamic_cast<GF*>(&rC);
	if(!pC) UG_THROW("GMG::apply: Expect Correction to be grid based.")
	const GF* pD = dynamic_cast<const GF*>(&rD);
	if(!pD) UG_THROW("GMG::apply: Expect Defect to be grid based.")
	GF& c = *pC;
	const GF& d = *pD;

	try{
// 	Check if surface level has been chosen correctly
//	Please note, that the approximation space returns the global number of levels,
//	i.e. the maximum of levels among all processes.
	if(m_topLev >= (int)m_spApproxSpace->num_levels())
		UG_THROW("GMG::apply: SurfaceLevel "<<m_topLev<<" does not exist.");

// 	Check if base level has been choose correctly
	if(m_baseLev > m_topLev)
		UG_THROW("GMG::apply: Base level must be smaller or equal to surface Level.");

//	debug output
	write_debug(d, "Defect_In");

//	project defect from surface to level
	GMG_PROFILE_BEGIN(GMG_ProjectDefectFromSurface);
	try{
		for(size_t i = 0; i < m_vSurfToLevelMap.size(); ++i){
			const int level = m_vSurfToLevelMap[i].level;
			const size_t index = m_vSurfToLevelMap[i].index;

			(*(m_vLevData[level]->sd))[index] = d[i];
		}
#ifdef UG_PARALLEL
		for(int lev = m_baseLev; lev <= m_topLev; ++lev)
			m_vLevData[lev]->sd->set_storage_type(d.get_storage_mask());
#endif
		for(int lev = m_baseLev; lev <= m_topLev; ++lev)
			m_vLevData[lev]->sc->set(0.0);
	}
	UG_CATCH_THROW("GMG::apply: Project d Surf -> Level failed.");
	GMG_PROFILE_END(); //GMGApply_ProjectDefectFromSurface

// 	Perform one multigrid cycle
	GMG_PROFILE_BEGIN(GMG_lmgc);
	UG_DLOG(LIB_DISC_MULTIGRID, 4, "gmg-apply lmgc (on level " << m_topLev << ")... \n");
	try{
		lmgc(m_topLev);
	}
	UG_CATCH_THROW("GMG: lmgc failed.");
	GMG_PROFILE_END(); //GMGApply_lmgc

//	project correction from level to surface
	GMG_PROFILE_BEGIN(GMG_ProjectCorrectionFromLevelToSurface);
	UG_DLOG(LIB_DISC_MULTIGRID, 4, "gmg-apply project_level_to_surface... \n");
	try{
		for(size_t i = 0; i < m_vSurfToLevelMap.size(); ++i){
			const int level = m_vSurfToLevelMap[i].level;
			const size_t index = m_vSurfToLevelMap[i].index;

			 c[i] = (*(m_vLevData[level]->sc))[index];
		}

		#ifdef UG_PARALLEL
		c.set_storage_type(PST_CONSISTENT);
		#endif
	}
	UG_CATCH_THROW("GMG::apply: Project c Level -> Surface failed.");
	GMG_PROFILE_END(); //GMGApply_ProjectCorrectionFromLevelToSurface

//	apply scaling
	const number kappa = this->damping()->damping(c, d, m_spSurfaceMat.template cast_dynamic<ILinearOperator<vector_type> >());
	if(kappa != 1.0) c *= kappa;

//	debug output
	write_debug(c, "Correction_Out");
	write_surface_debug(*m_spSurfaceMat, "SurfaceStiffness");

//	increase dbg counter
	if(m_spDebugWriter.valid()) m_dbgIterCnt++;

	} UG_CATCH_THROW("GMG::apply: Application failed.");

	UG_DLOG(LIB_DISC_MULTIGRID, 4, "gmg-apply done. \n");
	return true;
}

template <typename TDomain, typename TAlgebra>
bool AssembledMultiGridCycle<TDomain, TAlgebra>::
apply_update_defect(vector_type &c, vector_type& rD)
{
	PROFILE_FUNC_GROUP("gmg");

//	NOTE: 	This is the implementation of a multiplicative Multigrid. Thus, the
//			defect is kept up to date when traversing the grid. At the end of
//			the iteration the updated defect is stored in the level defects and
//			could be projected back to the surface in order to get an updated
//			surface defect. This is, however, not done. For these reasons:
//			a) A Matrix-Vector operation is not very expensive
//			b) In a 2d adaptive case, the update is difficult, but can be
//				performed. But if the implementation is incorrect, this will
//				hardly be detected.
//			c) In a 3d adaptive case, it is impossible to ensure, that assembled
//				level matrices and the surface matrix have the same couplings
//				(not even to inner points). This is due to the fact, that e.g.
//				finite volume geometries are computed using different
//				integration points. (Hanging fv used triangles as scvf in 3d,
//				while normal fv use quads). Therefore, the updated defect is
//				only approximately the correct defect. In order to return the
//				correct defect, we must recompute the defect in the
//				adaptive case.
//			d) If scaling of the correction is performed, the defect must be
//				recomputed anyway. Thus, only scale=1.0 allows optimizing.
//			e) Updated defects are only needed in LinearIterativeSolvers. In
//				Krylov-Methods (CG, BiCGStab, ...) only the correction is
//				needed. We optimize for that case.

//	compute correction
	if(!apply(c, rD)) return false;

//	update defect: d = d - A*c
	m_spSurfaceMat->matmul_minus(rD, c);

//	write for debugging
	const GF* pD = dynamic_cast<const GF*>(&rD);
	if(!pD) UG_THROW("GMG::apply: Expect Defect to be grid based.")
	const GF& d = *pD;
	write_debug(d, "Defect_Out");

	return true;
}

template <typename TDomain, typename TAlgebra>
bool AssembledMultiGridCycle<TDomain, TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - init(J, u)\n");

	// try to extract assembling routine
	SmartPtr<AssembledLinearOperator<TAlgebra> > spALO =
			J.template cast_dynamic<AssembledLinearOperator<TAlgebra> >();
	if(spALO.valid()){
		m_spAss = spALO->discretization();
	}

	// Store Surface Matrix
	m_spSurfaceMat = J.template cast_dynamic<matrix_type>();

	// Store Surface Solution
	m_pSurfaceSol = &u;

	// call init
	init();

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - init(J, u)\n");
	return true;
}

template <typename TDomain, typename TAlgebra>
bool AssembledMultiGridCycle<TDomain, TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > L)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - init(L)\n");

	// try to extract assembling routine
	SmartPtr<AssembledLinearOperator<TAlgebra> > spALO =
			L.template cast_dynamic<AssembledLinearOperator<TAlgebra> >();
	if(spALO.valid()){
		m_spAss = spALO->discretization();
	}

	// Store Surface Matrix
	m_spSurfaceMat = L.template cast_dynamic<matrix_type>();

	// Store Surface Solution
	m_pSurfaceSol = NULL;

	// call init
	init();

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - init(L)\n");
	return true;
}



template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init()
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start init_common\n");

	try{

// 	Cast Operator
	if(m_spSurfaceMat.invalid())
		UG_THROW("GMG:init: Can not cast Operator to Matrix.");

//	Check Approx Space
	if(m_spApproxSpace.invalid())
		UG_THROW("GMG::init: Approximation Space not set.");

//	check that grid given
	if(m_spApproxSpace->num_levels() == 0)
		UG_THROW("GMG:init: No grid levels");

	if(m_spAss.invalid())
		UG_THROW("GMG::init: Discretization not set.");

	if(m_spBaseSolver.invalid())
		UG_THROW("GMG::init: Base Solver not set.");

	if(m_spPreSmootherPrototype.invalid())
		UG_THROW("GMG::init: PreSmoother not set.");

	if(m_spPostSmootherPrototype.invalid())
		UG_THROW("GMG::init: PostSmoother not set.");

	if(m_spProlongationPrototype.invalid())
		UG_THROW("GMG::init: Prolongation not set.");

	if(m_spRestrictionPrototype.invalid())
		UG_THROW("GMG::init: Restriction not set.");

//	get current toplevel
	const GF* pSol = dynamic_cast<const GF*>(m_pSurfaceSol);
	if(pSol){
		m_surfaceLev = pSol->dof_distribution()->grid_level().level();
	}

	if(m_surfaceLev != GridLevel::TOP) m_topLev = m_surfaceLev;
	else m_topLev = m_spApproxSpace->num_levels() - 1;

	if(m_baseLev > m_topLev)
		UG_THROW("GMG::init: Base Level greater than Surface level.");

	if(m_ApproxSpaceRevision != m_spApproxSpace->revision())
	{
	//	Allocate memory for given top level
		try{
			init_level_memory(m_baseLev, m_topLev);
		}
		UG_CATCH_THROW("GMG::init: Cannot allocate memory.");

	//	init mapping from surface level to top level in case of full refinement
		GMG_PROFILE_BEGIN(GMG_InitSurfToLevelMapping);
		try{
			init_surface_to_level_mapping(m_vSurfToLevelMap, true);

			if(m_bUseRAP && (m_LocalFullRefLevel != m_topLev)){
				init_surface_to_level_mapping(m_vSurfToLevelMapLowerLev, false);
			} else {
				m_vSurfToLevelMapLowerLev.resize(0);
			}
		}
		UG_CATCH_THROW("GMG: Cannot create SurfaceToLevelMap.")
		GMG_PROFILE_END();

	//	init mapping from surface level to top level in case of full refinement
		GMG_PROFILE_BEGIN(GMG_CollectShadowing);
		try{
			for(int lev = m_baseLev; lev <= m_topLev; ++lev)
				collect_shadowing_indices(lev);
		}
		UG_CATCH_THROW("GMG: Cannot create SurfaceToLevelMap.")
		GMG_PROFILE_END();

	// 	Create Interpolation
		GMG_PROFILE_BEGIN(GMG_InitProlongation);
		try{
			init_transfer();
		}
		UG_CATCH_THROW("GMG:init: Cannot init Transfer (Prolongation/Restriction).");
		GMG_PROFILE_END();

	//	remember revision counter of approx space
		m_ApproxSpaceRevision = m_spApproxSpace->revision();
	}

//	Assemble coarse grid operators
	GMG_PROFILE_BEGIN(GMG_AssembleLevelGridOperator);
	try{
		if(m_bUseRAP){
			init_rap_operator();
		} else {
			assemble_level_operator();
		}
	}
	UG_CATCH_THROW("GMG: Initialization of Level Operator failed.");
	GMG_PROFILE_END();

//	Init smoother for coarse grid operators
	GMG_PROFILE_BEGIN(GMG_InitSmoother);
	try{
		init_smoother();
	}
	UG_CATCH_THROW("GMG:init: Cannot init Smoother.");
	GMG_PROFILE_END();

//	Init base solver
	GMG_PROFILE_BEGIN(GMG_InitBaseSolver);
	try{
		init_base_solver();
	}
	UG_CATCH_THROW("GMG:init: Cannot init Base Solver.");
	GMG_PROFILE_END();

	} UG_CATCH_THROW("GMG: Init failure for init(u)");

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop init_common\n");
}


template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
assemble_level_operator()
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start assemble_level_operator\n");

//	Create Projection
	try{
		if(m_pSurfaceSol) {
			UG_DLOG(LIB_DISC_MULTIGRID, 4, "  start assemble_level_operator: project sol\n");
			GMG_PROFILE_BEGIN(GMG_ProjectSurfaceSolution);
			#ifdef UG_PARALLEL
			if(!m_pSurfaceSol->has_storage_type(PST_CONSISTENT))
				UG_THROW("GMG::init: Can only project "
						"a consistent solution. Make sure to pass a consistent on.");
			#endif

			init_projection();

			for(size_t i = 0; i < m_vSurfToLevelMap.size(); ++i){
				const int level = m_vSurfToLevelMap[i].level;
				const size_t index = m_vSurfToLevelMap[i].index;

				(*(m_vLevData[level]->st))[index] = (*m_pSurfaceSol)[i];
			}
			GMG_PROFILE_END();
			UG_DLOG(LIB_DISC_MULTIGRID, 4, "  end   assemble_level_operator: project sol\n");
		}
	}
	UG_CATCH_THROW("GMG::init: Projection of Surface Solution failed.");

// 	Create coarse level operators
	for(int lev = m_topLev; lev >= m_baseLev; --lev)
	{
		LevData& ld = *m_vLevData[lev];

	//	In Full-Ref case we can copy the Matrix from the surface
		bool bCpyFromSurface = ((lev == m_topLev) && (lev <= m_LocalFullRefLevel));

		if(!bCpyFromSurface)
		{
			UG_DLOG(LIB_DISC_MULTIGRID, 4, "  start assemble_level_operator: assemble on lev "<<lev<<"\n");
			GMG_PROFILE_BEGIN(GMG_AssLevelMat);
			try{
			if(m_GridLevelType == GridLevel::LEVEL)
				m_spAss->ass_tuner()->set_force_regular_grid(true);
			m_spAss->assemble_jacobian(*ld.A, *ld.st, GridLevel(lev, m_GridLevelType, false));
			m_spAss->ass_tuner()->set_force_regular_grid(false);
			}
			UG_CATCH_THROW("GMG:init: Cannot init operator for level "<<lev);
			GMG_PROFILE_END();
			UG_DLOG(LIB_DISC_MULTIGRID, 4, "  end   assemble_level_operator: assemble on lev "<<lev<<"\n");
		}
		else
		{
		//	in case of full refinement we simply copy the matrix (with correct numbering)
			UG_DLOG(LIB_DISC_MULTIGRID, 4, "  start assemble_level_operator: copy mat on lev "<<lev<<"\n");
			GMG_PROFILE_BEGIN(GMG_CopySurfMat);

		//	loop all mapped indices
			UG_ASSERT(m_spSurfaceMat->num_rows() == m_vSurfToLevelMap.size(),
			          "Surface Matrix rows != Surf Level Indices")
			ld.A->resize_and_clear(m_spSurfaceMat->num_rows(), m_spSurfaceMat->num_cols());
			for(size_t surfFrom = 0; surfFrom < m_vSurfToLevelMap.size(); ++surfFrom)
			{
			//	get mapped level index
				UG_ASSERT(m_vSurfToLevelMap[surfFrom].level == m_topLev,
				          "All surface Indices must be on top level for full-ref.")
				const size_t lvlFrom = m_vSurfToLevelMap[surfFrom].index;

			//	loop all connections of the surface dof to other surface dofs
			//	and copy the matrix coupling into the level matrix
				typedef typename matrix_type::const_row_iterator const_row_iterator;
				const_row_iterator conn = m_spSurfaceMat->begin_row(surfFrom);
				const_row_iterator connEnd = m_spSurfaceMat->end_row(surfFrom);
				for( ;conn != connEnd; ++conn){
				//	get corresponding level connection index
					UG_ASSERT(m_vSurfToLevelMap[conn.index()].level == m_topLev,
					          "All surface Indices must be on top level for full-ref.")
					const size_t lvlTo = m_vSurfToLevelMap[conn.index()].index;

				//	copy connection to level matrix
					(*ld.A)(lvlFrom, lvlTo) = conn.value();
				}
			}

			#ifdef UG_PARALLEL
			ld.A->set_storage_type(m_spSurfaceMat->get_storage_mask());
			ld.A->set_layouts(ld.st->layouts());
			#endif
			GMG_PROFILE_END();
			UG_DLOG(LIB_DISC_MULTIGRID, 4, "  end   assemble_level_operator: copy mat on lev "<<lev<<"\n");
		}

		if(m_pSurfaceSol && lev > m_baseLev)
		{
			GMG_PROFILE_BEGIN(GMG_ProjectSolutionDown);
			LevData& lc = *m_vLevData[lev-1];

			copy_noghost_to_ghost(ld.t, ld.st, ld.vMapPatchToGlobal);
			copy_to_vertical_masters(*ld.t);

			try{
				ld.Projection->do_restrict(*lc.st, *ld.t);
			} UG_CATCH_THROW("GMG::init: Cannot project "
						"solution to coarse grid function of level "<<lev-1<<".\n");

			#ifdef UG_PARALLEL
			lc.st->set_storage_type(m_pSurfaceSol->get_storage_mask());
			#endif
			GMG_PROFILE_END();
		}
	}

//	if no ghosts are present we can simply use the whole grid. If the base
//	solver is carried out in serial (gathering to some processes), we have
//	to assemble the assemble the coarse grid matrix on the whole grid as
//	well
	if(m_bGatheredBaseUsed)
	{
		UG_DLOG(LIB_DISC_MULTIGRID, 4, "  start assemble_level_operator: ass gathered on lev "<<m_baseLev<<"\n");
		GMG_PROFILE_BEGIN(GMG_AssGatheredLevMat);
		LevData& ld = *m_vLevData[m_baseLev];

		if(m_pSurfaceSol){
			copy_noghost_to_ghost(ld.t, ld.st, ld.vMapPatchToGlobal);
			copy_to_vertical_masters(*ld.t);
		}

		// only assemble on gathering proc
		if(gathered_base_master())
		{
			try{
			if(m_GridLevelType == GridLevel::LEVEL)
				m_spAss->ass_tuner()->set_force_regular_grid(true);
			m_spAss->assemble_jacobian(*spGatheredBaseMat, *ld.t, GridLevel(m_baseLev, m_GridLevelType, true));
			m_spAss->ass_tuner()->set_force_regular_grid(false);
			}
			UG_CATCH_THROW("GMG:init: Cannot init operator base level operator");
		}
		GMG_PROFILE_END();
		UG_DLOG(LIB_DISC_MULTIGRID, 4, "  end   assemble_level_operator: ass gathered on lev "<<m_baseLev<<"\n");
	}

//	write computed level matrices for debug purpose
	for(int lev = m_baseLev; lev <= m_topLev; ++lev){
		write_smooth_level_debug(*m_vLevData[lev]->A, "LevelMatrix", lev);
	}

//	assemble missing coarse grid matrix contribution (only in adaptive case)
	try{
		assemble_missing_coarse_grid_coupling(m_pSurfaceSol);
	}
	UG_CATCH_THROW("GMG:init: Cannot init missing coarse grid coupling.");

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop assemble_level_operator\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
assemble_missing_coarse_grid_coupling(const vector_type* u)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - assemble_missing_coarse_grid_coupling " << "\n");

//	clear matrices
	for(int lev = m_baseLev; lev <= m_topLev; ++lev)
		m_vLevData[lev]->CoarseGridContribution.resize_and_clear(0,0);

//	if the grid is fully refined, nothing to do
	if(m_topLev <= m_LocalFullRefLevel){
		UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - assemble_missing_coarse_grid_coupling (non-adaptive)" << "\n");
		return;
	}

//	get the surface view
	const SurfaceView& surfView = *m_spApproxSpace->surface_view();

//	loop all levels to compute the missing contribution
	BoolMarker sel(*m_spApproxSpace->domain()->grid());
	for(int lev = m_LocalFullRefLevel; lev < m_topLev; ++lev)
	{
		LevData& lc = *m_vLevData[lev];
		LevData& lf = *m_vLevData[lev+1];

	//	select all elements, that have a shadow as a subelement, but are not
	//	itself a shadow -  we assemble on those elements only
		sel.clear();
		SelectNonShadowsAdjacentToShadowsOnLevel(sel, surfView, lev);

	//	create a surface matrix
		matrix_type surfMat;
		GridLevel surfLevel(m_surfaceLev, GridLevel::SURFACE);

	//	assemble the surface jacobian only for selected elements
		m_spAss->ass_tuner()->set_marker(&sel);
		if(u){
			m_spAss->assemble_jacobian(surfMat, *u, surfLevel);
		}else{
		//	\todo: not use tmp vector here
			vector_type tmpVec; tmpVec.resize(m_spApproxSpace->dof_distribution(surfLevel)->num_indices());
			m_spAss->assemble_jacobian(surfMat, tmpVec, surfLevel);
		}
		m_spAss->ass_tuner()->set_marker(NULL);

		UG_ASSERT(surfMat.num_rows() == m_vSurfToLevelMap.size(), "Size mismatch")
		UG_ASSERT(surfMat.num_cols() == m_vSurfToLevelMap.size(), "Size mismatch")

		lc.CoarseGridContribution.resize_and_clear(lf.sd->size(), lc.sc->size());

		for(size_t i = 0; i< lf.vShadowing.size(); ++i)
		{
			const size_t lvlFrom = lf.vShadowing[i];
			const size_t surfFrom = lf.vSurfShadowing[i];

			typedef typename matrix_type::row_iterator row_iterator;
			row_iterator conn = surfMat.begin_row(surfFrom);
			row_iterator connEnd = surfMat.end_row(surfFrom);
			for( ;conn != connEnd; ++conn)
			{
				const size_t surfTo = conn.index();
				if(m_vSurfToLevelMap[surfTo].level != lev) continue;

				const size_t lvlTo = m_vSurfToLevelMap[surfTo].index;

				(lc.CoarseGridContribution)(lvlFrom, lvlTo) = conn.value();
			}
		}

#ifdef UG_PARALLEL
		lc.CoarseGridContribution.set_storage_type(surfMat.get_storage_mask());
#endif
		//write_level_debug(m_vLevData[lev]->CoarseGridContribution, "MissingLevelMat", lev);
		write_surface_debug(surfMat, std::string("MissingSurfMat_").append(ToString(lev)));
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - assemble_missing_coarse_grid_coupling " << "\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_rap_operator()
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start init_rap_operator\n");

	for(int lev = m_topLev; lev >= m_baseLev; --lev)
	{
		LevData& ld = *m_vLevData[lev];
		ld.A->resize_and_clear(ld.st->size(), ld.st->size());
		#ifdef UG_PARALLEL
		ld.A->set_storage_type(m_spSurfaceMat->get_storage_mask());
		ld.A->set_layouts(ld.st->layouts());
		#endif
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 4, "  start init_rap_operator: copy from surface\n");
	GMG_PROFILE_BEGIN(GMG_CopyMatFromSurface);
	for(size_t surfFrom = 0; surfFrom < m_vSurfToLevelMap.size(); ++surfFrom)
	{
	//	loop all connections of the surface dof to other surface dofs
	//	and copy the matrix coupling into the level matrix
		typedef typename matrix_type::const_row_iterator const_row_iterator;
		const_row_iterator conn = m_spSurfaceMat->begin_row(surfFrom);
		const_row_iterator connEnd = m_spSurfaceMat->end_row(surfFrom);
		for( ;conn != connEnd; ++conn)
		{
		//	get mapped level index
			int lvlFrom = m_vSurfToLevelMap[surfFrom].level;
			size_t indexFrom = m_vSurfToLevelMap[surfFrom].index;

		//	get corresponding level connection index
			const size_t surfTo = conn.index();
			int lvlTo = m_vSurfToLevelMap[surfTo].level;
			size_t indexTo = m_vSurfToLevelMap[surfTo].index;

			if(lvlTo > lvlFrom){
				lvlTo = m_vSurfToLevelMapLowerLev[surfTo].level;
				indexTo = m_vSurfToLevelMapLowerLev[surfTo].index;
			}

			if(lvlFrom > lvlTo){
				lvlFrom = m_vSurfToLevelMapLowerLev[surfFrom].level;
				indexFrom = m_vSurfToLevelMapLowerLev[surfFrom].index;
			}

			if(lvlTo != lvlFrom)
				continue;

		//	copy connection to level matrix
			(*(m_vLevData[lvlTo]->A))(indexFrom, indexTo) = conn.value();
		}
	}
	GMG_PROFILE_END();
	UG_DLOG(LIB_DISC_MULTIGRID, 4, "  end   init_rap_operator: copy from surface\n");

//	write computed level matrices for debug purpose
	for(int lev = m_baseLev; lev <= m_topLev; ++lev){
		write_smooth_level_debug(*m_vLevData[lev]->A, "LevelMatrixFromSurf", lev);
	}

// 	Create coarse level operators
	UG_DLOG(LIB_DISC_MULTIGRID, 4, "  start init_rap_operator: build rap\n");
	GMG_PROFILE_BEGIN(GMG_BuildAllRAP);
	for(int lev = m_topLev; lev > m_baseLev; --lev)
	{
		UG_DLOG(LIB_DISC_MULTIGRID, 4, "  start init_rap_operator: build rap on lev "<<lev<<"\n");
		LevData& lf = *m_vLevData[lev];
		LevData& lc = *m_vLevData[lev-1];

		GMG_PROFILE_BEGIN(GMG_CreateRestrictionAndProlongation);
		SmartPtr<matrix_type> R = lf.Restriction->restriction();
		SmartPtr<matrix_type> P = lf.Prolongation->prolongation();
		GMG_PROFILE_END();

		GMG_PROFILE_BEGIN(GMG_CopyGhosts);
		SmartPtr<matrix_type> spA = lf.A;
#ifdef UG_PARALLEL
		if( !lf.t->layouts()->vertical_master().empty() ||
			!lf.t->layouts()->vertical_slave().empty())
		{
			SmartPtr<matrix_type> spGhostA = SmartPtr<matrix_type>(new matrix_type);
			spGhostA->resize_and_clear(lf.t->size(), lf.t->size());
			spGhostA->set_layouts(lf.t->layouts());

			if(!lf.t->layouts()->vertical_master().empty()){
				copy_noghost_to_ghost(spGhostA, lf.A, lf.vMapPatchToGlobal);
			} else {
				*spGhostA = *lf.A;
			}
			divide_vertical_slave_rows_by_number_of_masters(*spGhostA);

			ComPol_MatAddInnerInterfaceCouplings<matrix_type> cpMatAdd(*spGhostA, true);
			if(!spGhostA->layouts()->vertical_slave().empty())
				m_Com.send_data(spGhostA->layouts()->vertical_slave(), cpMatAdd);
			if(!spGhostA->layouts()->vertical_master().empty())
				m_Com.receive_data(spGhostA->layouts()->vertical_master(), cpMatAdd);
			m_Com.communicate();
			spA = spGhostA;
		}
#endif
		GMG_PROFILE_END();

		GMG_PROFILE_BEGIN(GMG_ComputeRAP);
		AddMultiplyOf(*lc.A, *R, *spA, *P);
		GMG_PROFILE_END();
		UG_DLOG(LIB_DISC_MULTIGRID, 4, "  end   init_rap_operator: build rap on lev "<<lev<<"\n");
	}
	GMG_PROFILE_END();
	UG_DLOG(LIB_DISC_MULTIGRID, 4, "  end   init_rap_operator: build rap\n");

	if(m_bGatheredBaseUsed)
	{
		GMG_PROFILE_BEGIN(GMG_GatherFromBaseSolver);
		LevData& ld = *m_vLevData[m_baseLev];

#ifdef UG_PARALLEL
		if(gathered_base_master()){
			spGatheredBaseMat->resize_and_clear(ld.t->size(), ld.t->size());
			copy_noghost_to_ghost(spGatheredBaseMat.template cast_dynamic<matrix_type>(), ld.A, ld.vMapPatchToGlobal);
		} else {
			spGatheredBaseMat = SmartPtr<MatrixOperator<matrix_type, vector_type> >(new MatrixOperator<matrix_type, vector_type>);
			*spGatheredBaseMat = *ld.A;
		}
		spGatheredBaseMat->set_storage_type(m_spSurfaceMat->get_storage_mask());
		spGatheredBaseMat->set_layouts(ld.t->layouts());

		divide_vertical_slave_rows_by_number_of_masters(*spGatheredBaseMat);
		ComPol_MatAddInnerInterfaceCouplings<matrix_type> cpMatAdd(*spGatheredBaseMat, true);
		if(!spGatheredBaseMat->layouts()->vertical_slave().empty())
			m_Com.send_data(spGatheredBaseMat->layouts()->vertical_slave(), cpMatAdd);
		if(gathered_base_master())
			m_Com.receive_data(spGatheredBaseMat->layouts()->vertical_master(), cpMatAdd);
		m_Com.communicate();

		if(!gathered_base_master()){
			spGatheredBaseMat = NULL;
		}
#endif
		GMG_PROFILE_END();
	}

//	write computed level matrices for debug purpose
	for(int lev = m_baseLev; lev <= m_topLev; ++lev){
		write_smooth_level_debug(*m_vLevData[lev]->A, "LevelMatrix", lev);
	}

//	coarse grid contributions
	try{
		init_rap_missing_coarse_grid_coupling();
	}
	UG_CATCH_THROW("GMG:init: Cannot init missing coarse grid coupling.");

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop init_rap_operator\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_rap_missing_coarse_grid_coupling()
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - init_rap_missing_coarse_grid_coupling " << "\n");

//	clear matrices
	for(int lev = m_baseLev; lev <= m_topLev; ++lev){
		LevData& ld = *m_vLevData[lev];
		ld.CoarseGridContribution.resize_and_clear(0,0);
	}

//	if the grid is fully refined, nothing to do
	if(m_topLev <= m_LocalFullRefLevel){
		UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - init_rap_missing_coarse_grid_coupling (non-adaptive)" << "\n");
		return;
	}

	for(int lev = m_LocalFullRefLevel; lev < m_topLev; ++lev)
	{
		LevData& lc = *m_vLevData[lev];
		LevData& lf = *m_vLevData[lev+1];
#ifdef UG_PARALLEL
		lc.CoarseGridContribution.set_storage_type(m_spSurfaceMat->get_storage_mask());
#endif
		lc.CoarseGridContribution.resize_and_clear(lf.sd->size(), lc.sc->size());

		for(size_t i = 0; i< lf.vShadowing.size(); ++i)
		{
			const size_t lvlFrom = lf.vShadowing[i];
			const size_t surfFrom = lf.vSurfShadowing[i];

			typedef typename matrix_type::const_row_iterator const_row_iterator;
			const_row_iterator conn = m_spSurfaceMat->begin_row(surfFrom);
			const_row_iterator connEnd = m_spSurfaceMat->end_row(surfFrom);
			for( ;conn != connEnd; ++conn)
			{
				const size_t surfTo = conn.index();
				if(m_vSurfToLevelMap[surfTo].level != lev) continue;

				const size_t lvlTo = m_vSurfToLevelMap[surfTo].index;

				(lc.CoarseGridContribution)(lvlFrom, lvlTo) = conn.value();
			}
		}
	}


	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - init_rap_missing_coarse_grid_coupling " << "\n");
}


template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_transfer()
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start init_transfer\n");

//	loop all levels
	for(int lev = m_baseLev+1; lev <= m_topLev; ++lev)
	{
		LevData& lf = *m_vLevData[lev];
		LevData& lc = *m_vLevData[lev-1];

	//	check if same operator for prolongation and restriction used
		bool bOneOperator = false;
		if(lf.Prolongation == lf.Restriction) bOneOperator = true;

		lf.Prolongation->set_levels(lc.st->grid_level(), lf.t->grid_level());
		lf.Prolongation->clear_constraints();
		for(size_t i = 0; i < m_spAss->num_constraints(); ++i){
			SmartPtr<IConstraint<TAlgebra> > pp = m_spAss->constraint(i);
			lf.Prolongation->add_constraint(pp);
		}
		lf.Prolongation->init();

		if(!bOneOperator){
			lf.Restriction->set_levels(lc.st->grid_level(), lf.t->grid_level());
			lf.Restriction->clear_constraints();
			for(size_t i = 0; i < m_spAss->num_constraints(); ++i){
				SmartPtr<IConstraint<TAlgebra> > pp = m_spAss->constraint(i);
				lf.Restriction->add_constraint(pp);
			}
			lf.Restriction->init();
		}
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop init_transfer\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_projection()
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start init_projection\n");

	for(int lev = m_baseLev+1; lev <= m_topLev; ++lev)
	{
		LevData& lf = *m_vLevData[lev];
		LevData& lc = *m_vLevData[lev-1];
		lf.Projection->set_levels(lc.st->grid_level(), lf.t->grid_level());
		lf.Projection->init();
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop init_projection\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_smoother()
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start init_smoother\n");

// 	Init smoother
	for(int lev = m_baseLev+1; lev <= m_topLev; ++lev)
	{
		LevData& ld = *m_vLevData[lev];

		UG_DLOG(LIB_DISC_MULTIGRID, 4, "  init_smoother: initializing pre-smoother on lev "<<lev<<"\n");
		if(!ld.PreSmoother->init(ld.A, *ld.st))
			UG_THROW("GMG::init: Cannot init pre-smoother for level "<<lev);

		UG_DLOG(LIB_DISC_MULTIGRID, 4, "  init_smoother: initializing post-smoother on lev "<<lev<<"\n");
		if(ld.PreSmoother != ld.PostSmoother)
		{
			if(!ld.PostSmoother->init(ld.A, *ld.st))
				UG_THROW("GMG::init: Cannot init post-smoother for level "<<lev);
		}
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop init_smoother\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_base_solver()
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start init_base_solver\n");
	LevData& ld = *m_vLevData[m_baseLev];

//	check, if a gathering base solver is required:
	if(m_bGatheredBaseUsed)
	{
	//  only init on gathering proc
		if(!gathered_base_master()) return;

	//	create layout with only v-master (v-slave do not exist on gathering proc)
	//	especially: remove h-layouts and set proc-comm to process-only
		#ifdef UG_PARALLEL
		SmartPtr<AlgebraLayouts> spOneProcLayout =
				SmartPtr<AlgebraLayouts>(new AlgebraLayouts);
		spOneProcLayout->clear();
		spOneProcLayout->vertical_master() = ld.t->layouts()->vertical_master();
		spOneProcLayout->proc_comm() = pcl::ProcessCommunicator(pcl::PCD_LOCAL);

		spGatheredBaseMat->set_layouts(spOneProcLayout);
		spGatheredBaseCorr->set_layouts(spOneProcLayout);
		ld.t->set_layouts(spOneProcLayout);
		#endif

	//	we init the base solver with the whole grid matrix
		if(!m_spBaseSolver->init(spGatheredBaseMat, *spGatheredBaseCorr))
			UG_THROW("GMG::init: Cannot init base solver on baselevel "<< m_baseLev);
	}
	else
	{
		if(ld.st->num_indices() == 0) return;

#ifdef UG_PARALLEL
		if(!ld.st->layouts()->master().empty() || !ld.st->layouts()->slave().empty())
			if(!m_spBaseSolver->supports_parallel())
				UG_THROW("GMG: Base level is distributed onto more than process, "
						"but the chosen base solver "<<m_spBaseSolver->name()<<
						" does not support parallel solving. Choose a parallel"
						" base solver or construct a situation, where a single"
						" process stores the whole base grid, to cure this issue.");
#endif

		if(!m_spBaseSolver->init(ld.A, *ld.st))
			UG_THROW("GMG::init: Cannot init base solver on baselevel "<< m_baseLev);
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop init_base_solver\n");
}

template <typename TDomain, typename TAlgebra>
template <typename TElem>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_surface_to_level_mapping(std::vector<LevelIndex>& vSurfToLevelMap,
                              bool bHigherLevWins)
{
/*	This Method is used to create a caching the transfer between surface and
 * 	level grid functions. Note, that in the adaptive case the surface grid is
 * 	distributed over several levels (all elements that do not have children). But
 * 	even in the full refinement case the surface level and the top level may
 * 	not have the same index pattern, since a different sorting may be used (however
 * 	the number of indices must be equal in the full-ref case).
 * 	In every case, each surface index has a corresponding index on some level
 * 	of the level grid functions. In this method this index is looked up and
 * 	stored in a vector (for fast access). I.e. for every surface index i, the
 * 	corresponding index in the level grid function is stored in m_vSurfToLevelMap[i]
 * 	together with the level of the grid function. Using this mapping copying
 * 	between surface <-> levels is implemented. Note: When Shadow-Copy are present
 * 	the dof manager usually assigns the same surface index to both, shadowing and
 * 	shadow-copy. Thus, there exist more than one corresponding level index for
 * 	such a surface index. In this case the shadowing index is used in the mapping
 * 	since this index will have the correct value at the end of the multigrid
 * 	cycle and at startup the projection to the level is necessary only to shadowing,
 * 	because shadows will get their value by transfer between the level. In order
 * 	to get the map this way, the surface loop below is only performed on
 * 	SURFACE_NONCOPY elements.
 */

	ConstSmartPtr<SurfaceView> spSurfView = m_spApproxSpace->surface_view();

	std::vector<ConstSmartPtr<DoFDistribution> > vLevelDD(m_topLev+1);
	for(int lev = m_baseLev; lev <= m_topLev; ++lev)
		vLevelDD[lev] = m_spApproxSpace->dof_distribution(GridLevel(lev, m_GridLevelType, false));

	ConstSmartPtr<DoFDistribution> surfDD =
			m_spApproxSpace->dof_distribution(GridLevel(m_surfaceLev, GridLevel::SURFACE));

//	iterators for subset
	// \todo: The loop below should only be on SURFACE_NONCOPY, since the
	//		 copy-shadows have the same indices as their shadowing. In such a
	//		 case the child index (shadowing) must win. This could be realized by
	//		 looping only non-copy elements. But it seems, that on a process
	//		 the shadow-copy may exist without its shadowing. In such a case
	//		 the mapping is invalid. To fix this, the loop is extended temporarily
	//		 below and doubling of appearance is checked.
	typedef typename DoFDistribution::traits<TElem>::const_iterator iter_type;
	iter_type iter = surfDD->begin<TElem>(SurfaceView::ALL);
	iter_type iterEnd = surfDD->end<TElem>(SurfaceView::ALL);

//	vector of indices
	std::vector<size_t> vSurfInd, vLevelInd;

//	loop all elements of type
	for( ; iter != iterEnd; ++iter)
	{
	//	get elem and its level
		TElem* elem = *iter;
		int level = m_spApproxSpace->domain()->grid()->get_level(elem);

		if (m_GridLevelType == GridLevel::SURFACE)
			level = m_topLev;

	//	check that coarse grid covers whole domain. If this is not the case,
	//	some surface indices are mapped to grid levels below baseLev. We
	//	do not allow that.
		if(level < m_baseLev)
			UG_THROW("GMG::init: Some index of the surface grid is located on "
					"level "<<level<<", which is below the chosen baseLev "<<
					m_baseLev<<". This is not allowed, since otherwise the "
					"gmg correction would only affect parts of the domain. Use"
					"gmg:set_base_level(lvl) to cure this issue.");

	//	remember minimal level, that contains a surface index on this proc
		m_LocalFullRefLevel = std::min(m_LocalFullRefLevel, level);

	//	extract all algebra indices for the element on surface and level
		vLevelDD[level]->inner_algebra_indices(elem, vLevelInd);
		surfDD->inner_algebra_indices(elem, vSurfInd);
		UG_ASSERT(vSurfInd.size() == vLevelInd.size(), "Number of indices does not match.");

	//	set mapping index
		for(size_t i = 0; i < vSurfInd.size(); ++i)
		{
		//  check if some shadowing has already set the value
			if(bHigherLevWins){
				if(vSurfToLevelMap[vSurfInd[i]].level > level) continue;
			} else {
				if(vSurfToLevelMap[vSurfInd[i]].level >= 0 &&
				   vSurfToLevelMap[vSurfInd[i]].level < level) continue;
			}

		//	store index + level
			vSurfToLevelMap[vSurfInd[i]] = LevelIndex(vLevelInd[i], level);
		}
	}
}


template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_surface_to_level_mapping(std::vector<LevelIndex>& vSurfToLevelMap,
                              bool bHigherLevWins)
{
	PROFILE_FUNC_GROUP("gmg");
	ConstSmartPtr<DoFDistribution> surfDD =
			m_spApproxSpace->dof_distribution(GridLevel(m_surfaceLev, GridLevel::SURFACE));

	vSurfToLevelMap.resize(0);
	vSurfToLevelMap.resize(surfDD->num_indices());
	m_LocalFullRefLevel = std::numeric_limits<int>::max();

	if(surfDD->max_dofs(VERTEX)) init_surface_to_level_mapping<VertexBase>(vSurfToLevelMap, bHigherLevWins);
	if(surfDD->max_dofs(EDGE))   init_surface_to_level_mapping<EdgeBase>(vSurfToLevelMap, bHigherLevWins);
	if(surfDD->max_dofs(FACE))   init_surface_to_level_mapping<Face>(vSurfToLevelMap, bHigherLevWins);
	if(surfDD->max_dofs(VOLUME)) init_surface_to_level_mapping<Volume>(vSurfToLevelMap, bHigherLevWins);

	if(m_baseLev > m_LocalFullRefLevel)
		UG_THROW("GMG: Base level "<<m_baseLev<<" does not cover whole grid. "
		         <<"Highest full-ref level is "<<m_LocalFullRefLevel);
}


template <typename TDomain, typename TAlgebra>
template <typename TElem>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_noghost_to_ghost_mapping(	std::vector<size_t>& vNoGhostToGhostMap,
								ConstSmartPtr<DoFDistribution> spNoGhostDD,
								ConstSmartPtr<DoFDistribution> spGhostDD)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator iter_type;
	iter_type iter = spNoGhostDD->begin<TElem>();
	iter_type iterEnd = spNoGhostDD->end<TElem>();

//	vector of indices
	std::vector<size_t> vGhostInd, vNoGhostInd;

//	loop all elements of type
	for( ; iter != iterEnd; ++iter){
	//	get elem
		TElem* elem = *iter;

	//	extract all algebra indices for the element
		  spGhostDD->inner_algebra_indices(elem, vGhostInd);
		spNoGhostDD->inner_algebra_indices(elem, vNoGhostInd);
		UG_ASSERT(vGhostInd.size() == vNoGhostInd.size(), "Number of indices does not match.");

	//	set mapping index
		for(size_t i = 0; i < vNoGhostInd.size(); ++i){
			vNoGhostToGhostMap[vNoGhostInd[i]] = vGhostInd[i];
		}
	}
}


template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_noghost_to_ghost_mapping(int lev)
{
	PROFILE_FUNC_GROUP("gmg");

	LevData& ld = *m_vLevData[lev];
	std::vector<size_t>& vMapPatchToGlobal = ld.vMapPatchToGlobal;
	vMapPatchToGlobal.resize(0);

	if(ld.st == ld.t) return;

	ConstSmartPtr<DoFDistribution> spNoGhostDD = ld.st->dd();
	ConstSmartPtr<DoFDistribution> spGhostDD = ld.t->dd();

	vMapPatchToGlobal.resize(spNoGhostDD->num_indices());

	if(spNoGhostDD->max_dofs(VERTEX)) init_noghost_to_ghost_mapping<VertexBase>(vMapPatchToGlobal, spNoGhostDD, spGhostDD);
	if(spNoGhostDD->max_dofs(EDGE))   init_noghost_to_ghost_mapping<EdgeBase>(vMapPatchToGlobal, spNoGhostDD, spGhostDD);
	if(spNoGhostDD->max_dofs(FACE))   init_noghost_to_ghost_mapping<Face>(vMapPatchToGlobal, spNoGhostDD, spGhostDD);
	if(spNoGhostDD->max_dofs(VOLUME)) init_noghost_to_ghost_mapping<Volume>(vMapPatchToGlobal, spNoGhostDD, spGhostDD);
}


template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
copy_ghost_to_noghost(SmartPtr<GF> spVecTo,
                      ConstSmartPtr<GF> spVecFrom,
                      const std::vector<size_t>& vMapPatchToGlobal)
{
	PROFILE_FUNC_GROUP("gmg");

	if(spVecTo == spVecFrom) return;

	if(spVecTo->grid_level() == spVecFrom->grid_level())
	{
		for(size_t i = 0; i < spVecTo->size(); ++i)
			(*spVecTo)[i] = (*spVecFrom)[i];
	}
	else
	{
		UG_ASSERT(vMapPatchToGlobal.size() == spVecTo->size(),
				  "Mapping range ("<<vMapPatchToGlobal.size()<<") != "
				  "To-Vec-Size ("<<spVecTo->size()<<")");

		for(size_t i = 0; i < vMapPatchToGlobal.size(); ++i)
			(*spVecTo)[i] = (*spVecFrom)[ vMapPatchToGlobal[i] ];
	}

	#ifdef UG_PARALLEL
	spVecTo->set_storage_type(spVecFrom->get_storage_mask());
	#endif
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
copy_noghost_to_ghost(SmartPtr<GF> spVecTo,
                      ConstSmartPtr<GF> spVecFrom,
                      const std::vector<size_t>& vMapPatchToGlobal)
{
	PROFILE_FUNC_GROUP("gmg");

	if(spVecTo == spVecFrom) return;

	if(spVecTo->grid_level() == spVecFrom->grid_level())
	{
		for(size_t i = 0; i < spVecTo->size(); ++i)
			(*spVecTo)[i] = (*spVecFrom)[i];
	}
	else
	{
		UG_ASSERT(vMapPatchToGlobal.size() == spVecFrom->size(),
				  "Mapping domain ("<<vMapPatchToGlobal.size()<<") != "
				  "From-Vec-Size ("<<spVecFrom->size()<<")");

		for(size_t i = 0; i < vMapPatchToGlobal.size(); ++i)
			(*spVecTo)[ vMapPatchToGlobal[i] ] = (*spVecFrom)[i];
	}
	#ifdef UG_PARALLEL
	spVecTo->set_storage_type(spVecFrom->get_storage_mask());
	#endif
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
copy_noghost_to_ghost(SmartPtr<matrix_type> spMatTo,
                      ConstSmartPtr<matrix_type> spMatFrom,
                      const std::vector<size_t>& vMapPatchToGlobal)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_ASSERT(vMapPatchToGlobal.size() == spMatFrom->num_rows(),
			  "Mapping domain ("<<vMapPatchToGlobal.size()<<") != "
			  "From-Mat-Num-Rows ("<<spMatFrom->num_rows()<<")");
	UG_ASSERT(vMapPatchToGlobal.size() == spMatFrom->num_cols(),
			  "Mapping domain ("<<vMapPatchToGlobal.size()<<") != "
			  "From-Mat-Num-Cols ("<<spMatFrom->num_cols()<<")");

	for(size_t noghostFrom = 0; noghostFrom < vMapPatchToGlobal.size(); ++noghostFrom)
	{
		typedef typename matrix_type::const_row_iterator row_iter;
		row_iter conn = spMatFrom->begin_row(noghostFrom);
		row_iter connEnd = spMatFrom->end_row(noghostFrom);
		for(; conn != connEnd; ++conn)
		{
			const typename matrix_type::value_type& block = conn.value();
			size_t noghostTo = conn.index();

			(*spMatTo)(vMapPatchToGlobal[noghostFrom], vMapPatchToGlobal[noghostTo])
					= block;
		}
	}

	#ifdef UG_PARALLEL
	spMatTo->set_storage_type(spMatFrom->get_storage_mask());
	#endif
}

template <typename TDomain, typename TAlgebra>
template <typename TElem>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
collect_shadowing_indices(std::vector<size_t>& vShadowing,
                          ConstSmartPtr<DoFDistribution> spDD,
                          std::vector<size_t>& vSurfShadowing,
                          ConstSmartPtr<DoFDistribution> spSurfDD)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator iter_type;
	iter_type iter = spDD->begin<TElem>();
	iter_type iterEnd = spDD->end<TElem>();
	ConstSmartPtr<SurfaceView> spSurfView = m_spApproxSpace->surface_view();

//	vector of indices
	std::vector<size_t> vInd, vSurfInd;

//	loop all elements of type
	for( ; iter != iterEnd; ++iter){
	//	get elem
		TElem* elem = *iter;

	//	check if shadowing
		if(!spSurfView->is_shadowing(elem)) continue;

	// 	get global indices
		spDD->inner_algebra_indices(elem, vInd);
		spSurfDD->inner_algebra_indices(elem, vSurfInd);

	//	add to list of shadowing indices
		for(size_t i = 0; i < vInd.size(); ++i){
			vShadowing.push_back(vInd[i]);
			vSurfShadowing.push_back(vSurfInd[i]);
		}
	}
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
collect_shadowing_indices(int lev)
{
	PROFILE_FUNC_GROUP("gmg");
	ConstSmartPtr<DoFDistribution> spDD =
			m_spApproxSpace->dof_distribution(GridLevel(lev, m_GridLevelType, false));
	ConstSmartPtr<DoFDistribution> spSurfDD =
			m_spApproxSpace->dof_distribution(GridLevel(m_surfaceLev, GridLevel::SURFACE));

	std::vector<size_t>& vShadowing = m_vLevData[lev]->vShadowing;
	std::vector<size_t>& vSurfShadowing = m_vLevData[lev]->vSurfShadowing;
	vShadowing.clear();
	vSurfShadowing.clear();

	if(lev >= m_LocalFullRefLevel && m_GridLevelType != GridLevel::SURFACE) {
		if(spDD->max_dofs(VERTEX)) collect_shadowing_indices<VertexBase>(vShadowing, spDD, vSurfShadowing, spSurfDD);
		if(spDD->max_dofs(EDGE))   collect_shadowing_indices<EdgeBase>(vShadowing, spDD, vSurfShadowing, spSurfDD);
		if(spDD->max_dofs(FACE))   collect_shadowing_indices<Face>(vShadowing, spDD, vSurfShadowing, spSurfDD);
		if(spDD->max_dofs(VOLUME)) collect_shadowing_indices<Volume>(vShadowing, spDD, vSurfShadowing, spSurfDD);
	}
}


////////////////////////////////////////////////////////////////////////////////
// Init Level Data
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
init_level_memory(int baseLev, int topLev)
{
	PROFILE_FUNC_GROUP("gmg");

	m_vLevData.resize(0);
	m_vLevData.resize(topLev+1);

	for(int lev = baseLev; lev <= topLev; ++lev)
	{
		m_vLevData[lev] = SmartPtr<LevData>(new LevData);
		LevData& ld = *m_vLevData[lev];

		GridLevel gl = GridLevel(lev, m_GridLevelType, false);
		ld.sc = SmartPtr<GF>(new GF(m_spApproxSpace, gl, false));
		ld.sd = SmartPtr<GF>(new GF(m_spApproxSpace, gl, false));
		ld.st = SmartPtr<GF>(new GF(m_spApproxSpace, gl, false));

		// TODO: all procs must call this, since a MPI-Group is created.
		//		 Think about optimizing this
		GridLevel glGhosts = GridLevel(lev, m_GridLevelType, true);
		ld.t = SmartPtr<GF>(new GF(m_spApproxSpace, glGhosts, false));

		if(!m_spApproxSpace->might_contain_ghosts(lev))
			ld.t = ld.st;

		ld.A = SmartPtr<MatrixOperator<matrix_type, vector_type> >(
				new MatrixOperator<matrix_type, vector_type>);

		ld.PreSmoother = m_spPreSmootherPrototype->clone();
		if(m_spPreSmootherPrototype == m_spPostSmootherPrototype)
			ld.PostSmoother = ld.PreSmoother;
		else
			ld.PostSmoother = m_spPostSmootherPrototype->clone();

		ld.Projection = m_spProjectionPrototype->clone();

		ld.Prolongation = m_spProlongationPrototype->clone();
		if(m_spProlongationPrototype == m_spRestrictionPrototype)
			ld.Restriction = ld.Prolongation;
		else
			ld.Restriction = m_spRestrictionPrototype->clone();

		init_noghost_to_ghost_mapping(lev);
	}

	LevData& ld = *m_vLevData[baseLev];
	bool bHasVertMaster = false;
	bool bHasVertSlave = false;
	bool bHasHorrMaster = false;
	bool bHasHorrSlave = false;
	#ifdef UG_PARALLEL
	if( !ld.t->layouts()->vertical_master().empty()) bHasVertMaster = true;
	if( !ld.t->layouts()->vertical_slave().empty()) bHasVertSlave = true;
	if( !ld.t->layouts()->master().empty()) bHasHorrMaster = true;
	if( !ld.t->layouts()->slave().empty()) bHasHorrSlave = true;
	#endif
	bool bHasVertConn = bHasVertMaster || bHasVertSlave;
	bool bHasHorrConn = bHasHorrMaster || bHasHorrSlave;

	m_bGatheredBaseUsed = m_bGatheredBaseIfAmbiguous;

//	if no vert-interfaces are present (especially in serial or on a proc that
//	does not have parts of the base level), we cannot gather. In this case we
//	overrule the choice of the user
//	Note: levels not containing any dof, are skipped from computation anyway
	if(!bHasVertConn) m_bGatheredBaseUsed = false;

//	check if parallel solver is available, if not, try to use gathered
	if(!m_bGatheredBaseUsed
		&& bHasHorrConn
		&& !m_spBaseSolver->supports_parallel())
	{
		if(!bHasVertConn)
			UG_THROW("GMG: base level distributed in parallel, without possibility"
					" to gather on one process, but chosen base solver is not"
					" parallel. Choose a parallel base solver.")

		m_bGatheredBaseUsed = true;
	}

//	check if gathering base solver possible: If some horizontal layouts are
//	given, we know, that still the grid is distributed. But, if no
//	vertical layouts are present in addition, we can not gather the vectors
//	to on proc.
	if(m_bGatheredBaseUsed){
		if(bHasVertMaster && bHasVertSlave)
			UG_THROW("GMG: Gathered base solver expects one proc to have all "
					"v-masters and no v-slaves. Use gmg:set_gathered_base_solver_if_ambiguous(false)"
					" to avoid this error.");

		if(!bHasVertConn && bHasHorrConn){
			UG_THROW("GMG: Base level distributed among processes and no possibility"
					" of gathering (vert. interfaces) present, but a gathering"
					" base solving is required. Use gmg:set_gathered_base_solver_if_ambiguous(false)"
					" to avoid this error.");
		}

#ifdef UG_PARALLEL
		if(bHasVertSlave && (ld.t->layouts()->vertical_slave().num_interfaces() != 1))
			UG_THROW("GMG: Gathered base solver expects a v-slave containing "
					"process to send all its values to exactly on master. But "
					"more then one master detected. Use gmg:set_gathered_base_solver_if_ambiguous(false)"
					" to avoid this error.");

		if(bHasVertSlave && (ld.t->layouts()->vertical_slave().num_interface_elements() != ld.t->size()))
			UG_THROW("GMG: Gathered base solver expects that all indices"
					"are sent to exactly on master. But not all indices are "
					"v-slaves. Use gmg:set_gathered_base_solver_if_ambiguous(false)"
					" to avoid this error.");
#endif
	}

	if(m_bGatheredBaseUsed && gathered_base_master()){
		// note: we can safely call this here, since the dd has already been
		//		created in the level loop
		GridLevel glGhosts = GridLevel(baseLev, m_GridLevelType, true);
		spGatheredBaseCorr = SmartPtr<GF>(new GF(m_spApproxSpace, glGhosts, false));

		spGatheredBaseMat = SmartPtr<MatrixOperator<matrix_type, vector_type> >(
								new MatrixOperator<matrix_type, vector_type>);
	} else {
		spGatheredBaseCorr = NULL;
		spGatheredBaseMat = NULL;
	}
}

template <typename TDomain, typename TAlgebra>
bool AssembledMultiGridCycle<TDomain, TAlgebra>::
gathered_base_master() const
{
	if(!m_bGatheredBaseUsed)
		UG_THROW("GMG: Should only be called when gather-base solver used.")

	#ifdef UG_PARALLEL
	if(!m_vLevData[m_baseLev]->t->layouts()->vertical_master().empty())
		return true;
	#endif

	return false;
}


////////////////////////////////////////////////////////////////////////////////
// Cycle - Methods
////////////////////////////////////////////////////////////////////////////////

// perform the smoothing
template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
smooth(SmartPtr<GF> sc, SmartPtr<GF> sd, SmartPtr<GF> st,
       MatrixOperator<matrix_type, vector_type>& A,
       ILinearIterator<vector_type>& S,
       int lev, int nu)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - smooth on level "<<lev<<"\n");

//	since no parallel communication takes place in this method, we may
//	return immediately, if d is empty.
	if(sd->size() == 0) return;

// 	smooth nu times
	for(int i = 0; i < nu; ++i)
	{
	// 	Compute Correction of one smoothing step, but do not update defect
	//	a)  Compute t = B*d with some iterator B
		if(!S.apply(*st, *sd))
			UG_THROW("GMG::smooth: Smoothing step "<<i+1<<" on level "<<lev<<" failed.");

	//	b) reset the correction to zero on the patch boundary.
		const std::vector<size_t>& vShadowing = m_vLevData[lev]->vShadowing;
		for(size_t i = 0; i < vShadowing.size(); ++i)
			(*st)[ vShadowing[i] ] = 0.0;

	//	c) update the defect with this correction ...
		A.apply_sub(*sd, *st);

	//	d) ... and add the correction to the overall correction
		(*sc) += (*st);
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - smooth on level "<<lev<<"\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
presmooth(int lev)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - presmooth on level "<<lev<<"\n");

	GMG_PROFILE_BEGIN(GMG_PreSmooth);
	try{
		LevData& ld = *m_vLevData[lev];
		smooth(ld.sc, ld.sd, ld.st, *ld.A, *ld.PreSmoother, lev, m_numPreSmooth);
	}
	UG_CATCH_THROW("GMG::lmgc: Pre-Smoothing on level "<<lev<<" failed.");
	GMG_PROFILE_END();

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - presmooth on level "<<lev<<"\n");
}


template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
restriction(int lev)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - restriction on level "<<lev<<"\n");

	LevData& lf = *m_vLevData[lev];
	LevData& lc = *m_vLevData[lev-1];

//	## PARALLEL CASE: gather vertical
//	v-slaves indicate, that part of the parent are stored on a different proc.
//	Thus the values of the defect must be transfered to the v-masters and will
//	have their parents on the proc of the master and must be restricted on
//	that process. Thus we have to transfer the defect (currently stored
//	additive in the noghost-vectors) to the v-masters. And have to make sure
//	that d is additive after this operation. We do that by ensuring that it
//	is additive-unique regarding v-masters and v-slaves (v-slaves will be set to 0).
//	We use a temporary vector including ghost, such that the no-ghost defect
//	remains valid and can be used when the cycle comes back to this level.

	SmartPtr<GF> spD = lf.sd;
	#ifdef UG_PARALLEL
	if( !lf.t->layouts()->vertical_slave().empty() ||
		!lf.t->layouts()->vertical_master().empty())
	{
		SetLayoutValues(&(*lf.t), lf.t->layouts()->vertical_master(), 0);
		copy_noghost_to_ghost(lf.t, lf.sd, lf.vMapPatchToGlobal);
		if(lf.t->size() > 0){
			divide_vertical_slaves_by_number_of_masters(*lf.t);
			add_to_vertical_masters_and_set_zero_vertical_slaves(*lf.t);
		}
		spD = lf.t;
	}
	#endif

//	Now we can restrict the defect from the fine level to the coarser level.
	if((lc.sd->size() > 0) && (spD->size() > 0)){
		GMG_PROFILE_BEGIN(GMG_RestrictDefect);
		try{
			m_vLevData[lev]->Restriction->do_restrict(*lc.sd, *spD);
		}
		UG_CATCH_THROW("GMG::lmgc: Restriction of Defect from level "<<lev<<
		               " to "<<lev-1<<" failed.");
		GMG_PROFILE_END();

	//	apply post processes
		for(size_t i = 0; i < m_vspRestrictionPostProcess.size(); ++i)
			m_vspRestrictionPostProcess[i]->post_process(lc.sd);
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - restriction on level "<<lev<<"\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
prolongation(int lev)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - prolongation on level "<<lev<<"\n");

	LevData& lf = *m_vLevData[lev];
	LevData& lc = *m_vLevData[lev-1];

//	## ADAPTIVE CASE
//	in the adaptive case there is a small part of the coarse coupling that
//	has not been used to update the defect. In order to ensure, that the
//	defect on this level still corresponds to the updated defect, we need
//	to add if here.
	if(lev > m_LocalFullRefLevel){
		GMG_PROFILE_BEGIN(GMG_AddCoarseGridContribution);
		lc.CoarseGridContribution.matmul_minus(*lf.sd, *lc.sc);
		GMG_PROFILE_END();
	}

//	## INTERPOLATE CORRECTION
	GMG_PROFILE_BEGIN(GMG_InterpolateCorr);
	try{
		lf.Prolongation->prolongate(*lf.t, *lc.sc);
	}
	UG_CATCH_THROW("GMG::lmgc: Prolongation from level "<<lev-1<<" to "<<lev<<" failed.");
	GMG_PROFILE_END();

//	PARALLEL CASE: Receive values of correction for vertical slaves
//	If there are vertical slaves/masters on the coarser level, we now copy
//	the correction values from the v-master DoFs to the v-slave	DoFs.
//	since dummies may exist, we'll copy the broadcasted correction to h-slave
//	interfaces (dummies are always h-slaves)
	copy_to_vertical_slaves(*lf.t);

//	## PROJECT COARSE GRID CORRECTION ONTO SMOOTH AREA
	copy_ghost_to_noghost(lf.st, lf.t, lf.vMapPatchToGlobal);

//	apply post processes
	for(size_t i = 0; i < m_vspProlongationPostProcess.size(); ++i)
		m_vspProlongationPostProcess[i]->post_process(lf.st);

// 	## ADD COARSE GRID CORRECTION
	GMG_PROFILE_BEGIN(GMG_AddCoarseGridCorr);
	if(lf.sc->size() > 0){
		(*lf.sc) += (*lf.st);
	}
	GMG_PROFILE_END(); // GMG_AddCoarseGridCorr

//	correction changed c := c + t. Thus, update the defect d := d - A*t
	GMG_PROFILE_BEGIN(GMG_UpdateDefectForCGCorr);
	if(lf.sd->size() > 0){
		lf.A->apply_sub(*lf.sd, *lf.st);
	}
	GMG_PROFILE_END(); // GMG_UpdateDefectForCGCorr

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - prolongation on level "<<lev<<"\n");
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
postsmooth(int lev)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - postsmooth on level "<<lev<<"\n");

	GMG_PROFILE_BEGIN(GMG_PostSmooth);
	try{
		LevData& ld = *m_vLevData[lev];
		smooth(ld.sc, ld.sd, ld.st, *ld.A, *ld.PostSmoother, lev, m_numPostSmooth);
	}
	UG_CATCH_THROW("GMG::lmgc: Post-Smoothing on level "<<lev<<" failed. ")
	GMG_PROFILE_END();

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - postsmooth on level "<<lev<<"\n");
}

// performs the base solving
template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
base_solve(int lev)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - base_solve on level "<<lev<<"\n");
	try{

	LevData& ld = *m_vLevData[lev];

//	SOLVE BASE PROBLEM
//	Here we distinguish two possibilities:
//	a) The coarse grid problem is solved in parallel, using a parallel solver
//	b) First all vectors are gathered to one process, solved on this one
//	   process and then again distributed

//	CASE a): We solve the problem in parallel (or normally for sequential code)
	if(!m_bGatheredBaseUsed)
	{
		UG_DLOG(LIB_DISC_MULTIGRID, 3, " GMG: entering distributed basesolver branch.\n");
		if(ld.sd->num_indices()){

			GMG_PROFILE_BEGIN(GMG_BaseSolver);
			ld.sc->set(0.0);
			try{
				if(!m_spBaseSolver->apply(*ld.sc, *ld.sd))
					UG_THROW("GMG::lmgc: Base solver on base level "<<lev<<" failed.");
			}
			UG_CATCH_THROW("GMG: BaseSolver::apply failed. (case: a).")

		//	UPDATE DEFECT
		//	*) if adaptive case, we also need to update the defect, such that on the
		//	   surface level the defect remains updated
		//	*) Only for full refinement and real coarser level, we can forget about
		//	   the defect on the base level, since only the correction is needed
		//	   on the higher level
			if(lev >= m_LocalFullRefLevel){
				ld.A->apply_sub(*ld.sd, *ld.sc);
			}

			GMG_PROFILE_END();
		}
		UG_DLOG(LIB_DISC_MULTIGRID, 3, " GMG: exiting distributed basesolver branch.\n");
	}

//	CASE b): We gather the processes, solve on one proc and distribute again
	else
	{
		UG_DLOG(LIB_DISC_MULTIGRID, 3, " GMG: entering gathered basesolver branch.\n");

	//	gather the defect
		SmartPtr<GF> spD = ld.sd;
		#ifdef UG_PARALLEL
		if( !ld.t->layouts()->vertical_slave().empty() ||
			!ld.t->layouts()->vertical_master().empty())
		{
			SetLayoutValues(&(*ld.t), ld.t->layouts()->vertical_master(), 0);
			copy_noghost_to_ghost(ld.t, ld.sd, ld.vMapPatchToGlobal);
			divide_vertical_slaves_by_number_of_masters(*ld.t);
			add_to_vertical_masters_and_set_zero_vertical_slaves(*ld.t);
			spD = ld.t;
		}
		#endif

		pcl::SynchronizeProcesses();
	//	only solve on gathering master
		if(gathered_base_master())
		{
			GMG_PROFILE_BEGIN(GMG_BaseSolver);
			UG_DLOG(LIB_DISC_MULTIGRID, 3, " GMG: Start serial base solver.\n");

		//	Reset correction
			spGatheredBaseCorr->set(0.0);

		//	compute coarse correction
			try{
				if(!m_spBaseSolver->apply(*spGatheredBaseCorr, *spD))
					UG_THROW("GMG::lmgc: Base solver on base level "<<lev<<" failed.");
			}
			UG_CATCH_THROW("GMG: BaseSolver::apply failed. (case: b).")

			GMG_PROFILE_END();
			UG_DLOG(LIB_DISC_MULTIGRID, 3, " GMG gathered base solver done.\n");
		}
		pcl::SynchronizeProcesses();

	//	broadcast the correction
		#ifdef UG_PARALLEL
		if(!ld.sc->layouts()->vertical_slave().empty()){
			ComPol_VecCopy<vector_type> cpVecCopy(&(*ld.sc));
			m_Com.receive_data(ld.sc->layouts()->vertical_slave(), cpVecCopy);
			m_Com.communicate();
		}
		#endif
		if(gathered_base_master()){
			#ifdef UG_PARALLEL
			ComPol_VecCopy<vector_type> cpVecCopy(&(*spGatheredBaseCorr));
			m_Com.send_data(spGatheredBaseCorr->layouts()->vertical_master(), cpVecCopy);
			m_Com.communicate();
			spGatheredBaseCorr->set_storage_type(PST_CONSISTENT);
			#endif
			copy_ghost_to_noghost(ld.sc, spGatheredBaseCorr, ld.vMapPatchToGlobal);
		}

		UG_DLOG(LIB_DISC_MULTIGRID, 3, " GMG: exiting gathered basesolver branch.\n");
	}

	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - base_solve on level "<<lev<<"\n");

	}
	UG_CATCH_THROW("GMG: Base Solver failed.");
}

// performs a  multi grid cycle on the level
template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
lmgc(int lev)
{
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - lmgc on level "<<lev<<"\n");

//	switch, if base level is reached. If so, call base Solver, else we will
//	perform smoothing, restrict the defect and call the lower level; then,
//	going up again in the level hierarchy the correction is interpolated and
//	used as coarse grid correction. Finally a post-smooth is performed.
	if((int)lev > m_baseLev)
	{
		for(int i = 0; i < m_cycleType; ++i)
		{
			log_debug_data(lev, "BeforePreSmooth");
			try{
				presmooth(lev);
			}
			UG_CATCH_THROW("GMG::lmgc: presmooth failed on level "<<lev);

			log_debug_data(lev, "BeforeRestrict");
			try{
				restriction(lev);
			}
			UG_CATCH_THROW("GMG::lmgc: restriction failed on level "<<lev);

			try{
				lmgc(lev-1);
			}
			UG_CATCH_THROW("GMG::lmgc: Linear multi grid "
							"cycle on level "<<lev-1<<" failed. (BaseLev="<<
							m_baseLev<<", TopLev="<<m_topLev<<").");

			log_debug_data(lev, "BeforeProlong");
			try{
				prolongation(lev);
			}
			UG_CATCH_THROW("GMG::lmgc: prolongation failed on level "<<lev);

			log_debug_data(lev, "BeforePostSmooth");
			try{
				postsmooth(lev);
			}
			UG_CATCH_THROW("GMG::lmgc: postsmooth failed on level "<<lev);
			log_debug_data(lev, "AfterPostSmooth");
		}

		UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - lmgc on level "<<lev<<"\n");
	}

//	if the base level has been reached, the coarse problem is solved exactly
	else if((int)lev == m_baseLev)
	{
		log_debug_data(lev, "BeforeBaseSolver");
		try{
			base_solve(lev);
		}
		UG_CATCH_THROW("GMG::lmgc: basesolver failed on level "<<lev);
		log_debug_data(lev, "AfterBaseSolver");

		UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - lmgc on level "<<lev<<" (base solver executed)\n");
	}
//	this case should never happen.
	else {
		UG_THROW("GMG::lmgc: Level index below baseLevel.");
	}
}

////////////////////////////////////////////////////////////////////////////////
// Parallel Communication
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
divide_vertical_slaves_by_number_of_masters(vector_type& d)
{
#ifdef UG_PARALLEL
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - divide_vertical_slaves_by_number_of_masters\n");
	PROFILE_FUNC_GROUP("gmg");

	if(d.layouts()->vertical_slave().empty()) return;

	GMG_PROFILE_BEGIN(GMG_DevideSlavesByNumberOfMasters);

//	there may be v-slaves with multiple v-masters. We only want to send
//	a fraction to each master, to keep d additive.
//	count number of occurrances in v-interfaces
	bool multiOccurance = false;
	std::vector<number> occurence;
	const IndexLayout& layout = d.layouts()->vertical_slave();

	if(layout.num_interfaces() > 1){
		occurence.resize(d.size(), 0);
		for(IndexLayout::const_iterator iiter = layout.begin();
			iiter != layout.end(); ++iiter)
		{
			const IndexLayout::Interface& itfc = layout.interface(iiter);
			for(IndexLayout::Interface::const_iterator iter = itfc.begin();
				iter != itfc.end(); ++iter)
			{
				const IndexLayout::Interface::Element& index = itfc.get_element(iter);

				occurence[index] += 1;
				if(occurence[index] > 1)
					multiOccurance = true;
			}
		}
	}

	if(multiOccurance){
		for(size_t i = 0; i < occurence.size(); ++i){
			if(occurence[i] > 1)
				d[i] *= (1./occurence[i]);
		}
	}

	GMG_PROFILE_END();
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - divide_vertical_slaves_by_number_of_masters\n");
#endif
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
divide_vertical_slave_rows_by_number_of_masters(matrix_type& mat)
{
#ifdef UG_PARALLEL
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - divide_vertical_slave_rows_by_number_of_masters\n");
	PROFILE_FUNC_GROUP("gmg");

	if(mat.layouts()->vertical_slave().empty()) return;

	GMG_PROFILE_BEGIN(GMG_DevideSlaveRowByNumberOfMasters);

//	there may be v-slaves with multiple v-masters. We only want to send
//	a fraction to each master, to keep d additive.
//	count number of occurrances in v-interfaces
	bool multiOccurance = false;
	std::vector<number> occurence;
	const IndexLayout& layout = mat.layouts()->vertical_slave();

	if(layout.num_interfaces() > 1){
		occurence.resize(mat.num_rows(), 0);
		for(IndexLayout::const_iterator iiter = layout.begin();
			iiter != layout.end(); ++iiter)
		{
			const IndexLayout::Interface& itfc = layout.interface(iiter);
			for(IndexLayout::Interface::const_iterator iter = itfc.begin();
				iter != itfc.end(); ++iter)
			{
				const IndexLayout::Interface::Element& index = itfc.get_element(iter);

				occurence[index] += 1;
				if(occurence[index] > 1)
					multiOccurance = true;
			}
		}
	}

	if(multiOccurance){
		for(size_t i = 0; i < occurence.size(); ++i){
			if(occurence[i] > 1)
				ScaleRow(mat, i, (1./occurence[i]));
		}
	}

	GMG_PROFILE_END();
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - divide_vertical_slave_rows_by_number_of_masters\n");
#endif
}


template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
add_to_vertical_masters_and_set_zero_vertical_slaves(vector_type& d)
{
#ifdef UG_PARALLEL
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - add_to_vertical_masters_and_set_zero_vertical_slaves\n");

//	send vertical-slaves -> vertical-masters
	GMG_PROFILE_BEGIN(GMG_BroadcastVerticalVector);
	ComPol_VecAddSetZero<vector_type> cpVecAdd(&d);

	if(!d.layouts()->vertical_slave().empty())
		m_Com.send_data(d.layouts()->vertical_slave(), cpVecAdd);

	if(!d.layouts()->vertical_master().empty())
		m_Com.receive_data(d.layouts()->vertical_master(), cpVecAdd);

//	perform communication
	m_Com.communicate();

	GMG_PROFILE_END();
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - add_to_vertical_masters_and_set_zero_vertical_slaves\n");
#endif
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
copy_to_vertical_slaves(vector_type& c)
{
#ifdef UG_PARALLEL
	PROFILE_FUNC_GROUP("gmg");
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - copy_to_vertical_slaves\n");

//	send vertical-masters -> vertical-slaves
	GMG_PROFILE_BEGIN(GMG_BroadcastVerticalVector);
	ComPol_VecCopy<vector_type> cpVecCopy(&c);

	if(!c.layouts()->vertical_slave().empty())
		m_Com.receive_data(c.layouts()->vertical_slave(), cpVecCopy);

	if(!c.layouts()->vertical_master().empty())
		m_Com.send_data(c.layouts()->vertical_master(), cpVecCopy);

//	communicate
	m_Com.communicate();

	GMG_PROFILE_END();
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - copy_to_vertical_slaves\n");
#endif
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
copy_to_vertical_masters(vector_type& c)
{
#ifdef UG_PARALLEL
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-start - copy_to_vertical_masters\n");
	PROFILE_FUNC_GROUP("gmg");

//	send vertical-slaves -> vertical-masters
	GMG_PROFILE_BEGIN(GMG_CopyToVerticalMasters);
	ComPol_VecCopy<vector_type> cpVecCopy(&c);

	if(!c.layouts()->vertical_master().empty())
		m_Com.receive_data(c.layouts()->vertical_master(), cpVecCopy);

	if(!c.layouts()->vertical_slave().empty())
		m_Com.send_data(c.layouts()->vertical_slave(), cpVecCopy);

//	communicate
	m_Com.communicate();

	GMG_PROFILE_END();
	UG_DLOG(LIB_DISC_MULTIGRID, 3, "gmg-stop - copy_to_vertical_masters\n");
#endif
}

////////////////////////////////////////////////////////////////////////////////
// Debug Methods
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
write_debug(ConstSmartPtr<GF> spGF, std::string name)
{
	write_debug(*spGF, name);
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
write_debug(const GF& rGF, std::string name)
{
	PROFILE_FUNC_GROUP("debug");

	if(m_spDebugWriter.invalid()) return;

//	build name
	GridLevel gl = rGF.grid_level();
	std::stringstream ss; ss << std::setfill('0') << "GMG_" << name << "_";

	if(gl.is_level()){
		if(gl.ghosts()) ss << "glev" << std::setw(3) << gl.level();
		else 			ss << "lev" << std::setw(3) << gl.level();
	} else if (gl.is_surface()){
		if(gl.ghosts()) ss << "gsurf" << std::setw(3) << gl.level();
		else 			ss << "surf" << std::setw(3) << gl.level();
	} else {
		UG_THROW("GMG: GridLevel not supported.")
	}

	ss << "_i" << std::setw(3) << m_dbgIterCnt << ".vec";

//	write
	GridLevel currGL = m_spDebugWriter->grid_level();
	m_spDebugWriter->set_grid_level(gl);
	m_spDebugWriter->write_vector(rGF, ss.str().c_str());
	m_spDebugWriter->set_grid_level(currGL);
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
write_smooth_level_debug(const matrix_type& mat, std::string name, int lev)
{
	PROFILE_FUNC_GROUP("debug");

	if(m_spDebugWriter.invalid()) return;

//	build name
	std::stringstream ss; ss << std::setfill('0') << "GMG_" << name << "_";
	ss << "l" << std::setw(3) << lev;
	ss << "_i" << std::setw(3) << m_dbgIterCnt << ".mat";

//	write
	GridLevel currGL = m_spDebugWriter->grid_level();
	m_spDebugWriter->set_grid_level(GridLevel(lev,m_GridLevelType, false));
	m_spDebugWriter->write_matrix(mat, ss.str().c_str());
	m_spDebugWriter->set_grid_level(currGL);
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
write_level_debug(const matrix_type& mat, std::string name, int lev)
{
	PROFILE_FUNC_GROUP("debug");

	if(m_spDebugWriter.invalid()) return;

//	build name
	std::stringstream ss; ss << std::setfill('0') << "GMG_" << name << "_";
	ss << "gl" << std::setw(3) << lev;
	ss << "_i" << std::setw(3) << m_dbgIterCnt << ".mat";

//	write
	GridLevel currGL = m_spDebugWriter->grid_level();
	m_spDebugWriter->set_grid_level(GridLevel(lev,m_GridLevelType));
	m_spDebugWriter->write_matrix(mat, ss.str().c_str());
	m_spDebugWriter->set_grid_level(currGL);
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
write_surface_debug(const matrix_type& mat, std::string name)
{
	PROFILE_FUNC_GROUP("debug");

	if(m_spDebugWriter.invalid()) return;

//	build name
	std::stringstream ss; ss << std::setfill('0') << "GMG_" << name << "_";
	ss << "surf";
	ss << "_i" << std::setw(3) << m_dbgIterCnt << ".mat";

//	write
	GridLevel currGL = m_spDebugWriter->grid_level();
	m_spDebugWriter->set_grid_level(GridLevel(GridLevel::TOP, GridLevel::SURFACE));
	m_spDebugWriter->write_matrix(mat, ss.str().c_str());
	m_spDebugWriter->set_grid_level(currGL);
}

template <typename TDomain, typename TAlgebra>
void AssembledMultiGridCycle<TDomain, TAlgebra>::
log_debug_data(int lvl, std::string name)
{
	if(m_spDebugWriter.valid()){
		std::string defName("Def_"); defName.append(name);
		std::string curName("Cor_"); curName.append(name);
		write_debug(m_vLevData[lvl]->sd, defName);
		write_debug(m_vLevData[lvl]->sc, curName);
	}

	const bool bEnableSerialNorm = false;
	const bool bEnableParallelNorm = false;

	if(!bEnableSerialNorm && !bEnableParallelNorm) return;

	std::string prefix;
	if(lvl < (int)m_vLevData.size())
		prefix.assign(2 + 2 * (m_vLevData.size() - lvl), ' ');
	prefix.append(name).append(" on lev ").append(ToString(lvl)).append(": ");

	LevData& ld = *m_vLevData[lvl];
	if(bEnableSerialNorm){
		UG_LOG(prefix << "local sd norm: " << sqrt(VecProd(*ld.sd, *ld.sd)) << std::endl);
		UG_LOG(prefix << "local sc norm: " << sqrt(VecProd(*ld.sc, *ld.sc)) << std::endl);
	}
	if(bEnableParallelNorm){
	#ifdef UG_PARALLEL
		uint oldStorageMask = ld.t->get_storage_mask();
		number norm = ld.t->norm();
		UG_LOG(prefix << " t norm: " << norm << "\n");
		if(oldStorageMask & PST_ADDITIVE)
			ld.t->change_storage_type(PST_ADDITIVE);
		else if(oldStorageMask & PST_CONSISTENT)
			ld.t->change_storage_type(PST_CONSISTENT);

		oldStorageMask = ld.sd->get_storage_mask();
		norm = ld.sd->norm();
		UG_LOG(prefix << "sd norm: " << norm << "\n");
		if(oldStorageMask & PST_ADDITIVE)
			ld.sd->change_storage_type(PST_ADDITIVE);
		else if(oldStorageMask & PST_CONSISTENT)
			ld.sd->change_storage_type(PST_CONSISTENT);

		oldStorageMask = ld.sc->get_storage_mask();
		norm = ld.sc->norm();
		UG_LOG(prefix << "sc norm: " << norm << "\n");
		if(oldStorageMask & PST_ADDITIVE)
			ld.sc->change_storage_type(PST_ADDITIVE);
		else if(oldStorageMask & PST_CONSISTENT)
			ld.sc->change_storage_type(PST_CONSISTENT);

		oldStorageMask = ld.st->get_storage_mask();
		norm = ld.st->norm();
		UG_LOG(prefix << "st norm: " << norm << "\n");
		if(oldStorageMask & PST_ADDITIVE)
			ld.st->change_storage_type(PST_ADDITIVE);
		else if(oldStorageMask & PST_CONSISTENT)
			ld.st->change_storage_type(PST_CONSISTENT);
	#endif
	}
}

template <typename TDomain, typename TAlgebra>
std::string
AssembledMultiGridCycle<TDomain, TAlgebra>::
config_string() const
{
	std::stringstream ss;
	ss << "GeometricMultigrid (";
	if(m_cycleType == 1) ss << "V-Cycle";
	else if(m_cycleType == 2) ss << "W-Cycle";
	else ss << " " << m_cycleType << "-Cycle";
	ss << ")\n";

	if(m_spPreSmootherPrototype==m_spPostSmootherPrototype)
		ss 	<< " Smoother (" << m_numPreSmooth << "x pre, " << m_numPostSmooth << "x post): "
			<< ConfigShift(m_spPreSmootherPrototype->config_string());
	else
	{
		ss << " Presmoother (" << m_numPreSmooth << "x): " << ConfigShift(m_spPreSmootherPrototype->config_string());
		ss << " Postsmoother ( " << m_numPostSmooth << "x): " << ConfigShift(m_spPostSmootherPrototype->config_string());
	}
	ss << "\n";
	ss << " Basesolver ( Baselevel = " << m_baseLev << ", gathered base = " << (m_bGatheredBaseIfAmbiguous ? "true" : "false") << "): ";
	ss << ConfigShift(m_spBaseSolver->config_string());
	return ss.str();

}

} // namespace ug


#endif /* __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__ */
