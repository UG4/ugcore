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
#include "mg_solver_util.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_disc/dof_manager/dof_manager_util.h"

#include "mg_solver.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"

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
	#define GMG_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define GMG_PROFILE_END()		PROFILE_END()
#else
	#define GMG_PROFILE_FUNC()
	#define GMG_PROFILE_BEGIN(name)
	#define GMG_PROFILE_END()
#endif

namespace ug{

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
apply(vector_type &c, const vector_type& d)
{
//	temporary vector for defect
	vector_type dTmp; dTmp.resize(d.size());

//	copy defect
	dTmp = d;

//	work on copy
	return apply_update_defect(c, dTmp);
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
apply_update_defect(vector_type &c, vector_type& d)
{
//	get MultiGrid for some checks
	typename TApproximationSpace::domain_type::grid_type& mg =
				m_pApproxSpace->get_domain().get_grid();

// 	Check if grid type is a Multigrid.
	if(dynamic_cast<MultiGrid*>(&mg) == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect':"
				"MultiGridSolver requires a Multigrid. Please use MultiGrid structure.\n");
		return false;
	}

// 	Check if surface level has been chosen correctly
//	Please note, that the approximation space returns the global number of levels,
//	i.e. the maximum of levels among all processes.
	if(m_topLev >= m_pApproxSpace->num_levels())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect':"
				" SurfaceLevel " << m_topLev << " does not exist.\n");
		return false;
	}

// 	Check if base level has been choose correctly
	if(m_baseLev > m_topLev)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect':"
				"Base level must be smaller or equal to surface Level.\n");
		return false;
	}

//	project defect from surface to level
	GMG_PROFILE_BEGIN(GMG_ProjectDefectFromSurface);
	if(!project_surface_to_level(level_defects(), d))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of defect to level failed.\n");
		return false;
	}
	GMG_PROFILE_END(); //GMGApply_ProjectDefectFromSurface

// 	Perform one multigrid cycle
//	At this point c, d are valid for m_vLevData[m_topLev]->c, m_vLevData[m_topLev]->d
	GMG_PROFILE_BEGIN(GMG_lmgc);
	if(!lmgc(m_topLev))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Cannot perform multi grid cycle on TopLevel "<<m_topLev<<".\n");
		return false;
	}
	GMG_PROFILE_END(); //GMGApply_lmgc

//	project defect from level to surface
	GMG_PROFILE_BEGIN(GMG_ProjectDefectFromLevelToSurface);
	if(!project_level_to_surface(d, const_level_defects()))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of defect to surface failed.\n");
		return false;
	}
	GMG_PROFILE_END(); //GMGApply_ProjectDefectFromLevelToSurface

//	project correction from level to surface
	GMG_PROFILE_BEGIN(GMG_ProjectCorrectionFromLevelToSurface);
	if(!project_level_to_surface(c, const_level_corrections()))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of correction to surface failed.\n");
		return false;
	}
	GMG_PROFILE_END(); //GMGApply_ProjectCorrectionFromLevelToSurface

//	increase dbg counter
	if(m_pDebugWriter) m_dbgIterCnt++;

//	we're done
	return true;
}

// perform the smoothing
template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
smooth(vector_type& c, vector_type& d, vector_type& tmp,
       MatrixOperator<vector_type, vector_type, matrix_type>& A,
       smoother_type& S,
       size_t lev, int nu)
{
// 	smooth nu times
	for(int i = 0; i < nu; ++i)
	{
	//	switch if adaptive case must be handled
		if(!m_bAdaptive)
		{
		// 	Compute Correction of one smoothing step:
		//	a)  Compute t = B*d with some iterator B
		//	b) 	Update Defect d := d - A * t
			if(!S.apply_update_defect(tmp, d))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Smoothing step "
						<< i+1 << " on level " << lev << " failed.\n");
				return false;
			}

		// 	add correction of smoothing step to level correction
		//	(Note: we do not work on c directly here, since we update the defect
		//	       after every smoothing step. The summed up correction corresponds
		//		   to the total correction of the whole smoothing.)
			c += tmp;
		}
		else
	//	This is the adaptive case. Here, we must ensure, that the added correction
	//	is zero on the adaptively refined patch boundary of this level
		{
		// 	Compute Correction of one smoothing step, but do not update defect
		//	a)  Compute t = B*d with some iterator B

			if(!S.apply(tmp, d))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Smoothing step "
						<< i+1 << " on level " << lev << " failed.\n");
				return false;
			}

		//	get surface view
			const SurfaceView& surfView = *m_pApproxSpace->get_surface_view();

		//	First we reset the correction to zero on the patch boundary.
			if(!SetZeroOnShadowing(tmp, *m_vLevData[lev]->pLevDD, surfView))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Could not "
						" reset the values to zero on patch boundary for correction "
						<< i+1 << " on level " << lev << ".\n");
				return false;
			}

		//	now, we can update the defect with this correction ...
			if(!A.apply_sub(d, tmp))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Could not "
						" update defect for patch correction in smooth step "
						<< i+1 << " on level " << lev << ".\n");
				return false;
			}

		//	... and add the correction to to overall correction
			c += tmp;
		}
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
presmooth(size_t lev)
{
//	Get all needed vectors and operators

//	get vectors used in smoothing operations. (This is needed if vertical
//	masters are present, since no smoothing is performed on those. In that case
//	only on a smaller part of the grid level - the smoothing patch - the
//	smoothing is performed)
	vector_type& sd = m_vLevData[lev]->get_smooth_defect();
	vector_type& sc = m_vLevData[lev]->get_smooth_correction();
	vector_type& sTmp = m_vLevData[lev]->get_smooth_tmp();

//	get smoother on this level and corresponding operator
	smoother_type& Smoother = m_vLevData[lev]->get_smoother();
	MatrixOperator<vector_type, vector_type, matrix_type>& SmoothMat =
		m_vLevData[lev]->get_smooth_mat();

// 	reset correction to zero on this level
	sc.set(0.0);

//	We start the multi grid cycle on this level by smoothing the defect. This
//	means that we compute a correction c, such that the defect is "smoother".
//	If ghosts are present in parallel, we only smooth on a patch. Thus we first
//	copy the values from the whole grid level to the smoothing patch.
	m_vLevData[lev]->copy_defect_to_smooth_patch();

// 	pre-smoothing
	GMG_PROFILE_BEGIN(GMG_PreSmooth);
	GMG_PARALLEL_DEBUG_BARRIER(sd.get_process_communicator());
	if(!smooth(sc, sd, sTmp, SmoothMat, Smoother, lev, m_numPreSmooth))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Pre-Smoothing on "
				"level " << lev << " failed. "
				"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
		return false;
	}
	GMG_PROFILE_END();

//	now copy the values of d back to the whole grid, since the restriction
//	acts on the whole grid. Since we will perform a addition of the vertical
//	slaves to the vertical masters, we now set the values on the ghosts to
//	zero.
	m_vLevData[lev]->copy_defect_from_smooth_patch(true);

	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
bool AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
restriction(size_t lev, bool* restrictionPerformedOut)
{
//	Get all needed vectors and operators

//	Get vectors defined on whole grid (including ghosts) on this level
//	denote by: c = Correction, d = Defect, tmp = Help vector
	vector_type& d = m_vLevData[lev]->d;

//	Lets get a reference to the coarser level correction, defect, help vector
	UG_ASSERT(lev > 0, "restriction can't be applied on level 0.");
	vector_type& cd = m_vLevData[lev-1]->d;

//	## PARALLEL CASE: gather vertical
//	Send vertical slave values to master and check resuming.
//	If there are vertical slaves/masters on the coarser level, we now copy
//	the restricted values of the defect from the slave DoFs to the master
//	DoFs.
	#ifdef UG_PARALLEL
//	Resume flag: Check if process should continue lmgc. This is the case if
//				 the process has still some grid level below this level. If
//				 there is no such level (e.g. if the process recieved an
//				 already refined grid level here; "vertical cut") the process
//				 stops execution until the other processes have performed the
//				 coarser levels and are back again at this level.
	if(!gather_vertical(d)){
	//	only continue if levels left
		*restrictionPerformedOut = false;
		return true;
	}
	#endif

//	Now we can restrict the defect from the fine level to the coarser level.
//	This is done using the transposed prolongation.
	GMG_PROFILE_BEGIN(GMG_RestrictDefect);
	if(!m_vLevData[lev]->Prolongation->apply_transposed(cd, d))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Restriction of "
				"Defect from level "<<lev<<" to "<<lev-1<<" failed. "
				"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
		return false;
	}
	GMG_PROFILE_END();

//	since we reached this point, the restriction was performed.
	*restrictionPerformedOut = true;

	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
prolongation(size_t lev, bool restrictionWasPerformed)
{
//	Get all needed vectors and operators

//	Get vectors defined on whole grid (including ghosts) on this level
//	denote by: c = Correction, d = Defect, tmp = Help vector
	vector_type& d = m_vLevData[lev]->d;
	vector_type& tmp = m_vLevData[lev]->t;

//	get vectors used in smoothing operations. (This is needed if vertical
//	masters are present, since no smoothing is performed on those. In that case
//	only on a smaller part of the grid level - the smoothing patch - the
//	smoothing is performed)
	vector_type& sd = m_vLevData[lev]->get_smooth_defect();
	vector_type& sc = m_vLevData[lev]->get_smooth_correction();
	vector_type& sTmp = m_vLevData[lev]->get_smooth_tmp();

//	Lets get a reference to the coarser level correction, defect, help vector
	UG_ASSERT(lev > 0, "prolongatoin can't be applied on level 0.");
	vector_type& cc = m_vLevData[lev-1]->c;
	vector_type& cTmp = m_vLevData[lev-1]->t;

//	get smoothing operator on this level
	MatrixOperator<vector_type, vector_type, matrix_type>& SmoothMat =
		m_vLevData[lev]->get_smooth_mat();

//	## INTERPOLATE CORRECTION
	if(restrictionWasPerformed){
	//	now we can interpolate the coarse grid correction from the coarse level
	//	to the fine level
		GMG_PROFILE_BEGIN(GMG_InterpolateCorr);
		if(!m_vLevData[lev]->Prolongation->apply(tmp, cc))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Prolongation from"
					" level " << lev-1 << " to " << lev << " failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");

			return false;
		}
		GMG_PROFILE_END();
	}

//	PARALLEL CASE: Receive values of correction for vertical slaves
//	If there are vertical slaves/masters on the coarser level, we now copy
//	the correction values from the master DoFs to the slave	DoFs.
	#ifdef UG_PARALLEL
	broadcast_vertical(tmp);
	#endif

//	## PROJECT COARSE GRID CORRECTION ONTO SMOOTH AREA
	m_vLevData[lev]->copy_tmp_to_smooth_patch();

// 	## ADD COARSE GRID CORRECTION
	GMG_PROFILE_BEGIN(GMG_AddCoarseGridCorr);
	sc += sTmp;
	GMG_PROFILE_END(); // GMG_AddCoarseGridCorr

//	## UPDATE DEFECT FOR COARSE GRID CORRECTION
//	the correction has changed c := c + t. Thus, we also have to update
//	the defect d := d - A*t
	GMG_PROFILE_BEGIN(GMG_UpdateDefectForCGCorr);
	if(!SmoothMat.apply_sub(sd, sTmp))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Updating of defect"
				" on level " << lev << " failed.\n");
		GMG_PROFILE_END(); // GMG_UpdateDefectForCGCorr
		return false;
	}
	GMG_PROFILE_END(); // GMG_UpdateDefectForCGCorr

//	## ADAPTIVE CASE
	if(m_bAdaptive)
	{
	//	in the adaptive case there is a small part of the coarse coupling that
	//	has not been used to update the defect. In order to ensure, that the
	//	defect on this level still corresponds to the updated defect, we need
	//	to add if here. This is done in three steps:
	//	a) Compute the coarse update of the defect induced by missing coupling
		cTmp.set(0.0);
		if(!m_vLevData[lev-1]->CoarseGridContribution.apply(cTmp, cc))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Could not compute"
					" missing update defect contribution on level "<<lev-1<<".\n");
			return false;
		}

	//	get surface view
		const SurfaceView& surfView = *m_pApproxSpace->get_surface_view();

	//	b) interpolate the coarse defect up
		if(!AddProjectionOfShadows(d, cTmp, -1.0,
		                           *m_vLevData[lev]->pLevDD, *m_vLevData[lev-1]->pLevDD,
		                           surfView))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Could not add"
					" missing update defect contribution to level "<<lev<<".\n");
			return false;
		}
	}
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
postsmooth(size_t lev)
{
//	get vectors used in smoothing operations. (This is needed if vertical
//	masters are present, since no smoothing is performed on those. In that case
//	only on a smaller part of the grid level - the smoothing patch - the
//	smoothing is performed)
	vector_type& sd = m_vLevData[lev]->get_smooth_defect();
	vector_type& sc = m_vLevData[lev]->get_smooth_correction();
	vector_type& sTmp = m_vLevData[lev]->get_smooth_tmp();

//	get smoother on this level and corresponding operator
	smoother_type& Smoother = m_vLevData[lev]->get_smoother();
	MatrixOperator<vector_type, vector_type, matrix_type>& SmoothMat =
		m_vLevData[lev]->get_smooth_mat();


// 	## POST-SMOOTHING
//	We smooth the updated defect again. This means that we compute a
//	correction c, such that the defect is "smoother".
	GMG_PROFILE_BEGIN(GMG_PostSmooth);
	GMG_PARALLEL_DEBUG_BARRIER(sd.get_process_communicator());
	if(!smooth(sc, sd, sTmp, SmoothMat, Smoother, lev, m_numPostSmooth))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Post-Smoothing on"
				" level " << lev << " failed. "
				"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
		return false;
	}
	GMG_PROFILE_END();

//	## PROJECT DEFECT, CORRECTION BACK TO WHOLE GRID FOR RESTRICTION
	m_vLevData[lev]->copy_defect_from_smooth_patch();
	m_vLevData[lev]->copy_correction_from_smooth_patch();

	return true;
}

// performs the base solving
template <typename TApproximationSpace, typename TAlgebra>
bool AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
base_solve(size_t lev)
{
//	get vectors used in smoothing operations. (This is needed if vertical
//	masters are present, since no smoothing is performed on those. In that case
//	only on a smaller part of the grid level - the smoothing patch - the
//	smoothing is performed)
	vector_type& sd = m_vLevData[lev]->get_smooth_defect();
	vector_type& sc = m_vLevData[lev]->get_smooth_correction();

//	SOLVE BASE PROBLEM
//	Here we distinguish two possibilities:
//	a) The coarse grid problem is solved in parallel, using a parallel solver
//	b) First all vectors are gathered to one process, solved on this one
//	   process and then again distributed

//	CASE a): We solve the problem in parallel (or normally for sequential code)
#ifdef UG_PARALLEL
//	vector defined on whole grid (including ghosts) on this level
	vector_type& d = m_vLevData[lev]->d;

	UG_DLOG(LIB_DISC_MULTIGRID, 2, "GMG: Start BaseSolver on level "<<lev<<".\n");
	if( m_bBaseParallel ||
	   (d.get_vertical_slave_layout().empty() &&
		d.get_vertical_master_layout().empty()))
	{
#endif

	//	LIFTING c TO SOLVING AREA
		m_vLevData[lev]->copy_defect_to_smooth_patch();

		GMG_PROFILE_BEGIN(GMG_BaseSolver);
		sc.set(0.0);
		if(!m_pBaseSolver->apply(sc, sd))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Base solver on"
					" base level " << lev << " failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");

			return false;
		}

	//	*) if baseLevel == surfaceLevel, we need also need the updated defect
	//	*) if adaptive case, we also need to update the defect, such that on the
	//	   surface level the defect remains updated
	//	*) Only for full refinement and real coarser level, we can forget about
	//	   the defect on the base level, since only the correction is needed
	//	   on the higher level
		if(m_baseLev == m_topLev || m_bAdaptive)
		{
		//	get smoothing matrix
			MatrixOperator<vector_type, vector_type, matrix_type>& SmoothMat
				= m_vLevData[lev]->get_smooth_mat();

		//	UPDATE DEFECT
			if(!SmoothMat.apply_sub(sd, sc))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Updating defect "
						" on base level " << lev << ". "
						"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
				return false;
			}

		//	copy back to whole grid
			m_vLevData[lev]->copy_defect_from_smooth_patch(true);
		}

	//	PROJECT CORRECTION BACK TO WHOLE GRID FOR PROLONGATION
		m_vLevData[lev]->copy_correction_from_smooth_patch(true);
		GMG_PROFILE_END();

#ifdef UG_PARALLEL
		UG_DLOG(LIB_DISC_MULTIGRID, 2, "GMG: BaseSolver done on level "<<lev<<".\n");
	}

//	CASE b): We gather the processes, solve on one proc and distribute again
	else
	{
	//	get whole grid correction
		vector_type& c = m_vLevData[lev]->c;

	//	gather the defect
		bool resume = gather_vertical(d);

	//	check, if this proc continues, else idle
		if(resume)
		{
			GMG_PROFILE_BEGIN(GMG_BaseSolver);
			UG_DLOG(LIB_DISC_MULTIGRID, 2, " GMG: Start BaseSolver on proc 1.\n");

		//	Reset correction
			c.set(0.0);

		//	compute coarse correction
			if(!m_pBaseSolver->apply(c, d))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Base solver on"
						" base level " << lev << " failed. "
						"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");

				return false;
			}

		//	update defect
			if(m_baseLev == m_topLev || m_bAdaptive)
				if(!m_vLevData[m_baseLev]->LevMat.apply_sub(d, c))
				{
					UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Updating defect "
							" on base level " << lev << ". "
							"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
					return false;
				}
			GMG_PROFILE_END();
			UG_DLOG(LIB_DISC_MULTIGRID, 2, " GMG Base solver done on 1 Proc.\n");
		}

	//	broadcast the correction
		broadcast_vertical(c);
		c.set_storage_type(PST_CONSISTENT);

	//	if baseLevel == surfaceLevel, we need also d
		if(m_baseLev == m_topLev || m_bAdaptive)
		{
			d.set_storage_type(PST_CONSISTENT);
			broadcast_vertical(d);
			d.change_storage_type(PST_ADDITIVE);
		}
	}
#endif

//	we're done for the solution of the base solver
	return true;

}

// performs a  multi grid cycle on the level
template <typename TApproximationSpace, typename TAlgebra>
bool AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
lmgc(size_t lev)
{
//	switch, if base level is reached. If so, call base Solver, else we will
//	perform smoothing, restrict the defect and call the lower level; then,
//	going up again in the level hierarchy the correction is interpolated and
//	used as coarse grid correction. Finally a post-smooth is performed.
	if(lev > m_baseLev)
	{
	//	check if dofs on that level. It not, skip level
		if(m_vLevData[lev]->num_indices() > 0){


			for(int i = 0; i < m_cycleType; ++i)
			{
			//	m_vLevData[lev]->c.set(0.0); // <<<< only for debug
			//	UG_LOG("Before presmooth:\n");	log_level_data(lev);
				if(!presmooth(lev))
					return false;

			//	store whether the restriction was resuming to the level below.
				bool restrictionPerformed = true;

			//	UG_LOG("Before restriction:\n");	log_level_data(lev);
				if(!restriction(lev, &restrictionPerformed))
					return false;

			//todo: It could make sense to call lmgc even if no restriction was
			//		performed. This could be required in a parallel environment
			//		where lmgc is applied to an adaptive redistributed grid.
			//		In that situation an empty level could be located between
			//		filled ones.
				if(restrictionPerformed){
					//UG_LOG("Before recursion:\n");	log_level_data(lev);
					if(!lmgc(lev-1))
					{
						UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Linear multi"
								" grid cycle on level " << lev-1 << " failed. "
								"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
						return false;
					}
				}

				//UG_LOG("Before prolongation:\n");	log_level_data(lev);
				if(!prolongation(lev, restrictionPerformed))
					return false;

				//UG_LOG("Before postsmooth:\n");	log_level_data(lev);
				if(!postsmooth(lev))
					return false;

				//UG_LOG("After postsmooth:\n");	log_level_data(lev);
			}
			return true;
		}
		else{
			for(int i = 0; i < m_cycleType; ++i){
				if(!lmgc(lev-1))
					return false;
			}
			return true;
		}
	}

//	if the base level has been reached, the coarse problem is solved exactly
	else if(lev == m_baseLev)
	{
		return base_solve(lev);
	}

//	this case should never happen.
	UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Level index below "
			" 'baseLevel' in lmgc. Aborting.\n");
	return false;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
{
// 	Cast Operator
	m_pSurfaceMat = dynamic_cast<matrix_type*>(&J);

//	Check that Operator type is correct
	if(m_pSurfaceMat == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Can not cast Operator to Matrix.\n");
		return false;
	}

//	check that grid given
	if(m_pApproxSpace->num_levels() == 0)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"No grid level in Approximation Space.\n");
		return false;
	}

//	get current toplevel
	m_topLev = m_pApproxSpace->num_levels() - 1;

//	Allocate memory for given top level
	if(!top_level_required(m_topLev))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init':"
				" Cannot allocate memory. Aborting.\n");
		return false;
	}

	if(m_pProjectionPrototype == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': "
				"Projection not set, although problem nonlinear.\n");
		return false;
	}

//	Create Projection
	GMG_PROFILE_BEGIN(GMG_InitProjection);
	if(!init_projection())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': "
				"Initialization of Projection failed.\n");
		return false;
	}
	GMG_PROFILE_END();

//	project
	GMG_PROFILE_BEGIN(GMG_ProjectSolutionFromSurface);
	if(!project_surface_to_level(level_solutions(), u))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': "
				"Projection of solution to level failed.\n");
		return false;
	}
	GMG_PROFILE_END();

// 	Project solution from surface grid to coarser grid levels
	GMG_PROFILE_BEGIN(GMG_ProjectSolutionDown);
	for(size_t lev = m_topLev; lev != m_baseLev; --lev)
	{
	//	skip void level
		if(m_vLevData[lev]->num_indices() == 0 ||
			m_vLevData[lev-1]->num_indices() == 0) continue;

		if(!m_vLevData[lev]->Projection->apply(m_vLevData[lev-1]->u, m_vLevData[lev]->u))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Cannot project "
					"solution to coarse grid function of level "<<lev-1<<".\n");
			return false;
		}
	}
	GMG_PROFILE_END();

//	init common
	if(!init_common(true))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Cannot init common part.\n");
		return false;
	}

//	assemble missing coarse grid matrix contribution (only in adaptive case)
	if(m_bAdaptive)
		if(!init_missing_coarse_grid_coupling(&u))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
					"Cannot init missing coarse grid coupling.\n");
			return false;
		}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init(ILinearOperator<vector_type, vector_type>& L)
{
// 	Cast Operator
	m_pSurfaceMat = dynamic_cast<matrix_type*>(&L);

//	Check that Operator type is correct
	if(m_pSurfaceMat == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Can not cast Operator to Matrix.\n");
		return false;
	}

//	check that grid given
	if(m_pApproxSpace->num_levels() == 0)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"No grid level in Approximation Space.\n");
		return false;
	}

//	get current toplevel
	m_topLev = m_pApproxSpace->num_levels() - 1;

//	Allocate memory for given top level
	GMG_PROFILE_BEGIN(GMG_CreateLevelStorage);
	if(!top_level_required(m_topLev))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init':"
				" Cannot allocate memory. Aborting.\n");
		return false;
	}
	GMG_PROFILE_END();

//	init common
	if(!init_common(false))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Cannot init common part.\n");
		return false;
	}

//	assemble missing coarse grid matrix contribution (only in adaptive case)
	GMG_PROFILE_BEGIN(GMG_AssMissingCoarseMat);
	if(m_bAdaptive)
		if(!init_missing_coarse_grid_coupling(NULL))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
					"Cannot init missing coarse grid coupling.\n");
			return false;
		}
	GMG_PROFILE_END();

//	we're done
	return true;
}



template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_common(bool nonlinear)
{
//	Perform some checks:
	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Discretization not set.\n");
		return false;
	}
	if(m_pApproxSpace == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Approximation Space not set.\n");
		return false;
	}
	if(m_pBaseSolver == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Base Solver not set.\n");
		return false;
	}
	if(m_pSmootherPrototype == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Smoother not set.\n");
		return false;
	}
	if(m_pProlongationPrototype == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Prolongation not set.\n");
		return false;
	}

	if(m_baseLev > m_topLev)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Base Level can not be greater than surface level.\n");
		return false;
	}

//	check, if grid is full-refined
//todo:	make sure that there are no vertical masters in topLevel. Otherwise
//		the grid can not be considered fully refined.
//todo: Even if there are vrtMasters and m_bFullRefined is false and the top
//		level matrix can't be copied, an injective SurfToTopLevMap might be useful...
	if(m_pApproxSpace->get_level_dof_distribution(m_topLev).num_indices() ==
		m_pApproxSpace->get_surface_dof_distribution().num_indices())
		m_bAdaptive = false;
	else
		m_bAdaptive = true;


//	init mapping from surface level to top level in case of full refinement
	if(!m_bAdaptive)
	{
		GMG_PROFILE_BEGIN(GMG_InitSurfToLevelMapping);
		if(!CreateSurfaceToToplevelMap(m_vSurfToTopMap,
									   m_pApproxSpace->get_surface_dof_distribution(),
									   m_pApproxSpace->get_level_dof_distribution(m_topLev)))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_common': "
					"Cannot init Mapping Surface -> TopLevel (full refinement case).\n");
			return false;
		}
		GMG_PROFILE_END();
	}

//	Assemble coarse grid operators
	GMG_PROFILE_BEGIN(GMG_AssembleLevelGridOperator);
	if(nonlinear){
		if(!init_non_linear_level_operator()){
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_common': "
					"Cannot init (nonlinear) Coarse Grid Operator.\n");
			return false;
		}
	}
	else {
		if(!init_linear_level_operator()){
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_common': "
					"Cannot init (linear) Coarse Grid Operator.\n");
			return false;
		}
	}
	GMG_PROFILE_END();

//	write computed level matrices for debug purpose
	for(size_t lev = m_baseLev; lev < m_vLevData.size(); ++lev)
		write_level_debug(m_vLevData[lev]->LevMat, "LevelMatrix", lev);

//	Init smoother for coarse grid operators
	GMG_PROFILE_BEGIN(GMG_InitSmoother);
	if(!init_smoother())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init_common': "
				"Cannot init Smoother.\n");
		return false;
	}
	GMG_PROFILE_END();

//	Init base solver
	GMG_PROFILE_BEGIN(GMG_InitBaseSolver);
	if(!init_base_solver())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init_common': "
				"Cannot init Base Solver.\n");
		return false;
	}
	GMG_PROFILE_END();

// 	Create Interpolation
	GMG_PROFILE_BEGIN(GMG_InitProlongation);
	if(!init_prolongation())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init_common': "
				"Cannot init Prolongation.\n");
		return false;
	}
	GMG_PROFILE_END();

//	we're done
	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_linear_level_operator()
{
// 	Create coarse level operators
	for(size_t lev = m_baseLev; lev < m_vLevData.size(); ++lev)
	{
	//	skip void level
		if(m_vLevData[lev]->num_indices() == 0) continue;

	//	in case of full refinement we simply copy the matrix (with correct numbering)
		if(lev == m_vLevData.size() - 1 && !m_bAdaptive)
		{
			GMG_PROFILE_BEGIN(GMG_CopySurfMat);
			matrix_type& levMat = m_vLevData[lev]->LevMat;
			matrix_type& surfMat = *m_pSurfaceMat;

			levMat.resize( surfMat.num_rows(), surfMat.num_cols());
			CopyMatrixByMapping(levMat, m_vSurfToTopMap, surfMat);

			GMG_PROFILE_END();
			continue;
		}

		GMG_PROFILE_BEGIN(GMG_AssLevelMat);
	//	if ghosts are present we have to assemble the matrix only on non-ghosts
	//	for the smoothing matrices
		if(m_vLevData[lev]->has_ghosts())
		{
		//	set this selector to the assembling, such that only those elements
		//	will be assembled and force grid to be considered as regular
			if(m_vLevData[lev]->has_ghosts()) m_pAss->set_selector(&m_vLevData[lev]->sel);
			else m_pAss->set_selector(NULL);
			m_pAss->force_regular_grid(true);

		//	init level operator
			if(!m_pAss->assemble_jacobian(m_vLevData[lev]->LevMat, m_vLevData[lev]->u, *m_vLevData[lev]->pLevDD))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator':"
						" Cannot init operator for level "<< lev << ".\n");
				return false;
			}

		//	remove force flag
			m_pAss->force_regular_grid(false);
			m_pAss->set_selector(NULL);

		//	copy the matrix into a new (smaller) one
			matrix_type& mat = m_vLevData[lev]->LevMat;
			matrix_type& smoothMat = m_vLevData[lev]->SmoothMat;

			const size_t numSmoothIndex = m_vLevData[lev]->num_smooth_indices();
			smoothMat.resize(numSmoothIndex, numSmoothIndex);
			CopyMatrixByMapping(smoothMat, m_vLevData[lev]->vMapFlag, mat);
		}

	//	if no ghosts are present we can simply use the whole grid. If the base
	//	solver is carried out in serial (gathering to some processes), we have
	//	to assemble the assemble the coarse grid matrix on the whole grid as
	//	well
		if(!m_vLevData[lev]->has_ghosts() ||
			(lev == m_baseLev && m_bBaseParallel == false))
		{
		//	init level operator
			m_pAss->force_regular_grid(true);
			if(!m_pAss->assemble_jacobian(m_vLevData[lev]->LevMat, m_vLevData[lev]->u, *m_vLevData[lev]->pLevDD))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator':"
						" Cannot init operator for level "<< lev << ".\n");
				return false;
			}
			m_pAss->force_regular_grid(false);
		}
	//	else we can forget about the whole-level matrix, since the needed
	//	smoothing matrix is stored in SmoothMat
		else
		{
			m_vLevData[lev]->LevMat.resize(0,0);
		}

		GMG_PROFILE_END();
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_non_linear_level_operator()
{
// 	Create coarse level operators
	for(size_t lev = m_baseLev; lev < m_vLevData.size(); ++lev)
	{
	//	skip void level
		if(m_vLevData[lev]->num_indices() == 0) continue;

	//	in case of full refinement we simply copy the matrix (with correct numbering)
		if(lev == m_vLevData.size() - 1 && !m_bAdaptive)
		{
			GMG_PROFILE_BEGIN(GMG_CopySurfMat);
			matrix_type& levMat = m_vLevData[lev]->LevMat;
			matrix_type& surfMat = *m_pSurfaceMat;

			levMat.resize( surfMat.num_rows(), surfMat.num_cols());
			CopyMatrixByMapping(levMat, m_vSurfToTopMap, surfMat);

			GMG_PROFILE_END();
			continue;
		}

		GMG_PROFILE_BEGIN(GMG_AssLevelMat);
	//	set this selector to the assembling, such that only those elements
	//	will be assembled and force grid to be considered as regular
		if(m_vLevData[lev]->has_ghosts()) m_pAss->set_selector(&m_vLevData[lev]->sel);
		else m_pAss->set_selector(NULL);
		m_pAss->force_regular_grid(true);

	//	init level operator
		if(!m_pAss->assemble_jacobian(m_vLevData[lev]->LevMat, m_vLevData[lev]->u, *m_vLevData[lev]->pLevDD))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator':"
					" Cannot init operator for level "<< lev << ".\n");
			return false;
		}

	//	remove force flag
		m_pAss->force_regular_grid(false);
		m_pAss->set_selector(NULL);
		GMG_PROFILE_END();

	//	if ghosts are present copy the matrix into a new (smaller) one
		if(m_vLevData[lev]->has_ghosts())
		{
			matrix_type& mat = m_vLevData[lev]->LevMat;
			matrix_type& smoothMat = m_vLevData[lev]->SmoothMat;

			const size_t numSmoothIndex = m_vLevData[lev]->num_smooth_indices();
			smoothMat.resize(numSmoothIndex, numSmoothIndex);
			CopyMatrixByMapping(smoothMat, m_vLevData[lev]->vMapFlag, mat);
		}
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_prolongation()
{
//	loop all levels
	for(size_t lev = m_baseLev+1; lev < m_vLevData.size(); ++lev)
	{
	//	skip void level
		if(m_vLevData[lev]->num_indices() == 0 ||
		   m_vLevData[lev-1]->num_indices() == 0) continue;

	//	set levels
		if(!m_vLevData[lev]->Prolongation->set_levels(lev-1, lev))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_prolongation':"
					" Cannot set level in interpolation matrices for level "
					<< lev << ", aborting.\n");
			return false;
		}

	//	init prolongation
		if(!m_vLevData[lev]->Prolongation->init())
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_prolongation':"
					" Cannot init interpolation operator for level "
					<< lev << ", aborting.\n");
			return false;
		}
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_projection()
{
//	loop all levels
	for(size_t lev = m_baseLev+1; lev < m_vLevData.size(); ++lev)
	{
	//	skip void level
		if(m_vLevData[lev]->num_indices() == 0 ||
		   m_vLevData[lev-1]->num_indices() == 0) continue;

	//	set levels
		if(!m_vLevData[lev]->Projection->set_levels(lev-1, lev))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_projection':"
					" Cannot set level in projection matrices for level "
					<< lev << ", aborting.\n");
			return false;
		}

	//	init projection
		if(!m_vLevData[lev]->Projection->init())
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_projection':"
					" Cannot init projection operator for level "
					<< lev << ", aborting.\n");
			return false;
		}
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_smoother()
{
// 	Init smoother
	for(size_t lev = m_baseLev; lev < m_vLevData.size(); ++lev)
	{
	//	skip void level
		if(m_vLevData[lev]->num_indices() == 0) continue;

	//	get smooth matrix and vector
		vector_type& u = m_vLevData[lev]->get_smooth_solution();
		MatrixOperator<vector_type, vector_type, matrix_type>& SmoothMat =
				m_vLevData[lev]->get_smooth_mat();

		if(!m_vLevData[lev]->Smoother->init(SmoothMat, u))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_smoother':"
					" Cannot init smoother for level "<< lev << ".\n");
			return false;
		}
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_base_solver()
{
//	skip void level
	if(m_vLevData[m_baseLev]->num_indices() == 0) return true;

#ifdef UG_PARALLEL
//	check, if a gathering base solver is required:
	if(!m_bBaseParallel)
	{
	//	check if gathering base solver possible: If some horizontal layouts are
	//	given, we know, that still the grid is distributed. But, if no
	//	vertical layouts are present in addition, we can not gather the vectors
	//	to on proc. Write a warning an switch to distributed coarse solver
		vector_type& d = m_vLevData[m_baseLev]->d;
		if((!d.get_master_layout().empty() || !d.get_slave_layout().empty()) &&
		   (d.get_vertical_slave_layout().empty() && d.get_vertical_master_layout().empty()))
		{
			UG_LOG("WARNING in 'AssembledMultiGridCycle::init_base_solver': "
					" Cannot init distributed base solver on level "<< m_baseLev << ":\n"
					" Base level distributed among processes and no possibility"
					" of gathering (vert. interfaces) present. But a gathering"
					" solving is required. Choose gmg:set_parallel_base_solver(true)"
					" to avoid this warning.\n");
			m_bBaseParallel = true;
		}
		else
		{
		//	we init the base solver with the whole grid matrix
			if(!m_pBaseSolver->init(m_vLevData[m_baseLev]->LevMat, m_vLevData[m_baseLev]->u))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::init_base_solver':"
						" Cannot init base solver on baselevel "<< m_baseLev << ".\n");
				return false;
			}
		}
	}

//	\todo: add a check if base solver can be run in parallel. This needs to
//		   introduce such a flag in the solver.
//	in Serial or in case of a distributed coarse grid solver, we can simply use
//	the smoothing matrices to set up the solver.
	if(m_bBaseParallel)
#endif
	{
	//	get smooth matrix and vector
		vector_type& u = m_vLevData[m_baseLev]->get_smooth_solution();
		MatrixOperator<vector_type, vector_type, matrix_type>& SmoothMat =
				m_vLevData[m_baseLev]->get_smooth_mat();

		if(!m_pBaseSolver->init(SmoothMat, u))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_base_solver':"
					" Cannot init base solver on baselevel "<< m_baseLev << ".\n");
			return false;
		}
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
project_level_to_surface(vector_type& surfVec,
                         std::vector<const vector_type*> vLevelVec)
{
//	level dof distributions
	const std::vector<const dof_distribution_type*>& vLevelDD =
								m_pApproxSpace->get_level_dof_distributions();

//	surface dof distribution
	const dof_distribution_type& surfDD =
								m_pApproxSpace->get_surface_dof_distribution();

//	surface view
	const SurfaceView* surfView = m_pApproxSpace->get_surface_view();

//	check that surface view exists
	if(surfView == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::project_level_to_surface':"
				" Surface View missing.\n");
		return false;
	}

//	Now we can project the surface vector to the levels
//	Note: even in case of full refinement this is necessary, since the ordering
//		  of DoFs may differ between surface grid and top level
	if(!m_bAdaptive){
		const vector_type& topVec = *(vLevelVec[m_topLev]);
		for(size_t surfIndex = 0; surfIndex < m_vSurfToTopMap.size(); ++surfIndex)
		{
		//	get corresponding level index
			const size_t levIndex = m_vSurfToTopMap[surfIndex];

		//	write value
			surfVec[surfIndex] = topVec[levIndex];

#ifdef UG_PARALLEL
		//	copy storage type into all vectors
			surfVec.copy_storage_type(topVec);
#endif
		}
	}
	else {
		if(!ProjectLevelToSurface(surfVec, surfDD, *surfView,
								  vLevelVec, vLevelDD, m_baseLev))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::project_level_to_surface': "
					"Projection of function from level to surface failed.\n");
			return false;
		}
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
project_surface_to_level(std::vector<vector_type*> vLevelVec,
                         const vector_type& surfVec)
{
//	level dof distributions
	const std::vector<const dof_distribution_type*>& vLevelDD =
								m_pApproxSpace->get_level_dof_distributions();

//	surface dof distribution
	const dof_distribution_type& surfDD =
								m_pApproxSpace->get_surface_dof_distribution();

//	surface view
	const SurfaceView* surfView = m_pApproxSpace->get_surface_view();

//	check that surface view exists
	if(surfView == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::project_surface_to_level':"
				" Surface View missing.\n");
		return false;
	}

//	reset vectors
//	\todo: Is this really necessary ?
	for(size_t lev = 0; lev < vLevelVec.size(); ++lev)
		if(vLevelVec[lev] != NULL)
			vLevelVec[lev]->set(0.0);

//	Now we can project the surface vector to the levels
//	Note: even in case of full refinement this is necessary, since the ordering
//		  of DoFs may differ between surface grid and top level
	if(!m_bAdaptive){
		vector_type& topVec = *(vLevelVec[m_topLev]);
		for(size_t surfIndex = 0; surfIndex < m_vSurfToTopMap.size(); ++surfIndex)
		{
		//	get corresponding level index
			const size_t levIndex = m_vSurfToTopMap[surfIndex];

		//	write value
			topVec[levIndex] = surfVec[surfIndex];

#ifdef UG_PARALLEL
		//	copy storage type into all vectors
			for(size_t lev = 0; lev < vLevelVec.size(); ++lev)
				if(vLevelVec[lev] != NULL)
					vLevelVec[lev]->copy_storage_type(surfVec);
#endif

		}
	}
	else{
		if(!ProjectSurfaceToLevel(vLevelVec, vLevelDD,
								  surfVec, surfDD, *surfView))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::project_surface_to_level': "
					"Projection of function from surface to level failed.\n");
			return false;
		}
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
write_level_debug(const vector_type& vec, const char* filename, size_t lev)
{
//	if no debug writer set, we're done
	if(m_pDebugWriter == NULL) return true;

//	typedef function type
	typedef typename approximation_space_type::function_type function_type;

//	create level function
	function_type* dbgFunc = m_pApproxSpace->create_level_function(lev);

//	cast dbg writer
	GridFunctionDebugWriter<function_type>* dbgWriter =
			dynamic_cast<GridFunctionDebugWriter<function_type>*>(m_pDebugWriter);

//	set grid function
	if(dbgWriter != NULL)
		dbgWriter->set_reference_grid_function(*dbgFunc);
	else
	{
		delete dbgFunc;
		UG_LOG("Cannot write debug on level "<< lev<<".\n");
		return false;
	}

//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_lev%03d_iter%03d", (int)lev, m_dbgIterCnt);
	name.append(ext);

//	write
	bool bRet = m_pDebugWriter->write_vector(vec, name.c_str());

//	remove dbgFunc
	delete dbgFunc;

	return bRet;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
write_level_debug(const matrix_type& mat, const char* filename, size_t lev)
{
//	if no debug writer set, we're done
	if(m_pDebugWriter == NULL) return true;

//	typedef function type
	typedef typename approximation_space_type::function_type function_type;

//	create level function
	function_type* dbgFunc = m_pApproxSpace->create_level_function(lev);

//	cast dbg writer
	GridFunctionDebugWriter<function_type>* dbgWriter =
			dynamic_cast<GridFunctionDebugWriter<function_type>*>(m_pDebugWriter);

//	set grid function
	if(dbgWriter != NULL)
		dbgWriter->set_reference_grid_function(*dbgFunc);
	else
	{
		delete dbgFunc;
		UG_LOG("Cannot write debug on level "<< lev<<".\n");
		return false;
	}

//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_lev%03d_iter%03d", (int)lev, m_dbgIterCnt);
	name.append(ext);

//	write
	bool bRet = m_pDebugWriter->write_matrix(mat, name.c_str());

//	remove dbgFunc
	delete dbgFunc;

	return bRet;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
write_surface_debug(const vector_type& vec, const char* filename)
{
//	if no debug writer set, we're done
	if(m_pDebugWriter == NULL) return true;

//	typedef function type
	typedef typename approximation_space_type::function_type function_type;

//	create level function
	function_type* dbgFunc = m_pApproxSpace->create_surface_function();

//	cast dbg writer
	GridFunctionDebugWriter<function_type>* dbgWriter =
			dynamic_cast<GridFunctionDebugWriter<function_type>*>(m_pDebugWriter);

//	set grid function
	if(dbgWriter != NULL)
		dbgWriter->set_reference_grid_function(*dbgFunc);
	else
	{
		delete dbgFunc;
		UG_LOG("Cannot write debug on surface.\n");
		return false;
	}

//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_surf_iter%03d", m_dbgIterCnt);
	name.append(ext);

//	write
	bool bRet = m_pDebugWriter->write_vector(vec, name.c_str());

//	remove dbgFunc
	delete dbgFunc;

	return bRet;
}


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
write_surface_debug(const matrix_type& mat, const char* filename)
{
//	if no debug writer set, we're done
	if(m_pDebugWriter == NULL) return true;

//	typedef function type
	typedef typename approximation_space_type::function_type function_type;

//	create level function
	function_type* dbgFunc = m_pApproxSpace->create_surface_function();

//	cast dbg writer
	GridFunctionDebugWriter<function_type>* dbgWriter =
			dynamic_cast<GridFunctionDebugWriter<function_type>*>(m_pDebugWriter);

//	set grid function
	if(dbgWriter != NULL)
		dbgWriter->set_reference_grid_function(*dbgFunc);
	else
	{
		delete dbgFunc;
		UG_LOG("Cannot write debug on surface.\n");
		return false;
	}

//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_surf_iter%03d", m_dbgIterCnt);
	name.append(ext);

//	write
	bool bRet = m_pDebugWriter->write_matrix(mat, name.c_str());

//	remove dbgFunc
	delete dbgFunc;

	return bRet;
}

template <typename TApproximationSpace, typename TAlgebra>
void
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
log_level_data(size_t lvl)
{
	std::string prefix;
	if(lvl < m_vLevData.size())
		prefix.assign(2 + 2 * (m_vLevData.size() - lvl), ' ');

	LevData& ld = *m_vLevData[lvl];
	UG_LOG(prefix << "local d norm: " << sqrt(VecProd(ld.d, ld.d)) << std::endl);
	UG_LOG(prefix << "local c norm: " << sqrt(VecProd(ld.c, ld.c)) << std::endl);
	UG_LOG(prefix << "local smooth_d norm: " << sqrt(VecProd(ld.get_smooth_defect(), ld.get_smooth_defect())) << std::endl);
	UG_LOG(prefix << "local smooth_c norm: " << sqrt(VecProd(ld.get_smooth_correction(), ld.get_smooth_correction())) << std::endl);

	#ifdef UG_PARALLEL
		uint oldStorageMask = ld.d.get_storage_mask();
		number norm = ld.d.two_norm();
		UG_LOG(prefix << "parallel d norm: " << norm << "\n");
		if(oldStorageMask & PST_ADDITIVE)
			ld.d.change_storage_type(PST_ADDITIVE);
		else if(oldStorageMask & PST_CONSISTENT)
			ld.d.change_storage_type(PST_CONSISTENT);

		oldStorageMask = ld.c.get_storage_mask();
		norm = ld.c.two_norm();
		UG_LOG(prefix << "parallel c norm: " << norm << "\n");
		if(oldStorageMask & PST_ADDITIVE)
			ld.c.change_storage_type(PST_ADDITIVE);
		else if(oldStorageMask & PST_CONSISTENT)
			ld.c.change_storage_type(PST_CONSISTENT);

		oldStorageMask = ld.get_smooth_defect().get_storage_mask();
		norm = ld.get_smooth_defect().two_norm();
		UG_LOG(prefix << "parallel smooth defect norm: " << norm << "\n");
		if(oldStorageMask & PST_ADDITIVE)
			ld.get_smooth_defect().change_storage_type(PST_ADDITIVE);
		else if(oldStorageMask & PST_CONSISTENT)
			ld.get_smooth_defect().change_storage_type(PST_CONSISTENT);

		oldStorageMask = ld.get_smooth_correction().get_storage_mask();
		norm = ld.get_smooth_correction().two_norm();
		UG_LOG(prefix << "parallel smooth correction norm: " << norm << "\n");
		if(oldStorageMask & PST_ADDITIVE)
			ld.get_smooth_correction().change_storage_type(PST_ADDITIVE);
		else if(oldStorageMask & PST_CONSISTENT)
			ld.get_smooth_correction().change_storage_type(PST_CONSISTENT);
	#endif
}

template <typename TApproximationSpace, typename TAlgebra>
typename AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::base_type*
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
clone()
{
	AssembledMultiGridCycle<TApproximationSpace, TAlgebra>* clone =
		new AssembledMultiGridCycle<TApproximationSpace, TAlgebra>(*m_pApproxSpace);

	clone->set_base_level(m_baseLev);
	clone->set_base_solver(*m_pBaseSolver);
	clone->set_cycle_type(m_cycleType);
	clone->set_debug(m_pDebugWriter);
	clone->set_discretization(*m_pAss);
	clone->set_num_postsmooth(m_numPostSmooth);
	clone->set_num_presmooth(m_numPreSmooth);
	clone->set_projection_operator(*m_pProjectionPrototype);
	clone->set_prolongation_operator(*m_pProlongationPrototype);
	clone->set_smoother(*m_pSmootherPrototype);

	return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
}



template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_missing_coarse_grid_coupling(const vector_type* u)
{
//	clear matrices
	for(size_t lev = 0; lev < m_vLevData.size(); ++lev)
		m_vLevData[lev]->CoarseGridContribution.resize(0,0);

//	if the grid is fully refined, nothing to do
	if(!m_bAdaptive) return true;

//	get the surface view
	const SurfaceView& surfView = *m_pApproxSpace->get_surface_view();

//	check that surface view exists
	if(&surfView == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_missing_coarse_grid_coupling':"
				" Surface View missing.\n");
		return false;
	}

//	create storage for matrices on the grid levels
	for(size_t lev = 0; lev < m_vLevData.size(); ++lev)
	{
	//	get dof distributions on levels
		const dof_distribution_type& dofDistr
							= m_pApproxSpace->get_level_dof_distribution(lev);

	//	resize the matrix
		m_vLevData[lev]->CoarseGridContribution.resize(dofDistr.num_indices(),
		                                               dofDistr.num_indices());
	}

///////////////////////////////////////
//	create surface -> level mappings
///////////////////////////////////////

//	level dof distributions
	const std::vector<const dof_distribution_type*>& vLevelDD =
								m_pApproxSpace->get_level_dof_distributions();

//	surface dof distribution
	const dof_distribution_type& surfDD =
								m_pApproxSpace->get_surface_dof_distribution();

//	check that surface dof distribution exists
	if(&surfDD == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_missing_coarse_grid_coupling':"
				"Surface DoF Distribution missing.\n");
		return false;
	}

//	create mappings
	std::vector<std::vector<int> > vSurfLevelMapping;
	if(!CreateSurfaceToLevelMapping(vSurfLevelMapping, vLevelDD, surfDD, surfView))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_missing_coarse_grid_coupling':"
				" Cannot build surface index to level index mappings.\n");
		return false;
	}

///////////////////////////////////////
//	assemble contribution for each level and project
///////////////////////////////////////

//	loop all levels to compute the missing contribution
//	\todo: this is implemented very resource consuming, re-think arrangement
	Selector sel(m_pApproxSpace->get_domain().get_grid());
	for(size_t lev = 0; lev < m_vLevData.size(); ++lev)
	{
	//	select all elements, that have a shadow as a subelement, but are not itself
	//	a shadow
		sel.clear();
		SelectNonShadowsAdjacentToShadowsOnLevel(sel, surfView, lev);

	//	now set this selector to the assembling, such that only those elements
	//	will be assembled
		m_pAss->set_selector(&sel);

	//	create a surface matrix
		matrix_type surfMat;

	//	assemble the surface jacobian only for selected elements
		if(u)
			m_pAss->assemble_jacobian(surfMat, *u, m_pApproxSpace->get_surface_dof_distribution());
		else
		{
		//	\todo: not use tmp vector here
			vector_type tmpVec; tmpVec.resize(m_pApproxSpace->get_surface_dof_distribution().num_indices());
			m_pAss->assemble_jacobian(surfMat, tmpVec, m_pApproxSpace->get_surface_dof_distribution());
		}

	//	write matrix for debug purpose
		std::stringstream ss; ss << "MissingSurfMat_" << lev;
		write_surface_debug(surfMat, ss.str().c_str());

	//	remove the selector from the assembling procedure
		m_pAss->set_selector(NULL);

	//	project
		if(!CopyMatrixSurfaceToLevel(m_vLevData[lev]->CoarseGridContribution,
		                             vSurfLevelMapping[lev],
		                             surfMat))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_missing_coarse_grid_coupling': "
					"Projection of matrix from surface to level failed.\n");
			return false;
		}
	}

//	write matrix for debug purpose
	for(size_t lev = 0; lev < m_vLevData.size(); ++lev)
		write_level_debug(m_vLevData[lev]->CoarseGridContribution, "MissingLevelMat", lev);

/////////////
// end project
/////////////

//	we're done
	return true;
}


#ifdef UG_PARALLEL
template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
gather_vertical(vector_type& d)
{
//	start with resume as true, i.e. process will continue computation
//	on the coarser level
	bool resume = true;

//	send vertical-slaves -> vertical-masters
//	one proc may not have both, a vertical-slave- and vertical-master-layout.
	GMG_PROFILE_BEGIN(GMG_GatherVerticalVector);
	ComPol_VecAdd<vector_type> cpVecAdd(&d);
	if(!d.get_vertical_slave_layout().empty()){
	//	do not resume if vertical slaves are present
		resume = false;
		UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
		  " Going down: SENDS vert. dofs.\n");

	//	schedule Sending of DoFs of vertical slaves
		m_Com.send_data(d.get_vertical_slave_layout(), cpVecAdd);
	}
	else if(!d.get_vertical_master_layout().empty()){
		UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
		 " Going down:  WAITS FOR RECIEVE of vert. dofs.\n");

	//	schedule Receive of DoFs on vertical masters
		m_Com.receive_data(d.get_vertical_master_layout(), cpVecAdd);
	}

//	perform communication
	m_Com.communicate();
	GMG_PROFILE_END();

	return resume;
}

template <typename TApproximationSpace, typename TAlgebra>
void
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
broadcast_vertical(vector_type& t)
{
//	send vertical-masters -> vertical-slaves
//	one proc may not have both, a vertical-slave- and vertical-master-layout.
	GMG_PROFILE_BEGIN(GMG_BroadcastVerticalVector);
	ComPol_VecCopy<vector_type> cpVecCopy(&t);
	if(!t.get_vertical_slave_layout().empty())
	{
		UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
		 " Going up: WAITS FOR RECIEVE of vert. dofs.\n");

	//	schedule slaves to receive correction
		m_Com.receive_data(t.get_vertical_slave_layout(), cpVecCopy);
	}
	else if(!t.get_vertical_master_layout().empty())
	{
		UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
		 " Going up: SENDS vert. dofs.\n");

	//	schedule masters to send correction
		m_Com.send_data(t.get_vertical_master_layout(), cpVecCopy);
	}

//	communicate
	m_Com.communicate();
	GMG_PROFILE_END();
}
#endif

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
top_level_required(size_t topLevel)
{
//	allocated level if needed
	while(num_levels() <= topLevel)
	{
		m_vLevData.push_back(new LevData);
	}

//	free level if needed
	while(num_levels() > topLevel+1)
	{
		delete m_vLevData.back();
		m_vLevData.pop_back();
	}

//	reinit all levels
	for(size_t lev = m_baseLev; lev < m_vLevData.size(); ++lev)
	{
		m_vLevData[lev]->update(lev,
		                       *m_pApproxSpace,
		                       *m_pAss,
		                       *m_pSmootherPrototype,
		                       *m_pProjectionPrototype,
		                       *m_pProlongationPrototype);
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
void
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
LevData::
update(size_t lev,
       approximation_space_type& approxSpace,
       assemble_type& ass,
       smoother_type& smoother,
       projection_operator_type& projection,
       prolongation_operator_type& prolongation)
{
//	get dof distribution
	pLevDD = &approxSpace.get_level_dof_distribution(lev);

//	resize vectors for operations on whole grid level
	const size_t numIndex = pLevDD->num_indices();
	u.resize(numIndex);
	c.resize(numIndex);
	d.resize(numIndex);
	t.resize(numIndex);

//	prepare level operator
#ifdef UG_PARALLEL
	CopyLayoutsAndCommunicatorIntoMatrix(LevMat, *pLevDD);
#endif

//	clone operators
	if(!Smoother) Smoother = smoother.clone();
	if(!Projection) Projection = projection.clone();
	if(!Prolongation) Prolongation = prolongation.clone();

//	IN PARALLEL:
//	In the parallel case one may have vertical slaves/masters. Those are needed
//	when the grid hierarchy ends at a certain point on one process but is still
//	continued on another. In this case the values are transfered from the
//	valishing processor and copied to another. In addition it may happen, that
//	at a certain level a process only starts to have a grid elements and the
//	values are copied from another process to this process.
//	While these parts of the grid is needed for the transfer of correction or
//	defect between the grid levels, on those elements no smoothing is performed.
//	Therefore, if vertical masters are present only on the part of the grid
//	without vertical masters is smoothed. To account for this, in addition
//	smoothing vectors and matrices of smaller size are created and assembled.
//	Please note that smoothing is performed on vertical slaves.
#ifdef UG_PARALLEL
//	copy the layouts into the level vectors
	CopyLayoutsAndCommunicatorIntoVector(u, *pLevDD);
	CopyLayoutsAndCommunicatorIntoVector(c, *pLevDD);
	CopyLayoutsAndCommunicatorIntoVector(d, *pLevDD);
	CopyLayoutsAndCommunicatorIntoVector(t, *pLevDD);

//	if no vertical masters, there can be no ghost and we're ready. By ghosts
//	we denote vertical masters, that are not horizontal master/slave
	if(pLevDD->get_vertical_master_layout().empty())
	{
		m_numSmoothIndices = numIndex;
		return;
	}

//	If ghosts are present we create the infrastructure for this. This includes
//	the creation of Smoother on a smaller patch and required vectors/matrices
//	Also a mapping between the index set on the whole grid and the index set
//	on the smoothing patch must be created

//	** 1. **: We create the mapping between the index sets
//	create a vector of size of the whole grid with 1 everywhere
	vMapFlag.clear(); vMapFlag.resize(numIndex, 1);

//	set the vector to -1 where vertical masters are present, the set all
//	indices back to 1 where the index is also a horizontal master/slave
	SetLayoutValues(&vMapFlag, pLevDD->get_vertical_master_layout(), -1);
	SetLayoutValues(&vMapFlag, pLevDD->get_master_layout(), 1);
	SetLayoutValues(&vMapFlag, pLevDD->get_slave_layout(), 1);

//	now we create the two mapping:
//	vMapFlag: mapping (whole grid index -> patch index): the non-ghost indices
//	are mapped to a patch index, while the ghosts are flagged by a -1 index
//	vMap: mapping (patch index -> whole grid index): For each patch index the
//	corresponding whole grid index is stored
	vMap.clear();
	for(size_t j = 0; j < vMapFlag.size(); ++j)
	{
	//	if the index is still negative (i.e. ghost, leave index at -1)
		if(vMapFlag[j] -1) continue;

	//	if the index is a non-ghost set the new index
		vMapFlag[j] = vMap.size();

	//	since in the patch we store the mapping index
		vMap.push_back(j);
	}

//	now we know the size of the smoothing patch index set and resize help vectors
//	by the preceeding 's' the relation to the smoothing is indicated
	const size_t numSmoothIndex = vMap.size();
	m_numSmoothIndices = numSmoothIndex;
	su.resize(numSmoothIndex);
	sc.resize(numSmoothIndex);
	sd.resize(numSmoothIndex);
	st.resize(numSmoothIndex);

//	** 2. **: We have to create new layouts for the smoothers since on the
//	smoothing patch the indices are labeled differently.
//	copy layouts
	SmoothMasterLayout = pLevDD->get_master_layout();
	SmoothSlaveLayout = pLevDD->get_slave_layout();

//	Replace indices in the layout with the smaller (smoothing patch) indices
	ReplaceIndicesInLayout(SmoothMasterLayout, vMapFlag);
	ReplaceIndicesInLayout(SmoothSlaveLayout, vMapFlag);

//	replace old layouts by new modified ones
	sc.set_layouts(SmoothMasterLayout, SmoothSlaveLayout);
	su.set_layouts(SmoothMasterLayout, SmoothSlaveLayout);
	sd.set_layouts(SmoothMasterLayout, SmoothSlaveLayout);
	st.set_layouts(SmoothMasterLayout, SmoothSlaveLayout);
	sc.set_communicator(pLevDD->get_communicator());
	su.set_communicator(pLevDD->get_communicator());
	sd.set_communicator(pLevDD->get_communicator());
	st.set_communicator(pLevDD->get_communicator());
	sc.set_process_communicator(pLevDD->get_process_communicator());
	su.set_process_communicator(pLevDD->get_process_communicator());
	sd.set_process_communicator(pLevDD->get_process_communicator());
	st.set_process_communicator(pLevDD->get_process_communicator());

//	set the layouts in the smooth matrix
	SmoothMat.set_master_layout(SmoothMasterLayout);
	SmoothMat.set_slave_layout(SmoothSlaveLayout);
	SmoothMat.set_communicator(pLevDD->get_communicator());
	SmoothMat.set_process_communicator(pLevDD->get_process_communicator());

//	** 3. **: Since smoothing is only performed on non-ghost elements, the
//	corresoding operator must be assembled only on those elements. So we
//	use a selector to mark all non-ghosts and assemble on those later
//	get distributed Grid manager
	DistributedGridManager* pDstGrdMgr
		= approxSpace.get_domain().get_distributed_grid_manager();

//	select all ghost geometric objects
	sel.clear();
	sel.assign_grid(approxSpace.get_domain().get_grid());
	for(int si = 0; si < pLevDD->num_subsets(); ++si)
	{
		SelectNonGhosts<VertexBase>(sel, *pDstGrdMgr,
								 pLevDD->template begin<VertexBase>(si),
								 pLevDD->template end<VertexBase>(si));
		SelectNonGhosts<EdgeBase>(sel, *pDstGrdMgr,
								 pLevDD->template begin<EdgeBase>(si),
								 pLevDD->template end<EdgeBase>(si));
		SelectNonGhosts<Face>(sel, *pDstGrdMgr,
								 pLevDD->template begin<Face>(si),
								 pLevDD->template end<Face>(si));
		SelectNonGhosts<Volume>(sel, *pDstGrdMgr,
								 pLevDD->template begin<Volume>(si),
								 pLevDD->template end<Volume>(si));
	}
	
#else //PARALLEL

//	We have to smooth on the entire level
	m_numSmoothIndices = numIndex;
	
#endif
}

template <typename TApproximationSpace, typename TAlgebra>
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
LevData::~LevData()
{
//	free operators if allocated
	if(Smoother) delete Smoother;
	if(Projection) delete Projection;
	if(Prolongation) delete Prolongation;
}





} // namespace ug


#endif /* __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__ */
