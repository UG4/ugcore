/*
 * mg_solver_impl.h
 *
 *  Created on: 04.01.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__

#include "common/profiler/profiler.h"
#include "mg_solver_util.h"
#include "lib_discretization/function_spaces/grid_function_util.h"
#include "lib_discretization/dof_manager/dof_manager_util.h"

#include "mg_solver.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

#include <iostream>
#include <sstream>


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
	vector_type d_copy; d_copy.resize(d.size());

//	copy defect
	d_copy = d;

//	work on copy
	return apply_update_defect(c, d_copy);
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

// 	Check if surface level has been choosen correctly
	if(m_topLev >= mg.num_levels())
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
	GMG_PROFILE_BEGIN(GMG_lmgc);
	if(!lmgc(*m_vLevData[m_topLev].c, *m_vLevData[m_topLev].d, m_topLev))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Cannot perform multi grid cycle.\n");
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
smooth(vector_type& c, vector_type& d, vector_type& t,
       IMatrixOperator<vector_type, vector_type, matrix_type>& A,
       smoother_type& S,
       size_t lev, int nu)
{
// 	smooth nu times
	for(int i = 0; i < nu; ++i)
	{
	//	switch if adaptive case must be handled
		if(m_bFullRefined)
		{
		// 	Compute Correction of one smoothing step:
		//	a)  Compute t = B*d with some iterator B
		//	b) 	Update Defect d := d - A * t
			if(!S.apply_update_defect(t, d))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Smoothing step "
						<< i+1 << " on level " << lev << " failed.\n");
				return false;
			}

		// 	add correction of smoothing step to level correction
		//	(Note: we do not work on c directly here, since we update the defect
		//	       after every smoothing step. The summed up correction corresponds
		//		   to the total correction of the whole smoothing.)
			c += t;
		}
		else
	//	This is the adaptive case. Here, we must ensure, that the added correction
	//	is zero on the adaptively refined patch boundary of this level
		{
		// 	Compute Correction of one smoothing step, but do not update defect
		//	a)  Compute t = B*d with some iterator B
			if(!S.apply(t, d))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Smoothing step "
						<< i+1 << " on level " << lev << " failed.\n");
				return false;
			}

		//	First we reset the correction to zero on the patch boundary.
			if(!SetZeroOnShadowingVertex(t, *m_pApproxSpace, lev))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Could not "
						" reset the values to zero on patch boundary for correction "
						<< i+1 << " on level " << lev << ".\n");
				return false;
			}

		//	now, we can update the defect with this correction ...
			if(!A.apply_sub(d, t))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Could not "
						" update defect for patch correction in smooth step "
						<< i+1 << " on level " << lev << ".\n");
				return false;
			}

		//	... and add the correction to to overall correction
			c += t;
		}
	}

//	we're done
	return true;
}

// performs a  multi grid cycle on the level
template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
lmgc(vector_type& c, vector_type& d, size_t lev)
{
// 	reset correction to zero on this level
  	c.set(0.0);

//	vector defined on whole grid (including ghosts) on this level
//	are given by c, d, t
	vector_type& t = *m_vLevData[lev].t;

//	Lets get a reference to the coarser level correction, defect, help vector
	vector_type& cc = *m_vLevData[lev-1].c;
	vector_type& dc = *m_vLevData[lev-1].d;
	vector_type& tc = *m_vLevData[lev-1].t;

//	if vertical masters present, we have to restrict the smoothing
	vector_type* sc = &c;
	vector_type* sd = &d;
	vector_type* st = m_vLevData[lev].t;

//	switch, if base level is reached. If so, call base Solver, else we will
//	perform smoothing, restrict the defect and call the lower level; then,
//	going up again in the level hierarchy the correction is interpolated and
//	used as coarse grid correction. Finally a post-smooth is performed.
	if(lev > m_baseLev)
	{
	//	get smomother on this level and corresponding operator
		smoother_type* S = m_vLevData[lev].Smoother;
		IMatrixOperator<vector_type, vector_type, matrix_type>* A = m_vLevData[lev].A;

	//	LIFTING c,d TO SMOOTHING AREA
		#ifdef UG_PARALLEL
		if(m_vLevData[lev].SmoothMat != NULL)
		{
			sc = m_vLevData[lev].sc;
			sd = m_vLevData[lev].sd;
			st = m_vLevData[lev].st;
			for(size_t i = 0; i < m_vLevData[lev].vMap.size(); ++i)
			{
				const size_t oldIndex =  m_vLevData[lev].vMap[i];
				(*sc)[i] = c[oldIndex];
				(*sd)[i] = d[oldIndex];
			}
			sc->copy_storage_type(c);
			sd->copy_storage_type(d);

			A = m_vLevData[lev].SmoothMat;
		}
		#endif

	// 	PRE-SMOOTHING
	//	We start the multi grid cycle on this level by smoothing the defect. This
	//	means that we compute a correction c, such that the defect is "smoother".
		GMG_PROFILE_BEGIN(GMG_PreSmooth);
		if(!smooth(*sc, *sd, *st, *A, *S, lev, m_numPreSmooth))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Pre-Smoothing on "
					"level " << lev << " failed. "
				"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
			return false;
		}
		GMG_PROFILE_END();

	//	PROJECT DEFECT BACK TO WHOLE GRID FOR RESTRICTION
		#ifdef UG_PARALLEL
		if(m_vLevData[lev].SmoothMat != NULL)
		{
			d.set(0.0);
			for(size_t i = 0; i < m_vLevData[lev].vMap.size(); ++i)
			{
				const size_t oldIndex =  m_vLevData[lev].vMap[i];
				d[oldIndex] = (*sd)[i];
			}
			d.copy_storage_type(*sd);
		}
		#endif

	//	PARALLEL CASE: Send vertical slave values to master and check resuming
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
		bool resume = gather_vertical(d);

	//	only continue if levels left
		if(resume) {
		#endif

	// 	RESTRICT DEFECT
	//	Now we can restrict the defect from the fine level to the coarser level.
	//	This is done using the transposed prolongation.
		GMG_PROFILE_BEGIN(GMG_RestrictDefect);
		if(!m_vLevData[lev].Prolongation->apply_transposed(dc, d))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Restriction of "
					"Defect from level "<<lev<<" to "<<lev-1<<" failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
			return false;
		}
		GMG_PROFILE_END();

	// 	COMPUTE COARSE GRID CORRECTION
	//	Now, we have to compute the coarse grid correction, i.e. the correction
	//	on the coarser level. This is done iteratively by calling lmgc again for
	//	the coarser level.
		for(int i = 0; i < m_cycleType; ++i)
		{
			if(!lmgc(cc, dc, lev-1))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Linear multi"
						" grid cycle on level " << lev-1 << " failed. "
						"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
				return false;
			}
		}

	//	INTERPOLATE CORRECTION
	//	now we can interpolate the coarse grid correction from the coarse level
	//	to the fine level
		GMG_PROFILE_BEGIN(GMG_InterpolateCorr);
		if(!m_vLevData[lev].Prolongation->apply(t, cc))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Prolongation from"
					" level " << lev-1 << " to " << lev << " failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");

			return false;
		}
		GMG_PROFILE_END();

	//	PARALLEL CASE: Receive values of correction for vertical slaves
	//	If there are vertical slaves/masters on the coarser level, we now copy
	//	the correction values from the master DoFs to the slave	DoFs.
		#ifdef UG_PARALLEL
		}
		broadcast_vertical(t);
		#endif

	//	PROJECT COARSE GRID CORRECTION ONTO SMOOTH AREA
		#ifdef UG_PARALLEL
		if(m_vLevData[lev].SmoothMat != NULL)
		{
			for(size_t i = 0; i < m_vLevData[lev].vMap.size(); ++i)
			{
				const size_t oldIndex =  m_vLevData[lev].vMap[i];
				(*st)[i] = t[oldIndex];
			}
			st->copy_storage_type(t);
		}
		#endif

	// 	ADD COARSE GRID CORRECTION
		GMG_PROFILE_BEGIN(GMG_AddCoarseGridCorr);
		(*sc) += (*st);
		GMG_PROFILE_END(); // GMG_AddCoarseGridCorr

	//	UPDATE DEFECT FOR COARSE GRID CORRECTION
	//	the correction has changed c := c + t. Thus, we also have to update
	//	the defect d := d - A*t
		GMG_PROFILE_BEGIN(GMG_UpdateDefectForCGCorr);
		if(!A->apply_sub(*sd, *st))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Updating of defect"
					" on level " << lev << " failed.\n");
			GMG_PROFILE_END(); // GMG_UpdateDefectForCGCorr
			return false;
		}
		GMG_PROFILE_END(); // GMG_UpdateDefectForCGCorr

	//	ADAPTIVE CASE
		if(!m_bFullRefined)
		{
		//	in the adaptive case there is a small part of the coarse coupling that
		//	has not been used to update the defect. In order to ensure, that the
		//	defect on this level still corresponds to the updated defect, we need
		//	to add if here. This is done in three steps:
		//	a) Compute the coarse update of the defect induced by missing coupling
			tc.set(0.0);
			if(!m_vLevData[lev-1].CoarseGridContribution->apply(tc, cc))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Could not compute"
						" missing update defect contribution on level "<<lev-1<<".\n");
				return false;
			}

			tc *= -1.0;

		//	b) interpolate the coarse defect up
			if(!AddProjectionOfVertexShadows(d, tc, *m_pApproxSpace, lev-1, lev))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Could not add"
						" missing update defect contribution to level "<<lev<<".\n");
				return false;
			}
		}

	// 	POST-SMOOTHING
	//	We smooth the updated defect again. This means that we compute a
	//	correction c, such that the defect is "smoother".
		GMG_PROFILE_BEGIN(GMG_PostSmooth);
		if(!smooth(*sc, *sd, *st, *A, *S, lev, m_numPostSmooth))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Post-Smoothing on"
					" level " << lev << " failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
			return false;
		}
		GMG_PROFILE_END();

	//	PROJECT DEFECT, CORRECTION BACK TO WHOLE GRID FOR RESTRICTION
		#ifdef UG_PARALLEL
		if(m_vLevData[lev].SmoothMat != NULL)
		{
			for(size_t i = 0; i < m_vLevData[lev].vMap.size(); ++i)
			{
				const size_t oldIndex =  m_vLevData[lev].vMap[i];
				d[oldIndex] = (*sd)[i];
				c[oldIndex] = (*sc)[i];
			}
			d.copy_storage_type(*sd);
			c.copy_storage_type(*sc);
		}
		#endif

	//	we are done on this level
		return true;
	}

//	if the base level has been reached, the coarse problem is solved exactly
	else if(lev == m_baseLev)
	{
	//	SOLVE BASE PROBLEM
	//	Here we distinguish two possibilities:
	//	a) The coarse grid problem is solved in parallel, using a parallel solver
	//	b) First all vectors are gathered to one process, solved on this one
	//		process and then again distributed

	//	CASE a): We solve the problem in parallel (or normally for sequential code)
		#ifdef UG_PARALLEL
		d.set_storage_type(PST_ADDITIVE);
		UG_DLOG(LIB_DISC_MULTIGRID, 2, " Starting Base solver on level "<<lev<< ".\n");

		if( m_bBaseParallel ||
		   (d.get_vertical_slave_layout().empty() &&
			d.get_vertical_master_layout().empty()))
		{
		#endif
			GMG_PROFILE_BEGIN(GMG_BaseSolver);
			if(!m_pBaseSolver->apply(c, d))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Base solver on"
						" base level " << lev << " failed. "
						"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");

				return false;
			}

		//	UPDATE DEFECT
			if(!m_vLevData[lev].A->apply_sub(d, c))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Updating defect "
						" on base level " << lev << ". "
						"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
				return false;
			}
			GMG_PROFILE_END();
		#ifdef UG_PARALLEL
		}

	//	CASE b): We gather the processes, solve on one proc and distribute again
		else
		{
		//	gather the defect
			bool resume = gather_vertical(d);

		//	check, if this proc continues, else idle
			if(resume){
				GMG_PROFILE_BEGIN(GMG_BaseSolver);
				UG_DLOG(LIB_DISC_MULTIGRID, 2, " Base solver starting on 1 Proc.\n");
				if(!m_pBaseSolver->apply(c, d))
				{
					UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Base solver on"
							" base level " << lev << " failed. "
							"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");

					return false;
				}

			//	UPDATE DEFECT
				if(m_vLevData[m_baseLev].sel == NULL){
					if(!m_vLevData[lev].A->apply_sub(d, c))
					{
						UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Updating defect "
								" on base level " << lev << ". "
								"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
						return false;
					}
				}
				else{
					if(!m_BaseOperator.apply_sub(d, c))
					{
						UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Updating defect "
								" on base level " << lev << ". "
								"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
						return false;
					}
				}
				GMG_PROFILE_END();
				UG_DLOG(LIB_DISC_MULTIGRID, 2, " Base solver done on 1 Proc.\n");
			}

		//	broadcast the correction
			broadcast_vertical(c);
			c.set_storage_type(PST_CONSISTENT);

		//	if baseLevel == surfaceLevel, we need also d
			if(m_baseLev == m_topLev){
				d.set_storage_type(PST_CONSISTENT);
				broadcast_vertical(d);
				d.change_storage_type(PST_ADDITIVE);
			}
		}
		#endif

	//	we're done for the solution of the base solver
		return true;
	}
//	this case should never happen.
	else
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Level index below "
				" 'baseLevel' in lmgc. Aborting.\n");
		return false;
	}
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
{
// 	Cast Operator
	m_pSurfaceOp = dynamic_cast<operator_type*>(&J);

//	Check that Operator type is correct
	if(m_pSurfaceOp == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Can not cast Operator to AssembledLinearizedOperator.\n");
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
		if(!m_vLevData[lev].Projection->apply(*m_vLevData[lev-1].u, *m_vLevData[lev].u))
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
	if(!m_bFullRefined)
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
	m_pSurfaceOp = dynamic_cast<operator_type*>(&L);

//	Check that Operator type is correct
	if(m_pSurfaceOp == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Can not cast Operator to AssembledLinearizedOperator.\n");
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
	if(!m_bFullRefined)
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
	if(m_pApproxSpace->get_level_dof_distribution(m_topLev).num_dofs() ==
		m_pApproxSpace->get_surface_dof_distribution().num_dofs())
		m_bFullRefined = true;
	else
		m_bFullRefined =false;

//	Assemble coarse grid operators
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
#endif
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
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
#endif
	GMG_PROFILE_END();

	for(size_t lev = m_baseLev; lev < m_vLevData.size(); ++lev)
		write_level_debug(m_vLevData[lev].A->get_matrix(), "LevelMatrix", lev);

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

//	init mapping from surface level to top level in case of full refinement
	if(m_bFullRefined)
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

//	we're done
	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_linear_level_operator()
{
	GMG_PROFILE_FUNC();
// 	Create coarse level operators
	for(size_t lev = m_baseLev; lev < m_vLevData.size(); ++lev)
	{
	//	get dof distribution
		const dof_distribution_type& levelDD =
							m_pApproxSpace->get_level_dof_distribution(lev);

	//	skip assembling if no dofs given
		if(levelDD.num_dofs() == 0)
			continue;

	//	now set this selector to the assembling, such that only those elements
	//	will be assembled and force grid to be considered as regular
		m_pAss->set_selector(m_vLevData[lev].sel);
		m_vLevData[lev].A->force_regular_grid(true);

	//	init level operator
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
#endif
		GMG_PROFILE_BEGIN(GMG_AssLevOp);
		if(!m_vLevData[lev].A->init())
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator':"
					" Cannot init operator for level "<< lev << ".\n");
			return false;
		}
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
#endif
		GMG_PROFILE_END();

	//	remove force flag
		m_vLevData[lev].A->force_regular_grid(false);
		m_pAss->set_selector(NULL);

	//	now we copy the matrix into a new (smaller) one
		if(m_vLevData[lev].sel != NULL)
		{
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
#endif
			GMG_PROFILE_BEGIN(GMG_CopySmoothMatrix);
			UG_ASSERT(m_vLevData[lev].SmoothMat != NULL, "SmoothMat missing");
			matrix_type& mat = m_vLevData[lev].A->get_matrix();
			matrix_type& smoothMat = m_vLevData[lev].SmoothMat->get_matrix();

			smoothMat.resize( m_vLevData[lev].vMap.size(), m_vLevData[lev].vMap.size());
			CopySmoothingMatrix(smoothMat, m_vLevData[lev].vMapMat, mat);
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
#endif
		GMG_PROFILE_END();
		}
	}

//	get dof distribution
	const dof_distribution_type& levelDD =
						m_pApproxSpace->get_level_dof_distribution(m_baseLev);

#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
#endif
//	assemble base operator
	if(levelDD.num_dofs() != 0)
	{
		GMG_PROFILE_BEGIN(GMG_AssBaseSolver);
		m_BaseOperator.set_discretization(*m_pAss);

	//	set dof distribution to level operator
		if(!m_BaseOperator.set_dof_distribution(levelDD))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator': "
					"Cannot set dof distribution on baselevel.\n");
			return false;
		}

	//	force grid to be considered as regular
		m_BaseOperator.force_regular_grid(true);

	//	init level operator
		if(!m_BaseOperator.init())
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator':"
					" Cannot init operator for baselevel.\n");
			return false;
		}

	//	remove force flag
		m_BaseOperator.force_regular_grid(false);
		GMG_PROFILE_END();
	}
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
#endif

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
	//	now set this selector to the assembling, such that only those elements
	//	will be assembled and force grid to be considered as regular
		m_pAss->set_selector(m_vLevData[lev].sel);
		m_vLevData[lev].A->force_regular_grid(true);

	//	init operator (i.e. assemble matrix)
		if(!m_vLevData[lev].A->init(*m_vLevData[lev].u))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_non_linear_level_operator':"
					" Cannot init coarse grid operator for level "<< lev << ".\n");
			return false;
		}

	//	remove force flag
		m_vLevData[lev].A->force_regular_grid(false);
		m_pAss->set_selector(NULL);

	//	now we copy the matrix into a new (smaller) one
		if(m_vLevData[lev].sel != NULL)
		{
			UG_ASSERT(m_vLevData[lev].SmoothMat != NULL, "SmoothMat missing");
			matrix_type& mat = m_vLevData[lev].A->get_matrix();
			matrix_type& smoothMat = m_vLevData[lev].SmoothMat->get_matrix();

			smoothMat.resize( m_vLevData[lev].vMap.size(), m_vLevData[lev].vMap.size());
			CopySmoothingMatrix(smoothMat, m_vLevData[lev].vMapMat, mat);
		}
	}

//	get dof distribution
	const dof_distribution_type& levelDD =
						m_pApproxSpace->get_level_dof_distribution(m_baseLev);

//	assemble base operator
	if(levelDD.num_dofs() != 0)
	{
		m_BaseOperator.set_discretization(*m_pAss);

	//	set dof distribution to level operator
		if(!m_BaseOperator.set_dof_distribution(levelDD))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator': "
					"Cannot set dof distribution on baselevel.\n");
			return false;
		}

	//	force grid to be considered as regular
		m_BaseOperator.force_regular_grid(true);

	//	init level operator
		if(!m_BaseOperator.init(*m_vLevData[m_baseLev].u))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator':"
					" Cannot init operator for baselevel.\n");
			return false;
		}

	//	remove force flag
		m_BaseOperator.force_regular_grid(false);
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
	//	set levels
		if(!m_vLevData[lev].Prolongation->set_levels(lev-1, lev))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_prolongation':"
					" Cannot set level in interpolation matrices for level "
					<< lev << ", aborting.\n");
			return false;
		}

	//	init prolongation
		if(!m_vLevData[lev].Prolongation->init())
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
	//	set levels
		if(!m_vLevData[lev].Projection->set_levels(lev-1, lev))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_projection':"
					" Cannot set level in projection matrices for level "
					<< lev << ", aborting.\n");
			return false;
		}

	//	init projection
		if(!m_vLevData[lev].Projection->init())
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
		if(m_vLevData[lev].SmoothMat == NULL)
		{
			if(!m_vLevData[lev].Smoother->init(*m_vLevData[lev].A, *m_vLevData[lev].su))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::init_smoother':"
						" Cannot init smoother for level "<< lev << ".\n");
				return false;
			}
		}
		else
		{
			if(!m_vLevData[lev].Smoother->init(*m_vLevData[lev].SmoothMat, *m_vLevData[lev].u))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::init_smoother':"
						" Cannot init smoother for level "<< lev << ".\n");
				return false;
			}
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
// 	Prepare base solver
	if(m_pApproxSpace->get_level_dof_distribution(m_baseLev).num_dofs() == 0)
		return true;

	if(m_bBaseParallel){
		if(!m_pBaseSolver->init(*m_vLevData[m_baseLev].A, *m_vLevData[m_baseLev].u))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_base_solver':"
					" Cannot init base solver on baselevel "<< m_baseLev << ".\n");
			return false;
		}
	}
	else{
		#ifdef UG_PARALLEL
		vector_type& d = *m_vLevData[m_baseLev].d;
		if((!d.get_master_layout().empty() || !d.get_slave_layout().empty()) &&
		   (d.get_vertical_slave_layout().empty() &&
			d.get_vertical_master_layout().empty()))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_base_solver': "
					" Cannot init base solver on baselevel "<< m_baseLev << ":\n"
					" Base level distributed among processes and no possibility"
					" of grouping (vert. interfaces) present. But a serial"
					" solving is required.\n");
			return false;
		}
		#endif

		if(!m_pBaseSolver->init(m_BaseOperator, *m_vLevData[m_baseLev].u))
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
	if(m_bFullRefined){
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
	if(m_bFullRefined){
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
top_level_required(size_t topLevel)
{
//	allocated level if needed
	while(num_levels() <= topLevel)
	{
		m_vLevData.resize(m_vLevData.size()+1);
	}

//	free level if needed
	while(num_levels() > topLevel+1)
	{
		m_vLevData.back().free();
		m_vLevData.pop_back();
	}

//	reinit all levels
	for(size_t lev = m_baseLev; lev < m_vLevData.size(); ++lev)
	{
		m_vLevData[lev].allocate(lev,
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
typename AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::base_type*
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
clone()
{
	AssembledMultiGridCycle<TApproximationSpace, TAlgebra>* clone =
		new AssembledMultiGridCycle<TApproximationSpace, TAlgebra>();

	clone->set_approximation_space(*m_pApproxSpace);
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
		m_vLevData[lev].CoarseGridContribution->resize(0,0);

//	if the grid is fully refined, nothing to do
	if(m_bFullRefined) return true;

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
		m_vLevData[lev].CoarseGridContribution->resize(dofDistr.num_dofs(),
		                                               dofDistr.num_dofs());
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
	for(size_t lev = 0; lev < m_vLevData.size(); ++lev)
	{
	//	select all elements, that have a shadow as a subelement, but are not itself
	//	a shadow
		Selector sel(m_pApproxSpace->get_domain().get_grid());
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
		// \todo: Currently assemble_linear needs a solution. Remove this and do
		//	not use tmp vector here
			vector_type tmpVec; tmpVec.resize(m_pApproxSpace->get_surface_dof_distribution().num_dofs());
			m_pAss->assemble_jacobian(surfMat, tmpVec, m_pApproxSpace->get_surface_dof_distribution());
		}

	//	write matrix for debug purpose
		std::stringstream ss; ss << "MissingSurfMat_" << lev;
		write_surface_debug(surfMat, ss.str().c_str());

	//	remove the selector from the assembling procedure
		m_pAss->set_selector(NULL);

	//	project
		if(!CopyMatrixSurfaceToLevel(*m_vLevData[lev].CoarseGridContribution,
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
		write_level_debug(*m_vLevData[lev].CoarseGridContribution, "MissingLevelMat", lev);

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
		m_Com.receive_data(t.get_vertical_slave_layout(),	cpVecCopy);
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
void
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
LevData::
allocate(size_t lev,
			  approximation_space_type& approxSpace,
              assemble_type& ass,
              smoother_type& smoother,
              projection_operator_type& projection,
              prolongation_operator_type& prolongation)
{
//	get dof distribution
	pLevDD = &approxSpace.get_level_dof_distribution(lev);

//	allocate memory if not already present
	if(!u) u = new vector_type;
	if(!c) c = new vector_type;
	if(!d) d = new vector_type;
	if(!t) t = new vector_type;

//	create vectors and matrix
	const size_t numDoFs = pLevDD->num_dofs();
	u->resize(numDoFs);
	c->resize(numDoFs);
	d->resize(numDoFs);
	t->resize(numDoFs);

//	IN PARALLEL: add parallel informations
#ifdef UG_PARALLEL
	CopyLayoutsAndCommunicatorIntoVector(*u, *pLevDD);
	CopyLayoutsAndCommunicatorIntoVector(*c, *pLevDD);
	CopyLayoutsAndCommunicatorIntoVector(*d, *pLevDD);
	CopyLayoutsAndCommunicatorIntoVector(*t, *pLevDD);

	if(!has_ghosts())
	{
		sc = NULL; sd = NULL; st = NULL; sel = NULL;
	}
	else
	{
		vMapMat.resize(pLevDD->num_dofs(), 1);

		SetLayoutValues(&vMapMat, pLevDD->get_vertical_master_layout(), -1);
		SetLayoutValues(&vMapMat, pLevDD->get_master_layout(), 1);
		SetLayoutValues(&vMapMat, pLevDD->get_slave_layout(), 1);

		vMap.clear();
		for(size_t j = 0; j < vMapMat.size(); ++j)
		{
			if(vMapMat[j] < 0) continue;

			vMapMat[j] = vMap.size();
			vMap.push_back(j);
		}

		if(!su) su = new vector_type;
		if(!sc) sc = new vector_type;
		if(!sd) sd = new vector_type;
		if(!st) st = new vector_type;

		const size_t numSmoothDoFs = vMap.size();
		su->resize(numSmoothDoFs);
		sc->resize(numSmoothDoFs);
		sd->resize(numSmoothDoFs);
		st->resize(numSmoothDoFs);

		if(masterLayout) delete masterLayout;
		if(slaveLayout) delete slaveLayout;

		masterLayout = new IndexLayout;
		slaveLayout = new IndexLayout;

	//	copy layouts
		*masterLayout = pLevDD->get_master_layout();
		*slaveLayout = pLevDD->get_slave_layout();

	//	Replace indices in the layout with the smaller (smoothing patch) indices
		ReplaceIndicesInLayout(*masterLayout, vMapMat);
		ReplaceIndicesInLayout(*slaveLayout, vMapMat);

	//	copy layouts and constructors
		CopyLayoutsAndCommunicatorIntoVector(*su, *pLevDD);
		CopyLayoutsAndCommunicatorIntoVector(*sc, *pLevDD);
		CopyLayoutsAndCommunicatorIntoVector(*sd, *pLevDD);
		CopyLayoutsAndCommunicatorIntoVector(*st, *pLevDD);

	//	replace old layouts by new modified ones
		sc->set_layouts(*masterLayout, *slaveLayout);
		su->set_layouts(*masterLayout, *slaveLayout);
		sd->set_layouts(*masterLayout, *slaveLayout);
		st->set_layouts(*masterLayout, *slaveLayout);


	//	create a new selector
		if(!sel) sel = new Selector;
		sel->assign_grid(approxSpace.get_domain().get_grid());

	//	get distributed Grid manager
		DistributedGridManager* pDstGrdMgr
			= approxSpace.get_domain().get_distributed_grid_manager();

	//	select all ghost geometric objects
		SelectNonGhosts<VertexBase>(*sel, *pDstGrdMgr,
		                         pLevDD->template begin<VertexBase>(),
		                         pLevDD->template end<VertexBase>());
		SelectNonGhosts<EdgeBase>(*sel, *pDstGrdMgr,
		                         pLevDD->template begin<EdgeBase>(),
		                         pLevDD->template end<EdgeBase>());
		SelectNonGhosts<Face>(*sel, *pDstGrdMgr,
		                         pLevDD->template begin<Face>(),
		                         pLevDD->template end<Face>());
		SelectNonGhosts<Volume>(*sel, *pDstGrdMgr,
		                         pLevDD->template begin<Volume>(),
		                         pLevDD->template end<Volume>());

		if(!SmoothMat) SmoothMat = new PureMatrixOperator<vector_type, vector_type, matrix_type>;
		SmoothMat->get_matrix().set_master_layout(*masterLayout);
		SmoothMat->get_matrix().set_slave_layout(*slaveLayout);
	}
#endif

	if(!CoarseGridContribution) CoarseGridContribution = new matrix_type;

//	allocate level operator
	if(!A) A = new operator_type;
	A->set_discretization(ass);
	A->set_dof_distribution(*pLevDD);

	if(!Smoother) Smoother = smoother.clone();
	if(!Projection) Projection = projection.clone();
	if(!Prolongation) Prolongation = prolongation.clone();
}

template <typename TApproximationSpace, typename TAlgebra>
void
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
LevData::
free()
{
//	free operators if allocated
	if(A) delete A;
	if(Smoother) delete Smoother;
	if(Projection) delete Projection;
	if(Prolongation) delete Prolongation;

//	free algebra
	if(u) delete u;
	if(c) delete c;
	if(d) delete d;
	if(t) delete t;
	if(CoarseGridContribution) delete CoarseGridContribution;
	if(SmoothMat) delete SmoothMat;

	if(su) delete su;
	if(sc) delete sc;
	if(sd) delete sd;
	if(st) delete st;

#ifdef UG_PARALLEL
	if(masterLayout) delete masterLayout;
	if(slaveLayout) delete slaveLayout;
#endif

	if(sel) delete sel;
}





} // namespace ug


#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__ */
