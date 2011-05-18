/*
 * mg_solver_impl.h
 *
 *  Created on: 04.01.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__

#include "common/profiler/profiler.h"
#include "projection_surface_level.h"
#include "lib_discretization/function_spaces/grid_function_util.h"
#include "lib_discretization/dof_manager/dof_manager_util.h"

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

// perform the smoothing
template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
smooth(function_type& d, function_type& c, size_t lev, int nu)
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
			if(!m_vSmoother[lev]->apply_update_defect(*m_t[lev], d))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Smoothing step "
						<< i+1 << " on level " << lev << " failed.\n");
				return false;
			}

		// 	add correction of smoothing step to level correction
		//	(Note: we do not work on c directly here, since we update the defect
		//	       after every smoothing step. The summed up correction corresponds
		//		   to the total correction of the whole smoothing.)
			c += *m_t[lev];
		}
		else
	//	This is the adaptive case. Here, we must ensure, that the added correction
	//	is zero on the adaptively refined patch boundary of this level
		{
		// 	Compute Correction of one smoothing step, but do not update defect
		//	a)  Compute t = B*d with some iterator B
			if(!m_vSmoother[lev]->apply(*m_t[lev], d))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Smoothing step "
						<< i+1 << " on level " << lev << " failed.\n");
				return false;
			}

		//	First we reset the correction to zero on the patch boundary.
			if(!SetZeroOnShadowingVertex(*m_t[lev], *m_pApproxSpace, lev))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Could not "
						" reset the values to zero on patch boundary for correction "
						<< i+1 << " on level " << lev << ".\n");
				return false;
			}

		//	now, we can update the defect with this correction ...
			if(!m_A[lev]->apply_sub(d, *m_t[lev]))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::smooth': Could not "
						" update defect for patch correction in smooth step "
						<< i+1 << " on level " << lev << ".\n");
				return false;
			}

		//	... and add the correction to to overall correction
			c += *m_t[lev];
		}
	}

//	we're done
	return true;
}

// performs a  multi grid cycle on the level
template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
lmgc(size_t lev)
{
  GMG_PROFILE_FUNC();
// 	reset correction to zero on this level
	m_c[lev]->set(0.0);

//	switch, if base level is reached. If so, call base Solver, else we will
//	perform smoothing, restrict the defect and call the lower level; then,
//	going up again in the level hierarchy the correction is interpolated and
//	used as coarse grid correction. Finally a post-smooth is performed.
	if(lev > m_baseLev)
	{
	// 	PRE-SMOOTHING
	//	We start the multi grid cycle on this level by smoothing the defect. This
	//	means that we compute a correction c, such that the defect is "smoother".
		GMG_PROFILE_BEGIN(GMG_PreSmooth);
		if(!smooth(*m_d[lev], *m_c[lev], lev, m_numPreSmooth))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Pre-Smoothing on "
					"level " << lev << " failed. "
				"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
			return false;
		}
		GMG_PROFILE_END();

	//	PARALLEL CASE: We have to take special care of the vertical interfaces.
	//	We set all DoFs to 0 on the vertical masters on the coarse level before
	//	we restrict the defect. This is important since we will add vertical
	//	slave values after restriction.
		#ifdef UG_PARALLEL
		GMG_PROFILE_BEGIN(GMG_PreSmConsistentToUnique);
		if(!m_d[lev-1]->get_vertical_master_layout().empty())
		{
			ConsistentToUnique(m_d[lev-1], m_d[lev-1]->get_vertical_master_layout());
		}
		GMG_PROFILE_END(); // GMG_PreSmConsistentToUnique
		#endif

	// 	RESTRICT DEFECT
	//	Now we can restrict the defect from the fine level to the coarser level.
	//	This is done using the transposed prolongation.
		GMG_PROFILE_BEGIN(GMG_RestrictDefect);
		if(!m_vProlongation[lev-1]->apply_transposed(*m_d[lev-1], *m_d[lev]))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Restriction of "
					"Defect from level "<<lev<<" to "<<lev-1<<" failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
			return false;
		}
		GMG_PROFILE_END();

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
		bool resume = true;
	//	send vertical-slaves -> vertical-masters
	//	one proc may not have both, a vertical-slave- and vertical-master-layout.
		GMG_PROFILE_BEGIN(GMG_SendVerticalDefect);
		ComPol_VecAdd<typename function_type::vector_type> cpVecAdd(m_d[lev-1]);
		if(!m_d[lev-1]->get_vertical_slave_layout().empty())
		{
		//	do not resume if vertical slaves are present
			resume = false;
			UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
			  " Going down: SENDS vertical dofs on level " << lev -1 << ".\n");

		//	schedule Sending of DoFs of vertical slaves
			m_Com.send_data(m_d[lev-1]->get_vertical_slave_layout(), cpVecAdd);
		}
		else if(!m_d[lev-1]->get_vertical_master_layout().empty())
		{
			UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
			 " Going down:  WAITS FOR RECIEVE of vertical dofs on level " << lev -1 << ".\n");

		//	schedule Receive of DoFs on vertical masters
			m_Com.receive_data(m_d[lev-1]->get_vertical_master_layout(), cpVecAdd);
		}

	//	perform communication
		m_Com.communicate();
		GMG_PROFILE_END();

	//	only continue if levels left
		if(resume) {
		#endif


	// 	COMPUTE COARSE GRID CORRECTION
	//	Now, we have to compute the coarse grid correction, i.e. the correction
	//	on the coarser level. This is done iteratively by calling lmgc again for
	//	the coarser level.
		for(int i = 0; i < m_cycleType; ++i)
		{
			if(!lmgc(lev-1))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Linear multi"
						" grid cycle on level " << lev-1 << " failed. "
						"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
				return false;
			}
		}

	//	PARALLEL CASE: Receive values of correction for vertical slaves
	//	If there are vertical slaves/masters on the coarser level, we now copy
	//	the correction values from the master DoFs to the slave	DoFs.
		#ifdef UG_PARALLEL
		}
	//	send vertical-masters -> vertical-slaves
	//	one proc may not have both, a vertical-slave- and vertical-master-layout.
		GMG_PROFILE_BEGIN(GMG_SendVerticalCorrection);
		ComPol_VecCopy<typename function_type::vector_type> cpVecCopy(m_c[lev-1]);
		if(!m_c[lev-1]->get_vertical_slave_layout().empty())
		{
			UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
			 " Going up: WAITS FOR RECIEVE of vertical dofs on level " << lev -1 << ".\n");

		//	schedule slaves to receive correction
			m_Com.receive_data(m_c[lev-1]->get_vertical_slave_layout(),	cpVecCopy);

		//	the correction is consistent then
			m_c[lev-1]->set_storage_type(PST_CONSISTENT);
		}
		else if(!m_c[lev-1]->get_vertical_master_layout().empty())
		{
			UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
			 " Going up: SENDS vertical dofs on level " << lev -1 << ".\n");

		//	schedule masters to send correction
			m_Com.send_data(m_c[lev-1]->get_vertical_master_layout(), cpVecCopy);

		//	the correction is consistent then
			m_c[lev-1]->set_storage_type(PST_CONSISTENT);
		}
	//	communicate
		m_Com.communicate();
		GMG_PROFILE_END();
		#endif

	//	INTERPOLATE CORRECTION
	//	now we can interpolate the coarse grid correction from the coarse level
	//	to the fine level
		GMG_PROFILE_BEGIN(GMG_InterpolateCorr);
		if(!m_vProlongation[lev-1]->apply(*m_t[lev], *m_c[lev-1]))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Prolongation from"
					" level " << lev-1 << " to " << lev << " failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");

			return false;
		}
		GMG_PROFILE_END();

	// 	ADD COARSE GRID CORRECTION
		GMG_PROFILE_BEGIN(GMG_AddCoarseGridCorr);
		*m_c[lev] += *m_t[lev];
		GMG_PROFILE_END(); // GMG_AddCoarseGridCorr

	//	UPDATE DEFECT FOR COARSE GRID CORRECTION
	//	the correction has changed c := c + t. Thus, we also have to update
	//	the defect d := d - A*t
		GMG_PROFILE_BEGIN(GMG_UpdateDefectForCGCorr);
		if(!m_A[lev]->apply_sub(*m_d[lev], *m_t[lev]))
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
			m_t[lev-1]->set(0.0);
			if(!m_vCoarseContributionMat[lev-1]->apply(*m_t[lev-1], *m_c[lev-1]))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Could not compute"
						" missing update defect contribution on level "<<lev-1<<".\n");
				return false;
			}

			*m_t[lev-1] *= -1.0;

		//	b) interpolate the coarse defect up
			if(!AddProjectionOfVertexShadows(*m_d[lev], *m_t[lev-1], *m_pApproxSpace,
			                                 lev-1, lev))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Could not add"
						" missing update defect contribution to level "<<lev<<".\n");
				return false;
			}
		}

	// 	POST-SMOOTHING
	//	We smooth the updated defect againt. This means that we compute a
	//	correction c, such that the defect is "smoother".
		GMG_PROFILE_BEGIN(GMG_PostSmooth);
		if(!smooth(*m_d[lev], *m_c[lev], lev, m_numPostSmooth))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Post-Smoothing on"
					" level " << lev << " failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
			return false;
		}
		GMG_PROFILE_END();

	//	we are done on this level
		return true;
	}

//	if the base level has been reached, the coarse problem is solved exactly
	else if(lev == m_baseLev)
	{
	// 	PARALLEL CASE: the d is additive
		#ifdef UG_PARALLEL
		GMG_PROFILE_BEGIN(GMG_SetStorageTypeBeforeBS);
		m_d[lev]->set_storage_type(PST_ADDITIVE);
		GMG_PROFILE_END(); // GMG_SetStorageTypeBeforeBS
		#endif

	//	begin profiling
		GMG_PROFILE_BEGIN(GMG_BaseSolver);

	//	some debug output
		#ifdef UG_PARALLEL
		UG_DLOG(LIB_DISC_MULTIGRID, 2, " Starting Base solver on level "<<lev<< ".\n");
		#endif

	//	SOLVE BASE PROBLEM
		if(!m_pBaseSolver->apply(*m_c[lev], *m_d[lev]))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Base solver on"
					" base level " << lev << " failed. "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");

			return false;
		}

	//	UPDATE DEFECT
		if(!m_A[lev]->apply_sub(*m_d[lev], *m_c[lev]))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::lmgc': Updating defect "
					" on base level " << lev << ". "
					"(BaseLev="<<m_baseLev<<", TopLev="<<m_topLev<<")\n");
			return false;
		}

	//	end profiling
		GMG_PROFILE_END();

	//	some debug output
		#ifdef UG_PARALLEL
		UG_DLOG(LIB_DISC_MULTIGRID, 2, " Base solver done.\n");
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
	if(!project_surface_to_level(m_u, u))
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
		if(!m_vProjection[lev-1]->apply(*m_u[lev-1], *m_u[lev]))
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
	if(!top_level_required(m_topLev))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init':"
				" Cannot allocate memory. Aborting.\n");
		return false;
	}

//	init common
	if(!init_common(false))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Cannot init common part.\n");
		return false;
	}

//	assemble missing coarse grid matrix contribution (only in adaptive case)
	if(!m_bFullRefined)
		if(!init_missing_coarse_grid_coupling(NULL))
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
	GMG_PROFILE_BEGIN(GMG_InitCoarseGridOperator);
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

	for(size_t lev = m_baseLev; lev < m_A.size(); ++lev)
		write_level_debug(m_A[lev]->get_matrix(), "LevelMatrix", lev);

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
	for(size_t lev = m_baseLev; lev < m_A.size(); ++lev)
	{
	//	skip if no operator needed
		if(m_A[lev] == NULL) continue;

	//	get dof distribution
		const dof_distribution_type& levelDD =
							m_pApproxSpace->get_level_dof_distribution(lev);

	//	skip assembling if no dofs given
		if(levelDD.num_dofs() == 0)
			continue;

	//	set dof distribution to level operator
		if(!m_A[lev]->set_dof_distribution(levelDD))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator': "
					"Cannot set dof distribution on level "<< lev << ".\n");
			return false;
		}

	//	force grid to be considered as regular
		m_A[lev]->force_regular_grid(true);

	//	init level operator
		if(!m_A[lev]->init())
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init_linear_level_operator':"
					" Cannot init operator for level "<< lev << ".\n");
			return false;
		}

	//	remove force flag
		m_A[lev]->force_regular_grid(false);
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
	for(size_t lev = m_baseLev; lev < m_A.size(); ++lev)
	{
	//	skip if no operator needed
		if(m_A[lev] == NULL) continue;

	//	get dof distribution
		const dof_distribution_type& levelDD =
							m_pApproxSpace->get_level_dof_distribution(lev);

	//	set correct dof distribution
		if(!m_A[lev]->set_dof_distribution(levelDD))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_non_linear_level_operator':"
					" Cannot set dof distribution for coarse operator on level "
					<<m_baseLev<<".\n");
			return false;
		}

	//	force grid to be considered as regular
		m_A[lev]->force_regular_grid(true);

	//	init operator (i.e. assemble matrix)
		if(!m_A[lev]->init(*m_u[lev]))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_non_linear_level_operator':"
					" Cannot init coarse grid operator for level "<< lev << ".\n");
			return false;
		}

	//	remove force flag
		m_A[lev]->force_regular_grid(false);
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
	for(size_t lev = m_baseLev; lev < m_topLev; ++lev)
	{
	//	skip if no prolongation needed
		if(m_vProlongation[lev] == NULL) continue;

	//	set levels
		if(!m_vProlongation[lev]->set_levels(lev, lev+1))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_prolongation':"
					" Cannot set level in interpolation matrices for level "
					<< lev << ", aborting.\n");
			return false;
		}

	//	init prolongation
		if(!m_vProlongation[lev]->init())
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
	for(size_t lev = m_baseLev; lev < m_topLev; ++lev)
	{
	//	skip if no projection needed
		if(m_vProjection[lev] == NULL) continue;

	//	set levels
		if(!m_vProjection[lev]->set_levels(lev, lev+1))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_projection':"
					" Cannot set level in projection matrices for level "
					<< lev << ", aborting.\n");
			return false;
		}

	//	init projection
		if(!m_vProjection[lev]->init())
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
	for(size_t lev = m_baseLev; lev < m_vSmoother.size(); ++lev)
	{
	//	skip if no smoother needed
		if(m_vSmoother[lev] == NULL) continue;

		if(!m_vSmoother[lev]->init(*m_A[lev], *m_u[lev]))
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
// 	Prepare base solver
	if(m_pApproxSpace->get_level_dof_distribution(m_baseLev).num_dofs() != 0)
		if(!m_pBaseSolver->init(*m_A[m_baseLev], *m_u[m_baseLev]))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_base_solver':"
					" Cannot init base solver on baselevel "<< m_baseLev << ".\n");
			return false;
		}

//	we're done
	return true;
}

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
	GMG_PROFILE_BEGIN(GMGApply_ProjectDefectFromSurface);
	if(!project_surface_to_level(m_d, d))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of defect to level failed.\n");
		return false;
	}
	GMG_PROFILE_END(); //GMGApply_ProjectDefectFromSurface

// 	Perform one multigrid cycle
	GMG_PROFILE_BEGIN(GMGApply_lmgc);
	if(!lmgc(m_topLev))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Cannot perform multi grid cycle.\n");
		return false;
	}
	GMG_PROFILE_END(); //GMGApply_lmgc

//	project defect from level to surface
	GMG_PROFILE_BEGIN(GMGApply_ProjectDefectFromLevelToSurface);
	if(!project_level_to_surface(d, m_d))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of defect to surface failed.\n");
		return false;
	}
	GMG_PROFILE_END(); //GMGApply_ProjectDefectFromLevelToSurface

//	project correction from level to surface
	GMG_PROFILE_BEGIN(GMGApply_ProjectCorrectionFromLevelToSurface);
	if(!project_level_to_surface(c, m_c))
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

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
project_level_to_surface(vector_type& surfFunc,
                         const std::vector<function_type*>& vLevelFunc)
{
  GMG_PROFILE_FUNC();
//	level dof distributions
	const std::vector<const dof_distribution_type*>& vLevelDD =
								m_pApproxSpace->get_level_dof_distributions();

//	surface dof distribution
	const dof_distribution_type& surfDD =
								m_pApproxSpace->get_surface_dof_distribution();

//	check that surface dof distribution exists
	if(&surfDD == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::project_level_to_surface':"
				"Surface DoF Distribution missing.\n");
		return false;
	}

//	surface view
	const SurfaceView* surfView = m_pApproxSpace->get_surface_view();

//	check that surface view exists
	if(surfView == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::project_level_to_surface':"
				" Surface View missing.\n");
		return false;
	}

//	std::vector of algebra level vectors
	std::vector<const vector_type*> cvLevelVec;

//	create std::vector of const level vectors
	ExtractVectorsFromGridFunction(cvLevelVec, vLevelFunc);

//	project
	if(!ProjectLevelToSurface(surfFunc, surfDD, *surfView, cvLevelVec, vLevelDD, m_baseLev))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::project_level_to_surface': "
				"Projection of function from level to surface failed.\n");
		return false;
	}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
project_surface_to_level(std::vector<function_type*>& vLevelFunc,
                         const vector_type& surfFunc)
{
  GMG_PROFILE_FUNC();
//	level dof distributions
	const std::vector<const dof_distribution_type*>& vLevelDD =
								m_pApproxSpace->get_level_dof_distributions();

//	surface dof distribution
	const dof_distribution_type& surfDD =
								m_pApproxSpace->get_surface_dof_distribution();

//	check that surface dof distribution exists
	if(&surfDD == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::project_surface_to_level':"
				"Surface DoF Distribution missing.\n");
		return false;
	}

//	surface view
	const SurfaceView* surfView = m_pApproxSpace->get_surface_view();

//	check that surface view exists
	if(surfView == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::project_surface_to_level':"
				" Surface View missing.\n");
		return false;
	}

//	std::vector of algebra level vectors
	std::vector<vector_type*> vLevelVec;

//	create std::vector of level vectors
	ExtractVectorsFromGridFunction(vLevelVec, vLevelFunc);

//	project
	if(!ProjectSurfaceToLevel(vLevelVec, vLevelDD, surfFunc, surfDD, *surfView))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::project_surface_to_level': "
				"Projection of function from surface to level failed.\n");
		return false;
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
		if(!allocate_level(num_levels()))
			return false;

//	free level if needed
	while(num_levels() > topLevel+1)
		free_level(num_levels());

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
allocate_level(size_t level)
{
//	allocate grid functions
	if(m_u.size() <= level) m_u.resize(level+1, NULL);
	if(m_u[level] == NULL)
		m_u[level] = m_pApproxSpace->create_level_function(level);

	if(m_c.size() <= level) m_c.resize(level+1, NULL);
	if(m_c[level] == NULL)
		m_c[level] = m_pApproxSpace->create_level_function(level);

	if(m_t.size() <= level) m_t.resize(level+1, NULL);
	if(m_t[level] == NULL)
		m_t[level] = m_pApproxSpace->create_level_function(level);

	if(m_d.size() <= level) m_d.resize(level+1, NULL);
	if(m_d[level] == NULL)
		m_d[level] = m_pApproxSpace->create_level_function(level);

//	allocate coarse grid operators
	if(m_A.size() <= level) m_A.resize(level+1, NULL);
	if(m_A[level] == NULL)
		m_A[level] = new operator_type(*m_pAss);

//	allocate prolongation
	if(m_vProlongation.size() <= level) m_vProlongation.resize(level+1, NULL);
	if(m_vProlongation[level] == NULL)
		m_vProlongation[level] = m_pProlongationPrototype->clone();

//	allocate projection
	if(m_vProjection.size() <= level) m_vProjection.resize(level+1, NULL);
	if(m_vProjection[level] == NULL)
		m_vProjection[level] = m_pProjectionPrototype->clone();

//	allocate smoother
	if(m_vSmoother.size() <= level) m_vSmoother.resize(level+1, NULL);
	if(m_vSmoother[level] == NULL)
		m_vSmoother[level] = m_pSmootherPrototype->clone();

//	check success
	if(m_A[level] == NULL)
	{UG_LOG("ERROR in 'AssembledMultiGridCycle::allocate_level': Cannot allocate"
			" memory for coarse grid operator on level "<<level<<".\n"); return false;}
	if(m_vProlongation[level] == NULL)
	{UG_LOG("ERROR in 'AssembledMultiGridCycle::allocate_level': Cannot allocate"
			" memory for prolongation operator on level "<<level<<".\n"); return false;}
	if(m_vProjection[level] == NULL)
	{UG_LOG("ERROR in 'AssembledMultiGridCycle::allocate_level': Cannot allocate"
			" memory for projection operator on level "<<level<<".\n"); return false;}
	if(m_vSmoother[level] == NULL)
	{UG_LOG("ERROR in 'AssembledMultiGridCycle::allocate_level': Cannot allocate"
			" memory for smoother operator on level "<<level<<".\n"); return false;}

	if(m_u[level] == NULL)
	{UG_LOG("ERROR in 'AssembledMultiGridCycle::allocate_level': Cannot allocate"
			" memory for solution on level "<<level<<".\n"); return false;}
	if(m_c[level] == NULL)
	{UG_LOG("ERROR in 'AssembledMultiGridCycle::allocate_level': Cannot allocate"
			" memory for correction on level "<<level<<".\n"); return false;}
	if(m_t[level] == NULL)
	{UG_LOG("ERROR in 'AssembledMultiGridCycle::allocate_level': Cannot allocate"
			" memory for help correction on level "<<level<<".\n"); return false;}
	if(m_d[level] == NULL)
	{UG_LOG("ERROR in 'AssembledMultiGridCycle::allocate_level': Cannot allocate"
			" memory for defect on level "<<level<<".\n"); return false;}

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
void
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
free_level(size_t level)
{
//	free operator
	if(m_A[level] != NULL) delete m_A[level];
	if(m_vProlongation[level] != NULL) delete m_vProlongation[level];
	if(m_vProjection[level] != NULL) delete m_vProjection[level];
	if(m_vSmoother[level] != NULL) delete m_vSmoother[level];

//	free grid functions
	if(m_u[level] != NULL) {delete m_u[level]; m_u[level] = NULL;}
	if(m_c[level] != NULL) {delete m_c[level]; m_c[level] = NULL;}
	if(m_t[level] != NULL) {delete m_t[level]; m_t[level] = NULL;}
	if(m_d[level] != NULL) {delete m_d[level]; m_d[level] = NULL;}

//	free unused elements from vectors
	while(m_A.back() == NULL) m_A.pop_back();
	while(m_vProlongation.back() == NULL) m_vProlongation.pop_back();
	while(m_vProjection.back() == NULL) m_vProjection.pop_back();
	while(m_vSmoother.back() == NULL) m_vSmoother.pop_back();

	while(m_u.back() == NULL) m_u.pop_back();
	while(m_c.back() == NULL) m_c.pop_back();
	while(m_t.back() == NULL) m_t.pop_back();
	while(m_d.back() == NULL) m_d.pop_back();
}


template <typename TApproximationSpace, typename TAlgebra>
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
~AssembledMultiGridCycle()
{
	top_level_required(0);
	free_level(0);
}


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
write_level_debug(const vector_type& vec, const char* filename, size_t lev)
{
//	if no debug writer set, we're done
	if(m_pDebugWriter == NULL) return true;

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
	clone->set_projection_operator(*(m_vProjection[0]));
	clone->set_prolongation_operator(*(m_vProlongation[0]));
	clone->set_smoother(*(m_vSmoother[0]));

	return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
}



template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_missing_coarse_grid_coupling(const vector_type* u)
{
//	create storage for matrix (one matrix for each level, except top level)
	size_t oldSize = m_vCoarseContributionMat.size();
//	a) delete not need matrices
	for(size_t i = m_topLev+1; i < m_vCoarseContributionMat.size(); ++i)
		if(m_vCoarseContributionMat[i])
			delete m_vCoarseContributionMat[i];
	m_vCoarseContributionMat.resize(m_topLev+1);

//	b) create needed storage
	for(size_t i = oldSize; i < m_vCoarseContributionMat.size(); ++i)
		m_vCoarseContributionMat[i] = new matrix_type;

//	clear matrices
	for(size_t lev = 0; lev < m_vCoarseContributionMat.size(); ++lev)
		m_vCoarseContributionMat[lev]->resize(0,0);

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
	for(size_t lev = 0; lev < m_vCoarseContributionMat.size(); ++lev)
	{
	//	get dof distributions on levels
		const dof_distribution_type& dofDistr = m_pApproxSpace->get_level_dof_distribution(lev);

	//	resize the matrix
		m_vCoarseContributionMat[lev]->resize(dofDistr.num_dofs(),
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
	for(size_t lev = 0; lev < m_vCoarseContributionMat.size(); ++lev)
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
		if(!CopyMatrixSurfaceToLevel(*m_vCoarseContributionMat[lev],
		                             vSurfLevelMapping[lev],
		                             surfMat))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_missing_coarse_grid_coupling': "
					"Projection of matrix from surface to level failed.\n");
			return false;
		}
	}

//	write matrix for debug purpose
	for(size_t lev = 0; lev < m_vCoarseContributionMat.size(); ++lev)
		write_level_debug(*m_vCoarseContributionMat[lev], "MissingLevelMat", lev);

/////////////
// end project
/////////////

//	we're done
	return true;
}


} // namespace ug


#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__ */
