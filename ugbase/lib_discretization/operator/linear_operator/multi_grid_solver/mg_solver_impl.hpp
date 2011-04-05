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

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

#include <iostream>
#include <sstream>


//#define PROFILE_GMG
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
smooth(function_type& d, function_type& c, size_t lev, int nu)
{

// 	smooth nu times
	for(int i = 0; i < nu; ++i)
	{
	//	For Full-Refinement, we can work on the defect directly
		if(m_bFullRefined)
		{
		// 	compute correction of one smoothing step (defect is updated d:= d - A m_t[l])
			if(!m_vSmoother[lev]->apply_update_defect(*m_t[lev], d))
			{
				UG_LOG("ERROR in smoothing step " << i+1 << " on level " << lev << ".\n");
				return false;
			}

		// 	add correction of smoothing step to level correction
			c += *m_t[lev];
		}
	//	For Adaptive-Refinement, we have to take care for inner-level nodes
	//  On these nodes, no smoothing is performed, i.e. the correction has to
	//	be zero
		else
		{
			static int cnt = 0;

		// 	compute correction of one smoothing step (defect is updated d:= d - A m_t[l])
			if(!m_vSmoother[lev]->apply(*m_t[lev], d))
			{
				UG_LOG("ERROR in smoothing step " << i+1 << " on level " << lev << ".\n");
				return false;
			}

			char ext[200]; sprintf(ext, "GMG_CorrSmoothOrig_run%03d", cnt);

			write_level_debug(*m_t[lev], ext, lev);

		//	reset correction to zero for inner patch-bnd nodes
			const std::vector<bool>& vSkip = m_vvSkipSmooth[lev];
			for(size_t i = 0; i < vSkip.size(); ++i)
			{
				if(vSkip[i])
					(*m_t[lev])[i] = 0.0;
			}

			sprintf(ext, "GMG_CorrSmoothUsed_run%03d", cnt++);
			write_level_debug(*m_t[lev], ext, lev);

		//	compute new defect
			m_A[lev]->apply_sub(d, *m_t[lev]);

		// 	add correction of smoothing step to level correction
			c += *m_t[lev];

		}
	}

//	we're done
	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
lmgc(size_t lev)
{
	// reset correction
	m_c[lev]->set(0.0);

	if(lev > m_baseLev)
	{
	// 	Presmooth
		GMG_PROFILE_BEGIN(GMG_PreSmooth);
		write_level_debug(*m_c[lev], "GMG_CorrBeforePreSmooth", lev);
		write_level_debug(*m_d[lev], "GMG_DefBeforePreSmooth", lev);
		if(!smooth(*m_d[lev], *m_c[lev], lev, m_numPreSmooth))
		{
			UG_LOG("ERROR in presmoothing on level " << lev << ".\n");
			return false;
		}
		write_level_debug(*m_d[lev], "GMG_DefAfterPreSmooth", lev);
		write_level_debug(*m_c[lev], "GMG_CorrAfterPreSmooth", lev);
		GMG_PROFILE_END();

		#ifdef UG_PARALLEL
			if(!m_d[lev-1]->get_vertical_master_layout().empty()){
			//	set all dofs to 0. This is important since we will add vertical slave values
			//	after restriction.
				ConsistentToUnique(m_d[lev-1],
						m_d[lev-1]->get_vertical_master_layout());
			}
		#endif

	// 	Restrict Defect
		GMG_PROFILE_BEGIN(GMG_RestrictDefect);
		if(!m_vProlongation[lev-1]->apply_transposed(*m_d[lev-1], *m_d[lev]))
		{
			UG_LOG("ERROR in restriction from level " << lev << " to " << lev-1 << ".\n");
			return false;
		}
		GMG_PROFILE_END();

		bool resume = true;

		#ifdef UG_PARALLEL
		//	send vertical-slaves -> vertical-masters
		//	one proc may not have both, a vertical-slave- and vertical-master-layout.
			ComPol_VecAdd<typename function_type::vector_type> cpVecAdd(m_d[lev-1]);
			if(!m_d[lev-1]->get_vertical_slave_layout().empty()){
				resume = false;
				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
				  " Going down: SENDS vertical dofs on level " << lev -1 << ".\n");
				m_Com.send_data(m_d[lev-1]->get_vertical_slave_layout(), cpVecAdd);
			}
			else if(!m_d[lev-1]->get_vertical_master_layout().empty()){

				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
				 " Going down:  WAITS FOR RECIEVE of vertical dofs on level " << lev -1 << ".\n");
				m_Com.receive_data(m_d[lev-1]->get_vertical_master_layout(), cpVecAdd);
			}
			m_Com.communicate();
		#endif

		if(resume)
		{
			// apply lmgc on coarser grid
			for(int i = 0; i < m_cycleType; ++i)
			{
				if(!lmgc(lev-1))
				{
					UG_LOG("ERROR in lmgc on level " << lev-1 << ".\n");
					return false;
				}
			}
		}

		#ifdef UG_PARALLEL
			//	send vertical-masters -> vertical-slaves
			//	one proc may not have both, a vertical-slave- and vertical-master-layout.
			ComPol_VecCopy<typename function_type::vector_type> cpVecCopy(m_c[lev-1]);
			if(!m_c[lev-1]->get_vertical_slave_layout().empty()){
				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
				 " Going up: WAITS FOR RECIEVE of vertical dofs on level " << lev -1 << ".\n");
				m_Com.receive_data(m_c[lev-1]->get_vertical_slave_layout(),	cpVecCopy);

				m_c[lev-1]->set_storage_type(PST_CONSISTENT);
			}
			else if(!m_c[lev-1]->get_vertical_master_layout().empty()){
				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
				 " Going up: SENDS vertical dofs on level " << lev -1 << ".\n");
				m_Com.send_data(m_c[lev-1]->get_vertical_master_layout(), cpVecCopy);
				m_c[lev-1]->set_storage_type(PST_CONSISTENT);
			}
		m_Com.communicate();
		#endif

	//	Interpolate Correction
		GMG_PROFILE_BEGIN(GMG_InterpolateCorr);
		if(!m_vProlongation[lev-1]->apply(*m_t[lev], *m_c[lev-1]))
		{
			UG_LOG("ERROR in prolongation from level " << lev-1 << " to " << lev << ".\n");
			return false;
		}
		GMG_PROFILE_END();

	// 	Add coarse grid correction to level correction
		*m_c[lev] += *m_t[lev];

		write_level_debug(*m_t[lev], "GMG_CorrIntpolCoarseCorr", lev);
		write_level_debug(*m_d[lev], "GMG_DefBeforeUpdateCoarseCorr", lev);

	//	\todo: This is a lousy quick-hack for test purpose. Replace!
	//	Update Defect for correction on coarse grid
		if(!m_bFullRefined)
		{
			clear_on_hidden_values(*m_c[lev-1], lev-1);
			write_level_debug(*m_c[lev-1], "AAA_ClearedCoarseCorr", lev-1);

			vector_type dTmpCoarse; dTmpCoarse.resize(m_d[lev-1]->num_dofs());
			vector_type dTmpFine; dTmpFine.resize(m_d[lev]->num_dofs());

			//m_A[lev-1]->apply(dTmpCoarse, *m_c[lev-1]);
			matrix_type mat;
			mat.set_as_copy_of(m_A[lev-1]->get_matrix());

			const std::vector<bool>& vSkipBnd = m_vvSkipHidden[lev-1];
			for(size_t i = 0; i < mat.num_rows(); ++i)
			{
//				mat(i,i) *= 0.5;
				if(!vSkipBnd[i]) continue;

				for(typename matrix_type::row_iterator conn = mat.begin_row(i);
						conn != mat.end_row(i); ++conn)
				{
					const size_t j = conn.index();
					if(!vSkipBnd[j] ) continue;

					conn.value() *= 0.5;
//					UG_LOG("Scal i="<<i<<", j="<<j<<"\n");

				}
			}

			vector_type* pVec = dynamic_cast<vector_type*>(m_c[lev-1]);
			if(pVec == NULL) UG_ASSERT(0, "Cannot cast.");
#ifdef UG_PARALLEL
			mat.set_storage_type(PST_ADDITIVE);
			pVec->set_storage_type(PST_CONSISTENT);
#endif
			mat.apply(dTmpCoarse, *pVec);
#ifdef UG_PARALLEL
			dTmpCoarse.set_storage_type(PST_ADDITIVE);
#endif

			write_level_debug(dTmpCoarse, "AAA_dCoarse", lev-1);

			m_vProjection[lev-1]->apply_transposed(dTmpFine, dTmpCoarse);
			write_level_debug(dTmpFine, "AAA_dFine", lev);

			//dTmpFine *= 0.5;

			const std::vector<bool>& vSkip = m_vvSkipSmooth[lev];
			for(size_t i = 0; i < vSkip.size(); ++i)
			{
				if(vSkip[i])
					(*m_d[lev])[i] -= dTmpFine[i];
			}
			write_level_debug(*m_d[lev], "AAA_DefAfterSubtrCoarseCorr", lev);
		}


	//	Update Defect for correction on this grid
		if(!m_A[lev]->apply_sub(*m_d[lev], *m_t[lev]))
		{
			UG_LOG("ERROR in updating defect on level " << lev << ".\n");
			return false;
		}
		write_level_debug(*m_d[lev], "GMG_DefAfterUpdateCoarseCorr", lev);


	// 	Postsmooth
		GMG_PROFILE_BEGIN(GMG_PostSmooth);
		write_level_debug(*m_c[lev], "GMG_CorrBeforePostSmooth", lev);
		write_level_debug(*m_d[lev], "GMG_DefBeforePostSmooth", lev);
		if(!smooth(*m_d[lev], *m_c[lev], lev, m_numPostSmooth))
		{
			UG_LOG("ERROR in postsmoothing on level " << lev << ".\n");
			return false;
		}
		write_level_debug(*m_c[lev], "GMG_CorrAfterPostSmooth", lev);
		write_level_debug(*m_d[lev], "GMG_DefAfterPostSmooth", lev);
		GMG_PROFILE_END();

		return true;
	}
	else if(lev == m_baseLev)
	{
		// set d to be additive
#ifdef UG_PARALLEL
		m_d[lev]->set_storage_type(PST_ADDITIVE);
#endif

		GMG_PROFILE_BEGIN(GMG_BaseSolver);
		// solve on base level
		m_c[lev]->set(0.0);

#ifdef UG_PARALLEL
		UG_DLOG(LIB_DISC_MULTIGRID, 2, " Starting Base solver on level " << lev << ".... \n");
#endif

		//PROFILE_BEGIN(baseSolver);
		write_level_debug(*m_d[lev], "GMG_DefBeforeBase", lev);
		if(!m_pBaseSolver->apply(*m_c[lev], *m_d[lev]))
			{UG_LOG("ERROR in base solver on level " << lev << ".\n"); return false;}
		//PROFILE_END();

	//update defect
		if(!m_A[lev]->apply_sub(*m_d[lev], *m_c[lev]))
			{UG_LOG("ERROR in updating defect on level " << lev << ".\n"); return false;}
		GMG_PROFILE_END();

		write_level_debug(*m_c[lev], "GMG_CorrAfterBase", lev);
		write_level_debug(*m_d[lev], "GMG_DefAfterBase", lev);


#ifdef UG_PARALLEL
		UG_DLOG(LIB_DISC_MULTIGRID, 2, " Base solver done.\n");
#endif
		return true;
	}
	else
	{
		UG_LOG("Level index below 'baseLevel' in lmgc. ERROR.\n");
		return false;
	}
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
{
// 	Cast Operator
	m_Op = dynamic_cast<AssembledLinearOperator<dof_distribution_impl_type, algebra_type>*>(&J);

//	Check that Operator type is correct
	if(m_Op == NULL)
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

//	check, if grid is full-refined
	if(m_pApproxSpace->get_level_dof_distribution(m_topLev).num_dofs() ==
		m_pApproxSpace->get_surface_dof_distribution().num_dofs())
		m_bFullRefined = true;
	else
		m_bFullRefined =false;


//	init common
	GMG_PROFILE_BEGIN(GMG_InitCommon);
	if(!init_common(true))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Can init common part.\n");
		return false;
	}
	GMG_PROFILE_END();

//	*****   Collect needs for Projection
//	level dof distributions
	std::vector<const dof_distribution_type*> vLevelDD =
			m_pApproxSpace->get_level_dof_distributions();

//	surface dof distribution
	const dof_distribution_type& surfDD = m_pApproxSpace->get_surface_dof_distribution();

//	check that surface dof distribution exists
	if(&surfDD == NULL)
	{
		UG_LOG("Surface DoF Distribution missing.\n");
		return false;
	}

//	surface view
	const SurfaceView* surfView = m_pApproxSpace->get_surface_view();

//	check that surface view exists
	if(surfView == NULL)
	{
		UG_LOG("Surface View missing.\n");
		return false;
	}

//	std::vector of algebra level vectors
	std::vector<vector_type*> vLevelVec;

// *****	Project current Solution from surface grid onto the level grids
//	create std::vector of level vectors
	ExtractVectorsFromGridFunction(vLevelVec, m_u);

//	project
	GMG_PROFILE_BEGIN(GMG_ProjectSolutionFromSurface);
	if(!ProjectSurfaceToLevel(vLevelVec, vLevelDD, u, surfDD, *surfView))
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
					"solution to coarse grid function of level "<< lev - 1<<".\n");
			return false;
		}
	}
	GMG_PROFILE_END();

// 	Create coarse level operators
	GMG_PROFILE_BEGIN(GMG_AssembleCoarseGridMatrices);
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
	{
	//	set correct dof distribution
		if(!m_A[lev]->set_dof_distribution(m_pApproxSpace->get_level_dof_distribution(lev)))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Cannot set dof "
					"distribution for coarse operator on level "<< m_baseLev << ".\n");
			return false;
		}

	//	force grid to be considered as regular
		m_A[lev]->force_regular_grid(true);

	//	init operator (i.e. assemble matrix)
		if(!m_A[lev]->init(*m_u[lev]))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Cannot init coarse "
					"grid operator for level "<< lev << ".\n");
			return false;
		}

	//	remove force flag
		m_A[lev]->force_regular_grid(false);
	}
	GMG_PROFILE_END();

//	Init smoother and base solver for coarse grid operators
	GMG_PROFILE_BEGIN(GMG_InitBaseAndSmoother);
	if(!init_smoother_and_base_solver())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Can init Smoother and Base Solver.\n");
		return false;
	}
	GMG_PROFILE_END();

	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init(ILinearOperator<vector_type, vector_type>& L)
{
// 	Cast Operator
	m_Op = dynamic_cast<AssembledLinearOperator<dof_distribution_impl_type, algebra_type>*>(&L);

//	Check that Operator type is correct
	if(m_Op == NULL)
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

//	check, if grid is full-refined
	if(m_pApproxSpace->get_level_dof_distribution(m_topLev).num_dofs() ==
		m_pApproxSpace->get_surface_dof_distribution().num_dofs())
		m_bFullRefined = true;
	else
		m_bFullRefined =false;

//	init common
	if(!init_common(!m_bFullRefined))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Can init common part.\n");
		return false;
	}

// 	Create coarse level operators
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
	{
	//	get dof distribution
		const dof_distribution_type& levelDD =
							m_pApproxSpace->get_level_dof_distribution(lev);

	//	skip assembling if no dofs given
		if(levelDD.num_dofs() == 0)
			continue;

	//	set dof distribution to level operator
		if(!m_A[lev]->set_dof_distribution(levelDD))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Cannot set dof "
					" distribution on level "<< lev << ".\n");
			return false;
		}

	//	force grid to be considered as regular
		m_A[lev]->force_regular_grid(true);

	//	init level operator
		if(!m_A[lev]->init())
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Cannot"
					" init operator for level "<< lev << ".\n");
			return false;
		}

	//	remove force flag
		m_A[lev]->force_regular_grid(false);
	}

//	Init smoother and base solver for coarse grid operators
	if(!init_smoother_and_base_solver())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': "
				"Can init Smoother and Base Solver.\n");
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
	if(m_vSmoother[0] == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Smoother not set.\n");
		return false;
	}
	if(m_vProlongation[0] == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Prolongation not set.\n");
		return false;
	}
	if(nonlinear)
		if(m_vProjection[0] == NULL)
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
					"Projection not set, although problem nonlinear.\n");
			return false;
		}

	if(m_baseLev > m_topLev)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common': "
				"Base Level can not be greater than surface level.\n");
		return false;
	}

// 	If grid may be changed, we reallocate all memory
	if(m_grid_changes && m_allocated)
		if(!free_memory())
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common':"
					" Cannot free memory. Aborting.\n");
			return false;
		}

//	Reallocate memory if needed
	bool reallocated = false;
	if(!(m_allocated))
	{
		if(!allocate_memory())
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common':"
					" Cannot allocate memory. Aborting.\n");
			return false;
		}
		reallocated = true;
	}

// 	Create Interpolation and Projection Operators
	if(reallocated)
	{
		for(size_t lev = m_baseLev; lev != m_topLev; ++lev)
		{
			if(!m_vProlongation[lev]->set_levels(lev, lev+1))
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common':"
						" Cannot set level in interpolation matrices for level "
						<< lev << ", aborting.\n");
				return false;
			}
			if(!m_vProlongation[lev]->init())
			{
				UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common':"
						" Cannot init interpolation operator for level "
						<< lev << ", aborting.\n");
				return false;
			}

			if(nonlinear)
			{
				if(!m_vProjection[lev]->set_levels(lev, lev+1))
				{
					UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common':"
							" Cannot set level for projection operator for level "
							<< lev << ", aborting.\n");
					return false;
				}
				if(!m_vProjection[lev]->init())
				{
					UG_LOG("ERROR in 'AssembledMultiGridCycle::init_common':"
							" Cannot init  projection operator for level "
							<< lev << ", aborting.\n");
					return false;
				}
			}
		}
	}

//	we're done
	return true;
}
/*
template<typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::write_connection_viewer_matrices()
{
	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		vector<MathVector<dim> > positionsFrom, positionsTo;
		ExtractPositions(*m_u[lev], positionsFrom);
		ExtractPositions(*m_u[lev+1], positionsTo);

		P1ProlongationOperator<TApproximationSpace, TAlgebra>* Op =
						dynamic_cast< P1ProlongationOperator<TApproximationSpace, TAlgebra> *>(m_vProlongation[lev]);
		if(Op)
		{
			stringstream s; s << "prolongation" << lev << ".mat";
			WriteMatrixToConnectionViewer (Op->get_matrix(), s.str().c_str(), positionsFrom, positionsTo, dim);
		}

		stringstream s; s << "A" << lev << ".mat";
		WriteMatrixToConnectionViewer (s.str().c_str(), m_A, &positionsFrom[0], dim);
	}

	stringstream s; s << "A" << m_surfaceLevel << ".mat";
	vector<MathVector<dim> > positions;
	ExtractPositions(*m_u[m_surfaceLevel], positions);
	WriteMatrixToConnectionViewer (s.str().c_str(), m_A, &positions[0], dim);

}
*/
template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_smoother_and_base_solver()
{
// 	Init smoother
//	Todo: maybe we should admit also smoothers depending explicitly on the current solution
//	todo: since currently m_u may be uninitiallized here if linear operator was passed
	GMG_PROFILE_BEGIN(GMG_InitSmoother);
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
	{
		if(!m_vSmoother[lev]->init(*m_A[lev], *m_u[lev]))
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init_smoother_and_base_solver':"
					" Cannot init smoother for level "<< lev << ".\n");
			return false;
		}

		write_level_debug(m_A[lev]->get_matrix(), "Operator", lev);
	}
	GMG_PROFILE_END();

	// init skip flags
	if(!m_bFullRefined)
	{
	//	clear skip flags
		m_vvSkipSmooth.clear();
		m_vvSkipHidden.clear();

	//	resize skip flags for level
		m_vvSkipSmooth.resize(m_topLev+1);
		m_vvSkipHidden.resize(m_topLev+1);

	//	create skip flags for each level
		for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
		{
			create_level_skip_flags(m_vvSkipSmooth[lev], lev);
			if(lev != m_topLev)
			create_hidden_flags(m_vvSkipHidden[lev], lev);
		}
	}
// 	TODO: For non-adaptive refinement this should use the passed (already assembled) operator

// 	Prepare base solver
	GMG_PROFILE_BEGIN(GMG_BaseSolver);
	if(!m_pBaseSolver->init(*m_A[m_baseLev], *m_u[m_baseLev]))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init_smoother_and_base_solver':"
				" Cannot init base solver on baselevel "<< m_baseLev << ".\n");
		return false;
	}
	GMG_PROFILE_END();

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

//	*****   Collect needs for Projection
//	level dof distributions
	std::vector<const dof_distribution_type*> vLevelDD =
			m_pApproxSpace->get_level_dof_distributions();

//	surface dof distribution
	const dof_distribution_type& surfDD = m_pApproxSpace->get_surface_dof_distribution();

//	check that surface dof distribution exists
	if(&surfDD == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect':"
				"Surface DoF Distribution missing.\n");
		return false;
	}

//	surface view
	const SurfaceView* surfView = m_pApproxSpace->get_surface_view();

//	check that surface view exists
	if(surfView == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect':"
				" Surface View missing.\n");
		return false;
	}

//	debug output
	write_surface_debug(d, "GMG_DefectIn");

//	std::vector of algebra level vectors
	std::vector<vector_type*> vLevelVec;

// *****	Project Defect from surface grid onto the level grids
//	create std::vector of level vectors
	ExtractVectorsFromGridFunction(vLevelVec, m_d);

//	project
	if(!ProjectSurfaceToLevel(vLevelVec, vLevelDD, d, surfDD, *surfView))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of defect to level failed.\n");
		return false;
	}

// *****	Project Correction from surface grid onto the level grids
// \todo: Do we really need the projected correction or is zero everywhere ok?
//	create std::vector of level vectors
	ExtractVectorsFromGridFunction(vLevelVec, m_c);

//	project
	if(!ProjectSurfaceToLevel(vLevelVec, vLevelDD, c, surfDD, *surfView))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of correction to level failed.\n");
		return false;
	}

//	debug output
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
		write_level_debug(*m_d[lev], "GMG_DefectInProject", lev);

// *****	Perform one multigrid cycle
	if(!lmgc(m_topLev))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Cannot perform multi grid cycle.\n");
		return false;
	}

//	std::vector of algebra level vectors
	std::vector<const vector_type*> cvLevelVec;

// *****	Project Defect from the level grids to surface grid
//	create std::vector of const level vectors
	ExtractVectorsFromGridFunction(cvLevelVec, m_d);

//	project
	if(!ProjectLevelToSurface(d, surfDD, *surfView, cvLevelVec, vLevelDD))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of defect to surface failed.\n");
		return false;
	}

// *****	Project Correction from the level grid to surface grid
//	create std::vector of level vectors
	ExtractVectorsFromGridFunction(cvLevelVec, m_c);

//	project
	if(!ProjectLevelToSurface(c, surfDD, *surfView, cvLevelVec, vLevelDD))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::apply_update_defect': "
				"Projection of correction to surface failed.\n");
		return false;
	}

//	debug output
	write_surface_debug(c, "GMG_CorrOut");
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
		write_level_debug(*m_c[lev], "GMG_CorrOutProject", lev);
	write_surface_debug(d, "GMG_DefectOut");
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
		write_level_debug(*m_d[lev], "GMG_DefectOutProject", lev);

//	increase dbg counter
	m_dbgIterCnt++;

//	we're done
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
allocate_memory()
{
	// dynamically created pointer for Coarse Grid Functions
	m_u.resize(m_topLev+1);
	m_c.resize(m_topLev+1);
	m_t.resize(m_topLev+1);
	m_d.resize(m_topLev+1);

	// create coarse level vectors
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
	{
		// create solution
		m_u[lev] = m_pApproxSpace->create_level_function(lev);

		// create correction
		m_c[lev] = m_pApproxSpace->create_level_function(lev);

		// create help vector
		m_t[lev] = m_pApproxSpace->create_level_function(lev);

		// create defect
		m_d[lev] = m_pApproxSpace->create_level_function(lev);
	}

	//	dynamically created pointer for Coarse Operators
	m_A.resize(m_topLev + 1);
	m_vProlongation.resize(m_topLev+1);
	m_vProjection.resize(m_topLev+1);

	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
	{
		// create prolongation operators
		if(lev != 0)
			m_vProlongation[lev] = m_vProlongation[0]->clone();

		// create prolongation operators
		if(m_vProjection[0] != NULL && lev != 0)
			m_vProjection[lev] = m_vProjection[0]->clone();

		// create coarse grid matrices
		m_A[lev] = new operator_type(*m_pAss);
	}

	// create smoother for all level
	m_vSmoother.resize(m_topLev+1);
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
	{
		if(lev == 0) continue;

		m_vSmoother[lev] = m_vSmoother[0]->clone();
	}

	m_allocated = true;
	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
free_memory()
{
	// do nothing, if no memory allocated
	if(m_allocated == false) return true;

	for(size_t j = m_baseLev; j <= m_topLev; ++j)
	{
		if(j < m_A.size())
			if(m_A[j] != NULL) {delete m_A[j]; m_A[j] = NULL;};

		if(j < m_vProlongation.size() && j != 0)
			if(m_vProlongation[j] != NULL)
			{
				delete m_vProlongation[j];
				m_vProlongation[j] = NULL;
			}

		if(j < m_vProjection.size() && j != 0)
			if(m_vProjection[j] != NULL)
			{
				delete m_vProjection[j];
				m_vProjection[j] = NULL;
			}

		if(j < m_u.size())
			if(m_u[j] != NULL) {delete m_u[j]; m_u[j] = NULL;};
		if(j < m_c.size())
			if(m_c[j] != NULL) {delete m_c[j]; m_c[j] = NULL;};
		if(j < m_t.size())
			if(m_t[j] != NULL) {delete m_t[j]; m_t[j] = NULL;};
		if(j < m_d.size())
			if(m_d[j] != NULL) {delete m_d[j]; m_d[j] = NULL;};
	}

	m_A.clear();
	m_vProlongation.resize(1);
	m_vProjection.resize(1);

	m_u.clear();
	m_c.clear();
	m_t.clear();
	m_d.clear();

	// delete smoother
	for(size_t lev = m_baseLev; lev <= m_topLev; ++lev)
	{
		if(lev == 0) continue;
		delete m_vSmoother[lev];
	}

	m_allocated = false;
	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
~AssembledMultiGridCycle()
{
	free_memory();
}


} // namespace ug


#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__ */
