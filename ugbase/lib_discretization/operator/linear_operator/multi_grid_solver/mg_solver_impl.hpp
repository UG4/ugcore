/*
 * mg_solver_impl.h
 *
 *  Created on: 04.01.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__

#include "common/profiler/profiler.h"
#include "lib_discretization/io/vtkoutput.h"
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

	// Presmooth
	for(int i = 0; i < nu; ++i)
	{
		// compute correction of one smoothing step (defect is updated d:= d - A m_t[l])
		if(!m_vSmoother[lev]->apply_update_defect(*m_t[lev], d))
		{
			UG_LOG("Error in smoothing step " << i+1 << " on level " << lev << ".\n");
			return false;
		}

		// add correction of smoothing step to level correction
		c += *m_t[lev];
	}

	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
lmgc(size_t lev)
{
	// reset correction
	m_c[lev]->set(0.0);

	if(lev > m_baseLevel)
	{
	// 	Presmooth
		GMG_PROFILE_BEGIN(GMG_PreSmooth);
		if(!smooth(*m_d[lev], *m_c[lev], lev, m_nu1))
		{
			UG_LOG("Error in premoothing on level " << lev << ".\n");
			return false;
		}
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
		if(!m_vProlongation[lev-1]->apply_transposed(*m_d[lev-1], *m_d[lev]))
		{
			UG_LOG("Error in restriction from level " << lev << " to " << lev-1 << ".\n");
			return false;
		}

		bool resume = true;

		#ifdef UG_PARALLEL
		//	send vertical-slaves -> vertical-masters
		//	one proc may not have both, a vertical-slave- and vertical-master-layout.
			ComPol_VecAdd<typename function_type::vector_type> cpVecAdd(m_d[lev-1]);
			if(!m_d[lev-1]->get_vertical_slave_layout().empty()){
				resume = false;
				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2, " Going down: SENDS vertical dofs on level " << lev -1 << ".\n");
				m_Com.send_data(m_d[lev-1]->get_vertical_slave_layout(), cpVecAdd);
			}
			else if(!m_d[lev-1]->get_vertical_master_layout().empty()){

				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2, " Going down:  WAITS FOR RECIEVE of vertical dofs on level " << lev -1 << ".\n");
				m_Com.receive_data(m_d[lev-1]->get_vertical_master_layout(), cpVecAdd);
			}
			m_Com.communicate();
		#endif

		if(resume)
		{
			// apply lmgc on coarser grid
			for(int i = 0; i < m_cycle_type; ++i)
			{
				if(!lmgc(lev-1))
				{
					UG_LOG("Error in lmgc on level " << lev-1 << ".\n");
					return false;
				}
			}
		}

		#ifdef UG_PARALLEL
			//	send vertical-masters -> vertical-slaves
			//	one proc may not have both, a vertical-slave- and vertical-master-layout.
			ComPol_VecCopy<typename function_type::vector_type> cpVecCopy(m_c[lev-1]);
			if(!m_c[lev-1]->get_vertical_slave_layout().empty()){
				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2, " Going up: WAITS FOR RECIEVE of vertical dofs on level " << lev -1 << ".\n");
				m_Com.receive_data(m_c[lev-1]->get_vertical_slave_layout(),	cpVecCopy);

				m_c[lev-1]->set_storage_type(PST_CONSISTENT);
			}
			else if(!m_c[lev-1]->get_vertical_master_layout().empty()){
				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2, " Going up: SENDS vertical dofs on level " << lev -1 << ".\n");
				m_Com.send_data(m_c[lev-1]->get_vertical_master_layout(), cpVecCopy);
				m_c[lev-1]->set_storage_type(PST_CONSISTENT);
			}
		m_Com.communicate();
		#endif

	//	Interpolate Correction
		if(!m_vProlongation[lev-1]->apply(*m_t[lev], *m_c[lev-1]))
		{
			UG_LOG("Error in prolongation from level " << lev-1 << " to " << lev << ".\n");
			return false;
		}

	// 	Add coarse grid correction to level correction
		*m_c[lev] += *m_t[lev];

	//	Update Defect
		if(!m_A[lev]->apply_sub(*m_d[lev], *m_t[lev]))
		{
			UG_LOG("Error in updating defect on level " << lev << ".\n");
			return false;
		}

	// 	Postsmooth
		GMG_PROFILE_BEGIN(GMG_PostSmooth);
		if(!smooth(*m_d[lev], *m_c[lev], lev, m_nu2))
		{
			UG_LOG("Error in postsmoothing on level " << lev << ".\n");
			return false;
		}
		GMG_PROFILE_END();

		return true;
	}
	else if(lev == m_baseLevel)
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
		if(!m_pBaseSolver->apply(*m_c[lev], *m_d[lev]))
			{UG_LOG("Error in base solver on level " << lev << ".\n"); return false;}
		//PROFILE_END();

	//update defect
		if(!m_A[lev]->apply_sub(*m_d[lev], *m_c[lev]))
			{UG_LOG("Error in updating defect on level " << lev << ".\n"); return false;}
		GMG_PROFILE_END();

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
	m_Op = dynamic_cast<AssembledLinearOperator<dof_distribution_type, algebra_type>*>(&J);

//	Check that Operator type is correct
	if(m_Op == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Can not cast Operator to AssembledLinearizedOperator.\n");
		return false;
	}

//	init common
	GMG_PROFILE_BEGIN(GMG_InitCommon);
	if(!init_common(true))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Can init common part.\n");
		return false;
	}
	GMG_PROFILE_END();

//	Project current Solution from surface grid onto the level grids
	GMG_PROFILE_BEGIN(GMG_ProjectSolutionFromSurface);
	if(!ProjectionSurfaceLevel<dof_distribution_type, vector_type>::
			surface_to_level(*m_u[m_surfaceLevel], m_pApproxSpace->get_level_dof_distribution(m_surfaceLevel),
								u, m_pApproxSpace->get_surface_dof_distribution()))
	{
		UG_LOG("Projection of solution failed in AssembledMultiGridCycle::prepare\n"); return false;
	}
	GMG_PROFILE_END();

// 	Project solution from surface grid to coarser grid levels
	GMG_PROFILE_BEGIN(GMG_ProjectSolutionDown);
	for(size_t lev = m_surfaceLevel; lev != m_baseLevel; --lev)
	{
		if(!m_vProjection[lev-1]->apply(*m_u[lev-1], *m_u[lev]))
			{UG_LOG("ERROR while projecting solution to coarse grid function of level "<< lev -1 << ", aborting.\n");return false;}
	}
	GMG_PROFILE_END();

// 	Create coarse level operators
	GMG_PROFILE_BEGIN(GMG_AssembleCoarseGridMatrices);
	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		if(!m_A[lev]->set_dof_distribution(m_pApproxSpace->get_level_dof_distribution(lev)))
			{UG_LOG("ERROR while setting dof distribution for coarse operator on level "<< m_baseLevel << ", aborting.\n");return false;}

		if(!m_A[lev]->init(*m_u[lev]))
			{UG_LOG("ERROR while constructing coarse grid matrices for level "<< lev << ", aborting.\n");return false;}
	}
	GMG_PROFILE_END();

//	Init smoother and base solver for coarse grid operators
	GMG_PROFILE_BEGIN(GMG_InitBaseAndSmoother);
	if(!init_smoother_and_base_solver())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Can init Smoother and Base Solver.\n");
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
	m_Op = dynamic_cast<AssembledLinearOperator<dof_distribution_type, algebra_type>*>(&L);

//	Check that Operator type is correct
	if(m_Op == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Can not cast Operator to AssembledLinearizedOperator.\n");
		return false;
	}

//	init common
	if(!init_common(false))
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Can init common part.\n");
		return false;
	}

// 	Create coarse level operators
	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		if(!m_A[lev]->set_dof_distribution(m_pApproxSpace->get_level_dof_distribution(lev)))
			{UG_LOG("ERROR while setting dof distribution for coarse operator on level "<< m_baseLevel << ", aborting.\n");return false;}

		// todo: we should not pass u here
		if(!m_A[lev]->init(*m_u[lev]))
			{UG_LOG("ERROR while constructing coarse grid matrices for level "<< lev << ", aborting.\n");return false;}
	}

//	Init smoother and base solver for coarse grid operators
	if(!init_smoother_and_base_solver())
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle:init': Can init Smoother and Base Solver.\n");
		return false;
	}

	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_common(bool nonlinear)
{
//	Preform some checks:
	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Discretization not set.\n");
		return false;
	}
	if(m_pApproxSpace == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Approximation Space not set.\n");
		return false;
	}
	if(m_pBaseSolver == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Base Solver not set.\n");
		return false;
	}
	if(m_vSmoother[0] == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Smoother not set.\n");
		return false;
	}
	if(m_vProlongation[0] == NULL)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Prolongation not set.\n");
		return false;
	}
	if(nonlinear)
		if(m_vProjection[0] == NULL)
		{
			UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Projection not set, although problem nonlinear.\n");
			return false;
		}

	if(m_baseLevel > m_surfaceLevel)
	{
		UG_LOG("ERROR in 'AssembledMultiGridCycle::init': Base Level can not be greater than surface level.\n");
		return false;
	}

// 	If grid may be changed, we reallocate all memory
	if(m_grid_changes && m_allocated)
		if(!free_memory())
		{
			UG_LOG("AssembledMultiGridCycle::prepare: Cannot free memory. Aborting.\n");
			return false;
		}

//	Reallocate memory if needed
	bool reallocated = false;
	if(!(m_allocated))
	{
		if(!allocate_memory())
		{
			UG_LOG("AssembledMultiGridCycle::prepare: Cannot allocate memory. Aborting.\n");
			return false;
		}
		reallocated = true;
	}

// 	Create Interpolation and Projection Operators
	if(reallocated)
	{
		for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
		{
			if(!m_vProlongation[lev]->set_levels(lev, lev+1))
				{UG_LOG("ERROR during setup for interpolation matrices for level "<< lev << ", aborting.\n");return false;}
			if(!m_vProlongation[lev]->init())
				{UG_LOG("ERROR while constructing interpolation matrices for level "<< lev << ", aborting.\n");return false;}

			if(nonlinear)
			{
				if(!m_vProjection[lev]->set_levels(lev, lev+1))
					{UG_LOG("ERROR during setup for projection matrices for level "<< lev << ", aborting.\n");return false;}
				if(!m_vProjection[lev]->init())
					{UG_LOG("ERROR while constructing projection matrices for level "<< lev << ", aborting.\n");return false;}
			}
		}
	}

	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
init_smoother_and_base_solver()
{
// 	Init smoother
//	Todo: maybe we should admit also smoothers depending explicitly on the current solution
//	todo: since currently m_u may be uninitiallized here if linear operator was passed
	GMG_PROFILE_BEGIN(GMG_InitSmoother);
	for(size_t lev = m_baseLevel; lev < m_surfaceLevel; ++lev)
	{
		if(!m_vSmoother[lev]->init(*m_A[lev], *m_u[lev]))
			{UG_LOG("ERROR while initializing smoother for level "<< lev << ", aborting.\n");return false;}
	}

// 	TODO: For non-adaptive refinement this should use the passed (already assembled) operator
	m_A[m_surfaceLevel] = m_Op;
	if(!m_vSmoother[m_surfaceLevel]->init(*m_A[m_surfaceLevel], *m_u[m_surfaceLevel]))
		{UG_LOG("ERROR while initializing smoother for level "<< m_surfaceLevel << ", aborting.\n");return false;}
	GMG_PROFILE_END();

// 	Prepare base solver
	GMG_PROFILE_BEGIN(GMG_BaseSolver);
	if(!m_pBaseSolver->init(*m_A[m_baseLevel], *m_u[m_baseLevel]))
		{UG_LOG("ERROR while initializing base solver on level "<< m_baseLevel << ", aborting.\n");return false;}
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
		UG_LOG("MultiGridSolver requires a Multigrid. Please use MultiGrid structure.\n");
		return false;
	}

// 	Check if surface level has been choosen correctly
	if(m_surfaceLevel >= mg.num_levels())
	{
		UG_LOG("SurfaceLevel " << m_surfaceLevel << " does not exist.\n");
		return false;
	}

// 	Check if base level has been choose correctly
	if(m_baseLevel > m_surfaceLevel)
	{
		UG_LOG("Base level must be smaller or equal to surface Level.\n");
		return false;
	}

//	Project Defect from surface grid onto the level grids
	if(!ProjectionSurfaceLevel<dof_distribution_type, vector_type>::
			surface_to_level(*m_d[m_surfaceLevel], m_pApproxSpace->get_level_dof_distribution(m_surfaceLevel),
								d, m_pApproxSpace->get_surface_dof_distribution()))
	{
		UG_LOG("Projection of defect failed in AssembledMultiGridCycle::apply\n"); return false;
	}

//	Project Correction from surface grid onto the level grids
	if(!ProjectionSurfaceLevel<dof_distribution_type, vector_type>::
			surface_to_level(*m_c[m_surfaceLevel], m_pApproxSpace->get_level_dof_distribution(m_surfaceLevel),
								c, m_pApproxSpace->get_surface_dof_distribution()))
	{
		UG_LOG("Projection of correction failed in AssembledMultiGridCycle::apply\n"); return false;
	}

// 	Perform one multigrid cycle
	if(!lmgc(m_surfaceLevel))
		{UG_LOG("MultiGridCycle: Error in step. Aborting.\n"); return false;}

//	Project Defect from level grids to surface grid
	if(!ProjectionSurfaceLevel<dof_distribution_type, vector_type>::
			level_to_surface(d, m_pApproxSpace->get_surface_dof_distribution(),
								*m_d[m_surfaceLevel], m_pApproxSpace->get_level_dof_distribution(m_surfaceLevel)))
	{
		UG_LOG("Projection of defect failed in AssembledMultiGridCycle::apply\n"); return false;
	}

//	Project Correction from level grids to surface grid
	if(!ProjectionSurfaceLevel<dof_distribution_type, vector_type>::
			level_to_surface(c, m_pApproxSpace->get_surface_dof_distribution(),
								*m_c[m_surfaceLevel], m_pApproxSpace->get_level_dof_distribution(m_surfaceLevel)))
	{
		UG_LOG("Projection of correction failed in AssembledMultiGridCycle::apply\n"); return false;
	}

	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
allocate_memory()
{
	// dynamically created pointer for Coarse Grid Functions
	m_u.resize(m_surfaceLevel+1);
	m_c.resize(m_surfaceLevel+1);
	m_t.resize(m_surfaceLevel+1);
	m_d.resize(m_surfaceLevel+1);

	// top level matrix and vectors
	m_u[m_surfaceLevel] = m_pApproxSpace->create_level_function(m_surfaceLevel);
	m_c[m_surfaceLevel] = m_pApproxSpace->create_level_function(m_surfaceLevel);
	m_t[m_surfaceLevel] = m_pApproxSpace->create_level_function(m_surfaceLevel);
	m_d[m_surfaceLevel] = m_pApproxSpace->create_level_function(m_surfaceLevel);

	// create coarse level vectors
	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
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
	m_A.resize(m_surfaceLevel + 1);
	m_vProlongation.resize(m_surfaceLevel);
	m_vProjection.resize(m_surfaceLevel);

	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
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
	m_vSmoother.resize(m_surfaceLevel+1);
	for(size_t lev = m_baseLevel; lev <= m_surfaceLevel; ++lev)
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

	for(size_t j = m_baseLevel; j != m_surfaceLevel; ++j)
	{
		if(j < m_A.size())
			if(m_A[j] != NULL) {delete m_A[j]; m_A[j] = NULL;};

		if(j < m_vProlongation.size() && j != 0)
			if(m_vProlongation[j] != NULL) {delete m_vProlongation[j]; m_vProlongation[j] = NULL;};

		if(j < m_vProjection.size() && j != 0)
			if(m_vProjection[j] != NULL) {delete m_vProjection[j]; m_vProjection[j] = NULL;};

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
	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
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
