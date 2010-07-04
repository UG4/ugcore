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
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

template <typename TApproximationSpace, typename TAlgebra>
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
AssembledMultiGridCycle(	IAssemble<level_function_type, algebra_type>& ass, approximation_space_type& approxSpace,
							size_t surfaceLevel, size_t baseLevel, int cycle_type,
							smoother_type& smoother, int nu1, int nu2, base_solver_type& baseSolver, bool grid_changes) :
				m_ass(ass), m_approxSpace(approxSpace), m_domain(approxSpace.get_domain()),
				m_surfaceLevel(surfaceLevel), m_baseLevel(baseLevel), m_cycle_type(cycle_type),
				m_nu1(nu1), m_nu2(nu2), m_baseSolver(baseSolver),
				m_grid_changes(grid_changes), m_allocated(false)

				{
					m_smoother.resize(1);
					m_smoother[0] = &smoother;
				};


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
smooth(level_function_type& d, level_function_type& c, size_t lev, int nu)
{

	// Presmooth
	for(int i = 0; i < nu; ++i)
	{
		// compute correction of one smoothing step (defect is updated d:= d - A m_t[l])
		if(!m_smoother[lev]->apply(d, *m_t[lev]))
			{UG_LOG("Error in smoothing step " << i << " on level " << lev << ".\n"); return false;}

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
		// presmooth
		if(!smooth(*m_d[lev], *m_c[lev], lev, m_nu1))
			{UG_LOG("Error in premoothing on level " << lev << ".\n"); return false;}

		#ifdef UG_PARALLEL
			typename level_function_type::dof_manager_type& dofMgr = m_d[lev-1]->get_dof_manager();
			if(!dofMgr.get_vertical_master_layout(lev-1).empty()){
			//	set all dofs to 0. This is important since we will add vertical slave values
			//	after restriction.
				ConsistentToUnique( &m_d[lev-1]->get_vector(),
									dofMgr.get_vertical_master_layout(lev-1));
			}
		#endif

		// restrict defect
		if(!m_I[lev-1]->apply_transposed(*m_d[lev-1], *m_d[lev]))
			{UG_LOG("Error in restriction from level " << lev << " to " << lev-1 << ".\n"); return false;}

		bool resume = true;

		#ifdef UG_PARALLEL
		//	send vertical-slaves -> vertical-masters
		//	one proc may not have both, a vertical-slave- and vertical-master-layout.
			ComPol_VecAdd<typename level_function_type::vector_type> cpVecAdd(&m_d[lev-1]->get_vector());
			if(!dofMgr.get_vertical_slave_layout(lev-1).empty()){
				resume = false;
				m_Com.send_data(dofMgr.get_vertical_slave_layout(lev-1), cpVecAdd);
			}
			else if(!dofMgr.get_vertical_master_layout(lev-1).empty()){

				m_Com.receive_data(dofMgr.get_vertical_master_layout(lev-1), cpVecAdd);
			}
			m_Com.communicate();
		#endif

		if(resume)
		{
			// apply lmgc on coarser grid
			for(int i = 0; i < m_cycle_type; ++i)
			{
				if(!lmgc(lev-1))
					{UG_LOG("Error in lmgc on level " << lev-1 << ".\n"); return false;}
			}
		}

		#ifdef UG_PARALLEL
			//	send vertical-masters -> vertical-slaves
			//	one proc may not have both, a vertical-slave- and vertical-master-layout.
			ComPol_VecCopy<typename level_function_type::vector_type> cpVecCopy(&m_c[lev-1]->get_vector());
			if(!dofMgr.get_vertical_slave_layout(lev-1).empty()){
				m_Com.receive_data(dofMgr.get_vertical_slave_layout(lev-1),	cpVecCopy);

				m_c[lev-1]->set_storage_type(PST_CONSISTENT);
			}
			else if(!dofMgr.get_vertical_master_layout(lev-1).empty()){
				m_Com.send_data(dofMgr.get_vertical_master_layout(lev-1), cpVecCopy);
				m_c[lev-1]->set_storage_type(PST_CONSISTENT);
			}
		m_Com.communicate();
		#endif

		//interpolate correction
		if(!m_I[lev-1]->apply(*m_t[lev], *m_c[lev-1]))
			{UG_LOG("Error in prolongation from level " << lev-1 << " to " << lev << ".\n"); return false;}

		// add coarse grid correction to level correction
		*m_c[lev] += *m_t[lev];

		//update defect
		if(!m_A[lev]->apply_sub(*m_t[lev], *m_d[lev]))
			{UG_LOG("Error in updating defect on level " << lev << ".\n"); return false;}

		// postsmooth
		if(!smooth(*m_d[lev], *m_c[lev], lev, m_nu2))
			{UG_LOG("Error in postsmoothing on level " << lev << ".\n"); return false;}

		return true;
	}
	else if(lev == m_baseLevel)
	{
		// set d to be additive
		m_d[lev]->set_storage_type(PST_ADDITIVE);

		// solve on base level
		m_c[lev]->set(0.0);

		PROFILE_BEGIN(baseSolver);
		if(!m_baseSolver.apply(*m_d[lev], *m_c[lev]))
			{UG_LOG("Error in base solver on level " << lev << ".\n"); return false;}

		PROFILE_END();

		// no update of the defect is needed, since the linear solver interface requires, that the updated defect is already returned.
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
init(ILinearizedOperator<surface_function_type,surface_function_type>& A)
{
	// TODO: This is full refined problem only
	m_Op = dynamic_cast<AssembledLinearizedOperator<surface_function_type>*>(&A);

	if(m_Op == NULL)
	{
		UG_LOG("AssembledMultiGridCycle.init: Can not cast Operator to matrix based operator.\n");
		return false;
	}

	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
apply(surface_function_type& d, surface_function_type &c)
{
	// TODO: currently the matrix on finest level exists twice

	// TODO: project defect onto grids
	m_d[m_surfaceLevel]->project_surface(d);
	m_c[m_surfaceLevel]->project_surface(c);

	// perform one multigrid cycle
	if(!lmgc(m_surfaceLevel))
		{UG_LOG("MultiGridCycle: Error in step. Aborting.\n"); return false;}

	// TODO: Project correction to surface correction
	m_d[m_surfaceLevel]->release_surface(d);
	m_c[m_surfaceLevel]->release_surface(c);

	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
prepare(surface_function_type &u, surface_function_type& d, surface_function_type &c)
{
	typename TApproximationSpace::domain_type::grid_type& mg = m_domain.get_grid();

	// check if grid type is a Multigrid.
	if(dynamic_cast<MultiGrid*>(&mg) == NULL)
	{
		UG_LOG("MultiGridSolver requires a Multigrid. Please use MultiGrid structure.\n");
		return false;
	}

	// check if surface level has been choosen correctly
	if(m_surfaceLevel >= mg.num_levels())
	{
		UG_LOG("SurfaceLevel " << m_surfaceLevel << " does not exist.\n");
		return false;
	}

	// Check if base level has been choose correctly
	if(m_baseLevel > m_surfaceLevel)
	{
		UG_LOG("Base level must be smaller than surface Level.\n");
		return false;
	}

	// if grid may be changed, we reallocate all memory
	if(m_grid_changes && m_allocated)
	{
		if(free_memory() != true)
		{
			UG_LOG("AssembledMultiGridCycle::prepare: Cannot free memory. Aborting.\n");
			return false;
		}
	}

	bool reallocated = false;

	if(!(m_allocated))
	{
		if(allocate_memory() != true)
		{
			UG_LOG("AssembledMultiGridCycle::prepare: Cannot allocate memory. Aborting.\n");
			return false;
		}
		reallocated = true;
	}

	// top level matrix and vectors
	m_u[m_surfaceLevel]->project_surface(u);

	// create Interpolation and Projection Operators
	if(reallocated)
	{
		for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
		{
			if(m_I[lev]->prepare(*m_d[lev], *m_d[lev+1]) != true)
				{UG_LOG("ERROR while constructing interpolation matrices for level "<< lev << ", aborting.\n");return false;}
			if(m_P[lev]->prepare(*m_u[lev], *m_u[lev+1]) != true)
				{UG_LOG("ERROR while constructing projection matrices for level "<< lev << ", aborting.\n");return false;}
		}
	}

	// project solution from surface grid to coarser grid levels
	for(size_t lev = m_surfaceLevel; lev != m_baseLevel; --lev)
	{
		if(m_P[lev-1]->apply_transposed(*m_u[lev-1], *m_u[lev]) != true)
			{UG_LOG("ERROR while projecting solution to coarse grid function of level "<< lev -1 << ", aborting.\n");return false;}
	}

	// create coarse level operators
	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		if(m_A[lev]->prepare(*m_u[lev], *m_c[lev], *m_d[lev]) != true)
			{UG_LOG("ERROR while constructing coarse grid matrices for level "<< lev << ", aborting.\n");return false;}
	}

	// init smoother
	for(size_t lev = m_baseLevel; lev < m_surfaceLevel; ++lev)
	{
		if(!m_smoother[lev]->init(*m_A[lev]))
			{UG_LOG("ERROR while initializing smoother for level "<< lev << ", aborting.\n");return false;}
		if(!m_smoother[lev]->prepare(*m_u[lev], *m_d[lev], *m_t[lev]))
			{UG_LOG("ERROR while preparing smoother for level "<< lev << ", aborting.\n");return false;}
	}

	// TODO: For non-adaptive refinement this should use the passed (already assembled) operator
	m_A[m_surfaceLevel] = m_Op;
	if(!m_smoother[m_surfaceLevel]->init(*m_A[m_surfaceLevel]))
		{UG_LOG("ERROR while initializing smoother for level "<< m_surfaceLevel << ", aborting.\n");return false;}
	if(!m_smoother[m_surfaceLevel]->prepare(*m_u[m_surfaceLevel], *m_d[m_surfaceLevel], *m_t[m_surfaceLevel]))
		{UG_LOG("ERROR while preparing smoother for level "<< m_surfaceLevel << ", aborting.\n");return false;}

	// prepare base solver
	if(!m_baseSolver.init(*m_A[m_baseLevel]))
		{UG_LOG("ERROR while initializing base solver on level "<< m_baseLevel << ", aborting.\n");return false;}
	if(!m_baseSolver.prepare(*m_u[m_baseLevel], *m_d[m_baseLevel], *m_c[m_baseLevel]))
		{UG_LOG("ERROR while preparing base solver on level "<< m_baseLevel << ", aborting.\n");return false;}

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
	m_u[m_surfaceLevel] = m_approxSpace.create_level_grid_function("u", m_surfaceLevel, false);
	m_c[m_surfaceLevel] = m_approxSpace.create_level_grid_function("c", m_surfaceLevel, false);
	m_t[m_surfaceLevel] = m_approxSpace.create_level_grid_function("t", m_surfaceLevel);
	m_d[m_surfaceLevel] = m_approxSpace.create_level_grid_function("d", m_surfaceLevel, false);

	// create coarse level vectors
	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		// create solution
		m_u[lev] = m_approxSpace.create_level_grid_function("u", lev);

		// create correction
		m_c[lev] = m_approxSpace.create_level_grid_function("c", lev);

		// create help vector
		m_t[lev] = m_approxSpace.create_level_grid_function("t", lev);

		// create defect
		m_d[lev] = m_approxSpace.create_level_grid_function("d", lev);
	}

	//	dynamically created pointer for Coarse Operators
	m_A.resize(m_surfaceLevel + 1);
	m_I.resize(m_surfaceLevel);
	m_P.resize(m_surfaceLevel);

	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		// create prolongation operators
		m_I[lev] = new prolongation_operator_type(m_approxSpace, m_ass, lev);
		// create prolongation operators
		m_P[lev] = new projection_operator_type(m_approxSpace, m_ass, lev);
		// create coarse grid matrices
		m_A[lev] = new level_operator_type(m_ass);
	}

	// create smoother for all level
	m_smoother.resize(m_surfaceLevel+1);
	for(size_t lev = m_baseLevel; lev <= m_surfaceLevel; ++lev)
	{
		if(lev == 0) continue;

		m_smoother[lev] = m_smoother[0]->clone();
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

		if(j < m_I.size())
			if(m_I[j] != NULL) {delete m_I[j]; m_I[j] = NULL;};

		if(j < m_P.size())
			if(m_P[j] != NULL) {delete m_P[j]; m_P[j] = NULL;};

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
	m_I.clear();

	m_u.clear();
	m_c.clear();
	m_t.clear();
	m_d.clear();

	// delete smoother
	for(size_t lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		if(lev == 0) continue;
		delete m_smoother[lev];
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
