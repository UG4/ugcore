/*
 * mg_solver_impl.h
 *
 *  Created on: 04.01.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__

namespace ug{

template <typename TApproximationSpace, typename TAlgebra>
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
AssembledMultiGridCycle(	IAssemble<algebra_type, level_function_type>& ass, approximation_space_type& approxSpace,
							uint surfaceLevel, uint baseLevel, int cycle_type,
							smoother_type& smoother, int nu1, int nu2, base_solver_type& baseSolver, bool grid_changes) :
				m_ass(ass), m_approxSpace(approxSpace), m_domain(approxSpace.get_domain()),
				m_surfaceLevel(surfaceLevel), m_baseLevel(baseLevel), m_cycle_type(cycle_type),
				m_smoother(smoother), m_nu1(nu1), m_nu2(nu2), m_baseSolver(baseSolver),
				m_A(NULL), m_P(NULL), m_I(NULL),
				m_grid_changes(grid_changes), m_allocated(false)

				{};


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
smooth(level_function_type& d, level_function_type& c, uint l, int nu)
{
	// init smoother iterators for Matrix
	UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::smooth on level " << l << ": Init smoother ... ");
	if(m_smoother.init(*m_A[l]) != true) return false;
	UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

	// Presmooth
	for(int i = 0; i < nu; ++i)
	{
		// prepare smoother
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Prepare smoother ... ");
		if(m_smoother.prepare(*m_u[l], d, *m_t[l]) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

		// compute correction of one smoothing step (defect is updated d:= d - A m_t[l])
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Apply smoother ... ");
		if(m_smoother.apply(d, *m_t[l]) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

		// add correction of smoothing step to level correction
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Update correction ... ");
		c += *m_t[l];
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");
	}

	return true;
}


template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
lmgc(uint l)
{
	// reset correction
	m_c[l]->set(0.0);

	if(l > m_baseLevel)
	{
		//UG_DLOG(LIB_DISC_MULTIGRID, 3, "\n --- Defect on level " << l << " after restriction of defect: " << m_d[l]->two_norm() << ".");
		UG_DLOG(LIB_DISC_MULTIGRID, 10, "\n ---------- BEGIN: Defect on level " << l << ":\n" << *m_d[l] << " \n---------- END: Defect on level " << l << ".\n");

		// presmooth
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Calling Smoother ... ");
		if(smooth(*m_d[l], *m_c[l], l, m_nu1) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

		UG_DLOG(LIB_DISC_MULTIGRID, 10, "\n ---------- BEGIN: Defect on level " << l << ":\n" << *m_d[l] << " \n---------- END: Defect on level " << l << ".\n");

		// restrict defect
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Restrict defect ... ");
		if(m_I[l-1]->apply_transposed(*m_d[l-1], *m_d[l]) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

		// apply lmgc on coarser grid
		for(int i = 0; i < m_cycle_type; ++i)
		{
			if(lmgc(l-1) == false)
			{
				UG_LOG("Error in lmgc on level " << l-1 << "." << std::endl);
				return false;
			}
		}

		//interpolate correction
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Interpolate correction ... ");
		if(m_I[l-1]->apply(*m_t[l], *m_c[l-1]) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

		// add coarse grid correction to level correction
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Add coarse grid correction ... ");
		*m_c[l] += *m_t[l];
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");
		UG_DLOG(LIB_DISC_MULTIGRID, 10, "\n ---------- BEGIN: Coarse grid correction on level " << l << ":\n" << *m_t[l] << " \n---------- END: Coarse grid correction on level " << l << ".\n");

		//update defect
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Update defect ... ");
		if(m_A[l]->apply_sub(*m_t[l], *m_d[l]) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

		UG_DLOG(LIB_DISC_MULTIGRID, 3, "\n --- Defect on level " << l << " after coarse grid correction: " << m_d[l]->two_norm() << ".");
		UG_DLOG(LIB_DISC_MULTIGRID, 10, "\n ---------- BEGIN: Defect on level " << l << ":\n" << *m_d[l] << " \n---------- END: Defect on level " << l << ".\n");

		// postsmooth
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Calling Smoother ... ");
		if(smooth(*m_d[l], *m_c[l], l, m_nu2) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

		UG_DLOG(LIB_DISC_MULTIGRID, 3, "\n --- Defect on level " << l << " after post-smoothing: " << m_d[l]->two_norm() << ".");
		UG_DLOG(LIB_DISC_MULTIGRID, 10, "\n ---------- BEGIN: Defect on level " << l << ":\n" << *m_d[l] << " \n---------- END: Defect on level " << l << ".\n");

		return true;
	}
	else if(l == m_baseLevel)
	{
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Init Coarse Grid solver ... ");
		if(m_baseSolver.init(*m_A[l]) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----\n");

		UG_DLOG(LIB_DISC_MULTIGRID, 10, "\n ---------- BEGIN: Defect on level " << l << ":\n" << *m_d[l] << " \n---------- END: Defect on level " << l << ".\n");

		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Prepare Coarse Grid solver ... ");
		if(m_baseSolver.prepare(*m_u[l], *m_d[l], *m_c[l]) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----\n");
		m_d[l]->set_storage_type(GFST_ADDITIVE);

		// solve on base level
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " ---- AssembledMultiGridCycle::lmgc on level " << l << ": Apply Coarse Grid solver ... ");
		m_c[l]->set(0.0);
//		m_c[l]->set_storage_type(GFST_CONSISTENT);
		if(m_baseSolver.apply(*m_d[l], *m_c[l]) != true) return false;
		UG_DLOG(LIB_DISC_MULTIGRID, 4, " done ----.\n");

		UG_DLOG(LIB_DISC_MULTIGRID, 3, "\n --- Defect on level " << l << " after exact solving: " << m_d[l]->two_norm() << ".");
		UG_DLOG(LIB_DISC_MULTIGRID, 10, "\n ---------- BEGIN: Defect on level " << l << ":\n" << *m_d[l] << " \n---------- END: Defect on level " << l << ".\n");

		// no update of the defect is needed, since the linear solver interface requires, that the updated defect is already returned.
		// TODO: needed ????
		//m_A[l]->apply_sub(*m_c[l], *m_d[l]);

		return true;
	}
	else
	{
		UG_LOG("Level index below 'baseLevel' in lmgc. ERROR." << std::endl);
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

	// TODO: For non-adaptive refinement this should use the passed (already assembled) operator
	m_A[m_surfaceLevel] = m_Op;

	/*
	if(m_A[m_surfaceLevel]->prepare(*m_c[m_surfaceLevel], *m_d[m_surfaceLevel]) != true)
	{
		UG_LOG("ERROR while constructing coarse grid matrices for level "<< m_surfaceLevel << ", aborting.\n");
		return false;
	}
	*/

	// perform one multigrid cycle
	if(lmgc(m_surfaceLevel) == false)
	{
		UG_LOG("MultiGridCycle: Error in step. Aborting.\n");
		return false;
	}

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

	UG_DLOG(LIB_DISC_MULTIGRID, 10, "\n ---------- BEGIN: Solution on level " << m_surfaceLevel << ":\n" << *m_u[m_surfaceLevel] << " \n---------- END: Solution on level " << m_surfaceLevel << ".\n");

	// create Interpolation and Projection Operators
	if(reallocated)
	{
		for(uint lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
		{
			if(m_I[lev]->prepare(*m_d[lev], *m_d[lev+1]) != true)
			{
				UG_LOG("ERROR while constructing coarse grid matrices for level "<< lev << ", aborting.\n");
				return false;
			}
			if(m_P[lev]->prepare(*m_u[lev], *m_u[lev+1]) != true)
			{
				UG_LOG("ERROR while constructing coarse grid matrices for level "<< lev << ", aborting.\n");
				return false;
			}
		}
	}

	// project solution from surface grid to coarser grid levels
	for(uint lev = m_surfaceLevel; lev != m_baseLevel; --lev)
	{
		if(m_P[lev-1]->apply_transposed(*m_u[lev-1], *m_u[lev]) != true)
		{
			UG_LOG("ERROR while projecting solution to coarse grid function of level "<< lev -1 << ", aborting.\n");
			return false;
		}
	}

	// create coarse level operators
	for(uint lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		if(m_A[lev]->prepare(*m_u[lev], *m_c[lev], *m_d[lev]) != true)
		{
			UG_LOG("ERROR while constructing coarse grid matrices for level "<< lev << ", aborting.\n");
			return false;
		}
	}

	return true;
}

template <typename TApproximationSpace, typename TAlgebra>
bool
AssembledMultiGridCycle<TApproximationSpace, TAlgebra>::
allocate_memory()
{
	//	dynamically created pointer for Coarse Operators
	m_A = new level_operator_type*[m_surfaceLevel + 1];
	m_I = new prolongation_operator_type*[m_surfaceLevel];
	m_P = new projection_operator_type*[m_surfaceLevel];

	// dynamically created pointer for Coarse Grid Functions
	m_u = new level_function_type*[m_surfaceLevel+1];
	m_c = new level_function_type*[m_surfaceLevel+1];
	m_t = new level_function_type*[m_surfaceLevel+1];
	m_d = new level_function_type*[m_surfaceLevel+1];

	// top level matrix and vectors
	m_u[m_surfaceLevel] = m_approxSpace.create_level_grid_function("u", m_surfaceLevel, false);
	m_c[m_surfaceLevel] = m_approxSpace.create_level_grid_function("c", m_surfaceLevel, false);
	m_t[m_surfaceLevel] = m_approxSpace.create_level_grid_function("t", m_surfaceLevel);
	m_d[m_surfaceLevel] = m_approxSpace.create_level_grid_function("d", m_surfaceLevel, false);

	// create coarse level vectors
	for(uint lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
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

	for(uint lev = m_baseLevel; lev != m_surfaceLevel; ++lev)
	{
		// create prolongation operators
		m_I[lev] = new prolongation_operator_type(m_approxSpace, m_ass, lev);
		// create prolongation operators
		m_P[lev] = new projection_operator_type(m_approxSpace, m_ass, lev);
		// create coarse grid matrices
		m_A[lev] = new level_operator_type(m_ass);
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

	for(uint j = m_baseLevel; j != m_surfaceLevel; ++j)
	{
		if(m_A != NULL)
			if(m_A[j] != NULL) {delete m_A[j]; m_A[j] = NULL;};

		if(m_I != NULL)
			if(m_I[j] != NULL) {delete m_I[j]; m_I[j] = NULL;};

		if(m_P != NULL)
			if(m_P[j] != NULL) {delete m_P[j]; m_P[j] = NULL;};

		if(m_u != NULL)
			if(m_u[j] != NULL) {delete m_u[j]; m_u[j] = NULL;};
		if(m_c != NULL)
			if(m_c[j] != NULL) {delete m_c[j]; m_c[j] = NULL;};
		if(m_t != NULL)
			if(m_t[j] != NULL) {delete m_t[j]; m_t[j] = NULL;};
		if(m_d != NULL)
			if(m_d[j] != NULL) {delete m_d[j]; m_d[j] = NULL;};

	}

	if(m_A != NULL) {delete[] m_A; m_A = NULL;};
	if(m_I != NULL) {delete[] m_I; m_I = NULL;};

	if(m_u != NULL) {delete[] m_u; m_u = NULL;};
	if(m_c != NULL) {delete[] m_c; m_c = NULL;};
	if(m_t != NULL) {delete[] m_t; m_t = NULL;};
	if(m_d != NULL) {delete[] m_d; m_d = NULL;};

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
