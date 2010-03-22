/*
 * mg_solver_impl.h
 *
 *  Created on: 04.01.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__

namespace ug{

template <typename TDomain, typename TAlgebra, typename TDiscreteFunction>
MultiGridCycle<TDomain, TAlgebra, TDiscreteFunction>::
MultiGridCycle(	IAssemble<algebra_type, discrete_function_type>& ass, domain_type& domain, discrete_function_type& u,
				uint surfaceLevel, uint baseLevel, int cycle_type,
				TransferOperator<TDomain, TAlgebra, typename TDiscreteFunction::dof_manager_type>& transferOperator,
				IIterativeStep<TAlgebra>& smoother, int nu1, int nu2, ILinearSolver<TAlgebra>& baseSolver) :
				m_ass(ass), m_domain(domain), m_u(u),
				m_surfaceLevel(surfaceLevel), m_baseLevel(baseLevel), m_cycle_type(cycle_type),
				m_smoother(smoother), m_nu1(nu1), m_nu2(nu2),
				m_baseSolver(baseSolver),
				m_c(u), m_d(u), m_t(u),
				m_trans(transferOperator)
				{};


template <typename TDomain, typename TAlgebra, typename TDiscreteFunction>
bool
MultiGridCycle<TDomain, TAlgebra, TDiscreteFunction>::
lmgc(matrix_type* A[], discrete_function_type& c, discrete_function_type& d, uint l)
{
	// reset correction
	c.get_vector(l).set(0.0);

	if(l > m_baseLevel)
	{
		// Presmooth
		for(int i = 0; i < m_nu1; ++i)
		{
			// compute correction of one smoothing step (defect is updated d:= d - A m_t[l])
			m_smoother.step(*A[l], m_t.get_vector(l), d.get_vector(l));

			// add correction of smoothing step to level correction
			c.get_vector(l) += m_t.get_vector(l);
		}

		// Restrict Defect
		m_trans.init(l-1);
		m_trans.applyTransposed(d, d);

		// apply lmgc on coarser grid
		for(int i = 0; i < m_cycle_type; ++i)
		{
			if(lmgc(A, c, d, l-1) == false)
			{
				UG_LOG("Error in lmgc on level " << l-1 << "." << std::endl);
				return false;
			}
		}

		//interpolate correction
		m_trans.init(l-1);
		m_trans.apply( m_t, c);

		// add coarse grid correction to level correction
		c.get_vector(l) += m_t.get_vector(l);

		//update defect
		A[l]->matmul_minus(d.get_vector(l), m_t.get_vector(l));

		// Postsmooth
		for(int i = 0; i < m_nu2; ++i)
		{
			// compute correction of one smoothing step (defect is updated d:= d - A m_t[l])
			m_smoother.step(*A[l], m_t.get_vector(l), d.get_vector(l));

			// add correction of smoothing step to level correction
			c.get_vector(l) += m_t.get_vector(l);
		}

		return true;
	}
	else if(l == m_baseLevel)
	{
		 // solve on base level
		 m_baseSolver.solve(*A[l], c.get_vector(l), d.get_vector(l));

		 // update defect
		 A[l]->matmul_minus(d.get_vector(l), c.get_vector(l));

		return true;
	}
	else
	{
		UG_LOG("Level index below 'baseLevel' in lmgc. ERROR." << std::endl);
		return false;
	}
}


template <typename TDomain, typename TAlgebra, typename TDiscreteFunction>
bool
MultiGridCycle<TDomain, TAlgebra, TDiscreteFunction>::
step(matrix_type& A, vector_type& c, vector_type& d)
{
	// set matrix_type for surface Level
	m_A[m_surfaceLevel] = &A;

	// perform one multigrid cycle
	if(lmgc(m_A, m_c, m_d, m_surfaceLevel) == false)
	{
		UG_LOG("MultiGridCycle: Error in step. Aborting.\n");
		return false;
	}

	return true;
}

template <typename TDomain, typename TAlgebra, typename TDiscreteFunction>
bool
MultiGridCycle<TDomain, TAlgebra, TDiscreteFunction>::
prepare()
{
	typename domain_type::grid_type& mg = m_domain.get_grid();

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

	// reset defect
	for(uint i = m_baseLevel; i < m_surfaceLevel; ++i)
	{
		m_d.set(0.0); //TODO: Needed ???
	}

	// assemble matrices on coarser grids
	for(uint i = m_baseLevel; i < m_surfaceLevel; ++i)
	{
		if(m_ass.assemble_jacobian(*m_A[i], m_u, i) != IAssemble_OK)
		{
			UG_LOG(" Error in assemble_stiffness matrix, aborting.\n");
			return false;
		}
	}

	// prepare transfer operator
	if(m_trans.prepare() != true) return false;

	// prepare smoother and base solver
	if(m_smoother.prepare() != true) return false;
	if(m_baseSolver.prepare() != true) return false;

	return true;
}

template <typename TDomain, typename TAlgebra, typename TDiscreteFunction>
bool
MultiGridCycle<TDomain, TAlgebra, TDiscreteFunction>::
finish()
{
	for(uint j = m_baseLevel; j < m_surfaceLevel; ++j)
	{
		delete m_A[j];
		delete m_I[j];
	}
	delete[] m_A;
	delete[] m_I;

	if(m_smoother.finish() != true) return false;
	if(m_baseSolver.finish() != true) return false;

	return true;
}


} // namespace ug


#endif /* __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_IMPL__ */
