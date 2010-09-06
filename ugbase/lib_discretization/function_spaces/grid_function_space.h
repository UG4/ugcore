/*
 * grid_function_space.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */


#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__

#ifdef UG_PARALLEL
	#include "lib_discretization/parallelization/parallelization.h"
#endif

#include "lib_discretization/dof_manager/mg_dof_manager.h"
#include "./grid_function.h"

namespace ug{

/**
 * This class describes a grid function space on a domain.
 * The Domain gives the partition of the Grid/Multigrid in terms of subsets.
 * Now the user can add fundamental Discrete functions on this subsets or unions of them.
 *
 * Once finalized, this function pattern is fixed. Internally dof indices are created.
 * Using this Approximation Space the user can create GridFunctions of the following types:
 *
 * - surface grid function = a grid function representing the space on the surface grid
 * - level grid function = a grid function representing the space on a level grid
 *  (NOTE: For a fully refined Multigrid a level grid covers the whole domain. However
 *         for a locally/adaptively refined MultiGrid the level grid solution is only
 *         living on a part of the domain)
 */
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
class ApproximationSpace{
	private:
		// to make it more readable
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> this_type;

	public:
		// domain type
		typedef TDomain domain_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::grid_type grid_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::subset_handler_type subset_handler_type;

		// algebra of this factory
		typedef TAlgebra algebra_type;

		// dof manager used
		#ifdef UG_PARALLEL
			typedef ParallelMGDoFManager<MGDoFManager<TDoFDistribution> > dof_manager_type;
		#else
			typedef MGDoFManager<TDoFDistribution> dof_manager_type;
		#endif

		// grid function type
		#ifdef UG_PARALLEL
			typedef ParallelGridFunction<GridFunction<TDomain, TDoFDistribution, TAlgebra> > function_type;
		#else
			typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> function_type;
		#endif

	public:
		ApproximationSpace(std::string name, domain_type& domain) :
			m_name(name), m_pDomain(NULL), m_pMGSubsetHandler(NULL), m_pMGDoFManager(NULL)
		{
			assign_domain(domain);
		};

		void assign_domain(domain_type& domain)
		{
			m_pDomain = &domain;
			m_pMGSubsetHandler = &domain.get_subset_handler();
		}

		bool assign_function_pattern(FunctionPattern& fp)
		{
			if(!fp.is_locked())
				{UG_LOG("Function pattern nor locked."); return false;}

			m_pFunctionPattern = &fp;
			m_pMGDoFManager = new dof_manager_type(*m_pMGSubsetHandler, fp);
			if(m_pMGDoFManager == NULL)
				{UG_LOG("Cannot allocate MGDoFManager."); return false;}

#ifdef UG_PARALLEL
			m_pMGDoFManager->set_distributed_grid_manager(
					*m_pDomain->get_distributed_grid_manager());
#endif
			if(!m_pMGDoFManager->distribute_dofs())
				{UG_LOG("Cannot distribute dofs.\n"); return false;}

			return true;
		}

		// return the domain
		const domain_type& get_domain() const {return *m_pDomain;}

		// return the domain
		domain_type& get_domain() {return *m_pDomain;}

		// create a new grid function of this approximation space
		function_type* create_level_function(std::string name, size_t level, bool allocate = true)
		{
			if(m_pMGDoFManager == NULL)
				{UG_LOG("Function pattern not set.\n"); return false;}

			function_type* gridFct = new function_type(name, *this, *m_pMGDoFManager->get_level_dof_distribution(level), allocate);
			return gridFct;
		}

		// create a new grid function of this approximation space
		function_type* create_surface_function(std::string name, bool allocate = true)
		{
			if(m_pMGDoFManager == NULL)
				{UG_LOG("Function pattern not set.\n"); return false;}

			function_type* gridFct = new function_type(name, *this, *m_pMGDoFManager->get_surface_dof_distribution(), allocate);
			return gridFct;
		}

		~ApproximationSpace()
		{
			if(m_pMGDoFManager != NULL)
				delete m_pMGDoFManager;
		}

	protected:
		// name of this Approximation Space
		std::string m_name;

		// Domain, where solution lives
		domain_type* m_pDomain;

		// grid or multigrid or subsethandler, where elements are stored
		subset_handler_type* m_pMGSubsetHandler;

		// dof manager used for this Approximation Space
		dof_manager_type* m_pMGDoFManager;

		// function pattern
		FunctionPattern* m_pFunctionPattern;
};


}


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__ */
