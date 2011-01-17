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

// \todo: Remove
// typedef for non-parallel case (to make it compileable in serial)
#ifndef UG_PARALLEL
	typedef boost::function<int (int)>		Callback_ProcessIDToSubdomainID;
#endif

template <typename TDomain>
class IApproximationSpace
{
	public:
	//	Domain type
		typedef TDomain domain_type;

	// 	Subset Handler type
		typedef typename domain_type::subset_handler_type subset_handler_type;

	public:
	// constructor
		IApproximationSpace()
			: m_pDomain(NULL), m_pMGSubsetHandler(NULL), m_pFunctionPattern(NULL)
		{};

		virtual ~IApproximationSpace()	{}

	//	Assign domain
		void assign_domain(domain_type& domain)
		{
			m_pDomain = &domain;
			m_pMGSubsetHandler = &domain.get_subset_handler();
		}

	//	Assign Function Pattern
		bool assign_function_pattern(FunctionPattern& fp)
		{
			if(!fp.is_locked())
			{
				UG_LOG("Function pattern not locked.");
				return false;
			}

			m_pFunctionPattern = &fp;
			return true;
		}

	// 	Return Function Pattern
		const FunctionPattern& get_function_pattern() const {return *m_pFunctionPattern;}

	// 	Return the domain
		const domain_type& get_domain() const {return *m_pDomain;}

	// 	Return the domain
		domain_type& get_domain() {return *m_pDomain;}

		virtual void enable_domain_decomposition(Callback_ProcessIDToSubdomainID)	{}
		virtual bool domain_decomposition_enabled()						{return false;}

	protected:
	// 	Domain, where solution lives
		domain_type* m_pDomain;

	// grid or multigrid or subsethandler, where elements are stored
		subset_handler_type* m_pMGSubsetHandler;

	// 	Function pattern
		FunctionPattern* m_pFunctionPattern;
};

struct UG_ERROR_DoFDistributionMissing{};

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
class ApproximationSpace : public IApproximationSpace<TDomain>{
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

		// Type of DoF Distribution
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

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
		ApproximationSpace() :
			m_bInit(false)
		{};

		~ApproximationSpace(){}

		bool init()
		{
			if(this->m_pFunctionPattern == NULL)
			{
				UG_LOG("No Function Pattern assigned to Approximation Space.\n");
				return false;
			}

			if(this->m_pMGSubsetHandler == NULL)
			{
				UG_LOG("No domain assigned to Approximation Space.\n");
				return false;
			}

			if(!m_MGDoFManager.assign_multi_grid_subset_handler(*(this->m_pMGSubsetHandler)))
			{
				UG_LOG("In 'ApproximationSpace::init':"
						" Cannot assign multi grid subset handler.\n");
				return false;
			}
			if(!m_MGDoFManager.assign_function_pattern(*(this->m_pFunctionPattern)))
			{
				UG_LOG("In 'ApproximationSpace::init':"
						" Cannot assign Function Pattern.\n");
				return false;
			}

#ifdef UG_PARALLEL
			m_MGDoFManager.set_distributed_grid_manager(
					*this->m_pDomain->get_distributed_grid_manager());
#endif

			if(!m_MGDoFManager.distribute_dofs())
			{
				UG_LOG("In 'ApproximationSpace::init':"
						" Cannot distribute dofs.\n");
				return false;
			}

			m_bInit = true;
			return true;
		}

		void print_statistic() const
		{
			m_MGDoFManager.print_statistic();
		}

		// create a new grid function of this approximation space
		function_type* create_level_function(std::string name, size_t level, bool allocate = true)
		{
			if(!m_bInit)
			{
				UG_LOG("Approximation Space not initialized.\n");
				return NULL;
			}

			const dof_distribution_type* dofDistr = m_MGDoFManager.get_level_dof_distribution(level);
			if(dofDistr == NULL)
			{
				throw(UG_ERROR_DoFDistributionMissing());
			}

			function_type* gridFct = new function_type(name, *this,
														*dofDistr,
														allocate);
			return gridFct;
		}

		// create a new grid function of this approximation space
		function_type* create_surface_function(std::string name, bool allocate = true)
		{
			if(!m_bInit)
			{
				UG_LOG("Approximation Space not initialized.\n");
				return NULL;
			}

			const dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();
			if(dofDistr == NULL)
			{
				throw(UG_ERROR_DoFDistributionMissing());
			}

			function_type* gridFct = new function_type(name, *this,
														*dofDistr,
														allocate);
			return gridFct;
		}

		const dof_distribution_type& get_surface_dof_distribution() const
		{
			if(!m_bInit)
				throw(UG_ERROR_DoFDistributionMissing());

			const dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();

			if(dofDistr == NULL)
				throw(UG_ERROR_DoFDistributionMissing());

			return *dofDistr;
		}

		const dof_distribution_type& get_level_dof_distribution(size_t level) const
		{
			if(!m_bInit)
				throw(UG_ERROR_DoFDistributionMissing());

			return *(m_MGDoFManager.get_level_dof_distribution(level));
		}

	///	Influences the way in which parallel interfaces are created from the grids interfaces.
	/**	If enabled, process and subdomain layouts are generated.*/
		virtual void enable_domain_decomposition(
						Callback_ProcessIDToSubdomainID cb_ProcIDToSubdomID)
		{
			#ifdef UG_PARALLEL
				m_MGDoFManager.enable_domain_decomposition(cb_ProcIDToSubdomID);
			#endif
		}

		virtual bool domain_decomposition_enabled()
		{
			#ifdef UG_PARALLEL
				return m_MGDoFManager.domain_decomposition_enabled();
			#endif
			return false;
		}
	protected:
	//	Init flag
		bool m_bInit;

	// 	Dof manager used for this Approximation Space
		dof_manager_type m_MGDoFManager;
};


}


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__ */
