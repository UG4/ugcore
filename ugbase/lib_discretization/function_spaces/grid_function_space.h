/*
 * grid_function_space.h
 *
 *  Created on: 19.02.2010
 *      Author: andreasvogel
 */

// TODO: remove 'Callback_ProcessIDToSubdomainID' stuff ... (27012011)

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__

#ifdef UG_PARALLEL
	#include "lib_discretization/parallelization/parallelization.h"
#endif

#include "lib_discretization/dof_manager/mg_dof_manager.h"
#include "./grid_function.h"

namespace ug{

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
			: m_pDomain(NULL), m_pMGSH(NULL), m_pFuncPattern(NULL)
		{};

		virtual ~IApproximationSpace()	{}

	//	Assign domain
		void assign_domain(domain_type& domain)
		{
			m_pDomain = &domain;
			m_pMGSH = &domain.get_subset_handler();
		}

	//	Assign Function Pattern
		bool assign_function_pattern(FunctionPattern& fp)
		{
			if(!fp.is_locked())
			{
				UG_LOG("Function pattern not locked.");
				return false;
			}

			m_pFuncPattern = &fp;
			return true;
		}

	// 	Return Function Pattern
		const FunctionPattern& get_function_pattern() const {return *m_pFuncPattern;}

	// 	Return the domain
		const domain_type& get_domain() const {return *m_pDomain;}

	// 	Return the domain
		domain_type& get_domain() {return *m_pDomain;}

	protected:
	// 	Domain, where solution lives
		domain_type* m_pDomain;

	// grid or multigrid or subsethandler, where elements are stored
		subset_handler_type* m_pMGSH;

	// 	Function pattern
		FunctionPattern* m_pFuncPattern;
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
	//	 to make it more readable
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> this_type;

	public:
	///	Type of Domain, where DoFs are defined
		typedef TDomain domain_type;

	///	Type of Grid, where DoFs are defined
		typedef typename domain_type::grid_type grid_type;

	///	Type of Subset Handler, where DoFs are defined
		typedef typename domain_type::subset_handler_type subset_handler_type;

	///	Type of Algebra used
		typedef TAlgebra algebra_type;

	///	Type of DoF Distribution
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Type of DoF Manager used
		#ifdef UG_PARALLEL
			typedef ParallelMGDoFManager<MGDoFManager<TDoFDistribution> > dof_manager_type;
		#else
			typedef MGDoFManager<TDoFDistribution> dof_manager_type;
		#endif

	///	Type of Grid function used
		#ifdef UG_PARALLEL
			typedef ParallelGridFunction<GridFunction<TDomain, TDoFDistribution, TAlgebra> > function_type;
		#else
			typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> function_type;
		#endif

	public:
	///	Constructor
		ApproximationSpace() :
			m_bInit(false)
		{};

	///	Destructor
		~ApproximationSpace(){}

	///	initializes the Approximation Space
		bool init()
		{
			if(this->m_pFuncPattern == NULL)
			{
				UG_LOG("No Function Pattern assigned to Approximation Space.\n");
				return false;
			}

			if(this->m_pMGSH == NULL)
			{
				UG_LOG("No domain assigned to Approximation Space.\n");
				return false;
			}

			if(!m_MGDoFManager.assign_multi_grid_subset_handler(*(this->m_pMGSH)))
			{
				UG_LOG("In 'ApproximationSpace::init':"
						" Cannot assign multi grid subset handler.\n");
				return false;
			}
			if(!m_MGDoFManager.assign_function_pattern(*(this->m_pFuncPattern)))
			{
				UG_LOG("In 'ApproximationSpace::init':"
						" Cannot assign Function Pattern.\n");
				return false;
			}

#ifdef UG_PARALLEL
			m_MGDoFManager.set_distributed_grid_manager(
					*this->m_pDomain->get_distributed_grid_manager());
#endif

			if(!m_MGDoFManager.enable_dofs())
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

		void print_layout_statistic() const
		{
			m_MGDoFManager.print_layout_statistic();
		}

		// create a new grid function of this approximation space
		function_type* create_level_function(size_t level)
		{
			if(!m_bInit)
			{
				UG_LOG("Approximation Space not initialized.\n");
				return NULL;
			}

			dof_distribution_type* dofDistr = m_MGDoFManager.get_level_dof_distribution(level);
			if(dofDistr == NULL)
			{
				throw(UG_ERROR_DoFDistributionMissing());
			}

			function_type* gridFct = new function_type(*this, *dofDistr);
			return gridFct;
		}

		// create a new grid function of this approximation space
		function_type* create_surface_function()
		{
			if(!m_bInit)
			{
				UG_LOG("Approximation Space not initialized.\n");
				return NULL;
			}

			dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();
			if(dofDistr == NULL)
			{
				throw(UG_ERROR_DoFDistributionMissing());
			}

			function_type* gridFct = new function_type(*this, *dofDistr);

			return gridFct;
		}

		dof_distribution_type& get_surface_dof_distribution()
		{
			if(!m_bInit)
				throw(UG_ERROR_DoFDistributionMissing());

			dof_distribution_type* dofDistr = m_MGDoFManager.get_surface_dof_distribution();

			if(dofDistr == NULL)
				throw(UG_ERROR_DoFDistributionMissing());

			return *dofDistr;
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

		const SurfaceView* get_surface_view() const {return m_MGDoFManager.get_surface_view();}

		std::vector<const dof_distribution_type*> get_level_dof_distributions() const
		{
			return m_MGDoFManager.get_level_dof_distributions();
		}

		dof_distribution_type& get_level_dof_distribution(size_t level)
		{
			if(!m_bInit)
				throw(UG_ERROR_DoFDistributionMissing());

			return *(m_MGDoFManager.get_level_dof_distribution(level));
		}

		const dof_distribution_type& get_level_dof_distribution(size_t level) const
		{
			if(!m_bInit)
				throw(UG_ERROR_DoFDistributionMissing());

			return *(m_MGDoFManager.get_level_dof_distribution(level));
		}

	protected:
	//	Init flag
		bool m_bInit;

	// 	Dof manager used for this Approximation Space
		dof_manager_type m_MGDoFManager;
};


}


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__ */
