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
 *         for a locally/adaptivly refined MultiGrid the level grid solution is only
 *         living on a part of the domain)
 */
struct nonadaptive{};
struct adaptive{};

template <	typename TDomain,
			typename TAlgebra,
			template <class TDomain> class TDoFManager,
			typename TCategory = typename TDomain::categroy_type>
class ApproximationSpace;


// non adaptive Approximation space
template <typename TDomain, typename TAlgebra, template <class TDomain> class TDoFManager>
class ApproximationSpace<TDomain, TAlgebra, TDoFManager, nonadaptive>{
	private:
		// to make it more readable
		typedef ApproximationSpace<TDomain, TAlgebra, TDoFManager, nonadaptive> approximation_space_type;

	public:
		// DOMAIN
		// global coordinate type
		typedef TDomain domain_type;

		// global coordinate type
		typedef typename domain_type::position_type position_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::grid_type grid_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::subset_handler_type level_subset_handler_type;
		typedef typename domain_type::subset_handler_type surface_subset_handler_type;


		// GRID FUNCTION of this factory
		typedef TAlgebra algebra_type;

		// dof manager used
		#ifdef UG_PARALLEL
			typedef ParallelDoFManager<TDoFManager<level_subset_handler_type> >
				level_dof_manager_type;
			typedef ParallelDoFManager<TDoFManager<surface_subset_handler_type> >
				surface_dof_manager_type;
		#else
			typedef TDoFManager<level_subset_handler_type> level_dof_manager_type;
			typedef TDoFManager<surface_subset_handler_type> surface_dof_manager_type;
		#endif

		typedef GridFunction<approximation_space_type, level_dof_manager_type, algebra_type> level_function_type;
		typedef GridFunction<approximation_space_type, surface_dof_manager_type, algebra_type> surface_function_type;
		typedef surface_function_type function_type;

	public:
		ApproximationSpace(std::string name, domain_type& domain) :
			m_name(name), m_domain(domain),
			m_subset_handler(domain.get_subset_handler()), m_level_dof_manager(m_subset_handler)
		{
			#ifdef UG_PARALLEL
				m_level_dof_manager.set_grid_layout_map(
					domain.get_distributed_grid_manager()->
					grid_layout_map());
			#endif
		};


		/// add a single solution of LocalShapeFunctionSetID to selected subsets in dimension 'dim'
		/**
		 * \param[in] name			Name of this Single Solution
		 * \param[in] TrialSpace	Trial Space for this function
		 * \param[in] SubsetIndecex	Std::Vector of subset indeces, where this solution lives
		 */
		bool add_fundamental_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim, std::vector<int>& SubsetIndices)
		{
			if(m_level_dof_manager.add_discrete_function(name, id, SubsetIndices, dim) != true)
			{
				UG_LOG("ERROR in add_fundamental_discrete_function: cannot add discret function to level dof manager.\n");
				return false;
			}
			return true;
		}

		/// add a single solution of LocalShapeFunctionSetID to selected subsets for the world dimension
		bool add_fundamental_discrete_function(std::string name, LocalShapeFunctionSetID id, std::vector<int>& SubsetIndices)
		{
			const int dim = TDomain::dim;
			return this->add_fundamental_discrete_function(name, id, dim, SubsetIndices);
		}

		// add function on all domains for the world dimension
		bool add_fundamental_discrete_function(std::string name, LocalShapeFunctionSetID id)
		{
			const int dim = TDomain::dim;
			if(m_level_dof_manager.add_discrete_function(name, id, dim) != true)
			{
				UG_LOG("ERROR in add_fundamental_discrete_function: cannot add discret function to level dof manager.\n");
				return false;
			}
			return true;
		}

		// finish function pattern
		bool finalize()
		{
			if(m_level_dof_manager.finalize() != true)
			{
				UG_LOG("ERROR in finalize(): cannot finalize level dof manager.\n");
				return false;
			}
			return true;
		}

		// adaptive domain used or not ?
		inline bool is_adaptive()
		{
			return true;
		}

		// number of fundamental discrete functions
		inline uint num_fct() const
		{
			return m_level_dof_manager.num_fct();
		}

		// number of subsets
		inline int num_subsets() const
		{
			return m_subset_handler.num_subsets();
		}

		// type of fundamental function
		LocalShapeFunctionSetID get_local_shape_function_set_id(uint fct) const
		{
			return m_level_dof_manager.get_local_shape_function_set_id(fct);
		}

		// where is function defined
		inline bool fct_is_def_in_subset(uint fct, int subsetIndex) const
		{
			return m_level_dof_manager.fct_is_def_in_subset(fct, subsetIndex);
		}

		// number of levels
		uint num_levels()
		{
			return m_level_dof_manager.num_levels();
		}

		// return the domain
		const domain_type& get_domain() const
		{
			return m_domain;
		}

		// return the domain
		domain_type& get_domain()
		{
			return m_domain;
		}

		// create a new grid function of this approximation space
		level_function_type* create_level_grid_function(std::string name, uint level, bool allocate = true)
		{
			if(m_level_dof_manager.is_locked() != true)
			{
				std::cout << "Approximation Space has not been locked. Can not create a grid function of an unlocked space.\n";
				return NULL;
			}

			level_function_type* gridFct = new level_function_type(name, *this, m_level_dof_manager, level, allocate);
			return gridFct;
		}

		// create a new grid function of this approximation space
		surface_function_type* create_surface_grid_function(std::string name, bool allocate = true)
		{
			if(m_level_dof_manager.is_locked() != true)
			{
				std::cout << "Approximation Space has not been locked. Can not create a grid function of an unlocked space.\n";
				return NULL;
			}

			surface_function_type* gridFct = new surface_function_type(name, *this, m_level_dof_manager, m_level_dof_manager.num_levels() - 1, allocate);
			return gridFct;
		}

	protected:
		// name of this Approximation Space
		std::string m_name;

		// Domain, where solution lives
		domain_type& m_domain;

		// grid or multigrid or subsethandler, where elements are stored
		typename domain_type::subset_handler_type& m_subset_handler;

		// dof manager used for this Approximation Space
		level_dof_manager_type m_level_dof_manager;
};


}


#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION_SPACE__ */
