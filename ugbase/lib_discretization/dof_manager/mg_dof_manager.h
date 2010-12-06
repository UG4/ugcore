/*
 * mg_dof_manager.h
 *
 *  Created on: 12.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__

#include <vector>

#include "lib_grid/lib_grid.h"
#include "./function_pattern.h"
#include "lib_discretization/dof_manager/dof_distribution.h"

namespace ug{

/**
 * A MultiGridDoFManager handles the distribution of degrees of freedom on a
 * MultiGrid. It distributes the dof on each grid level and for the surface grid.
 * Thus, it creates num_level + 1 DoFDistributions
 *
 * \tparam 	TDoFDistribution	Type of DoF Distribution
 */
template <typename TDoFDistribution>
class MGDoFManager
{
	public:
	///	DoF Distribution type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	public:
	///	Default Constructor
		MGDoFManager()
			: m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL), m_pSurfaceView(NULL),
				m_pFunctionPattern(NULL), m_pSurfaceDoFDistribution(NULL)
		{
			m_vLevelDoFDistribution.clear();
		};

	///	Constructor setting Function Pattern and Multi Grid Subset Handler
		MGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp)
			: m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL), m_pSurfaceView(NULL),
				m_pFunctionPattern(NULL), m_pSurfaceDoFDistribution(NULL)
		{
			m_vLevelDoFDistribution.clear();
			assign_multi_grid_subset_handler(mgsh);
			assign_function_pattern(dp);
		};

	/// set multi grid subset handler
		bool assign_multi_grid_subset_handler(MultiGridSubsetHandler& mgsh);

	/// set function pattern
		bool assign_function_pattern(FunctionPattern& dp);

	/// number of levels
		inline size_t num_levels() const
		{
		//	without SubsetHandler, we have no level information
			if(m_pMGSubsetHandler == NULL) return 0;

		//	forward request
			return m_pMGSubsetHandler->num_levels();
		}

	/// distribute dofs on all levels + surface level
		bool distribute_dofs();

	/// distribute dofs on all levels
		bool distribute_level_dofs();

	///	distribute dofs on surface grid
		bool distribute_surface_dofs();

	///	returns Surface DoF Distribution
		const dof_distribution_type* get_surface_dof_distribution() const
		{
			//return m_pSurfaceDoFDistribution;

			if(num_levels() == m_vLevelDoFDistribution.size())
				return m_vLevelDoFDistribution[num_levels()-1];
			else
				return NULL;
		}

	///	print a statistic on dof distribution
		void print_statistic() const;

	///	returns Level DoF Distribution
		const dof_distribution_type* get_level_dof_distribution(size_t level) const
		{
			if(level < m_vLevelDoFDistribution.size())
				return m_vLevelDoFDistribution[level];
			else
				return NULL;
		}

	///	Destructor
		virtual ~MGDoFManager()
		{
			m_levelStorageManager.clear_subset_handler();
			m_surfaceStorageManager.clear_subset_handler();

			delete_distributions();
		}

	protected:
	///	creates the surface view
		virtual bool surface_view_required();

	///	creates the surface distribution
		bool surface_distribution_required();

	///	creates level DoF Distributions iff needed
		bool level_distribution_required(size_t level);

	///	deletes all distributions
		void delete_distributions();

	/// print statistic for a DoFDistribution
		void print_statistic(const dof_distribution_type& dd) const;

	protected:
	// 	MultiGridSubsetHandler this DofManager works on
		MultiGridSubsetHandler* m_pMGSubsetHandler;

	// 	MultiGrid associated to the SubsetHandler
		MultiGrid* m_pMultiGrid;

	// 	Surface View
		SurfaceView* m_pSurfaceView;

	// 	DoF Pattern
		FunctionPattern* m_pFunctionPattern;

	// 	Level DoF Distributors
		std::vector<dof_distribution_type*> m_vLevelDoFDistribution;

	// 	Surface Grid DoF Distributor
		dof_distribution_type* m_pSurfaceDoFDistribution;

	// 	Storage manager
		typename TDoFDistribution::storage_manager_type	m_levelStorageManager;
		typename TDoFDistribution::storage_manager_type	m_surfaceStorageManager;
};

} // end namespace ug

// include implementation
#include "mg_dof_manager_impl.h"

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__ */
