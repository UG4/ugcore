/*
 * parallel_dof_manager.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel
 */

// TODO: remove 'Callback_ProcessIDToSubdomainID' stuff ... (27012011)

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER__

#include "pcl/pcl.h"
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/parallelization/util/parallel_subset_util.h"
#include "parallelization_util.h"

namespace ug
{

/**
 * A ParallelMGDoFManager is a wrapper class to use a MGDoFManager in parallel.
 * It reimplements the function that need special care in parallel. Especially,
 * IndexLayouts are created during dof distribution.
 */
template <typename TMGDoFManager>
class ParallelMGDoFManager : public TMGDoFManager
{
	public:
	///	Default Constructor
		ParallelMGDoFManager()
		: TMGDoFManager(), m_pDistGridManager(NULL), m_pLayoutMap(NULL)
		{}

	///	Constructor setting MultiGrid and FunctionPattern
		ParallelMGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp)
		: TMGDoFManager(mgsh, dp), m_pDistGridManager(NULL), m_pLayoutMap(NULL)
		{}

	///	Constructor setting MultiGrid, FunctionPattern and DistributedGridManager
		ParallelMGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp,
		                     DistributedGridManager& distGridManager)
		: TMGDoFManager(mgsh, dp), m_pDistGridManager(&distGridManager),
		  m_pLayoutMap(&distGridManager.grid_layout_map())
		{}

	///	assign Distributed Grid Manager
		void set_distributed_grid_manager(DistributedGridManager& distGridManager)
		{
			m_pDistGridManager = &distGridManager;
			m_pLayoutMap = &distGridManager.grid_layout_map();
		}

	///	distribute dofs on levels and surface
		bool enable_dofs();

	///	distribute dofs on levels
		bool enable_level_dofs();

	///	distribute dofs on surface
		bool enable_surface_dofs();

	///	print a statistic on dof distribution
		void print_statistic() const;

	///	defragments the index set
		void defragment();

	protected:
	///	create the layouts on all levels
		bool create_level_index_layouts(size_t numGlobalLevels);

	///	a helper, which will create the layouts for a given geometric object
		template <class TElem> bool
		create_level_index_layouts(size_t numGlobalLevels);

	///	create the layouts on the surface level
		bool create_surface_index_layouts();

	///	creates the surface view iff needed
		virtual bool surface_view_required();

	/// print statistic for a DoFDistribution
		void print_statistic(typename TMGDoFManager::dof_distribution_type& dd) const;

	private:
	/// Distributed Grid Manager
		DistributedGridManager* m_pDistGridManager;

	/// Layout map of grid
		GridLayoutMap* m_pLayoutMap;
};

} // end namespace ug

// include implementation
#include "parallel_dof_manager_impl.h"

#endif
