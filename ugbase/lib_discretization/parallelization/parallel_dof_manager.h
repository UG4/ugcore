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
		: TMGDoFManager(), m_pDistGridManager(NULL), m_pLayoutMap(NULL),
		  m_bDomainDecompositionEnabled(false)
		{}

	///	Constructor setting MultiGrid and FunctionPattern
		ParallelMGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp)
		: TMGDoFManager(mgsh, dp), m_pDistGridManager(NULL), m_pLayoutMap(NULL),
		  m_bDomainDecompositionEnabled(false)
		{}

	///	Constructor setting MultiGrid, FunctionPattern and DistributedGridManager
		ParallelMGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp,
		                     DistributedGridManager& distGridManager)
		: TMGDoFManager(mgsh, dp), m_pDistGridManager(&distGridManager),
		  m_pLayoutMap(&distGridManager.grid_layout_map()),
		  m_bDomainDecompositionEnabled(false)
		{}

	///	assign Distributed Grid Manager
		void set_distributed_grid_manager(DistributedGridManager& distGridManager)
		{
			m_pDistGridManager = &distGridManager;
			m_pLayoutMap = &distGridManager.grid_layout_map();
		}

	///	distribute dofs on levels and surface
		bool distribute_dofs();

	///	distribute dofs on levels
		bool distribute_level_dofs();

	///	distribute dofs on surface
		bool distribute_surface_dofs();

	///	print a statistic on dof distribution
		void print_statistic() const;

	///	if enabled parallel interfaces are build with domain decomposition in mind.
		void enable_domain_decomposition(pcl::IDomainDecompositionInfo& ddInfo) /*(Callback_ProcessIDToSubdomainID cb_ProcIDToSubdomID)*/
		{
			m_bDomainDecompositionEnabled = true;
			//m_cbProcIDToSubdomID = cb_ProcIDToSubdomID;
			m_pDDInfo = &ddInfo;
		}

		bool domain_decomposition_enabled()				{return m_bDomainDecompositionEnabled;}

	protected:
	///	creates the surface view iff needed
		virtual bool surface_view_required();

	/// print statistic for a DoFDistribution
		void print_statistic(typename TMGDoFManager::dof_distribution_type& dd) const;

	private:
	/// Distributed Grid Manager
		DistributedGridManager* m_pDistGridManager;

	/// Layout map of grid
		GridLayoutMap* m_pLayoutMap;

		// pointer to Domain decomposition info object
		pcl::IDomainDecompositionInfo* m_pDDInfo; // instead of: Callback_ProcessIDToSubdomainID m_cbProcIDToSubdomID;
		bool m_bDomainDecompositionEnabled;
};

} // end namespace ug

// include implementation
#include "parallel_dof_manager_impl.h"

#endif
