/*
 * parallel_dof_manager.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER__

#include "pcl/pcl.h"
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/parallelization/util/parallel_subset_util.h"
#include "parallelization_util.h"

namespace ug
{

template <typename TMGDoFManager>
class ParallelMGDoFManager : public TMGDoFManager
{
	public:
		ParallelMGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp)
		: TMGDoFManager(mgsh, dp), m_DistGridManager(NULL), m_pLayoutMap(NULL)
		{}

		ParallelMGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp, DistributedGridManager& distGridManager)
		: TMGDoFManager(mgsh, dp), m_DistGridManager(&distGridManager), m_pLayoutMap(&distGridManager.grid_layout_map())
		{}

		void set_distributed_grid_manager(DistributedGridManager& distGridManager)
		{
			m_DistGridManager = &distGridManager;
			m_pLayoutMap = &distGridManager.grid_layout_map();
		}

		//bool create_index_layouts()
		bool distribute_dofs()
		{
			if(!m_pLayoutMap)
				{UG_LOG("  no layout map specified. aborting.\n"); return false;}

			// distribute dofs in sequential
			if(!TMGDoFManager::distribute_dofs()) return false;

			//TODO:	this communicator should be specified from the application
			pcl::ProcessCommunicator commWorld;

			bool bRetVal = true;
			// if no cut level has appeared
			bool no_cut = true;
			for(size_t l = TMGDoFManager::num_levels() - 1; ; --l)
			{
				typename TMGDoFManager::dof_distribution_type& distr = *TMGDoFManager::get_level_dof_distribution(l);

				bRetVal &= CreateIndexLayout(distr.get_master_layout(), distr, *m_pLayoutMap, INT_MASTER,l);
				bRetVal &= CreateIndexLayout(distr.get_slave_layout(), distr, *m_pLayoutMap, INT_SLAVE,l);
				bRetVal &= CreateIndexLayout(distr.get_vertical_master_layout(), distr, *m_pLayoutMap, INT_VERTICAL_MASTER,l);
				bRetVal &= CreateIndexLayout(distr.get_vertical_slave_layout(), distr, *m_pLayoutMap, INT_VERTICAL_SLAVE,l);

			//	create local process communicator
			//	if a process has only vertical slaves, it is not involved in process communication.
			//TODO: perform a more precise check
				if(!distr.get_vertical_slave_layout().empty()) no_cut = false;
				bool participate = no_cut
								   && !commWorld.empty()
								   && (distr.num_dofs() > 0);

				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2, "  Says: Participate = "<< participate << " for level " << l << ""
															" (num_dofs = " << distr.num_dofs() << ").\n");

				distr.get_process_communicator() = commWorld.create_sub_communicator(participate);

				if(l==0) break;
			}

			return bRetVal;
		}

		bool distribute_surface_dofs()
		{
			if(!m_pLayoutMap){
				UG_LOG("  no layout map specified. aborting.\n");
				return false;
			}

			TMGDoFManager::distribute_surface_dofs();

			// TODO: Implement surface layouts
		}

	protected:
		virtual SurfaceView* create_surface_view()
		{
			// parallel version
			SurfaceView* pSurfaceView = new SurfaceView(*this->m_pMultiGrid);
			if(pSurfaceView == NULL)
				{UG_LOG("Allocation of Surface View failed.\n"); return false;}

			// Create surface view for all elements
			CreateSurfaceView(*pSurfaceView, *m_DistGridManager, *this->m_pMGSubsetHandler,
							this->m_pMultiGrid->vertices_begin(), this->m_pMultiGrid->vertices_end());
			CreateSurfaceView(*pSurfaceView, *m_DistGridManager, *this->m_pMGSubsetHandler,
							this->m_pMultiGrid->edges_begin(), this->m_pMultiGrid->edges_end());
			CreateSurfaceView(*pSurfaceView, *m_DistGridManager, *this->m_pMGSubsetHandler,
							this->m_pMultiGrid->faces_begin(), this->m_pMultiGrid->faces_end());
			CreateSurfaceView(*pSurfaceView, *m_DistGridManager, *this->m_pMGSubsetHandler,
							this->m_pMultiGrid->volumes_begin(), this->m_pMultiGrid->volumes_end());

			// set storage manager
			this->m_surfaceStorageManager.set_subset_handler(*pSurfaceView);

			return pSurfaceView;
		}

	private:
		// Distributed Grid Manager
		DistributedGridManager* m_DistGridManager;

		// Layout map of grid
		GridLayoutMap* m_pLayoutMap;
};

} // end namespace ug

#endif
