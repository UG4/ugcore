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
		ParallelMGDoFManager()
		: TMGDoFManager(), m_pDistGridManager(NULL), m_pLayoutMap(NULL)
		{}

		ParallelMGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp)
		: TMGDoFManager(mgsh, dp), m_pDistGridManager(NULL), m_pLayoutMap(NULL)
		{}

		ParallelMGDoFManager(MultiGridSubsetHandler& mgsh, FunctionPattern& dp,
		                     DistributedGridManager& distGridManager)
		: TMGDoFManager(mgsh, dp), m_pDistGridManager(&distGridManager),
		  m_pLayoutMap(&distGridManager.grid_layout_map())
		{}

		void set_distributed_grid_manager(DistributedGridManager& distGridManager)
		{
			m_pDistGridManager = &distGridManager;
			m_pLayoutMap = &distGridManager.grid_layout_map();
		}

		//bool create_index_layouts()
		bool distribute_dofs()
		{
		//	check that layout map has been set
			if(!m_pLayoutMap)
			{
				UG_LOG("  no layout map specified. aborting.\n");
				return false;
			}

		// 	Distribute dofs in sequential
			if(!TMGDoFManager::distribute_dofs()) return false;

		//	TODO:	this communicator should be specified from the application
			pcl::ProcessCommunicator commWorld;

		//	return flag
			bool bRet = true;

		// 	if no cut level has appeared
			bool no_cut = true;

		//	if no levels given, we're done
			if(TMGDoFManager::num_levels() > 0) return true;

			for(size_t l = TMGDoFManager::num_levels() - 1; ; --l)
			{
			//	get dof distribution
				typename TMGDoFManager::dof_distribution_type& distr =
					*const_cast<typename TMGDoFManager::dof_distribution_type*>
						(TMGDoFManager::get_level_dof_distribution(l));

			//	create index layouts
				bRet &= CreateIndexLayout(distr.get_master_layout(),
										  distr, *m_pLayoutMap, INT_MASTER,l);
				bRet &= CreateIndexLayout(distr.get_slave_layout(),
										  distr, *m_pLayoutMap, INT_SLAVE,l);
				bRet &= CreateIndexLayout(distr.get_vertical_master_layout(),
										  distr, *m_pLayoutMap, INT_VERTICAL_MASTER,l);
				bRet &= CreateIndexLayout(distr.get_vertical_slave_layout(),
										  distr, *m_pLayoutMap, INT_VERTICAL_SLAVE,l);

			//	create local process communicator
			//	The idea  of local processes is to exclude processes form
			//	e.g. norm computation that does not have a grid on a given
			//	level. If a grid exist, but all dofs are vertical slaces
			//	this process is also reguarded as empty for norm computations.
			//	In those cases the process votes false for the subcommunicator.

			//TODO: perform a more precise check
				if(!distr.get_vertical_slave_layout().empty())
					no_cut = false;

			// 	choose if this process participates
				bool participate = no_cut
								   && !commWorld.empty()
								   && (distr.num_dofs() > 0);

				UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
								  ": Participate = "<< participate <<
								  " for level " << l << " (num_dofs=" <<
								  distr.num_dofs() << ",no_cut=" << no_cut
								  << ",!empty=" << !commWorld.empty() << ").\n");

			//	create process communicator
				distr.get_process_communicator()
						= commWorld.create_sub_communicator(participate);

				if(l==0) break;
			}

		//	we're done
			return bRet;
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
		//	Check that Distr. Grid Manager has been set
			if(!m_pDistGridManager)
			{
				UG_LOG("No Distr. Grid Manager specified. aborting.\n");
				return false;
			}

		// 	Parallel version
			SurfaceView* pSurfaceView = new SurfaceView(*this->m_pMultiGrid);
			if(pSurfaceView == NULL)
			{
				UG_LOG("Allocation of Surface View failed.\n");
				return false;
			}

		// 	Create surface view for all elements
			CreateSurfaceView(*pSurfaceView, *m_pDistGridManager,
			                  *this->m_pMGSubsetHandler,
			                  this->m_pMultiGrid->vertices_begin(),
			                  this->m_pMultiGrid->vertices_end());
			CreateSurfaceView(*pSurfaceView, *m_pDistGridManager,
			                  *this->m_pMGSubsetHandler,
			                  this->m_pMultiGrid->edges_begin(),
			                  this->m_pMultiGrid->edges_end());
			CreateSurfaceView(*pSurfaceView, *m_pDistGridManager,
			                  *this->m_pMGSubsetHandler,
			                  this->m_pMultiGrid->faces_begin(),
			                  this->m_pMultiGrid->faces_end());
			CreateSurfaceView(*pSurfaceView, *m_pDistGridManager,
			                  *this->m_pMGSubsetHandler,
			                  this->m_pMultiGrid->volumes_begin(),
			                  this->m_pMultiGrid->volumes_end());

		// 	Set storage manager
			this->m_surfaceStorageManager.set_subset_handler(*pSurfaceView);

		//	return view
			return pSurfaceView;
		}

	private:
	/// Distributed Grid Manager
		DistributedGridManager* m_pDistGridManager;

	/// Layout map of grid
		GridLayoutMap* m_pLayoutMap;
};

} // end namespace ug

#endif
