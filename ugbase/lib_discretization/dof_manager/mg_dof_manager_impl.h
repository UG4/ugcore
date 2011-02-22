/*
 * mg_dof_manager_impl.h
 *
 *  Created on: 03.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER_IMPL__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER_IMPL__

#include "mg_dof_manager.h"
#include "lib_grid/algorithms/multi_grid_util.h"

namespace ug{

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
assign_multi_grid_subset_handler(MultiGridSubsetHandler& mgsh)
{
// 	Remember SubsetHandler and MultiGrid
	m_pMGSubsetHandler = &mgsh;
	m_pMultiGrid = m_pMGSubsetHandler->get_assigned_multi_grid();

// 	Get StorageManager for levels
	m_levelStorageManager.set_subset_handler(mgsh);

// 	Set Function pattern if already assigned
	if(m_pFuncPattern != NULL)
		if(!assign_function_pattern(*m_pFuncPattern))
			return false;

	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
assign_function_pattern(FunctionPattern& dp)
{
	m_pFuncPattern = &dp;

//	 if already subsethandler set
	if(m_pMGSubsetHandler != NULL)
	{
	// 	remember current levels
		size_t num_level = m_vLevelDD.size();

	//	update level dofs
		if(level_dofs_enabled())
		{
		// free memory
			disable_level_dofs();

		// reallocate for new pattern
		if(!level_distribution_required(num_level))
			return false;
		}

	//	update surface dofs
		if(surface_dofs_enabled())
		{
		// free memory
			disable_surface_dofs();

		// reallocate for new pattern
			if(!surface_distribution_required())
				return false;
		}
	}

//	we're done
	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
enable_dofs()
{
//	distribute level dofs
	if(!enable_level_dofs()) return false;

// 	distribute surface dofs
	if(!enable_surface_dofs()) return false;

//	we're done
	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
enable_level_dofs()
{
//	Checks
	if(m_pMGSubsetHandler == NULL)
	{
		UG_LOG("No Subset Handler set to MultiGrid DoF Manager.\n");
		return false;
	}

	if(m_pFuncPattern == NULL)
	{
		UG_LOG("No Function Pattern set to MultiGrid DoF Manager.\n");
		return false;
	}

// 	require distributions on all levels
	if(!level_distribution_required(num_levels()))
	{
		UG_LOG("Cannot access distribution of level.\n");
		return false;
	}

// 	distribute on level grids
	for(size_t l = 0; l < num_levels(); ++l)
	{
		if(!m_vLevelDD[l]->distribute_dofs())
		{
			UG_LOG("Cannot distribute dofs on level "<<l<<".\n");
			return false;
		}
	}

	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
enable_surface_dofs()
{
	if(m_pMGSubsetHandler == NULL)
	{
		UG_LOG("No Subset Handler set to MultiGrid DoF Manager.\n");
		return false;
	}
	if(m_pFuncPattern == NULL)
	{
		UG_LOG("No Function Pattern set to MultiGrid DoF Manager.\n");
		return false;
	}

// 	update surface distribution
	if(!surface_distribution_required())
	{
		UG_LOG("Cannot update surface distribution.\n");
		return false;
	}

// 	distribute on surface grid
	if(!m_pSurfDD->distribute_dofs())
	{
		UG_LOG("Cannot distribute dofs on surface.\n");
		return false;
	}

//	we're done
	return true;
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_statistic(const dof_distribution_type& dd) const
{
//	Total number of DoFs
	UG_LOG(std::setw(10) << dd.num_dofs() <<" | ");

//	Overall block size
	if(dd.blocksize() != -1){
		UG_LOG(std::setw(8)  << dd.blocksize());}
	else{
		UG_LOG("variable");};
	UG_LOG("  | " );

//	Subset informations
	for(int si = 0; si < dd.num_subsets(); ++si)
	{
		UG_LOG( " (" << si << ",");
		UG_LOG(dd.blocksize(si) <<",");
		UG_LOG(std::setw(8) << dd.num_dofs(si) << ") ");
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_statistic() const
{
//	Write info
	UG_LOG("DoFDistribution");
#ifdef UG_PARALLEL
	UG_LOG(" on Process " << pcl::GetProcRank());
#endif
	UG_LOG(":\n");

//	Write header line
	UG_LOG("  Level  |   Total   | BlockSize | "
			"(SubsetIndex, BlockSize, DoFs per Subset) \n");

//	Write Infos for Levels
	for(size_t l = 0; l < m_vLevelDD.size(); ++l)
	{
		UG_LOG("    " << l << "    |");
		print_statistic(*m_vLevelDD[l]);
		UG_LOG(std::endl);
	}

//	Write Infos for Surface Grid
	if(m_pSurfDD != NULL)
	{
		UG_LOG("  surf   |");
		print_statistic(*m_pSurfDD);
		UG_LOG(std::endl);
	}

}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_layout_statistic(const dof_distribution_type& dd) const
{
#ifdef UG_PARALLEL
//	Total number of DoFs
	UG_LOG(std::setw(8) << dd.num_master_dofs() <<" | ");

	UG_LOG(std::setw(8) << dd.num_slave_dofs() <<" | ");

	UG_LOG(std::setw(12) << dd.num_vertical_master_dofs() <<" | ");

	UG_LOG(std::setw(12) << dd.num_vertical_slave_dofs());
#endif
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_layout_statistic() const
{
//	Write info
#ifndef UG_PARALLEL
	UG_LOG(" No Layouts in sequential code.\n");
#else
	UG_LOG("Layouts on Process " <<  pcl::GetOutputProcRank() << ":\n");

//	Write header line
	UG_LOG(" Level |  Master  |  Slave   | vert. Master | vert. Slave\n");
	UG_LOG("----------------------------------------------------------\n");

//	Write Infos for Levels
	for(size_t l = 0; l < m_vLevelDD.size(); ++l)
	{
		UG_LOG(" " << std::setw(5)<< l << " | ");
		print_layout_statistic(*m_vLevelDD[l]);
		UG_LOG("\n");
	}

//	Write Infos for Surface Grid
	if(m_pSurfDD != NULL)
	{
		UG_LOG("  surf | ");
		print_layout_statistic(*m_pSurfDD);
		UG_LOG(std::endl);
	}
#endif
}



template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
surface_view_required()
{
// 	serial version
	if(m_pSurfaceView != NULL)
		return true;

//	Create Surface view if not already created
	if(m_pSurfaceView == NULL)
		m_pSurfaceView = new SurfaceView(*m_pMultiGrid);

//	Check success
	if(m_pSurfaceView == NULL)
	{
		UG_LOG("Allocation of Surface View failed.\n");
		return false;
	}

// 	Create surface view for all elements
	CreateSurfaceView(*m_pSurfaceView, *m_pMultiGrid, *m_pMGSubsetHandler);

// 	set storage manager
	m_surfaceStorageManager.set_subset_handler(*m_pSurfaceView);

//	we're done
	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
surface_distribution_required()
{
//	Create surface view iff needed
	if(!surface_view_required())
		return false;

// 	Create surface dof distributions
	if(m_pSurfDD == NULL)
	{
	//	create dof distribution on surface
		m_pSurfDD =
			new TDoFDistribution(m_pSurfaceView->get_geometric_object_collection(),
								 *m_pSurfaceView,
								 m_surfaceStorageManager,
								 *m_pFuncPattern,
								 *m_pSurfaceView);
	}

//	Check success
	if(m_pSurfDD == NULL)
	{
		UG_LOG("Cannot allocate Surface DoF Distribution.\n");
		return false;
	}

	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
level_distribution_required(size_t numLevel)
{
	if(numLevel > m_pMGSubsetHandler->num_levels())
	{
		UG_LOG("Level DoF Distribution required for "<< numLevel << " Level"
			   ", but MGSubsetHandler has only " <<
			   m_pMGSubsetHandler->num_levels()<< ".\n");
		return false;
	}

// 	Create level dof distributions
	for(size_t l = m_vLevelDD.size(); l < numLevel; ++l)
	{
		m_vLevelDD.push_back(
				new TDoFDistribution(m_pMGSubsetHandler->get_goc_by_level(l),
									 *m_pMGSubsetHandler,
									 m_levelStorageManager,
									 *m_pFuncPattern));

		if(m_vLevelDD[l] == NULL)
		{
			UG_LOG("Cannot allocate Level DoF Distribution on Level " << l << ".\n");
			return false;
		}
	}
	return true;
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
disable_surface_dofs()
{
// 	Delete surface dof distributions
	m_surfaceStorageManager.clear();
	if(m_pSurfDD != NULL)
		delete m_pSurfDD;
	m_pSurfDD = NULL;

// 	delete surface view
	if(m_pSurfaceView != NULL)
		delete m_pSurfaceView;
	m_pSurfaceView = NULL;
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
disable_level_dofs()
{
// 	Delete level dof distributions
	m_levelStorageManager.clear();
	for(size_t l = 0; l < m_vLevelDD.size(); ++l)
	{
		delete m_vLevelDD[l];
		m_vLevelDD[l] = NULL;
	}
	m_vLevelDD.clear();
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER_IMPL__ */
