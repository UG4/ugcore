/*
 * mg_dof_manager_impl.h
 *
 *  Created on: 03.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER_IMPL__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER_IMPL__

#include "mg_dof_manager.h"
#include "lib_grid/algorithms/subset_util.h"

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
	if(m_pFunctionPattern != NULL)
		if(!assign_function_pattern(*m_pFunctionPattern))
			return false;

	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
assign_function_pattern(FunctionPattern& dp)
{
	m_pFunctionPattern = &dp;

	// if already subsethandler set
	if(m_pMGSubsetHandler != NULL)
	{
		// remember current levels
		size_t num_level = m_vLevelDoFDistribution.size();

		// free memory
		delete_distributions();

		// reallocate for new pattern
		if(num_level > 0)
			if(!level_distribution_required(num_level-1))
				return false;
	}
	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
distribute_dofs()
{
// 	no levels -> nothing to do
	if(num_levels() == 0) return true;

//	distribute level dofs
	if(!distribute_level_dofs()) return false;

// 	distribute surface dofs
	if(!distribute_surface_dofs()) return false;

//	we're done
	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
distribute_level_dofs()
{
//	Checks
	if(m_pMGSubsetHandler == NULL)
	{
		UG_LOG("No Subset Handler set to MultiGrid DoF Manager.\n");
		return false;
	}

	if(m_pFunctionPattern == NULL)
	{
		UG_LOG("No Function Pattern set to MultiGrid DoF Manager.\n");
		return false;
	}

// 	no levels -> nothing to do
	if(num_levels() == 0) return true;

// 	require distributions on all levels
	if(!level_distribution_required(num_levels()-1))
	{
		UG_LOG("Cannot access distribution of level.\n");
		return false;
	}

// 	distribute on level grids
	for(size_t l = 0; l < num_levels(); ++l)
	{
		if(!m_vLevelDoFDistribution[l]->distribute_dofs())
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
distribute_surface_dofs()
{
	if(m_pMGSubsetHandler == NULL)
	{
		UG_LOG("No Subset Handler set to MultiGrid DoF Manager.\n");
		return false;
	}
	if(m_pFunctionPattern == NULL)
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
	if(!m_pSurfaceDoFDistribution->distribute_dofs())
	{
		UG_LOG("Cannot distribute dofs on surface.\n");
		return false;
	}

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
	for(size_t l = 0; l < m_vLevelDoFDistribution.size(); ++l)
	{
		UG_LOG("    " << l << "    |");
		print_statistic(*m_vLevelDoFDistribution[l]);
		UG_LOG(std::endl);
	}

//	Write Infos for Surface Grid
	if(m_pSurfaceDoFDistribution != NULL)
	{
		UG_LOG("  surf   |");
		print_statistic(*m_pSurfaceDoFDistribution);
		UG_LOG(std::endl);
	}

}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_layout_statistic(const dof_distribution_type& dd, size_t ddlev) const
{
#ifdef UG_PARALLEL
//	Total number of DoFs
	UG_LOG(std::setw(8) << dd.num_master_dofs(ddlev) <<" | ");

	UG_LOG(std::setw(8) << dd.num_slave_dofs(ddlev) <<" | ");

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
	UG_LOG(" No Layouts in sequentiel code.\n");
#else
	UG_LOG("Layouts on Process " << 0 << ":\n");

//	Write header line
	UG_LOG(" Level | DD Level |  Master  |  Slave   | vert. Master | vert. Slave\n");
	UG_LOG("---------------------------------------------------------------------\n");

//	Write Infos for Levels
	for(size_t l = 0; l < m_vLevelDoFDistribution.size(); ++l)
	{
		for(size_t ddlev = 0; ddlev <
			m_vLevelDoFDistribution[l]->num_domain_decomposition_level(); ++ddlev)
		{
			if(ddlev == 0) {UG_LOG("   " << l << "   |");}
			else {UG_LOG("       |");}

			UG_LOG(std::setw(9) << ddlev << " | ");
			print_layout_statistic(*m_vLevelDoFDistribution[l], ddlev);
			UG_LOG("\n");
		}
		UG_LOG("---------------------------------------------------------------------\n");
	}

//	Write Infos for Surface Grid
	if(m_pSurfaceDoFDistribution != NULL)
	{
		for(size_t ddlev = 0; ddlev <
			m_pSurfaceDoFDistribution->num_domain_decomposition_level(); ++ddlev)
		{
			UG_LOG("  " << ddlev << "  | ");
			print_layout_statistic(*m_pSurfaceDoFDistribution, ddlev);
			UG_LOG(std::endl);
		}
	}
#endif
}



template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
surface_view_required()
{
// 	serial version

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
	CreateSurfaceView(*m_pSurfaceView, *m_pMultiGrid, *m_pMGSubsetHandler,
				m_pMultiGrid->vertices_begin(), m_pMultiGrid->vertices_end());
	CreateSurfaceView(*m_pSurfaceView, *m_pMultiGrid, *m_pMGSubsetHandler,
				m_pMultiGrid->edges_begin(), m_pMultiGrid->edges_end());
	CreateSurfaceView(*m_pSurfaceView, *m_pMultiGrid, *m_pMGSubsetHandler,
				m_pMultiGrid->faces_begin(), m_pMultiGrid->faces_end());
	CreateSurfaceView(*m_pSurfaceView, *m_pMultiGrid, *m_pMGSubsetHandler,
				m_pMultiGrid->volumes_begin(), m_pMultiGrid->volumes_end());

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
	if(m_pSurfaceDoFDistribution == NULL)
		m_pSurfaceDoFDistribution =
			new TDoFDistribution(m_pSurfaceView->get_geometric_object_collection(),
								 *m_pSurfaceView,
								 m_surfaceStorageManager,
								 *m_pFunctionPattern);

//	Check success
	if(m_pSurfaceDoFDistribution == NULL)
	{
		UG_LOG("Cannot allocate Surface DoF Distribution.\n");
		return false;
	}

	return true;
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
level_distribution_required(size_t level)
{
	if(level >= m_pMGSubsetHandler->num_levels())
	{
		UG_LOG("Level DoF Distribution required on Level " << level <<
			   ", but MGSubsetHandler has only " <<
			   m_pMGSubsetHandler->num_levels()<< ".\n");
		return false;
	}

// 	Create level dof distributions
	for(size_t l = m_vLevelDoFDistribution.size(); l <= level; ++l)
	{
		m_vLevelDoFDistribution.push_back(
				new TDoFDistribution(m_pMGSubsetHandler->get_goc_by_level(l),
									 *m_pMGSubsetHandler,
									 m_levelStorageManager,
									 *m_pFunctionPattern));

		if(m_vLevelDoFDistribution[l] == NULL)
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
delete_distributions()
{
// 	Delete surface dof distributions
	m_surfaceStorageManager.clear();
	if(m_pSurfaceDoFDistribution != NULL)
		delete m_pSurfaceDoFDistribution;
	m_pSurfaceDoFDistribution = NULL;

// 	delete surface view
	if(m_pSurfaceView != NULL)
		delete m_pSurfaceView;
	m_pSurfaceView = NULL;

// 	Delete level dof distributions
	m_levelStorageManager.clear();
	for(size_t l = 0; l < m_vLevelDoFDistribution.size(); ++l)
	{
		delete m_vLevelDoFDistribution[l];
		m_vLevelDoFDistribution[l] = NULL;
	}
	m_vLevelDoFDistribution.clear();
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER_IMPL__ */
