/*
 * parallel_dof_manager_impl.h
 *
 *  Created on: 03.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER_IMPL__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER_IMPL__

#include "parallel_dof_manager.h"

namespace ug{

template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
distribute_dofs()
{
// 	no levels -> nothing to do
	if(TMGDoFManager::num_levels() == 0) return true;

//	distribute level dofs
	if(!distribute_level_dofs()) return false;

// 	distribute surface dofs
	return distribute_surface_dofs();
}

template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
distribute_level_dofs()
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
	if(TMGDoFManager::num_levels() == 0) return true;

	for(size_t l = TMGDoFManager::num_levels() - 1; ; --l)
	{
	//	get dof distribution
		typename TMGDoFManager::dof_distribution_type& distr =
			*const_cast<typename TMGDoFManager::dof_distribution_type*>
				(TMGDoFManager::get_level_dof_distribution(l));

		if(domain_decomposition_enabled()){
		//	create process and subdomain layouts
			bRet &= CreateIndexLayouts_DomainDecomposition(
									  distr.get_master_layout(0),
									  distr.get_master_layout(1),
									  distr, *m_pLayoutMap, INT_MASTER,l,
									  m_cbProcIDToSubdomID);
			bRet &= CreateIndexLayouts_DomainDecomposition(
									distr.get_slave_layout(0),
									distr.get_slave_layout(1),
									distr, *m_pLayoutMap, INT_SLAVE,l,
									m_cbProcIDToSubdomID);
		}
		else{
		//	create index layouts
			bRet &= CreateIndexLayout(distr.get_master_layout(0),
									  distr, *m_pLayoutMap, INT_MASTER,l);
			bRet &= CreateIndexLayout(distr.get_slave_layout(0),
									  distr, *m_pLayoutMap, INT_SLAVE,l);
		}

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

template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
distribute_surface_dofs()
{
	if(!m_pLayoutMap){
		UG_LOG("  no layout map specified. aborting.\n");
		return false;
	}

	TMGDoFManager::distribute_surface_dofs();

	// TODO: Implement surface layouts

	return true;
}


template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
surface_view_required()
{
//	Check that Distr. Grid Manager has been set
	if(!m_pDistGridManager)
	{
		UG_LOG("No Distr. Grid Manager specified. aborting.\n");
		return false;
	}

// 	Parallel version

//	Create Surface View if not already created
	if(this->m_pSurfaceView == NULL)
		this->m_pSurfaceView = new SurfaceView(*this->m_pMultiGrid);

//	Check Success
	if(this->m_pSurfaceView == NULL)
	{
		UG_LOG("Allocation of Surface View failed.\n");
		return false;
	}

// 	Create surface view for all elements
	CreateSurfaceView(*(this->m_pSurfaceView), *m_pDistGridManager,
	                  *this->m_pMGSubsetHandler,
	                  this->m_pMultiGrid->vertices_begin(),
	                  this->m_pMultiGrid->vertices_end());
	CreateSurfaceView(*(this->m_pSurfaceView), *m_pDistGridManager,
	                  *this->m_pMGSubsetHandler,
	                  this->m_pMultiGrid->edges_begin(),
	                  this->m_pMultiGrid->edges_end());
	CreateSurfaceView(*(this->m_pSurfaceView), *m_pDistGridManager,
	                  *this->m_pMGSubsetHandler,
	                  this->m_pMultiGrid->faces_begin(),
	                  this->m_pMultiGrid->faces_end());
	CreateSurfaceView(*(this->m_pSurfaceView), *m_pDistGridManager,
	                  *this->m_pMGSubsetHandler,
	                  this->m_pMultiGrid->volumes_begin(),
	                  this->m_pMultiGrid->volumes_end());

// 	Set storage manager
	this->m_surfaceStorageManager.
			set_subset_handler(*(this->m_pSurfaceView));

//	return view
	return true;
}

template <typename TMGDoFManager>
void
ParallelMGDoFManager<TMGDoFManager>::
print_statistic(typename TMGDoFManager::dof_distribution_type& dd) const
{
//	Get Process communciator;
	pcl::ProcessCommunicator pCom = dd.get_process_communicator();

//	global and local values
	std::vector<int> tNumGlobal, tNumLocal;

//	write local dof numbers
	tNumLocal.push_back(dd.num_dofs() - dd.num_slave_dofs());
	for(int si = 0; si < dd.num_subsets(); ++si)
		tNumLocal.push_back(dd.num_dofs(si));

//	resize recieve array
	tNumGlobal.resize(tNumLocal.size());

//	sum up over processes
	if(!pCom.empty())
	{
		pCom.allreduce(&tNumLocal[0], &tNumGlobal[0], tNumGlobal.size(),
	               	   PCL_DT_INT, PCL_RO_SUM);
	}
	else if (pcl::GetNumProcesses() == 1)
	{
		for(size_t i = 0; i < tNumGlobal.size(); ++i)
			tNumGlobal[i] = tNumLocal[i];
	}
	else
	{
		UG_LOG(" Unable to compute informations.");
		return;
	}

//	Total number of DoFs
	UG_LOG(std::setw(10) << tNumGlobal[0] <<" | ");

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
		UG_LOG(std::setw(8) << tNumGlobal[si+1] << ") ");
	}
}

template <typename TMGDoFManager>
void
ParallelMGDoFManager<TMGDoFManager>::
print_statistic() const
{
//	Write info
	UG_LOG("DoFDistribution on all Processes (m= master, s=slave):\n");

//	Write header line
	UG_LOG("  Level  | Total (m) | BlockSize | "
			"(SubsetIndex (m+s), BlockSize, DoFs per Subset) \n");

//	Write Infos for Levels
	for(size_t l = 0; l < this->m_vLevelDoFDistribution.size(); ++l)
	{
		UG_LOG("    " << l << "    |");
		print_statistic(*this->m_vLevelDoFDistribution[l]);
		UG_LOG(std::endl);
	}

//	Write Infos for Surface Grid
	if(this->m_pSurfaceDoFDistribution != NULL)
	{
		UG_LOG("  surf   |");
		print_statistic(*this->m_pSurfaceDoFDistribution);
		UG_LOG(std::endl);
	}

	TMGDoFManager::print_statistic();
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER_IMPL__ */
