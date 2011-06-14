/*
 * parallel_dof_manager_impl.h
 *
 *  Created on: 03.12.2010
 *      Author: andreasvogel
 */


#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER_IMPL__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER_IMPL__

#include "parallel_dof_manager.h"
#include "pcl/pcl_communication_structs.h"
#include <time.h>

namespace ug{

template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
enable_dofs()
{
//	distribute level dofs
	if(!enable_level_dofs()) return false;

// 	distribute surface dofs
	return enable_surface_dofs();
}

template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
enable_level_dofs()
{
//	check that layout map has been set
	if(!m_pLayoutMap)
	{
		UG_LOG("  no layout map specified. aborting.\n");
		return false;
	}

// 	distribute dofs in sequential
	if(!TMGDoFManager::enable_dofs()) return false;

//	proc local number of level
	int numLevLocal = 0;

//	without SubsetHandler, we have no level information
	if(this->m_pMGSubsetHandler != NULL)
		numLevLocal = this->m_pMGSubsetHandler->num_levels();

//	storage for global number of levels
	int numLevGlobal;

//	\todo: This is MPI World, should only be a subgroup if not all
//			processes carry the grid.
	pcl::ProcessCommunicator pCom;

	pCom.allreduce(&numLevLocal, &numLevGlobal, 1,
									PCL_DT_INT, PCL_RO_MAX);

//	in addition create the index layouts
	return create_level_index_layouts(numLevGlobal);
}

template <typename TMGDoFManager>
template <class TElem>
bool
ParallelMGDoFManager<TMGDoFManager>::
create_level_index_layouts(size_t numGlobalLevels)
{
	typedef typename GridLayoutMap::Types<TElem>::Layout	TLayout;

	bool bRet = true;

//	in order to correctly build horizontal interfaces at vertical master
//	nodes, we have to exchange whether a horizontal interface connects
//	a vertical master with another vertical master, or no.
//	this data is then passed to CreateLevelIndexLayout, so that only necessary
//	interfaces will be build.
/*
	ComPol_InterfaceStatus<TLayout>
		comPol_MasterStatus(m_pDistGridManager, INT_V_SLAVE);
	ComPol_InterfaceStatus<TLayout>
		comPol_SlaveStatus(m_pDistGridManager, INT_V_SLAVE);
	pcl::ParallelCommunicator<TLayout> com;
	com.exchange_data(m_pDistGridManager->grid_layout_map(),
					  INT_H_SLAVE, INT_H_MASTER, comPol_MasterStatus);
	com.exchange_data(m_pDistGridManager->grid_layout_map(),
					  INT_H_MASTER, INT_H_SLAVE, comPol_SlaveStatus);
	com.communicate();
*/

	GridLayoutMap& layoutMap = m_pDistGridManager->grid_layout_map();

	for(size_t l = numGlobalLevels - 1; ; --l)
	{
	//	get dof distribution
		typename TMGDoFManager::dof_distribution_type& distr =
			*const_cast<typename TMGDoFManager::dof_distribution_type*>
				(TMGDoFManager::get_level_dof_distribution(l));

		UG_ASSERT(&distr != NULL, "No level dof distribution.");

		if(l < TMGDoFManager::num_levels())
		{
		//	create horizontal index layouts
			bRet &= AddEntriesToLevelIndexLayout(distr.get_master_layout(), distr,
						  	  	  layoutMap.get_layout<TElem>(INT_H_MASTER).
						  	  	  	  layout_on_level(l)/*,
						  	  	  &comPol_MasterStatus.get_result_map(l)*/);

			bRet &= AddEntriesToLevelIndexLayout(distr.get_slave_layout(), distr,
						  	  	  layoutMap.get_layout<TElem>(INT_H_SLAVE).
						  	  	  	  layout_on_level(l)/*,
						  	  	  &comPol_SlaveStatus.get_result_map(l)*/);

		//	create vertical layouts
			bRet &= AddEntriesToLevelIndexLayout(distr.get_vertical_master_layout(),
								distr, layoutMap.get_layout<TElem>(INT_V_MASTER).
						  	  				layout_on_level(l));

			bRet &= AddEntriesToLevelIndexLayout(distr.get_vertical_slave_layout(),
								distr, layoutMap.get_layout<TElem>(INT_V_SLAVE).
						  	  				layout_on_level(l));
		}
		else
		{
			distr.get_master_layout().clear();
			distr.get_slave_layout().clear();
			distr.get_vertical_master_layout().clear();
			distr.get_vertical_slave_layout().clear();
		}

		if(l==0) break;
	}

	return bRet;
}

template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
create_level_index_layouts(size_t numGlobalLevels)
{
//	TODO:	this communicator should be specified from the application
	pcl::ProcessCommunicator commWorld;

//	return flag
	bool bRet = true;

// 	if no cut level has appeared
	bool no_cut = true;

//	if no levels given, we're done
	if(numGlobalLevels == 0) return true;

	bRet &= create_level_index_layouts<VertexBase>(numGlobalLevels);
	bRet &= create_level_index_layouts<EdgeBase>(numGlobalLevels);
	bRet &= create_level_index_layouts<Face>(numGlobalLevels);
	bRet &= create_level_index_layouts<Volume>(numGlobalLevels);

	for(size_t l = numGlobalLevels - 1; ; --l)
	{
	//	get dof distribution
		typename TMGDoFManager::dof_distribution_type& distr =
			*const_cast<typename TMGDoFManager::dof_distribution_type*>
				(TMGDoFManager::get_level_dof_distribution(l));

		UG_ASSERT(&distr != NULL, "No level dof distribution.");

	//	create local process communicator
	//	The idea  of local processes is to exclude processes from
	//	e.g. norm computation that does not have a grid on a given
	//	level. If a grid exist, but all dofs are vertical slaves
	//	this process is also reguarded as empty for norm computations.
	//	In those cases the process votes false for the subcommunicator.

	//TODO: perform a more precise check
//		if(!distr.get_vertical_slave_layout().empty())
//			no_cut = false;

	// 	choose if this process participates
		bool participate = no_cut
						   && !commWorld.empty()
						   && (distr.num_dofs() > 0);

		UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
						  ": Participate = "<< participate <<
						  " for level " << l << " (num_dofs=" <<
						  distr.num_dofs() << ",no_cut=" << no_cut
						  << ",!empty=" << !commWorld.empty() << ").\n");

	//	create process communicator for interprocess layouts
		distr.get_process_communicator()
				= commWorld.create_sub_communicator(participate);

		if(l==0) break;
	}

//	we're done
	return bRet;
}

template <typename TMGDoFManager>
void
ParallelMGDoFManager<TMGDoFManager>::
defragment()
{
//	defragment dofs on each process
	TMGDoFManager::defragment();

//	proc local number of level
	int numLevLocal = 0;

//	without SubsetHandler, we have no level information
	if(this->m_pMGSubsetHandler != NULL)
		numLevLocal = this->m_pMGSubsetHandler->num_levels();

//	storage for global number of levels
	int numLevGlobal;

//	\todo: This is MPI World, should only be a subgroup if not all
//			processes carry the grid.
	pcl::ProcessCommunicator pCom;

	pCom.allreduce(&numLevLocal, &numLevGlobal, 1,
									PCL_DT_INT, PCL_RO_MAX);

//	request global number of levels
	if(this->level_dofs_enabled())
		this->level_distribution_required(numLevGlobal);

//	build up interfaces
	bool bRet = true;

	if(this->level_dofs_enabled())
	{
		bRet = create_level_index_layouts(numLevGlobal);
	}

	if(this->surface_dofs_enabled())
	{
		bRet &= create_surface_index_layouts();
	}

	if(!bRet) throw(UGFatalError("Cannot build interfaces.\n"));
}


template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
enable_surface_dofs()
{
	if(!m_pLayoutMap){
		UG_LOG("  no layout map specified. aborting.\n");
		return false;
	}

//	create dofs on each process
	TMGDoFManager::enable_surface_dofs();

//	build up interfaces
	return create_surface_index_layouts();
}


template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
create_surface_index_layouts()
{
//	get dof distribution
	typename TMGDoFManager::dof_distribution_type& distr =
		*const_cast<typename TMGDoFManager::dof_distribution_type*>
			(TMGDoFManager::get_surface_dof_distribution());

//	if no levels given, error
	if(TMGDoFManager::num_levels() == 0) return false;

//	return flag
	bool bRet = true;

//	create index layouts
	distr.get_master_layout().clear();
	distr.get_slave_layout().clear();

//	create surface index layouts
	bRet &= CreateSurfaceIndexLayout(distr.get_master_layout(),
	                                 distr, *m_pLayoutMap, INT_H_MASTER,
									 *this->m_pMultiGrid,
									 *m_pDistGridManager);
	bRet &= CreateSurfaceIndexLayout(distr.get_slave_layout(),
	                                 distr, *m_pLayoutMap, INT_H_SLAVE,
	                                 *this->m_pMultiGrid,
	                                 *m_pDistGridManager);
// 	INFO: THIS MUST NOT BE USED AND WILL SOON BE DELETED AFTER SOME TESTS
/*	bRet &= CreateSurfaceIndexLayout(distr.get_master_layout(),
	                                 distr, *m_pLayoutMap, INT_V_MASTER,
									 *this->m_pMultiGrid,
									 *m_pDistGridManager);
	bRet &= CreateSurfaceIndexLayout(distr.get_slave_layout(),
	                                 distr, *m_pLayoutMap, INT_V_SLAVE,
	                                 *this->m_pMultiGrid,
	                                 *m_pDistGridManager);
*/

//	check
//	pcl::ParallelCommunicator<IndexLayout> comTmp;
//	pcl::PrintLayout(comTmp, distr.get_master_layout(), distr.get_slave_layout());

//	on the surface view, there are no vertical layouts.
	distr.get_vertical_master_layout().clear();
	distr.get_vertical_slave_layout().clear();

//	TODO:	this communicator should be specified from the application
	pcl::ProcessCommunicator commWorld;

//	create communicator
//	\todo: Here we assume, that all procs are involved in the surface grid,
//			but this must not be true
	distr.get_process_communicator()
			= commWorld.create_sub_communicator(true);

//	we're done
	return bRet;
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
	if(this->m_pSurfaceView != NULL)
		return true;

//	Create Surface View if not already created
	if(this->m_pSurfaceView == NULL)
		this->m_pSurfaceView = new SurfaceView(*this->m_pMultiGrid);

//	Check Success
	if(this->m_pSurfaceView == NULL)
	{
		UG_LOG("Allocation of Surface View failed.\n");
		return false;
	}

//  \todo: Use m_pDistGridManager to take care of ghost elements when creating
//			the surface view (e.g. exclude vertical master/slave)
// 	Create surface view for all elements
	CreateSurfaceView(*this->m_pSurfaceView,
	                  *m_pDistGridManager,
	                  *this->m_pMGSubsetHandler);

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
//	Get Process communicator;
	pcl::ProcessCommunicator pCom = dd.get_process_communicator();

//	global and local values
	std::vector<int> tNumGlobal, tNumLocal;

//	write local dof numbers of all masters; this the number of all dofs
//	minus the number of all slave dofs (double counting of slaves can not
//	occur, since each slave is only slave of one master
	tNumLocal.push_back(dd.num_dofs() - dd.num_slave_dofs());

//	write local number of dof in a subset for all subsets. For simplicity, we
//	only communicate the total number of dofs (i.e. master+slave)
//	\todo: count slaves in subset and subtract them to get only masters
	for(int si = 0; si < dd.num_subsets(); ++si)
		tNumLocal.push_back(dd.num_dofs(si));

//	resize receive array
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
	for(size_t l = 0; l < this->m_vLevelDD.size(); ++l)
	{
		UG_LOG("  " << std::setw(5) << l << "  |");
		print_statistic(*this->m_vLevelDD[l]);
		UG_LOG(std::endl);
	}

//	Write Infos for Surface Grid
	if(this->m_pSurfDD != NULL)
	{
		UG_LOG("  surf   |");
		print_statistic(*this->m_pSurfDD);
		UG_LOG(std::endl);
	}

	TMGDoFManager::print_statistic();
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_DOF_MANAGER_IMPL__ */
