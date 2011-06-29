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
	bool bRet = create_level_index_layouts(numLevGlobal);

//	return if all procs succeeded
	return pcl::AllProcsTrue(bRet, pCom);
}

template <typename TMGDoFManager>
template <class TElem>
bool
ParallelMGDoFManager<TMGDoFManager>::
create_level_index_layouts(serial_dd_type& dd, size_t lev)
{
//	type of layout for element type
	typedef typename GridLayoutMap::Types<TElem>::Layout TLayout;

//	get the grid layout map
	GridLayoutMap& layoutMap = m_pDistGridManager->grid_layout_map();

//	set the return flag to true (initially)
	bool bRet = true;

//	only in case that this proc really has a grid on the requested level
//	we have to build up the interfaces
	if(lev < TMGDoFManager::num_levels())
	{
	//	create horizontal index layouts
		bRet &= AddEntriesToLevelIndexLayout(dd.get_master_layout(), dd,
							  layoutMap.get_layout<TElem>(INT_H_MASTER).
								  layout_on_level(lev));

		bRet &= AddEntriesToLevelIndexLayout(dd.get_slave_layout(), dd,
							  layoutMap.get_layout<TElem>(INT_H_SLAVE).
								  layout_on_level(lev));

	//	create vertical layouts
		bRet &= AddEntriesToLevelIndexLayout(dd.get_vertical_master_layout(),
							dd, layoutMap.get_layout<TElem>(INT_V_MASTER).
										layout_on_level(lev));

		bRet &= AddEntriesToLevelIndexLayout(dd.get_vertical_slave_layout(),
							dd, layoutMap.get_layout<TElem>(INT_V_SLAVE).
										layout_on_level(lev));
	}

//	return the success
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

//	loop all levels to create the communicator
	for(size_t l = 0; l < numGlobalLevels; ++l)
	{
	//	get serial dof distribution
		serial_dd_type* pDD = TMGDoFManager::get_level_dof_distribution(l);

	//	check serial DoF Distribution
		if(pDD == NULL)
		{
			UG_LOG("ERROR in 'ParallelMGDoFManager::create_level_index_layouts':"
					" Level DoF Distribution is missing on level "<<l<<".\n");
			bRet = false;
			continue;
		}

	//	get a reference for convenience
		serial_dd_type& dd = *pDD;

	//  -----------------------------------
	//	CREATE PROCESS COMMUNICATOR
	//  -----------------------------------
	//	The idea  of local processes is to exclude processes from
	//	e.g. norm computation that does not have a grid on a given
	//	level. If no DoFs exist on the level, that level is excluded from
	//	norm computations. In those cases the process votes false for
	//	the subcommunicator.

	// 	choose if this process participates
		bool participate = !commWorld.empty() && (dd.num_dofs() > 0);

		UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
						  ": Participate = "<< participate <<
						  " for level "<<l<<" (num_dofs="<<dd.num_dofs()<<
						  ",!empty=" << !commWorld.empty() << ").\n");

	//	create process communicator for interprocess layouts
		dd.get_process_communicator()
				= commWorld.create_sub_communicator(participate);

	//  -----------------------------------
	//	CREATE INDEX LAYOUTS ON LEVEL
	//  -----------------------------------
	//	clear all Index Layouts
		dd.get_master_layout().clear();
		dd.get_slave_layout().clear();
		dd.get_vertical_master_layout().clear();
		dd.get_vertical_slave_layout().clear();

	//	create the index layouts
		bRet &= create_level_index_layouts<VertexBase>(dd, l);
		bRet &= create_level_index_layouts<EdgeBase>(dd, l);
		bRet &= create_level_index_layouts<Face>(dd, l);
		bRet &= create_level_index_layouts<Volume>(dd, l);
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
	bool bRet = create_surface_index_layouts();

//	\todo: This is MPI World, should only be a subgroup if not all
//			processes carry the grid.
	pcl::ProcessCommunicator pCom;

//	return if all procs succeeded
	return pcl::AllProcsTrue(bRet, pCom);
}


template <typename TMGDoFManager>
bool
ParallelMGDoFManager<TMGDoFManager>::
create_surface_index_layouts()
{
//	get serial dof distribution
	serial_dd_type* pDD = TMGDoFManager::get_surface_dof_distribution();

//	check serial DoF Distribution
	if(pDD == NULL)
	{
		UG_LOG("ERROR in 'ParallelMGDoFManager::create_level_index_layouts':"
				" Surface DoF Distribution is missing.\n");
		return false;
	}

//	get a reference for convenience
	serial_dd_type& dd = *pDD;

//	return flag
	bool bRet = true;

//	create index layouts
	dd.get_master_layout().clear();
	dd.get_slave_layout().clear();

//	create surface index layouts
	bRet &= CreateSurfaceIndexLayout(dd.get_master_layout(),
	                                 dd, *m_pLayoutMap, INT_H_MASTER,
									 *this->m_pMultiGrid,
									 *m_pDistGridManager);
	bRet &= CreateSurfaceIndexLayout(dd.get_slave_layout(),
	                                 dd, *m_pLayoutMap, INT_H_SLAVE,
	                                 *this->m_pMultiGrid,
	                                 *m_pDistGridManager);

//	check
//	pcl::ParallelCommunicator<IndexLayout> comTmp;
//	pcl::PrintLayout(comTmp, distr.get_master_layout(), distr.get_slave_layout());

//	on the surface view, there are no vertical layouts.
	dd.get_vertical_master_layout().clear();
	dd.get_vertical_slave_layout().clear();

//	TODO:	this communicator should be specified from the application
	pcl::ProcessCommunicator commWorld;

//	only those procs are included, that have DoFs on the surface level
	bool participate = true;
	if(dd.num_dofs() == 0) participate = false;

//	create communicator
	dd.get_process_communicator() = commWorld.create_sub_communicator(participate);

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

//	compute local dof numbers of all masters; this the number of all dofs
//	minus the number of all slave dofs (double counting of slaves can not
//	occur, since each slave is only slave of one master), minus the number of
//	all vertical master dofs (double counting can occure, since a vertical master
//	can be in a horizontal slave interface.
//	Therefore: 	if there are no vert. master dofs, we can simply calculate the
//				number of DoFs by counting. If there are vert. master dofs
//				we need to remove doubles when counting.
	int numMasterDoF;
	if(dd.num_vertical_master_dofs() > 0)
	{
	//	create vector of vertical masters and horiz. slaves, since those are
	//	not reguarded as masters on the dd. All others are masters. (Especially,
	//	also vertical slaves are masters w.r.t. computing like smoothing !!)
		std::vector<IndexLayout::Element> vIndex;
		CollectElements(vIndex, dd.get_vertical_master_layout());
		CollectElements(vIndex, dd.get_slave_layout(), false);

	//	sort vector and remove doubles
		std::sort(vIndex.begin(), vIndex.end());
		vIndex.erase(std::unique(vIndex.begin(), vIndex.end()),
		                         vIndex.end());

	//	now remove all horizontal masters, since they are still master though
	//	in the vert master interface
		std::vector<IndexLayout::Element> vIndex2;
		CollectElements(vIndex2, dd.get_master_layout());
		for(size_t i = 0; i < vIndex2.size(); ++i)
			vIndex.erase(std::remove(vIndex.begin(), vIndex.end(), vIndex2[i]),
			             vIndex.end());

	//	the remaining DoFs are the "slaves"
		numMasterDoF = dd.num_dofs() - vIndex.size();
	}
	else
	{
	//	easy case: only subtract masters from slaves
		numMasterDoF = dd.num_dofs() - dd.num_slave_dofs();
	}

//	global and local values
	std::vector<int> tNumGlobal, tNumLocal;

//	write number of Masters on this process
	tNumLocal.push_back(numMasterDoF);

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
