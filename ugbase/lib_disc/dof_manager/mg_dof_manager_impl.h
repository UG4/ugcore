/*
 * mg_dof_manager_impl.h
 *
 *  Created on: 03.12.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_MANAGER_IMPL__
#define __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_MANAGER_IMPL__

#include "mg_dof_manager.h"
#include "lib_grid/algorithms/multi_grid_util.h"

namespace ug{

template <typename TDoFDistribution>
void MGDoFManager<TDoFDistribution>::enable_indices()
{
//	distribute level dofs
	try{
		enable_level_indices();
	}
	UG_CATCH_THROW("Cannot enable Level Indices.");

// 	distribute surface dofs
	try{
		enable_surface_indices();
	}
	UG_CATCH_THROW("Cannot enable Surface Indices.");
}

template <typename TDoFDistribution>
void MGDoFManager<TDoFDistribution>::enable_level_indices()
{
//	Checks
	if(m_pMGSubsetHandler == NULL)
		UG_THROW_FATAL("No Subset Handler set to MultiGrid DoF Manager.");

	if(m_pFuncPattern == NULL)
		UG_THROW_FATAL("No Function Pattern set to MultiGrid DoF Manager.");

//	check, that all geom objects are assigned to a subset
	if(	m_pMGSubsetHandler->num<VertexBase>() != m_pMultiGrid->num<VertexBase>())
		UG_THROW_FATAL("All Vertices "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_pMGSubsetHandler->num<EdgeBase>() != m_pMultiGrid->num<EdgeBase>())
		UG_THROW_FATAL("All Edges "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_pMGSubsetHandler->num<Face>() != m_pMultiGrid->num<Face>())
		UG_THROW_FATAL("All Faces "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_pMGSubsetHandler->num<Volume>() != m_pMultiGrid->num<Volume>())
		UG_THROW_FATAL("All Volumes "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

// 	require distributions on all levels
	try{
		level_distribution_required(num_levels());
	}
	UG_CATCH_THROW("Cannot access distribution of level.");

// 	distribute on level grids
	for(size_t l = 0; l < num_levels(); ++l)
	{
		if(!m_vLevelDD[l]->distribute_indices())
			UG_THROW_FATAL("Cannot distribute dofs on level "<<l);
	}

//	register DoFManager as observer
	try{
	m_pMultiGrid->register_observer(this,	OT_GRID_OBSERVER | OT_VERTEX_OBSERVER |
	                       	   			    OT_EDGE_OBSERVER | OT_FACE_OBSERVER |
	                       	   		        OT_VOLUME_OBSERVER);
	}
	UG_CATCH_THROW("Cannot register grid observer.");
}

template <typename TDoFDistribution>
void MGDoFManager<TDoFDistribution>::enable_surface_indices()
{
	if(m_pMGSubsetHandler == NULL)
		UG_THROW_FATAL("No Subset Handler set to MultiGrid DoF Manager.");

	if(m_pFuncPattern == NULL)
		UG_THROW_FATAL("No Function Pattern set to MultiGrid DoF Manager.");

//	check, that all geom objects are assigned to a subset
	if(	m_pMGSubsetHandler->num<VertexBase>() != m_pMultiGrid->num<VertexBase>())
		UG_THROW_FATAL("ERROR in 'MGDoFManager::enable_level_indices': All Vertices "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_pMGSubsetHandler->num<EdgeBase>() != m_pMultiGrid->num<EdgeBase>())
		UG_THROW_FATAL("ERROR in 'MGDoFManager::enable_level_indices': All Edges "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_pMGSubsetHandler->num<Face>() != m_pMultiGrid->num<Face>())
		UG_THROW_FATAL("ERROR in 'MGDoFManager::enable_level_indices': All Faces "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

	if(	m_pMGSubsetHandler->num<Volume>() != m_pMultiGrid->num<Volume>())
		UG_THROW_FATAL("ERROR in 'MGDoFManager::enable_level_indices': All Volumes "
			   " must be assigned to a subset. The passed subset handler "
			   " contains non-assigned elements, thus the dof distribution"
			   " is not possible, aborting.");

// 	update surface distribution
	try{
		surface_distribution_required();
	}
	UG_CATCH_THROW("Cannot update surface distribution.");

// 	distribute on surface grid
	if(!m_pSurfDD->distribute_indices())
		UG_THROW_FATAL("Cannot distribute dofs on surface.");

	try{
	//	unregister DoFManager
		m_pMultiGrid->unregister_observer(this);

	//	register SurfaceView first
		m_pMultiGrid->register_observer(m_pSurfaceView,
											  OT_GRID_OBSERVER | OT_VERTEX_OBSERVER |
											  OT_EDGE_OBSERVER | OT_FACE_OBSERVER |
											  OT_VOLUME_OBSERVER);

	//	register DoFManager again
		m_pMultiGrid->register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER |
											  OT_EDGE_OBSERVER | OT_FACE_OBSERVER |
											  OT_VOLUME_OBSERVER);
	}
	UG_CATCH_THROW("Cannot register observer.");
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_statistic(const dof_distribution_type& dd, int verboseLev) const
{
//	Total number of DoFs
	//UG_LOG(std::setw(10) << dd.num_indices() <<" | "); // orig
	UG_LOG(std::setw(10) << ConvertNumber(dd.num_indices(),10,6) << " | ");

//	Overall block size
	if(dd.blocksize() != -1){
		UG_LOG(std::setw(8)  << dd.blocksize());}
	else{
		UG_LOG("variable");};
	UG_LOG("  | " );

//	Subset informations
	if(verboseLev >= 1)
	{
		for(int si = 0; si < dd.num_subsets(); ++si)
		{
			UG_LOG( " (" << dd.subset_name(si) << ",");
			UG_LOG(dd.blocksize(si) <<",");
			//UG_LOG(std::setw(8) << dd.num_indices(si) << ") "); // orig
			UG_LOG(std::setw(8) << ConvertNumber(dd.num_indices(si),8,4) << ") ");

		}
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_statistic(int verboseLev) const
{
//	Write info
	UG_LOG("DoFDistribution");
#ifdef UG_PARALLEL
	UG_LOG(" on Process " << pcl::GetProcRank());
#endif
	UG_LOG(":\n");

//	Write header line
	UG_LOG("  Level  |   Total   | BlockSize | ");
	if(verboseLev >= 1) UG_LOG("(Subset, BlockSize, DoFs) ");
	UG_LOG("\n");

//	Write Infos for Levels
	for(size_t l = 0; l < m_vLevelDD.size(); ++l)
	{
		UG_LOG("  " << std::setw(5) << l << "  |");
		print_statistic(*m_vLevelDD[l], verboseLev);
		UG_LOG(std::endl);
	}

//	Write Infos for Surface Grid
	if(m_pSurfDD != NULL)
	{
		UG_LOG("  surf   |");
		print_statistic(*m_pSurfDD, verboseLev);
		UG_LOG(std::endl);
	}

}


template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_local_dof_statistic(const dof_distribution_type& dd, int verboseLev) const
{
//	Subset informations
	UG_LOG(dd.num_subsets() << " Subset(s) used (Subset Name, dim): ");
	for(int si = 0; si < dd.num_subsets(); ++si)
	{
		if(si > 0) UG_LOG(", ");
		UG_LOG("(" << dd.subset_name(si) << ", " << dd.dim_subset(si) << ")");
	}
	UG_LOG("\n");

//	Function informations
	UG_LOG(dd.num_fct() << " Function(s) defined (Symbolic Name, dim): ");
	for(size_t fct = 0; fct < dd.num_fct(); ++fct)
	{
		if(fct > 0) UG_LOG(", ");
		UG_LOG("(" << dd.name(fct) << ", " << dd.dim(fct) << ")");
	}
	UG_LOG("\n");

//	print subsets of functions
	if(verboseLev >= 2)
	{
		UG_LOG("Function definition on subsets: \n");
		for(size_t fct = 0; fct < dd.num_fct(); ++fct)
		{
			UG_LOG("   "<<dd.name(fct) << ": ");
			if(dd.is_def_everywhere(fct)) UG_LOG("[everywhere] ");
			bool bFirst = true;
			for(int si = 0; si < dd.num_subsets(); ++si)
			{
				if(bFirst) bFirst = false; else UG_LOG(", ");
				if(!dd.is_def_in_subset(fct, si)) continue;
				UG_LOG(dd.subset_name(si));
			}
		}
		UG_LOG("\n");
	}

//	print info of dofdistribution
	dd.print_local_dof_statistic(verboseLev);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_local_dof_statistic(int verboseLev) const
{
	if(surface_indices_enabled())
	{
		UG_LOG("\nLocal DoF Statistic for Surface DoFDistribution:");
		UG_LOG("\n-----------------------------------------------\n");
		print_local_dof_statistic(*m_pSurfDD, verboseLev);
	}

	if(level_indices_enabled())
	{
		const int topLev = (int)num_levels() - 1;
		if(topLev < 0) return;
		UG_LOG("\nLocal DoF Statistic for Level DoFDistribution on Level "<<topLev<<":");
		UG_LOG("\n----------------------------------------------------------\n");
		print_local_dof_statistic(*m_vLevelDD[topLev], verboseLev);
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
print_layout_statistic(int verboseLev) const
{
//	Write info
	UG_LOG(" No Layouts in sequential code.\n");
}

template <typename TDoFDistribution>
void MGDoFManager<TDoFDistribution>::surface_view_required()
{
// 	serial version
	if(m_pSurfaceView != NULL) return;

//	Create Surface view if not already created
	try{
		m_pSurfaceView = new SurfaceView(*m_pMultiGrid);
	}
	UG_CATCH_THROW("Cannot allocate SurfaceView.");

// 	Create surface view for all elements
	try{
		CreateSurfaceView(*m_pSurfaceView, *m_pMultiGrid, *m_pMGSubsetHandler);
	}
	UG_CATCH_THROW("Cannot create surface view.");
}

template <typename TDoFDistribution>
void MGDoFManager<TDoFDistribution>::surface_distribution_required()
{
//	Create surface view iff needed
	try{
		surface_view_required();
	}
	UG_CATCH_THROW("Cannot create surface view.");

// 	set storage manager
	try{
		m_surfaceStorageManager.set_subset_handler(m_pSurfaceView);
	}
	UG_CATCH_THROW("Cannot set surface view for storage manager.");

// 	Create surface dof distributions
	if(m_pSurfDD == NULL)
	{
	//	create dof distribution on surface
		try{
		m_pSurfDD =
			new TDoFDistribution(m_pSurfaceView->get_geometric_objects(),
								 *m_pSurfaceView,
								 m_surfaceStorageManager,
								 *m_pFuncPattern,
								 *m_pSurfaceView);
		}
		UG_CATCH_THROW("Cannot create Surface DoFDistribution.");

		m_pSurfDD->set_grouping(m_bGrouped);
	}
}

template <typename TDoFDistribution>
void MGDoFManager<TDoFDistribution>::level_distribution_required(size_t numLevel)
{
// 	set StorageManager for levels
	try{
		m_levelStorageManager.set_subset_handler(m_pMGSubsetHandler);
	}
	UG_CATCH_THROW("Cannot assign SubsetHandler to Storage Manager.");

// 	Create level dof distributions
	for(size_t l = m_vLevelDD.size(); l < numLevel; ++l)
	{
		try{
			m_vLevelDD.push_back(
				new TDoFDistribution(m_pMGSubsetHandler->get_geometric_objects_in_level(l),
									 *m_pMGSubsetHandler,
									 m_levelStorageManager,
									 *m_pFuncPattern));
		}
		UG_CATCH_THROW("Cannot create DoFDistribution on level "<<l);

		m_vLevelDD[l]->set_grouping(m_bGrouped);
	}
}

template <typename TDoFDistribution>
void MGDoFManager<TDoFDistribution>::disable_surface_indices()
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
disable_level_indices()
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 	Handling of surface view
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
is_in_surface_view(GeometricObject* obj)
{
	uint type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:
			return is_in_surface_view(reinterpret_cast<VertexBase*>(obj));
		case EDGE:
			return is_in_surface_view(reinterpret_cast<EdgeBase*>(obj));
		case FACE:
			return is_in_surface_view(reinterpret_cast<Face*>(obj));
		case VOLUME:
			return is_in_surface_view(reinterpret_cast<Volume*>(obj));
	}

	throw(UGFatalError("GeomObject type not known."));
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
is_in_surface_view(VertexBase* vrt)
{
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");
	return (m_pSurfaceView->get_subset_index(vrt) >= 0);
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
is_in_surface_view(EdgeBase* edge)
{
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");
	return (m_pSurfaceView->get_subset_index(edge) >= 0);
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
is_in_surface_view(Face* face)
{
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");
	return (m_pSurfaceView->get_subset_index(face) >= 0);
}

template <typename TDoFDistribution>
bool
MGDoFManager<TDoFDistribution>::
is_in_surface_view(Volume* vol)
{
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");
	return (m_pSurfaceView->get_subset_index(vol) >= 0);
}


template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_surface_view(GeometricObject* obj)
{
	uint type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:
			add_to_surface_view(reinterpret_cast<VertexBase*>(obj)); return;
		case EDGE:
			add_to_surface_view(reinterpret_cast<EdgeBase*>(obj)); return;
		case FACE:
			add_to_surface_view(reinterpret_cast<Face*>(obj)); return;
		case VOLUME:
			add_to_surface_view(reinterpret_cast<Volume*>(obj)); return;
	}

	throw(UGFatalError("GeomObject type not known."));
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_surface_view(VertexBase* vrt)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");

	const int si = m_pMGSubsetHandler->get_subset_index(vrt);
	m_pSurfaceView->assign_subset(vrt, si);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_surface_view(EdgeBase* edge)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");

	const int si = m_pMGSubsetHandler->get_subset_index(edge);
	m_pSurfaceView->assign_subset(edge, si);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_surface_view(Face* face)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");

	const int si = m_pMGSubsetHandler->get_subset_index(face);
	m_pSurfaceView->assign_subset(face, si);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_surface_view(Volume* vol)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");

	const int si = m_pMGSubsetHandler->get_subset_index(vol);
	m_pSurfaceView->assign_subset(vol, si);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_surface_view(GeometricObject* obj)
{
	uint type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:
			remove_from_surface_view(reinterpret_cast<VertexBase*>(obj)); return;
		case EDGE:
			remove_from_surface_view(reinterpret_cast<EdgeBase*>(obj)); return;
		case FACE:
			remove_from_surface_view(reinterpret_cast<Face*>(obj)); return;
		case VOLUME:
			remove_from_surface_view(reinterpret_cast<Volume*>(obj)); return;
	}

	throw(UGFatalError("GeomObject type not known."));
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_surface_view(VertexBase* vrt)
{
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");

	m_pSurfaceView->assign_subset(vrt, -1);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_surface_view(EdgeBase* edge)
{
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");

	m_pSurfaceView->assign_subset(edge, -1);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_surface_view(Face* face)
{
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");

	m_pSurfaceView->assign_subset(face, -1);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_surface_view(Volume* vol)
{
	UG_ASSERT(m_pSurfaceView != NULL, "Missing SurfaceView.");

	m_pSurfaceView->assign_subset(vol, -1);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 	Handling of level dof distribution
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_level_dof_distribution(GeometricObject* obj)
{
	uint type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:
			add_to_level_dof_distribution(reinterpret_cast<VertexBase*>(obj)); return;
		case EDGE:
			add_to_level_dof_distribution(reinterpret_cast<EdgeBase*>(obj)); return;
		case FACE:
			add_to_level_dof_distribution(reinterpret_cast<Face*>(obj)); return;
		case VOLUME:
			add_to_level_dof_distribution(reinterpret_cast<Volume*>(obj)); return;
	}

	throw(UGFatalError("GeomObject type not known."));
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_level_dof_distribution(VertexBase* vrt)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");

//	get level of Vertex
	const int level = m_pMGSubsetHandler->get_level(vrt);

//	check if level dof distribution must be created
	try{
		level_distribution_required(level+1);
	}
	UG_CATCH_THROW("Cannot create level dof distribution");

//	announce element to level dof distribution
	level_dof_distribution(level)->grid_obj_added(vrt);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_level_dof_distribution(EdgeBase* edge)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");

//	get level of Vertex
	const int level = m_pMGSubsetHandler->get_level(edge);

//	check if level dof distribution must be created
	try{
		level_distribution_required(level+1);
	}
	UG_CATCH_THROW("Cannot create level dof distribution");

//	announce element to level dof distribution
	level_dof_distribution(level)->grid_obj_added(edge);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_level_dof_distribution(Face* face)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");

//	get level of Vertex
	const int level = m_pMGSubsetHandler->get_level(face);

//	check if level dof distribution must be created
	try{
		level_distribution_required(level+1);
	}
	UG_CATCH_THROW("Cannot create level dof distribution");

//	announce element to level dof distribution
	level_dof_distribution(level)->grid_obj_added(face);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
add_to_level_dof_distribution(Volume* vol)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");

//	get level of Vertex
	const int level = m_pMGSubsetHandler->get_level(vol);

//	check if level dof distribution must be created
	try{
		level_distribution_required(level+1);
	}
	UG_CATCH_THROW("Cannot create level dof distribution");

//	announce element to level dof distribution
	level_dof_distribution(level)->grid_obj_added(vol);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_level_dof_distribution(GeometricObject* obj)
{
	uint type = obj->base_object_type_id();
	switch(type)
	{
		case VERTEX:
			remove_from_level_dof_distribution(reinterpret_cast<VertexBase*>(obj)); return;
		case EDGE:
			remove_from_level_dof_distribution(reinterpret_cast<EdgeBase*>(obj)); return;
		case FACE:
			remove_from_level_dof_distribution(reinterpret_cast<Face*>(obj)); return;
		case VOLUME:
			remove_from_level_dof_distribution(reinterpret_cast<Volume*>(obj)); return;
	}

	throw(UGFatalError("GeomObject type not known."));
}


template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_level_dof_distribution(VertexBase* vrt)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");

//	get level of Vertex
	const int level = m_pMGSubsetHandler->get_level(vrt);

//	check if level dof distribution exists
	try{level_distribution_required(level+1);}
	UG_CATCH_THROW("Level dof distribution does not exist.");

//	announce removal of element to level dof distribution
	level_dof_distribution(level)->grid_obj_to_be_removed(vrt);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_level_dof_distribution(EdgeBase* edge)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");

//	get level of Vertex
	const int level = m_pMGSubsetHandler->get_level(edge);

//	check if level dof distribution exists
	try{level_distribution_required(level+1);}
	UG_CATCH_THROW("Level dof distribution does not exist.");

//	announce removal of element to level dof distribution
	level_dof_distribution(level)->grid_obj_to_be_removed(edge);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_level_dof_distribution(Face* face)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");

//	get level of Vertex
	const int level = m_pMGSubsetHandler->get_level(face);

//	check if level dof distribution exists
	try{level_distribution_required(level+1);}
	UG_CATCH_THROW("Level dof distribution does not exist.");

//	announce removal of element to level dof distribution
	level_dof_distribution(level)->grid_obj_to_be_removed(face);
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
remove_from_level_dof_distribution(Volume* vol)
{
	UG_ASSERT(m_pMGSubsetHandler != NULL, "Missing MGSubsetHandler.");

//	get level of Vertex
	const int level = m_pMGSubsetHandler->get_level(vol);

//	check if level dof distribution exists
	try{level_distribution_required(level+1);}
	UG_CATCH_THROW("Level dof distribution does not exist.");

//	announce removal of element to level dof distribution
	level_dof_distribution(level)->grid_obj_to_be_removed(vol);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 	Grid Observer Callbacks
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
vertex_created(Grid* grid, VertexBase* vrt,	GeometricObject* pParent, bool replacesParent)
{
//	if level dofs enabled, add the element to level dofs
	if(level_indices_enabled())
		add_to_level_dof_distribution(vrt);

//	if surface dofs enabled, add the element to surface dofs
	if(surface_indices_enabled())
	{
	//	handle: we have two possibilities here, what the creation could be:
	//		1. a regular refinement, thus we exclude the parent from SurfaceView
	//			and add the created obj the the Surface View
	//		2. a replaced obj; we also exclude the parent from the SurfaceView
	//			since it will be erased and add the created obj depending if
	//			the created obj has already children

	//	1. case: regular element creation
		if(!replacesParent)
		{
		//	first, handle the parent: remove it from surface view
			if(pParent != NULL)
				if(is_in_surface_view(pParent))
				{
				//	a) Release index for parent (which may not be part of surface
				//     view after adding the child; created shadows are not part of
				//	   the surface view at this stage.)
					surface_dof_distribution()->grid_obj_to_be_removed(pParent);

				//	b) Remove parent from surface view
					remove_from_surface_view(pParent);
				}

		//	now add the obj to the surface view:
		//	a) Add the created element to the surface view
			add_to_surface_view(vrt);

		//	b) Add index
			surface_dof_distribution()->grid_obj_added(vrt);
		}
	//	2. case: object is replaced
		else
		{
		//	replaced object
			VertexBase* pReplaced = reinterpret_cast<VertexBase*>(pParent);

		//	check if created element has no children. In this case, we copy
		//	the dofs from the replaced object ...
			if(!m_pMultiGrid->has_children(vrt))
			{
			//	a) Add the created element to the surface view
				add_to_surface_view(vrt);

			//	b) Add index
				surface_dof_distribution()->grid_obj_replaced(vrt, pReplaced);
			}

		//	in any case the old obj is dropped. The index is not removed, since
		//	the index remains valid on the shadowing obj.
			remove_from_surface_view(pReplaced);
		}

	//	NOTE: at this point the parent element may be a shadow. But we
	//		  do not add the shadow here, but on a later call of defragment
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
edge_created(Grid* grid, EdgeBase* edge,	GeometricObject* pParent, bool replacesParent)
{
//	if level dofs enabled, add the element to level dofs
	if(level_indices_enabled())
		add_to_level_dof_distribution(edge);

//	if surface dofs enabled, add the element to surface dofs
	if(surface_indices_enabled())
	{
	//	handle: we have two possibilities here, what the creation could be:
	//		1. a regular refinement, thus we exclude the parent from SurfaceView
	//			and add the created obj the the Surface View
	//		2. a replaced obj; we also exclude the parent from the SurfaceView
	//			since it will be erased and add the created obj depending if
	//			the created obj has already children

	//	1. case: regular element creation
		if(!replacesParent)
		{
		//	first, handle the parent: remove it from surface view
			if(pParent != NULL)
				if(is_in_surface_view(pParent))
				{
				//	a) Release index for parent (which may not be part of surface
				//     view after adding the child; created shadows are not part of
				//	   the surface view at this stage.)
					surface_dof_distribution()->grid_obj_to_be_removed(pParent);

				//	b) Remove parent from surface view
					remove_from_surface_view(pParent);
				}

		//	now add the obj to the surface view:
		//	a) Add the created element to the surface view
			add_to_surface_view(edge);

		//	b) Add index
			surface_dof_distribution()->grid_obj_added(edge);
		}
	//	2. case: object is replaced
		else
		{
		//	replaced object
			EdgeBase* pReplaced = reinterpret_cast<EdgeBase*>(pParent);

		//	check if created element has no children. In this case, we copy
		//	the dofs from the replaced object ...
			if(!m_pMultiGrid->has_children(edge))
			{
			//	a) Add the created element to the surface view
				add_to_surface_view(edge);

			//	b) Add index
				surface_dof_distribution()->grid_obj_replaced(edge, pReplaced);
			}
			else
			{
			//	remove also subelements from surface view
				std::vector<VertexBase*> vVertex;
				CollectAssociated(vVertex, *grid, pReplaced);
				for(size_t i = 0; i < vVertex.size(); ++i){
					remove_from_surface_view(vVertex[i]);
				}
			}

		//	in any case the old obj is dropped. The index is not removed, since
		//	the index remains valid on the shadowing obj.
			remove_from_surface_view(pReplaced);
		}

	//	NOTE: at this point the parent element may be a shadow. But we
	//		  do not add the shadow here, but on a later call of defragment
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
face_created(Grid* grid, Face* face,	GeometricObject* pParent, bool replacesParent)
{
//	if level dofs enabled, add the element to level dofs
	if(level_indices_enabled())
		add_to_level_dof_distribution(face);

//	if surface dofs enabled, add the element to surface dofs
	if(surface_indices_enabled())
	{
	//	handle: we have two possibilities here, what the creation could be:
	//		1. a regular refinement, thus we exclude the parent from SurfaceView
	//			and add the created obj the the Surface View
	//		2. a replaced obj; we also exclude the parent from the SurfaceView
	//			since it will be erased and add the created obj depending if
	//			the created obj has already children

	//	1. case: regular element creation
		if(!replacesParent)
		{
		//	first, handle the parent: remove it from surface view
			if(pParent != NULL)
				if(is_in_surface_view(pParent))
				{
				//	a) Release index for parent (which may not be part of surface
				//     view after adding the child; created shadows are not part of
				//	   the surface view at this stage.)
					surface_dof_distribution()->grid_obj_to_be_removed(pParent);

				//	b) Remove parent from surface view
					remove_from_surface_view(pParent);
				}

		//	now add the obj to the surface view:
		//	a) Add the created element to the surface view
			add_to_surface_view(face);

		//	b) Add index
			surface_dof_distribution()->grid_obj_added(face);
		}
	//	2. case: object is replaced
		else
		{
		//	replaced object
			Face* pReplaced = reinterpret_cast<Face*>(pParent);

		//	check if created element has no children. In this case, we copy
		//	the dofs from the replaced object ...
			if(!m_pMultiGrid->has_children(face))
			{
			//	a) Add the created element to the surface view
				add_to_surface_view(face);

			//	b) Add index
				surface_dof_distribution()->grid_obj_replaced(face, pReplaced);
			}
			else
			{
			//	remove also subelements from surface view
				std::vector<VertexBase*> vVertex;
				CollectAssociated(vVertex, *grid, pReplaced);
				for(size_t i = 0; i < vVertex.size(); ++i){
					remove_from_surface_view(vVertex[i]);
				}
				std::vector<EdgeBase*> vEdge;
				CollectAssociated(vEdge, *grid, pReplaced);
				for(size_t i = 0; i < vEdge.size(); ++i){
					remove_from_surface_view(vEdge[i]);
				}
			}

		//	in any case the old obj is dropped. The index is not removed, since
		//	the index remains valid on the shadowing obj.
			remove_from_surface_view(pReplaced);
		}

	//	NOTE: at this point the parent element may be a shadow. But we
	//		  do not add the shadow here, but on a later call of defragment
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
volume_created(Grid* grid, Volume* vol,	GeometricObject* pParent, bool replacesParent)
{
//	if level dofs enabled, add the element to level dofs
	if(level_indices_enabled())
		add_to_level_dof_distribution(vol);

//	if surface dofs enabled, add the element to surface dofs
	if(surface_indices_enabled())
	{
	//	1. case: regular element creation
		if(!replacesParent)
		{
		//	first, handle the parent: remove it from surface view
			if(pParent != NULL)
				if(is_in_surface_view(pParent))
				{
				//	a) Release index for parent (which may not be part of surface
				//     view after adding the child; created shadows are not part of
				//	   the surface view at this stage.)
					surface_dof_distribution()->grid_obj_to_be_removed(pParent);

				//	b) Remove parent from surface view
					remove_from_surface_view(pParent);
				}

		//	now add the obj to the surface view:
		//	a) Add the created element to the surface view
			add_to_surface_view(vol);

		//	b) Add index
			surface_dof_distribution()->grid_obj_added(vol);
		}
	//	2. case: object is replaced
		else
		{
			throw(UGFatalError("There are no Volume replacements in 3D."));
		}

	//	NOTE: at this point the parent element may be a shadow. But we
	//		  do not add the shadow here, but on a later call of defragment
	}
}


template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
vertex_to_be_erased(Grid* grid, VertexBase* vrt, VertexBase* replacedBy)
{
	UG_ASSERT(m_pMultiGrid != NULL, "No MultiGrid.");

//	if level dofs enabled, remove element from level dofs
	if(level_indices_enabled())
		remove_from_level_dof_distribution(vrt);

//	if surface dofs enabled, remove element from surface dofs
	if(surface_indices_enabled())
	{
	//	The case of replaced objects has been handled in create part
		if(replacedBy != NULL) return;

	//	1. Remove index for erased element
		surface_dof_distribution()->grid_obj_to_be_removed(vrt);

	//	get parent
		GeometricObject* pParent = m_pMultiGrid->get_parent(vrt);
		if(pParent != NULL && !is_in_surface_view(pParent))
		{
		//	2. Add parent to surface view
			add_to_surface_view(pParent);

		//	3. Add index for parent
			surface_dof_distribution()->grid_obj_added(pParent);
		}
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
edge_to_be_erased(Grid* grid, EdgeBase* edge, EdgeBase* replacedBy)
{
	UG_ASSERT(m_pMultiGrid != NULL, "No MultiGrid.");

//	if level dofs enabled, remove element from level dofs
	if(level_indices_enabled())
		remove_from_level_dof_distribution(edge);

//	if surface dofs enabled, remove element from surface dofs
	if(surface_indices_enabled())
	{
	//	The case of replaced objects has been handled in create part
		if(replacedBy != NULL) return;

	//	1. Remove index for erased element
		surface_dof_distribution()->grid_obj_to_be_removed(edge);

	//	get parent
		GeometricObject* pParent = m_pMultiGrid->get_parent(edge);
		if(pParent != NULL && !is_in_surface_view(pParent))
		{
		//	2. Add parent to surface view
			add_to_surface_view(pParent);

		//	3. Add index for parent
			surface_dof_distribution()->grid_obj_added(pParent);
		}
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
face_to_be_erased(Grid* grid, Face* face, Face* replacedBy)
{
	UG_ASSERT(m_pMultiGrid != NULL, "No MultiGrid.");

//	if level dofs enabled, remove element from level dofs
	if(level_indices_enabled())
		remove_from_level_dof_distribution(face);

//	if surface dofs enabled, remove element from surface dofs
	if(surface_indices_enabled())
	{
	//	The case of replaced objects has been handled in create part
		if(replacedBy != NULL) return;

	//	1. Remove index for erased element
		surface_dof_distribution()->grid_obj_to_be_removed(face);

	//	get parent
		GeometricObject* pParent = m_pMultiGrid->get_parent(face);
		if(pParent != NULL && !is_in_surface_view(pParent))
		{
		//	2. Add parent to surface view
			add_to_surface_view(pParent);

		//	3. Add index for parent
			surface_dof_distribution()->grid_obj_added(pParent);
		}
	}
}

template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
volume_to_be_erased(Grid* grid, Volume* vol, Volume* replacedBy)
{
	UG_ASSERT(m_pMultiGrid != NULL, "No MultiGrid.");

//	if level dofs enabled, remove element from level dofs
	if(level_indices_enabled())
		remove_from_level_dof_distribution(vol);

//	if surface dofs enabled, remove element from surface dofs
	if(surface_indices_enabled())
	{
	//	The case of replaced objects has been handled in create part
		if(replacedBy == NULL) return;

	//	1. Remove index for erased element
		surface_dof_distribution()->grid_obj_to_be_removed(vol);

	//	get parent
		GeometricObject* pParent = m_pMultiGrid->get_parent(vol);
		if(pParent != NULL && !is_in_surface_view(pParent))
		{
		//	2. Add parent to surface view
			add_to_surface_view(pParent);

		//	3. Add index for parent
			surface_dof_distribution()->grid_obj_added(pParent);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// 	Defragment
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


template <typename TDoFDistribution>
void
MGDoFManager<TDoFDistribution>::
defragment()
{
//	we add the shadows to the surface view, that may have been created
//	due to grid adaption

//	check grid options
//	\todo: all needed ?
	if(!m_pMultiGrid->option_is_enabled(VOLOPT_AUTOGENERATE_FACES)){
	  UG_LOG("WARNING: Auto-enabling grid option VOLOPT_AUTOGENERATE_FACES.\n");
	  m_pMultiGrid->enable_options(VOLOPT_AUTOGENERATE_FACES);
	}
/*	if(!m_pMultiGrid->option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)){
	  UG_LOG("WARNING: Auto-enabling grid option VOLOPT_AUTOGENERATE_EDGES.\n");
	  m_pMultiGrid->enable_options(VOLOPT_AUTOGENERATE_EDGES);
	}
*/	if(!m_pMultiGrid->option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
	  UG_LOG("WARNING: Auto-enabling grid option FACEOPT_AUTOGENERATE_EDGES.\n");
	  m_pMultiGrid->enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

//	add missing shadows to surface view
	if(surface_indices_enabled())
	{
		add_associated_sides_to_surface_view<Volume>();
		add_associated_sides_to_surface_view<Face>();
		add_associated_sides_to_surface_view<EdgeBase>();
	}

//	defragment dof distributions
	if(surface_indices_enabled())
		surface_dof_distribution()->defragment();
	if(level_indices_enabled())
		for(size_t lev = 0; lev < num_levels(); ++lev)
			level_dof_distribution(lev)->defragment();

}

template <typename TDoFDistribution>
template <class TElem>
void
MGDoFManager<TDoFDistribution>::
add_associated_sides_to_surface_view()
{
	typedef typename geometry_traits<TElem>::const_iterator iterator;
	typedef typename TElem::lower_dim_base_object Side;
	std::vector<Side*> vSides;
	Grid& grid = *m_pSurfaceView->get_assigned_grid();

	for(size_t l  = 0; l < m_pSurfaceView->num_levels(); ++l){
		for(int si = 0; si < m_pSurfaceView->num_subsets(); ++si){
			for(iterator iter = m_pSurfaceView->begin<TElem>(si, l);
				iter != m_pSurfaceView->end<TElem>(si, l); ++iter)
			{
				TElem* e = *iter;
				CollectAssociated(vSides, grid, e);

				for(size_t i = 0; i < vSides.size(); ++i)
				{
					Side* s = vSides[i];

					if(!is_in_surface_view(s))
					{
						add_to_surface_view(s);
						surface_dof_distribution()->grid_obj_added(s);
					}
				}
			}
		}
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__MG_DOF_MANAGER_IMPL__ */
