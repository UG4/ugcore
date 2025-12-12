/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "grid.h"

#include <cassert>
#include <algorithm>

//ø #include "grid_util.h"
#include "common/common.h"
#include "lib_grid/attachments/attached_list.h"
#include "lib_grid/tools/periodic_boundary_manager.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#endif


using namespace std;

namespace ug {

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	implementation of Grid

////////////////////////////////////////////////////////////////////////
//	constructors
Grid::Grid() :
	m_aVertexContainer("Grid_VertexContainer", false),
	m_aEdgeContainer("Grid_EdgeContainer", false),
	m_aFaceContainer("Grid_FaceContainer", false),
	m_aVolumeContainer("Grid_VolumeContainer", false),
	m_bMarking(false),
	m_aMark("Grid_Mark", false),
	m_distGridMgr(nullptr),
	m_periodicBndMgr(nullptr)
{
	m_hashCounter = 0;
	m_currentMark = 0;
	m_options = GridOptions::GRIDOPT_NONE;
	m_messageHub = SPMessageHub(new MessageHub());

	change_options(GridOptions::GRIDOPT_DEFAULT);
}

Grid::Grid(uint options) :
	m_aVertexContainer("Grid_VertexContainer", false),
	m_aEdgeContainer("Grid_EdgeContainer", false),
	m_aFaceContainer("Grid_FaceContainer", false),
	m_aVolumeContainer("Grid_VolumeContainer", false),
	m_bMarking(false),
	m_aMark("Grid_Mark", false),
	m_distGridMgr(nullptr),
	m_periodicBndMgr(nullptr)
{
	m_hashCounter = 0;
	m_currentMark = 0;
	m_options = GridOptions::GRIDOPT_NONE;
	m_messageHub = SPMessageHub(new MessageHub());

	change_options(options);
}

Grid::Grid(const Grid& grid) :
	m_aVertexContainer("Grid_VertexContainer", false),
	m_aEdgeContainer("Grid_EdgeContainer", false),
	m_aFaceContainer("Grid_FaceContainer", false),
	m_aVolumeContainer("Grid_VolumeContainer", false),
	m_bMarking(false),
	m_aMark("Grid_Mark", false),
	m_distGridMgr(nullptr),
	m_periodicBndMgr(nullptr)
{
	m_hashCounter = 0;
	m_currentMark = 0;
	m_options = GridOptions::GRIDOPT_NONE;
	m_messageHub = SPMessageHub(new MessageHub());

	assign_grid(grid);
}

Grid::~Grid()
{
	notify_and_clear_observers_on_grid_destruction();

//	erase all elements
	clear_geometry();

//	remove marks - would be done anyway...
	remove_marks();

//	erase any internal managers and handlers
	#ifdef UG_PARALLEL
		if(m_distGridMgr) delete m_distGridMgr;
	#endif

	if(m_periodicBndMgr) delete m_periodicBndMgr;
}

void Grid::notify_and_clear_observers_on_grid_destruction(GridObserver* initiator)
{
//	tell registered grid-observers that the grid is to be destroyed.
//	do this in reverse order, so that the danger of accessing invalid observers
//	is minimized.
	for(auto iter = m_gridObservers.rbegin();
	    iter != m_gridObservers.rend(); ++iter)
	{
		if(*iter != initiator)
			(*iter)->grid_to_be_destroyed(this);
	}

//	unregister all observers	
	while(!m_gridObservers.empty())
		unregister_observer(m_gridObservers.back());
	
	while(!m_vertexObservers.empty())
		unregister_observer(m_vertexObservers.back());

	while(!m_edgeObservers.empty())
		unregister_observer(m_edgeObservers.back());

	while(!m_faceObservers.empty())
		unregister_observer(m_faceObservers.back());

	while(!m_volumeObservers.empty())
		unregister_observer(m_volumeObservers.back());
}

void Grid::
set_parallel(bool parallel)
{
	if(parallel){
		#ifdef UG_PARALLEL
			if(!is_parallel()){
			//	we currently only support parallel mutli-grids, sadly...
				auto mg = dynamic_cast<MultiGrid*>(this);
				if(mg){
					if(m_distGridMgr) delete m_distGridMgr;
					m_distGridMgr = new DistributedGridManager;
					m_distGridMgr->assign(*mg);
				}
				else{
					UG_THROW("Error during Grid::set_parallel: "
							"The DistributedGridManager can currently only be used to "
							"parallelize ug::MultiGrid. ug::Grid can currently not be used "
							"for parallel computations. Sorry.");
				}
			}
		#else
			UG_THROW("Parallelism can only be activated, if ug was compiled with "
					"PARALLEL=ON. This is not the case for this application.");
		#endif
	}
	else if(is_parallel()){
		#ifdef UG_PARALLEL
			if(m_distGridMgr) delete m_distGridMgr;
			m_distGridMgr = nullptr;
		#endif
	}
}

void Grid::set_periodic_boundaries(bool is_periodic)
{
	if(is_periodic)
	{
		if(m_periodicBndMgr) delete m_periodicBndMgr;
		m_periodicBndMgr = new PeriodicBoundaryManager();
		m_periodicBndMgr->set_grid(this);
	}
	else if(m_periodicBndMgr){
		delete m_periodicBndMgr;
		m_periodicBndMgr = nullptr;
	}
}

bool Grid::has_periodic_boundaries() const
{
	return m_periodicBndMgr != nullptr;
}

PeriodicBoundaryManager* Grid::periodic_boundary_manager()
{
	return m_periodicBndMgr;
}

const PeriodicBoundaryManager* Grid::periodic_boundary_manager() const
{
	return m_periodicBndMgr;
}


void Grid::clear()
{
	clear_geometry();
	clear_attachments();
}

void Grid::clear_geometry()
{
//	disable all options to speed it up
	uint opts = get_options();
	set_options(GridOptions::GRIDOPT_NONE);
	
	clear<Volume>();
	clear<Face>();
	clear<Edge>();
	clear<Vertex>();
	
//	reset options
	set_options(opts);
}

template <typename TElem>
void Grid::clear_attachments()
{
	using AttachmentPipe = typename traits<TElem>::AttachmentPipe;

	vector<AttachmentEntry>	vEntries;

//	iterate through all attachment pipes
	AttachmentPipe& ap = get_attachment_pipe<TElem>();

//	collect all attachment entries
	for(typename AttachmentPipe::ConstAttachmentEntryIterator iter = ap.attachments_begin();
		iter != ap.attachments_end(); ++iter)
	{
			vEntries.push_back(*iter);
	}

//	iterate through the entries in the vector and delete the ones
//	that have an enabled pass-on behaviour
	for(size_t j = 0; j < vEntries.size(); ++j)
	{
		const AttachmentEntry& ae = vEntries[j];
		
		if(ae.m_userData == 1){
			ap.detach(*ae.m_pAttachment);
		}
	}
}

void Grid::clear_attachments()
{
	clear_attachments<Vertex>();
	clear_attachments<Edge>();
	clear_attachments<Face>();
	clear_attachments<Volume>();
}

Grid& Grid::operator = (const Grid& grid)
{
//	clears the grid and calls assign_grid afterwards.
//	we're disabling any options, since new options will
//	be set during assign_grid anyway. This might speed
//	things up a little
	set_options(GridOptions::GRIDOPT_NONE);
	clear_geometry();
	assign_grid(grid);
	
	return *this;
}

template <typename TAttachmentPipe>
void Grid::copy_user_attachments(const TAttachmentPipe& apSrc, TAttachmentPipe& apDest,
								vector<int>& srcDataIndices)
{
	for(typename TAttachmentPipe::ConstAttachmentEntryIterator iter = apSrc.attachments_begin();
		iter != apSrc.attachments_end(); ++iter)
	{
		const AttachmentEntry& ae = *iter;
		if(ae.m_userData == 1){
		//	attach the attachment to this grid
			apDest.attach(*ae.m_pAttachment, ae.m_userData);
			const IAttachmentDataContainer& conSrc = *ae.m_pContainer;
			IAttachmentDataContainer& conDest = *apDest.get_data_container(*ae.m_pAttachment);

		//	we use the containers copy-method
			conSrc.copy_to_container(&conDest, &srcDataIndices.front(),
									 (int)srcDataIndices.size());
		}
	}
}

void Grid::assign_grid(const Grid& grid)
{
//TODO: notify a grid observer that copying has started

//	we need a vertex-map that allows us to find a vertex in the new grid
//	given a vertex in the old one.
	vector<Vertex*>	vrtMap(grid.attachment_container_size<Vertex>(), nullptr);

//	we need index-lists that allow us to copy attachments later on
	vector<int> vSrcDataIndex[GridBaseObjectId::NUM_GEOMETRIC_BASE_OBJECTS];
	vector<int>& vSrcDataIndexVRT = vSrcDataIndex[GridBaseObjectId::VERTEX];
	vector<int>& vSrcDataIndexEDGE = vSrcDataIndex[GridBaseObjectId::EDGE];
	vector<int>& vSrcDataIndexFACE = vSrcDataIndex[GridBaseObjectId::FACE];
	vector<int>& vSrcDataIndexVOL = vSrcDataIndex[GridBaseObjectId::VOLUME];

//	copy all vertices
	vSrcDataIndexVRT.resize(grid.num<Vertex>());
	ConstVertexIterator vrtsEnd = grid.end<Vertex>();
	for(ConstVertexIterator iter = grid.begin<Vertex>(); iter != vrtsEnd; ++iter)
	{
		Vertex* vrt = *iter;
		Vertex* nVrt = *create_by_cloning(vrt);
		vrtMap[grid.get_attachment_data_index(vrt)] = nVrt;
		vSrcDataIndexVRT[get_attachment_data_index(nVrt)] = grid.get_attachment_data_index(vrt);
	}

//	copy all edges
	vSrcDataIndexEDGE.resize(grid.num<Edge>());
	ConstEdgeIterator edgesEnd = grid.end<Edge>();
	for(ConstEdgeIterator iter = grid.begin<Edge>(); iter != edgesEnd; ++iter)
	{
		Edge* e = *iter;
		Edge* nE = *create_by_cloning(e, EdgeDescriptor(
											vrtMap[grid.get_attachment_data_index(e->vertex(0))],
											vrtMap[grid.get_attachment_data_index(e->vertex(1))]));
		vSrcDataIndexEDGE[get_attachment_data_index(nE)] = grid.get_attachment_data_index(e);
	}

//	copy all faces
	vSrcDataIndexFACE.resize(grid.num<Face>());
	FaceDescriptor fd;
	ConstFaceIterator facesEnd = grid.end<Face>();
	for(ConstFaceIterator iter = grid.begin<Face>(); iter != facesEnd; ++iter)
	{
		Face* f = *iter;
		uint numVrts = f->num_vertices();
		Face::ConstVertexArray vrts = f->vertices();

	//	fill the face descriptor
		if(numVrts != fd.num_vertices())
			fd.set_num_vertices(numVrts);
		
		for(uint i = 0; i < numVrts; ++i)
			fd.set_vertex(i, vrtMap[grid.get_attachment_data_index(vrts[i])]);

	//	create the new face
		Face* nF = *create_by_cloning(f, fd);

		vSrcDataIndexFACE[get_attachment_data_index(nF)] = grid.get_attachment_data_index(f);
	}

//	copy all volumes
	vSrcDataIndexVOL.resize(grid.num<Volume>());
	VolumeDescriptor vd;
	ConstVolumeIterator volsEnd = grid.end<Volume>();
	for(ConstVolumeIterator iter = grid.begin<Volume>(); iter != volsEnd; ++iter)
	{
		Volume* v = *iter;
		uint numVrts = v->num_vertices();
		Volume::ConstVertexArray vrts = v->vertices();

	//	fill the volume descriptor
		if(numVrts != vd.num_vertices())
			vd.set_num_vertices(numVrts);

		for(uint i = 0; i < numVrts; ++i)
			vd.set_vertex(i, vrtMap[grid.get_attachment_data_index(vrts[i])]);

	//	create the volume
		Volume* nV = *create_by_cloning(v, vd);

		vSrcDataIndexVOL[get_attachment_data_index(nV)] = grid.get_attachment_data_index(v);
	}

//	enable options
	enable_options(grid.get_options());

//	copy attachments that may be passed on
	copy_user_attachments(grid.m_vertexElementStorage.m_attachmentPipe,
					m_vertexElementStorage.m_attachmentPipe, vSrcDataIndexVRT);
	copy_user_attachments(grid.m_edgeElementStorage.m_attachmentPipe,
					m_edgeElementStorage.m_attachmentPipe, vSrcDataIndexEDGE);
	copy_user_attachments(grid.m_faceElementStorage.m_attachmentPipe,
					m_faceElementStorage.m_attachmentPipe, vSrcDataIndexFACE);
	copy_user_attachments(grid.m_volumeElementStorage.m_attachmentPipe,
					m_volumeElementStorage.m_attachmentPipe, vSrcDataIndexVOL);

//	parallelism
	if(grid.is_parallel()){
		set_parallel(true);
	//todo:	copy interfaces from grid.
	}

//TODO: notify a grid observer that copying has ended
}


VertexIterator Grid::create_by_cloning(Vertex* pCloneMe, GridObject* pParent)
{
	auto pNew = reinterpret_cast<Vertex*>(pCloneMe->create_empty_instance());
	register_vertex(pNew, pParent);
	return iterator_cast<VertexIterator>(get_iterator(pNew));
}

EdgeIterator Grid::create_by_cloning(Edge* pCloneMe, const IVertexGroup& ev, GridObject* pParent)
{
	auto pNew = reinterpret_cast<Edge*>(pCloneMe->create_empty_instance());
	pNew->set_vertex(0, ev.vertex(0));
	pNew->set_vertex(1, ev.vertex(1));
	register_edge(pNew, pParent);
	return iterator_cast<EdgeIterator>(get_iterator(pNew));
}

FaceIterator Grid::create_by_cloning(Face* pCloneMe, const IVertexGroup& fv, GridObject* pParent)
{
	auto pNew = reinterpret_cast<Face*>(pCloneMe->create_empty_instance());
	uint numVrts = fv.num_vertices();
	Face::ConstVertexArray vrts = fv.vertices();
	for(uint i = 0; i < numVrts; ++i)
		pNew->set_vertex(i, vrts[i]);
	register_face(pNew, pParent);
	return iterator_cast<FaceIterator>(get_iterator(pNew));
}

VolumeIterator Grid::create_by_cloning(Volume* pCloneMe, const IVertexGroup& vv, GridObject* pParent)
{
	auto pNew = reinterpret_cast<Volume*>(pCloneMe->create_empty_instance());
	uint numVrts = vv.num_vertices();
	Volume::ConstVertexArray vrts = vv.vertices();
	for(uint i = 0; i < numVrts; ++i)
		pNew->set_vertex(i, vrts[i]);
	register_volume(pNew, pParent);
	return iterator_cast<VolumeIterator>(get_iterator(pNew));
}

////////////////////////////////////////////////////////////////////////
//	erase functions
void Grid::erase(GridObject* geomObj)
{
	assert(geomObj->container_section() != -1
			&& "ERROR in Grid::erase(Vertex*). Invalid pipe section!");

	uint objType = geomObj->base_object_id();
	switch(objType)
	{
		case GridBaseObjectId::VERTEX:
			erase(dynamic_cast<Vertex*>(geomObj));
			break;
		case GridBaseObjectId::EDGE:
			erase(dynamic_cast<Edge*>(geomObj));
			break;
		case GridBaseObjectId::FACE:
			erase(dynamic_cast<Face*>(geomObj));
			break;
		case GridBaseObjectId::VOLUME:
			erase(dynamic_cast<Volume*>(geomObj));
			break;
	};
}

void Grid::erase(Vertex* vrt)
{
	assert((vrt != nullptr) && "ERROR in Grid::erase(Vertex*): invalid pointer)");
	assert(vrt->container_section() != -1
			&& "ERROR in Grid::erase(Vertex*). Invalid pipe section!");

	unregister_vertex(vrt);

	delete vrt;
}

void Grid::erase(Edge* edge)
{
	assert((edge != nullptr) && "ERROR in Grid::erase(Edge*): invalid pointer)");
	assert(edge->container_section() != -1
			&& "ERROR in Grid::erase(Edge*). Invalid pipe section!");

	unregister_edge(edge);

	delete edge;
}

void Grid::erase(Face* face)
{
	assert((face != nullptr) && "ERROR in Grid::erase(Face*): invalid pointer)");
	assert(face->container_section() != -1
			&& "ERROR in Grid::erase(Face*). Invalid pipe section!");

	unregister_face(face);

	delete face;
}

void Grid::erase(Volume* vol)
{
	assert((vol != nullptr) && "ERROR in Grid::erase(Volume*): invalid pointer)");
	assert(vol->container_section() != -1
			&& "ERROR in Grid::erase(Volume*). Invalid pipe section!");

	unregister_volume(vol);

	delete vol;
}

//	the geometric-object-collection:
GridObjectCollection Grid::get_grid_objects()
{
	return GridObjectCollection(&m_vertexElementStorage.m_sectionContainer,
									 &m_edgeElementStorage.m_sectionContainer,
									 &m_faceElementStorage.m_sectionContainer,
									 &m_volumeElementStorage.m_sectionContainer);
}

void Grid::flip_orientation(Edge* e)
{
	swap(e->m_vertices[0], e->m_vertices[1]);
}

void Grid::flip_orientation(Face* f)
{
//	inverts the order of vertices.
	uint numVrts = (int)f->num_vertices();
	vector<Vertex*> vVrts(numVrts);
	
	uint i;
	for(i = 0; i < numVrts; ++i)
		vVrts[i] = f->vertex(i);
		
	for(i = 0; i < numVrts; ++i)
		f->set_vertex(i, vVrts[numVrts - 1 - i]);

//	update associated edge list
	if(option_is_enabled(FaceOptions::FACEOPT_STORE_ASSOCIATED_EDGES)){
		m_aaEdgeContainerFACE[f].clear();
		EdgeDescriptor ed;
		for(size_t ind = 0; ind < f->num_edges(); ++ind){
		//	get the descriptor of the i-th edge
			f->edge_desc(ind, ed);
		//	find the edge by checking vertices.
			Edge* e = find_edge_in_associated_edges(ed.vertex(0), ed);
			if(e)
				m_aaEdgeContainerFACE[f].push_back(e);
		}
	}
}

void Grid::flip_orientation(Volume* vol)
{
//	flips the orientation of volumes
//	get the descriptor for the flipped volume
	VolumeDescriptor vd;
	vol->get_flipped_orientation(vd);
	
//	change vertex order of the original volume
	size_t numVrts = vol->num_vertices();
	for(size_t i = 0; i < numVrts; ++i)
		vol->set_vertex(i, vd.vertex(i));

//	update associated edge list
	if(option_is_enabled(VolumeOptions::VOLOPT_STORE_ASSOCIATED_EDGES)){
		m_aaEdgeContainerVOLUME[vol].clear();
		EdgeDescriptor ed;
		for(size_t ind = 0; ind < vol->num_edges(); ++ind){
		//	get the descriptor of the i-th edge
			vol->edge_desc(ind, ed);
		//	find the edge by checking vertices.
			Edge* e = find_edge_in_associated_edges(ed.vertex(0), ed);
			if(e)
				m_aaEdgeContainerVOLUME[vol].push_back(e);
		}
	}

//	update associated face list
	if(option_is_enabled(VolumeOptions::VOLOPT_STORE_ASSOCIATED_FACES)){
		m_aaFaceContainerVOLUME[vol].clear();
		FaceDescriptor fd;
		for(size_t ind = 0; ind < vol->num_faces(); ++ind){
		//	get the descriptor of the i-th face
			vol->face_desc(ind, fd);
		//	find the face by checking vertices.
			Face* f = find_face_in_associated_faces(fd.vertex(0), fd);
			if(f)
				m_aaFaceContainerVOLUME[vol].push_back(f);
		}
	}
}

size_t Grid::vertex_fragmentation()
{
	return m_vertexElementStorage.m_attachmentPipe.num_data_entries()
			- m_vertexElementStorage.m_attachmentPipe.num_elements();
}

size_t Grid::edge_fragmentation()
{
	return m_edgeElementStorage.m_attachmentPipe.num_data_entries()
			- m_edgeElementStorage.m_attachmentPipe.num_elements();
}

size_t Grid::face_fragmentation()
{
	return m_faceElementStorage.m_attachmentPipe.num_data_entries()
			- m_faceElementStorage.m_attachmentPipe.num_elements();
}

size_t Grid::volume_fragmentation()
{
	return m_volumeElementStorage.m_attachmentPipe.num_data_entries()
			- m_volumeElementStorage.m_attachmentPipe.num_elements();
}


GridObject* Grid::
get_opposing_object(Vertex* vrt, Face* elem)
{
	std::pair<GridBaseObjectId, int> id = elem->get_opposing_object(vrt);
	switch(id.first){
		case GridBaseObjectId::VERTEX:
			return elem->vertex(id.second);
		case GridBaseObjectId::EDGE:
			return get_edge(elem, id.second);
		default:
			UG_THROW("Unsupported geometric base object type returned by "
					"Face::get_opposing_object(vrt)");
	}
}

GridObject* Grid::
get_opposing_object(Vertex* vrt, Volume* elem)
{
	std::pair<GridBaseObjectId, int> id = elem->get_opposing_object(vrt);
	switch(id.first){
		case GridBaseObjectId::VERTEX:
			return elem->vertex(id.second);
		case GridBaseObjectId::EDGE:
			return get_edge(elem, id.second);
		case GridBaseObjectId::FACE:
			return get_face(elem, id.second);
		default:
			UG_THROW("Unsupported geometric base object type returned by "
					"Volume::get_opposing_object(vrt)");
	}
}


////////////////////////////////////////////////////////////////////////
//	pass_on_values
template <typename TAttachmentPipe, typename TElem>
void Grid::pass_on_values(TAttachmentPipe& attachmentPipe,
							TElem* pSrc, TElem* pDest)
{
	for(auto iter = attachmentPipe.attachments_begin(); iter != attachmentPipe.attachments_end(); iter++)
	{
		if((*iter).m_userData == 1)
			(*iter).m_pContainer->copy_data(get_attachment_data_index(pSrc),
											get_attachment_data_index(pDest));
	}
}

void Grid::pass_on_values(Vertex* objSrc, Vertex* objDest)
{
	pass_on_values(m_vertexElementStorage.m_attachmentPipe, objSrc, objDest);
}

void Grid::pass_on_values(Edge* objSrc, Edge* objDest)
{
	pass_on_values(m_edgeElementStorage.m_attachmentPipe, objSrc, objDest);
}

void Grid::pass_on_values(Face* objSrc, Face* objDest)
{
	pass_on_values(m_faceElementStorage.m_attachmentPipe, objSrc, objDest);
}

void Grid::pass_on_values(Volume* objSrc, Volume* objDest)
{
	pass_on_values(m_volumeElementStorage.m_attachmentPipe, objSrc, objDest);
}

////////////////////////////////////////////////////////////////////////
//	options
void Grid::set_options(uint options)
{
	change_options(options);
}

uint Grid::get_options() const
{
	return m_options;
}

void Grid::enable_options(uint options)
{
	change_options(m_options | options);
}

void Grid::disable_options(uint options)
{
	change_options(m_options & (~options));
}

bool Grid::option_is_enabled(uint option) const
{
	return (m_options & option) == option;
}

void Grid::change_options(uint optsNew)
{
	change_vertex_options(optsNew &	0x000000FF);
	change_edge_options(optsNew & 	0x0000FF00);
	change_face_options(optsNew & 	0x00FF0000);
	change_volume_options(optsNew &	0xFF000000);
	assert((m_options == optsNew) && "Grid::change_options failed");
}
/*
void Grid::register_observer(GridObserver* observer, uint observerType)
{
//	check which elements have to be observed and store pointers to the observers.
//	avoid double-registration!
	ObserverContainer* observerContainers[] = {&m_gridObservers, &m_vertexObservers,
												&m_edgeObservers, & m_faceObservers, &m_volumeObservers};

	uint observerTypes[] = {OT_GRID_OBSERVER, OT_VERTEX_OBSERVER, OT_EDGE_OBSERVER, OT_FACE_OBSERVER, OT_VOLUME_OBSERVER};
	for(int i = 0; i < 5; ++i)
	{
		if((observerType & observerTypes[i]) == observerTypes[i])
		{
			ObserverContainer::iterator iter = find(observerContainers[i]->begin(), observerContainers[i]->end(), observer);
			if(iter == observerContainers[i]->end())
				observerContainers[i]->push_back(observer);
		}
	}

//	if the observer is a grid observer, notify him about the registration
	if((observerType & OT_GRID_OBSERVER) == OT_GRID_OBSERVER)
		observer->registered_at_grid(this);
}

void Grid::unregister_observer(GridObserver* observer)
{
//	check where the observer has been registered and erase the corresponding entries.
	ObserverContainer* observerContainers[] = {&m_gridObservers, &m_vertexObservers,
												&m_edgeObservers, & m_faceObservers, &m_volumeObservers};

	bool unregisterdFromGridObservers = false;
	for(int i = 0; i < 5; ++i)
	{
		ObserverContainer::iterator iter = find(observerContainers[i]->begin(), observerContainers[i]->end(), observer);
		if(iter != observerContainers[i]->end())
		{
			if(i == 0)
				unregisterdFromGridObservers = true;
			observerContainers[i]->erase(iter);
		}
	}

//	if the observer is a grid observer, notify him about the unregistration
	if(unregisterdFromGridObservers)
		observer->unregistered_from_grid(this);
}
*/
void Grid::register_observer(GridObserver* observer, ObserverType observerType)
{
//	check which elements have to be observed and store pointers to the observers.
//	avoid double-registration!
	if((observerType & ObserverType::OT_GRID_OBSERVER) == ObserverType::OT_GRID_OBSERVER)
	{
		auto iter = find(m_gridObservers.begin(),
		                 m_gridObservers.end(), observer);
		if(iter == m_gridObservers.end())
			m_gridObservers.push_back(observer);
	}

	if((observerType & ObserverType::OT_VERTEX_OBSERVER) == ObserverType::OT_VERTEX_OBSERVER)
	{
		auto iter = find(m_vertexObservers.begin(),
		                 m_vertexObservers.end(), observer);
		if(iter == m_vertexObservers.end())
			m_vertexObservers.push_back(observer);
	}

	if((observerType & ObserverType::OT_EDGE_OBSERVER) == ObserverType::OT_EDGE_OBSERVER)
	{
		auto iter = find(m_edgeObservers.begin(),
												m_edgeObservers.end(), observer);
		if(iter == m_edgeObservers.end())
			m_edgeObservers.push_back(observer);
	}

	if((observerType & ObserverType::OT_FACE_OBSERVER) == ObserverType::OT_FACE_OBSERVER)
	{
		auto iter = find(m_faceObservers.begin(),
												m_faceObservers.end(), observer);
		if(iter == m_faceObservers.end())
			m_faceObservers.push_back(observer);
	}

	if((observerType & ObserverType::OT_VOLUME_OBSERVER) == ObserverType::OT_VOLUME_OBSERVER)
	{
		auto iter = find(m_volumeObservers.begin(),
												m_volumeObservers.end(), observer);
		if(iter == m_volumeObservers.end())
			m_volumeObservers.push_back(observer);
	}

//	if the observer is a grid observer, notify him about the registration
//	if((observerType & ObserverType::OT_GRID_OBSERVER) == ObserverType::OT_GRID_OBSERVER)
//		observer->registered_at_grid(this);
}

void Grid::unregister_observer(GridObserver* observer)
{
//	check where the observer has been registered and erase the corresponding entries.
	//bool unregisterdFromGridObservers = false;

	{
		auto iter = find(m_gridObservers.begin(), m_gridObservers.end(), observer);
		if(iter != m_gridObservers.end()){
			m_gridObservers.erase(iter);
		}

//		unregisterdFromGridObservers = true;
	}

	{
		auto iter = find(m_vertexObservers.begin(), m_vertexObservers.end(), observer);
		if(iter != m_vertexObservers.end())
			m_vertexObservers.erase(iter);
	}

	{
		auto iter = find(m_edgeObservers.begin(), m_edgeObservers.end(), observer);
		if(iter != m_edgeObservers.end())
			m_edgeObservers.erase(iter);
	}

	{
		auto iter = find(m_faceObservers.begin(), m_faceObservers.end(), observer);
		if(iter != m_faceObservers.end())
			m_faceObservers.erase(iter);
	}

	{
		auto iter = find(m_volumeObservers.begin(), m_volumeObservers.end(), observer);
		if(iter != m_volumeObservers.end())
			m_volumeObservers.erase(iter);
	}

//	if the observer is a grid observer, notify him about the unregistration
//	if(unregisterdFromGridObservers)
//		observer->unregistered_from_grid(this);

}


////////////////////////////////////////////////////////////////////////
//	associated edge access
Grid::AssociatedEdgeIterator Grid::associated_edges_begin(Vertex* vrt)
{
	if(!option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_EDGES))
	{
		UG_LOG("WARNING in associated_edges_begin(vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES." << endl);
		vertex_store_associated_edges(true);
	}
	return m_aaEdgeContainerVERTEX[vrt].begin();
}

Grid::AssociatedEdgeIterator Grid::associated_edges_end(Vertex* vrt)
{
	if(!option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_EDGES))
	{
		UG_LOG("WARNING in associated_edges_end(vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES." << endl);
		vertex_store_associated_edges(true);
	}
	return m_aaEdgeContainerVERTEX[vrt].end();
}

Grid::AssociatedEdgeIterator Grid::associated_edges_begin(Face* face)
{
	if(!option_is_enabled(FaceOptions::FACEOPT_STORE_ASSOCIATED_EDGES))
	{
		UG_LOG("WARNING in associated_edges_begin(face): auto-enabling FACEOPT_STORE_ASSOCIATED_EDGES." << endl);
		face_store_associated_edges(true);
	}
	return m_aaEdgeContainerFACE[face].begin();
}

Grid::AssociatedEdgeIterator Grid::associated_edges_end(Face* face)
{
	if(!option_is_enabled(FaceOptions::FACEOPT_STORE_ASSOCIATED_EDGES))
	{
		UG_LOG("WARNING in associated_edges_end(face): auto-enabling FACEOPT_STORE_ASSOCIATED_EDGES." << endl);
		face_store_associated_edges(true);
	}
	return m_aaEdgeContainerFACE[face].end();
}

Grid::AssociatedEdgeIterator Grid::associated_edges_begin(Volume* vol)
{
	if(!option_is_enabled(VolumeOptions::VOLOPT_STORE_ASSOCIATED_EDGES))
	{
		UG_LOG("WARNING in associated_edges_begin(vol): auto-enabling VOLOPT_STORE_ASSOCIATED_EDGES." << endl);
		volume_store_associated_edges(true);
	}
	return m_aaEdgeContainerVOLUME[vol].begin();
}

Grid::AssociatedEdgeIterator Grid::associated_edges_end(Volume* vol)
{
	if(!option_is_enabled(VolumeOptions::VOLOPT_STORE_ASSOCIATED_EDGES))
	{
		UG_LOG("WARNING in associated_edges_end(vol): auto-enabling VOLOPT_STORE_ASSOCIATED_EDGES." << endl);
		volume_store_associated_edges(true);
	}
	return m_aaEdgeContainerVOLUME[vol].end();
}

////////////////////////////////////////////////////////////////////////
//	associated face access
Grid::AssociatedFaceIterator Grid::associated_faces_begin(Vertex* vrt)
{
	if(!option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_FACES))
	{
		UG_LOG("WARNING in associated_faces_begin(vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_FACES." << endl);
		vertex_store_associated_faces(true);
	}
	return m_aaFaceContainerVERTEX[vrt].begin();
}

Grid::AssociatedFaceIterator Grid::associated_faces_end(Vertex* vrt)
{
	if(!option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_FACES))
	{
		UG_LOG("WARNING in associated_faces_end(vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_FACES." << endl);
		vertex_store_associated_faces(true);
	}
	return m_aaFaceContainerVERTEX[vrt].end();
}

Grid::AssociatedFaceIterator Grid::associated_faces_begin(Edge* edge)
{
	if(!option_is_enabled(EdgeOptions::EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		UG_LOG("WARNING in associated_faces_begin(edge): auto-enabling EDGEOPT_STORE_ASSOCIATED_FACES." << endl);
		edge_store_associated_faces(true);
	}
	return m_aaFaceContainerEDGE[edge].begin();
}

Grid::AssociatedFaceIterator Grid::associated_faces_end(Edge* edge)
{
	if(!option_is_enabled(EdgeOptions::EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		UG_LOG("WARNING in associated_faces_end(edge): auto-enabling EDGEOPT_STORE_ASSOCIATED_FACES." << endl);
		edge_store_associated_faces(true);
	}
	return m_aaFaceContainerEDGE[edge].end();
}

Grid::AssociatedFaceIterator Grid::associated_faces_begin(Volume* vol)
{
	if(!option_is_enabled(VolumeOptions::VOLOPT_STORE_ASSOCIATED_FACES))
	{
		UG_LOG("WARNING in associated_faces_begin(vol): auto-enabling VOLOPT_STORE_ASSOCIATED_FACES." << endl);
		volume_store_associated_faces(true);
	}
	return m_aaFaceContainerVOLUME[vol].begin();
}

Grid::AssociatedFaceIterator Grid::associated_faces_end(Volume* vol)
{
	if(!option_is_enabled(VolumeOptions::VOLOPT_STORE_ASSOCIATED_FACES))
	{
		UG_LOG("WARNING in associated_faces_end(vol): auto-enabling VOLOPT_STORE_ASSOCIATED_FACES." << endl);
		volume_store_associated_faces(true);
	}
	return m_aaFaceContainerVOLUME[vol].end();
}

////////////////////////////////////////////////////////////////////////
//	associated volume access
Grid::AssociatedVolumeIterator Grid::associated_volumes_begin(Vertex* vrt)
{
	if(!option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
		UG_LOG("WARNING in associated_volumes_begin(vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_VOLUMES." << endl);
		vertex_store_associated_volumes(true);
	}
	return m_aaVolumeContainerVERTEX[vrt].begin();
}

Grid::AssociatedVolumeIterator Grid::associated_volumes_end(Vertex* vrt)
{
	if(!option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
		UG_LOG("WARNING in associated_volumes_end(vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_VOLUMES." << endl);
		vertex_store_associated_volumes(true);
	}
	return m_aaVolumeContainerVERTEX[vrt].end();
}

Grid::AssociatedVolumeIterator Grid::associated_volumes_begin(Edge* edge)
{
	if(!option_is_enabled(EdgeOptions::EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		UG_LOG("WARNING in associated_volumes_begin(edge): auto-enabling EDGEOPT_STORE_ASSOCIATED_VOLUMES." << endl);
		edge_store_associated_volumes(true);
	}
	return m_aaVolumeContainerEDGE[edge].begin();
}

Grid::AssociatedVolumeIterator Grid::associated_volumes_end(Edge* edge)
{
	if(!option_is_enabled(EdgeOptions::EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		UG_LOG("WARNING in associated_volumes_end(edge): auto-enabling EDGEOPT_STORE_ASSOCIATED_VOLUMES." << endl);
		edge_store_associated_volumes(true);
	}
	return m_aaVolumeContainerEDGE[edge].end();
}

Grid::AssociatedVolumeIterator Grid::associated_volumes_begin(Face* face)
{
	if(!option_is_enabled(FaceOptions::FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		UG_LOG("WARNING in associated_volumes_begin(face): auto-enabling FACEOPT_STORE_ASSOCIATED_VOLUMES." << endl);
		face_store_associated_volumes(true);
	}
	return m_aaVolumeContainerFACE[face].begin();
}

Grid::AssociatedVolumeIterator Grid::associated_volumes_end(Face* face)
{
	if(!option_is_enabled(FaceOptions::FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		UG_LOG("WARNING in associated_volumes_end(face): auto-enabling FACEOPT_STORE_ASSOCIATED_VOLUMES." << endl);
		face_store_associated_volumes(true);
	}
	return m_aaVolumeContainerFACE[face].end();
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	neighbourhood access
Edge* Grid::get_edge(Vertex* v1, Vertex* v2)
{
	EdgeDescriptor ed(v1, v2);
	return find_edge_in_associated_edges(v1, ed);
}

Edge* Grid::get_edge(const EdgeVertices& ev)
{
	return find_edge_in_associated_edges(ev.vertex(0), ev);
}

Edge* Grid::get_edge(Face* f, int ind)
{
//	check whether the face stores associated edges
	if(option_is_enabled(FaceOptions::FACEOPT_STORE_ASSOCIATED_EDGES))
	{
		if(option_is_enabled(FaceOptions::FACEOPT_AUTOGENERATE_EDGES))
			return m_aaEdgeContainerFACE[f][ind];
		else{
			EdgeDescriptor ed;
			f->edge_desc(ind, ed);
			return find_edge_in_associated_edges(f, ed);
		}
	}
	else
	{
	//	get the descriptor of the i-th edge
		EdgeDescriptor ed;
		f->edge_desc(ind, ed);
	//	it doesn't. find the edge by checking vertices.
		return find_edge_in_associated_edges(ed.vertex(0), ed);
	}
	
	return nullptr;
}

Edge* Grid::get_edge(Volume* v, int ind)
{
//	check whether the face stores associated edges
	if(option_is_enabled(VolumeOptions::VOLOPT_STORE_ASSOCIATED_EDGES))
	{
	//	if autogenerate is enabeld, edges are sorted.
		if(option_is_enabled(VolumeOptions::VOLOPT_AUTOGENERATE_EDGES)
			|| option_is_enabled(VolumeOptions::VOLOPT_AUTOGENERATE_FACES
								| FaceOptions::FACEOPT_AUTOGENERATE_EDGES))
		{
			return m_aaEdgeContainerVOLUME[v][ind];
		}
		else{
			EdgeDescriptor ed;
			v->edge_desc(ind, ed);
			return find_edge_in_associated_edges(v, ed);
		}
	}
	else
	{
	//	get the descriptor of the i-th edge
		EdgeDescriptor ed;
		v->edge_desc(ind, ed);
	//	it doesn't. find the edge by checking vertices.
		return find_edge_in_associated_edges(ed.vertex(0), ed);
	}

	return nullptr;
}

Face* Grid::get_face(const FaceVertices& fv)
{
	return find_face_in_associated_faces(fv.vertex(0), fv);
}

Face* Grid::get_face(Volume* v, int ind)
{
//	check whether the volume stores associated faces
	if(option_is_enabled(VolumeOptions::VOLOPT_STORE_ASSOCIATED_FACES))
	{
	//	if autogenerate is enabeld, faces are sorted.
		if(option_is_enabled(VolumeOptions::VOLOPT_AUTOGENERATE_FACES))
			return m_aaFaceContainerVOLUME[v][ind];
		else{
			FaceDescriptor fd;
			v->face_desc(ind, fd);
			return find_face_in_associated_faces(v, fd);
		}
	}
	else {
		FaceDescriptor fd;
		v->face_desc(ind, fd);
	//	it does not. check associated faces of the first vertex of fd.
		return find_face_in_associated_faces(fd.vertex(0), fd);
	}
	return nullptr;
}

Volume* Grid::get_volume(const VolumeVertices& vv)
{
	return find_volume_in_associated_volumes(vv.vertex(0), vv);
}

////////////////////////////////////////////////////////////////////////
//	sides
Vertex::side* Grid::get_side(Vertex* obj, size_t side)
{
	GRID_PROFILE_FUNC();
	assert(!"ERROR in Grid::get_side(Vertex*, ...): A vertex doesn't have sides!");
	return nullptr;
}

Edge::side* Grid::get_side(Edge* obj, size_t side)
{
	GRID_PROFILE_FUNC();
	assert(side < 2 && "ERROR in Grid::get_side(Edge*, ...): Bad side index!");
	return obj->vertex(side);
}

Face::side* Grid::get_side(Face* obj, size_t side)
{
	GRID_PROFILE_FUNC();
	assert(side < obj->num_edges() && "ERROR in Grid::get_side(Face*, ...): Bad side index!");
	return get_edge(obj, side);
}

Volume::side* Grid::get_side(Volume* obj, size_t side)
{
	GRID_PROFILE_FUNC();
	assert(side < obj->num_faces() && "ERROR in Grid::get_side(Volume*, ...): Bad side index!");
	return get_face(obj, side);
}

////////////////////////////////////////////////////////////////////////
//	marks
void Grid::init_marks()
{
//	attach marks to the elements
	if(m_currentMark == 0)
	{
	//	marks have not yet been initialized - do that now
	//	attach m_currentMark with default value 0
	//	(0 is never the currentMark while marks are active).
		m_currentMark = 1;
		attach_to_vertices_dv(m_aMark, 0);
		attach_to_edges_dv(m_aMark, 0);
		attach_to_faces_dv(m_aMark, 0);
		attach_to_volumes_dv(m_aMark, 0);
		
		m_aaMarkVRT.access(*this, m_aMark);
		m_aaMarkEDGE.access(*this, m_aMark);
		m_aaMarkFACE.access(*this, m_aMark);
		m_aaMarkVOL.access(*this, m_aMark);

		m_bMarking = false;
	}
}

void Grid::reset_marks()
{
//	set all marks to 0 and m_currentMark to 1
	m_currentMark = 1;

//	reset vertex marks
	AMark::ContainerType *pContainer = get_attachment_data_container<Vertex>(m_aMark);
	for(uint i = 0; i < pContainer->size(); ++i)
		pContainer->get_elem(i) = 0;

//	reset edge marks
	pContainer = get_attachment_data_container<Edge>(m_aMark);
	for(uint i = 0; i < pContainer->size(); ++i)
		pContainer->get_elem(i) = 0;

//	reset face marks
	pContainer = get_attachment_data_container<Face>(m_aMark);
	for(uint i = 0; i < pContainer->size(); ++i)
		pContainer->get_elem(i) = 0;

//	reset volume marks
	pContainer = get_attachment_data_container<Volume>(m_aMark);
	for(uint i = 0; i < pContainer->size(); ++i)
		pContainer->get_elem(i) = 0;	
}

void Grid::remove_marks()
{
	if(m_currentMark != 0)
	{
		m_currentMark = 0;
		detach_from_vertices(m_aMark);
		detach_from_edges(m_aMark);
		detach_from_faces(m_aMark);
		detach_from_volumes(m_aMark);
	}
}

void Grid::begin_marking()
{
	if(m_currentMark == 0)
	{
	//	marks are disabled. we have to activate them
		init_marks();
	}
	
	if(m_bMarking){
		throw(UGError("ERROR in Grid::begin_marking(): marking is already active. Don't forget to call end_marking when you're done with marking."));
	}
	
//	increase currentMark
	++m_currentMark;
	
//	check whether we have to reset-marks
	if(m_currentMark == -1)
		reset_marks();
		
//	set m_bMarking to true
	m_bMarking = true;
}

void Grid::end_marking()
{
	m_bMarking = false;
}

void Grid::clear_marks()
{
	if(m_bMarking){
		end_marking();
		begin_marking();
	}
	else{
		begin_marking();
		end_marking();
	}
}

void Grid::test_attached_linked_lists()
{
	UG_LOG("empty\n");
}

}//	end of namespace
