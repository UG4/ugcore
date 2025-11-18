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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/*
 * In this file the the grids connectivity-management is implemented.
 * The methods change_vertex_options(...), change_edge_options, ...
 * attach and detach Containers to vertices, edges, ... in order to
 * store associated elements, as specified in the passed options.
 * The generated data can be accessed through Grids
 * get_associated_edges, get_associated_faces, ... methods later on.
 *
 * New edges, faces and volumes have to be registered with other objects
 * directly after creation. For this task the methods
 * register_edge, register_face, ... are supplied.
 * On the other hand Grid has to update connectivity information, as
 * soon as grid-elements are removed. for these tasks the methods
 * unregister_edge, unregister_face, ... are supplied.
 */
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


#include <algorithm>
#include "grid.h"
#include "grid_util.h"
#include "common/common.h"
#include "common/profiler/profiler.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
//define PROFILE_GRID_CONNECTION_MANAGMENT if you want to profile
//the code below.
//#define PROFILE_GRID_CONNECTION_MANAGMENT
#ifdef PROFILE_GRID_CONNECTION_MANAGMENT
	#define GCM_PROFILE_FUNC()	PROFILE_FUNC()
	#define GCM_PROFILE(name)	PROFILE_BEGIN(name)
	#define GCM_PROFILE_END()	PROFILE_END()
#else
	#define GCM_PROFILE_FUNC()
	#define GCM_PROFILE(name)
	#define GCM_PROFILE_END()
#endif


////////////////////////////////////////////////////////////////////////
///	this macro helps calling callbacks of different observers.
/**
 * Be sure that callback is a complete function call - including parameters.
 */
#define NOTIFY_OBSERVERS(observerContainer, callback)	{for(Grid::ObserverContainer::iterator iter = observerContainer.begin(); iter != observerContainer.end(); iter++) (*iter)->callback;}
#define NOTIFY_OBSERVERS_REVERSE(observerContainer, callback)	{for(Grid::ObserverContainer::reverse_iterator iter = observerContainer.rbegin(); iter != observerContainer.rend(); iter++) (*iter)->callback;}

////////////////////////////////////////////////////////////////////////
///	a useful macro that checks if a set of options contains the specified option.
#define OPTIONS_CONTAIN_OPTION(options, option) (((options) & (option)) == (option))


namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	VERTICES
///	creates and removes connectivity data, as specified in optsNew.
void Grid::register_vertex(Vertex* v, GridObject* pParent)
{
	GCM_PROFILE_FUNC();

//	store the element and register it at the pipe.
	m_vertexElementStorage.m_attachmentPipe.register_element(v);
	m_vertexElementStorage.m_sectionContainer.insert(v, v->container_section());

//	assign the hash-value
	assign_hash_value(v);

	GCM_PROFILE(GCM_notify_vertex_observers);
//	inform observers about the creation
	NOTIFY_OBSERVERS(m_vertexObservers, vertex_created(this, v, pParent));
	GCM_PROFILE_END();
}

void Grid::register_and_replace_element(Vertex* v, Vertex* pReplaceMe)
{
	m_vertexElementStorage.m_attachmentPipe.register_element(v);
	m_vertexElementStorage.m_sectionContainer.insert(v, v->container_section());

//	assign the hash-value
	assign_hash_value(v);

//	pass on values
	pass_on_values(pReplaceMe, v);

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_vertexObservers, vertex_created(this, v, pReplaceMe, true));
//	inform observers about the deletion
	NOTIFY_OBSERVERS_REVERSE(m_vertexObservers,
							 vertex_to_be_erased(this, pReplaceMe, v));

//	all edges, faces and volumes associated with pReplaceMe have to be updated.
//	the options in GRIDOPT_VERTEXCENTRIC_INTERCONNECTION have to be enabled.
//	we do this per element only, to avoid unnecessary memory overhead.

//	update edges
	if(num_edges()){
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES)){
			LOG("WARNING in Grid::register_and_replace_element(...) - Vertex: autoenabling grid-option VRTOPT_STORE_ASSOCIATED_EDGES.");
			enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
		}

		for(AssociatedEdgeIterator iter = associated_edges_begin(pReplaceMe);
			iter != associated_edges_end(pReplaceMe); ++iter)
		{
			Edge* e = *iter;
		//	replace the vertex and push e into v's associated edges.
			for(uint i = 0; i < 2; ++i)
			{
				if(e->vertex(i) == pReplaceMe)
					e->set_vertex(i, v);
			}

			m_aaEdgeContainerVERTEX[v].push_back(e);
		}
	}

//	update faces
	if(num_faces()){
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES)){
			LOG("WARNING in Grid::register_and_replace_element(...) - Vertex: autoenabling grid-option VRTOPT_STORE_ASSOCIATED_FACES.");
			enable_options(VRTOPT_STORE_ASSOCIATED_FACES);
		}

		for(AssociatedFaceIterator iter = associated_faces_begin(pReplaceMe);
			iter != associated_faces_end(pReplaceMe); ++iter)
		{
			Face* f = *iter;
		//	replace the vertex and push f into v's associated faces.
			uint numVrts = f->num_vertices();
			Face::ConstVertexArray vrts = f->vertices();
			for(uint i = 0; i < numVrts; ++i)
			{
				if(vrts[i] == pReplaceMe)
					f->set_vertex(i, v);
			}

			m_aaFaceContainerVERTEX[v].push_back(f);
		}
	}

//	update volumes
	if(num_volumes()){
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES)){
			LOG("WARNING in Grid::register_and_replace_element(...) - Vertex: autoenabling grid-option VRTOPT_STORE_ASSOCIATED_VOLUMES.");
			enable_options(VRTOPT_STORE_ASSOCIATED_VOLUMES);
		}

		for(AssociatedVolumeIterator iter = associated_volumes_begin(pReplaceMe);
			iter != associated_volumes_end(pReplaceMe); ++iter)
		{
			Volume* vol = *iter;
		//	replace the vertex and push vol into v's associated volumes.
			uint numVrts = vol->num_vertices();
			Volume::ConstVertexArray vrts = vol->vertices();
			for(uint i = 0; i < numVrts; ++i)
			{
				if(vrts[i] == pReplaceMe)
					vol->set_vertex(i, v);
			}

			m_aaVolumeContainerVERTEX[v].push_back(vol);
		}
	}

//	clear the containers of pReplaceMe
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
		m_aaEdgeContainerVERTEX[pReplaceMe].clear();
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
		m_aaFaceContainerVERTEX[pReplaceMe].clear();
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
		m_aaVolumeContainerVERTEX[pReplaceMe].clear();

//	remove pReplaceMe
	m_vertexElementStorage.m_attachmentPipe.unregister_element(pReplaceMe);
	m_vertexElementStorage.m_sectionContainer.erase(get_iterator(pReplaceMe), pReplaceMe->container_section());
	delete pReplaceMe;
}

void Grid::unregister_vertex(Vertex* v)
{
//	notify observers that the vertex is being erased
	NOTIFY_OBSERVERS_REVERSE(m_vertexObservers, vertex_to_be_erased(this, v));

//	perform some checks in order to assert grid consistency.
//	all edges, faces and volume referencing this vertex have to be erased.

	if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
	//	if there are edges this option has to be enabled!
		if(num_edges() > 0)
		{
			LOG("WARNING in Grid::unregister_vertex(...): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES." << std::endl);
			vertex_store_associated_edges(true);
		}
	}

	if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
	{
	//	if there are faces we have to consider the following
		if(num_faces() > 0)
		{
		//	if edges store associated faces and if faces auto-generate edges, then
		//	nothing has to be performed here.
			if(!(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES)
				 &&	option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)))
			{
			//	since adjacent faces have to be removed, we have to enable VRTOPT_STORE_ASSOCIATED_FACES
				LOG("WARNING in Grid::unregister_vertex(...): auto-enabling VRTOPT_STORE_ASSOCIATED_FACES." << std::endl);
				vertex_store_associated_faces(true);
			}
		}
	}

	if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
	//	if there are volumes we have to consider the following
		if(num_volumes() > 0)
		{
		//	if edges store associated volumes and volumes auto-generate edges or
		//	if faces store associated volumes and volumes auto-generate faces then
		//	nothing has to be performed here.
			bool cond1 =	option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES) &&
							option_is_enabled(VOLOPT_AUTOGENERATE_EDGES);
			bool cond2 =	option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES) &&
							option_is_enabled(VOLOPT_AUTOGENERATE_FACES);

			if(!(cond1 || cond2))
			{
			//	since we have to remove adjacent volumes we have to enable VRTOPT_STORE_ASSOCIATED_VOLUMES
				LOG("WARNING in Grid::unregister_vertex(...): auto-enabling VRTOPT_STORE_ASSOCIATED_VOLUMES." << std::endl);
				vertex_store_associated_volumes(true);
			}
		}
	}

//	remove associated volumes
	if(num_volumes() > 0)
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
		{
		//	remove associated volumes. Make sure that there are no problems
		//	when volumes try to unregister from the vertex.
			VolumeContainer vols;
			vols.swap(m_aaVolumeContainerVERTEX[v]);

			for(auto iter = vols.begin(); iter != vols.end();)
			{
				Volume* eraseVol = *iter;
				++iter;
				erase(eraseVol);
			}
		}
	}

//	remove associated faces
	if(num_faces() > 0)
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
		{
		//	remove all associated faces. Make sure that there are no problems
		//	when faces try to unregister from the vertex.
			FaceContainer faces;
			faces.swap(m_aaFaceContainerVERTEX[v]);
			for(auto iter = faces.begin(); iter != faces.end();)
			{
				Face* eraseFace = *iter;
				++iter;
				erase(eraseFace);
			}
		}
	}

//	remove associated edges
	if(num_edges() > 0)
	{
		assert(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES) && "unexpected internal error in Grid::unregister_vertex - rae.");
	//	remove all associated edges. Make sure that there are no problems
	//	when edges try to unregister from the vertex.
		EdgeContainer edges;
		edges.swap(m_aaEdgeContainerVERTEX[v]);
		for(auto iter = edges.begin(); iter != edges.end();)
		{
			Edge* eraseEdge = *iter;
			++iter;
			erase(eraseEdge);
		}
	}

//	remove the element from the storage
	m_vertexElementStorage.m_sectionContainer.erase(get_iterator(v), v->container_section());
	m_vertexElementStorage.m_attachmentPipe.unregister_element(v);
}

void Grid::change_vertex_options(uint optsNew)
{
//	check if associated edge information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, VRTOPT_STORE_ASSOCIATED_EDGES))
	{
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES)){
			vertex_store_associated_edges(true);
		}
	}
	else
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES)){
			vertex_store_associated_edges(false);
		}
	}

//	check if associated face information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, VRTOPT_STORE_ASSOCIATED_FACES))
	{
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
			vertex_store_associated_faces(true);
	}
	else
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
			vertex_store_associated_faces(false);
	}

//	check if associated volume information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
			vertex_store_associated_volumes(true);
	}
	else
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
			vertex_store_associated_volumes(false);
	}
}

void Grid::vertex_store_associated_edges(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
		{
		//	store associated edges
			attach_to_vertices(m_aEdgeContainer);
			m_aaEdgeContainerVERTEX.access(*this, m_aEdgeContainer);

		//	iterate through all edges and store them in the associated vertices
			for(EdgeIterator iter = this->edges_begin(); iter != this->edges_end(); iter++)
			{
				Edge* e = *iter;
				m_aaEdgeContainerVERTEX[e->vertex(0)].push_back(e);
				m_aaEdgeContainerVERTEX[e->vertex(1)].push_back(e);
			}

			m_options |= VRTOPT_STORE_ASSOCIATED_EDGES;
		}
	}
	else
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
		{
			detach_from_vertices(m_aEdgeContainer);
			m_options &= (~VRTOPT_STORE_ASSOCIATED_EDGES);
		}
	}
}

void Grid::vertex_store_associated_faces(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
		{
		//	store associated faces
			attach_to_vertices(m_aFaceContainer);
			m_aaFaceContainerVERTEX.access(*this, m_aFaceContainer);

		//	iterator through all faces and store them in the referenced vertices
			for(FaceIterator iter = this->faces_begin(); iter != this->faces_end(); iter++)
			{
				Face* f = *iter;
				int numVrts = f->num_vertices();
				Face::ConstVertexArray vrts = f->vertices();
				for(int i = 0; i < numVrts; i++)
					m_aaFaceContainerVERTEX[vrts[i]].push_back(f);
			}

		//	store the option
			m_options |= VRTOPT_STORE_ASSOCIATED_FACES;
		}
	}
	else
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
		{
			detach_from_vertices(m_aFaceContainer);
			m_options &= (~VRTOPT_STORE_ASSOCIATED_FACES);
		}
	}
}

void Grid::vertex_store_associated_volumes(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
		{
		//	store associated volumes
			attach_to_vertices(m_aVolumeContainer);
			m_aaVolumeContainerVERTEX.access(*this, m_aVolumeContainer);

		//	iterator through all faces and store them in the referenced vertices
			for(VolumeIterator iter = this->volumes_begin(); iter != this->volumes_end(); iter++)
			{
				Volume* v = *iter;
				int numVrts = v->num_vertices();
				Volume::ConstVertexArray vrts = v->vertices();
				for(int i = 0; i < numVrts; i++)
					m_aaVolumeContainerVERTEX[vrts[i]].push_back(v);
			}

		//	store the option
			m_options |= VRTOPT_STORE_ASSOCIATED_VOLUMES;
		}
	}
	else
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
		{
			detach_from_vertices(m_aVolumeContainer);
			m_options &= (~VRTOPT_STORE_ASSOCIATED_VOLUMES);
		}
	}
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	EDGES
///	creates and removes connectivity data, as specified in optsNew.
void Grid::register_edge(Edge* e, GridObject* pParent,
						 Face* createdByFace, Volume* createdByVol)
{
	GCM_PROFILE_FUNC();

//	store the element and register it at the pipe.
	m_edgeElementStorage.m_attachmentPipe.register_element(e);
	m_edgeElementStorage.m_sectionContainer.insert(e, e->container_section());

//	register edge at vertices, faces and volumes, if the according options are enabled.
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
		m_aaEdgeContainerVERTEX[e->vertex(0)].push_back(e);
		m_aaEdgeContainerVERTEX[e->vertex(1)].push_back(e);
	}

//	register edges at faces and vice versa
//	this is only necessary, if faces do not autogenerate edges
	if(createdByFace != nullptr){
	//	e has been autogenerated by f. We thus have don't have to check
	//	for other faces. e is already contained in f's edge-container.
		if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
			m_aaFaceContainerEDGE[e].push_back(createdByFace);
	}
	else
	{
		if(!option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
			GCM_PROFILE(GCM_reg_edge_processing_faces);
			int switchVar = 0;	// this var will be used to determine what to register where.
			if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
				switchVar = 1;
			if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
				switchVar += 2;

			if(switchVar > 0)
			{
			//	find the faces that contain the edge. search the associated vertices for connected faces.
				auto iterStart = associated_faces_begin(e->vertex(0));
				auto iterEnd = associated_faces_end(e->vertex(0));

				EdgeDescriptor ed;
				for(auto iter = iterStart; iter != iterEnd; iter++)
				{
					Face* f = *iter;
					uint numEdges = f->num_edges();
					for(uint i = 0; i < numEdges; ++i)
					{
						f->edge_desc(i, ed);
						if(CompareVertices(e, &ed))
						{
						//	f contains e
							switch(switchVar)
							{
								case 1: m_aaFaceContainerEDGE[e].push_back(f);
										break;
								case 2:	m_aaEdgeContainerFACE[f].push_back(e);
										break;
								case 3:	m_aaFaceContainerEDGE[e].push_back(f);
										m_aaEdgeContainerFACE[f].push_back(e);
										break;
							}
						}
					}
				}
			}
		}
	}

//	register edges at volumes and vice versa
	if(createdByVol != nullptr){
	//	e was at least indirectly created by a volume.
	//	We thus do not have to check for other volumes.
	//	e is already registered in v's edge-container.
		if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
			m_aaVolumeContainerEDGE[e].push_back(createdByVol);
	}
	else
	{
		bool ignoreVolumes = option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)
						|| (option_is_enabled(VOLOPT_AUTOGENERATE_FACES)
							&& option_is_enabled(FACEOPT_AUTOGENERATE_EDGES));

		if(!ignoreVolumes){
			GCM_PROFILE(GCM_reg_edge_processing_volumes);
			int switchVar = 0;	// this var will be used to determine what to register where.
			if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
				switchVar = 1;
			if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
				switchVar += 2;

			if(switchVar > 0)
			{
			//	find the volumes that contain the edge.
				auto iterStart = associated_volumes_begin(e->vertex(0));
				auto iterEnd = associated_volumes_end(e->vertex(0));

				for(auto iter = iterStart; iter != iterEnd; iter++)
				{
					Volume* v = *iter;
					uint numEdges = v->num_edges();
					for(uint i = 0; i < numEdges; ++i)
					{
						EdgeDescriptor ed;
						v->edge_desc(i, ed);
						if(CompareVertices(e, &ed))
						{
						//	f contains e
							switch(switchVar)
							{
								case 1: m_aaVolumeContainerEDGE[e].push_back(v);
										break;
								case 2:	m_aaEdgeContainerVOLUME[v].push_back(e);
										break;
								case 3:	m_aaVolumeContainerEDGE[e].push_back(v);
										m_aaEdgeContainerVOLUME[v].push_back(e);
										break;
							}
						}
					}
				}
			}
		}
	}

	GCM_PROFILE(GCM_notify_edge_observers);
//	inform observers about the creation
	NOTIFY_OBSERVERS(m_edgeObservers, edge_created(this, e, pParent));
	GCM_PROFILE_END();
}

void Grid::register_and_replace_element(Edge* e, Edge* pReplaceMe)
{
//	store the element and register it at the pipe.
	m_edgeElementStorage.m_attachmentPipe.register_element(e);
	m_edgeElementStorage.m_sectionContainer.insert(e, e->container_section());

	pass_on_values(pReplaceMe, e);

//	assign vertices
	e->set_vertex(0, pReplaceMe->vertex(0));
	e->set_vertex(1, pReplaceMe->vertex(1));

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_edgeObservers, edge_created(this, e, pReplaceMe, true));
//	inform observers about the deletion
	NOTIFY_OBSERVERS_REVERSE(m_edgeObservers, edge_to_be_erased(this, pReplaceMe, e));

//	check if vertices, faces and volumes reference pReplaceMe.
//	if so, correct those references.
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
		for(uint i = 0; i < 2; ++i)
			replace(associated_edges_begin(e->vertex(i)),
					associated_edges_end(e->vertex(i)),
					pReplaceMe, e);
	}

	if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
	{
	//	collect all faces that are associated with pReplaceMe
		vector<Face*> vFaces;
		CollectFaces(vFaces, *this, pReplaceMe);
		for(uint i = 0; i < vFaces.size(); ++i)
			replace(associated_edges_begin(vFaces[i]),
					associated_edges_end(vFaces[i]),
					pReplaceMe, e);
	}

	if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
	{
	//	collect all volumes that are associated with pReplaceMe
		vector<Volume*> vVolumes;
		CollectVolumes(vVolumes, *this, pReplaceMe);
		for(uint i = 0; i < vVolumes.size(); ++i)
			replace(associated_edges_begin(vVolumes[i]),
					associated_edges_end(vVolumes[i]),
					pReplaceMe, e);
	}

//	now we have to copy all associated elements of pReplaceMe to e
	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		m_aaFaceContainerEDGE[e].assign(associated_faces_begin(pReplaceMe),
										associated_faces_end(pReplaceMe));
		m_aaFaceContainerEDGE[pReplaceMe].clear();
	}

	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		m_aaVolumeContainerEDGE[e].assign(associated_volumes_begin(pReplaceMe),
										associated_volumes_end(pReplaceMe));
		m_aaVolumeContainerEDGE[pReplaceMe].clear();
	}

//	remove the element from the storage and delete it.
	m_edgeElementStorage.m_sectionContainer.erase(get_iterator(pReplaceMe), pReplaceMe->container_section());
	m_edgeElementStorage.m_attachmentPipe.unregister_element(pReplaceMe);
	delete pReplaceMe;
}

void Grid::unregister_edge(Edge* e)
{
//	notify observers that the edge is being erased
	NOTIFY_OBSERVERS_REVERSE(m_edgeObservers, edge_to_be_erased(this, e));

//	delete associated faces or unregister from associated faces
	if(num_volumes() > 0)
	{
		if(option_is_enabled(VOLOPT_AUTOGENERATE_EDGES) ||
			option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
		{
		//	check if anything has to be performed at all.
			if(!(option_is_enabled(VOLOPT_AUTOGENERATE_FACES) &&
				 option_is_enabled(FACEOPT_AUTOGENERATE_EDGES) &&
				 option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES) &&
				 option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES)))
			{
			//	we have to perform something...
				vector<Volume*> vVolumes;
				CollectVolumes(vVolumes, *this, e);
				auto vIter = vVolumes.begin();

				if(option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
				{
					if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
						m_aaVolumeContainerEDGE[e].clear();

				//	erase the collected volumes
					while(vIter != vVolumes.end())
					{
						Volume* eraseVol = *vIter;
						++vIter;
						erase(eraseVol);
					}
				}
				else
				{
				//	unlink from the collected volumes
					for(; vIter != vVolumes.end(); ++vIter)
					{
						Volume* v = *vIter;
					//	find the edge and remove the corresponding entry
						auto iter = find(m_aaEdgeContainerVOLUME[v].begin(),
															m_aaEdgeContainerVOLUME[v].end(), e);
						if(iter != m_aaEdgeContainerVOLUME[v].end())
							m_aaEdgeContainerVOLUME[v].erase(iter);
					}
				}
			}
		}
	}

//	delete associated faces or unregister from associated faces
	if(num_faces() > 0)
	{
		if(option_is_enabled(FACEOPT_AUTOGENERATE_EDGES) ||
			option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
		{
			vector<Face*> vFaces;
			CollectFaces(vFaces, *this, e);
			auto fIter = vFaces.begin();

			if(option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
			{
				if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
					m_aaFaceContainerEDGE[e].clear();

			//	erase the collected faces.
				while(fIter != vFaces.end())
				{
					Face* eraseFace = *fIter;
					++fIter;
					erase(eraseFace);
				}
			}
			else
			{
			//	unlink
				for(; fIter != vFaces.end(); ++fIter)
				{
					Face* f = *fIter;
				//	find the matching entry and remove it
					auto iter = find(m_aaEdgeContainerFACE[f].begin(),
														m_aaEdgeContainerFACE[f].end(), e);
					if(iter != m_aaEdgeContainerFACE[f].end())
						m_aaEdgeContainerFACE[f].erase(iter);
				}
			}
		}
	}

//	unregister from vertices
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
	//	iterate through the associated vertices and remove the edge from their edge-list.
		for(int i = 0; i < 2; ++i)
		{
			Vertex* vrt = e->vertex(i);
			auto iter = find(m_aaEdgeContainerVERTEX[vrt].begin(),
												m_aaEdgeContainerVERTEX[vrt].end(), e);
			if(iter != m_aaEdgeContainerVERTEX[vrt].end())
				m_aaEdgeContainerVERTEX[e->vertex(i)].erase(iter);
		}
	}

//	remove the element from the storage
	m_edgeElementStorage.m_sectionContainer.erase(get_iterator(e), e->container_section());
	m_edgeElementStorage.m_attachmentPipe.unregister_element(e);
}

void Grid::change_edge_options(uint optsNew)
{
//	check if associated face information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		if(!option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
			edge_store_associated_faces(true);
	}
	else
	{
		if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
			edge_store_associated_faces(false);
	}

//	check if associated volume information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		if(!option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
			edge_store_associated_volumes(true);
	}
	else
	{
		if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
			edge_store_associated_volumes(false);
	}
}

void Grid::edge_store_associated_faces(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
		{
		//	store associated faces
			attach_to_edges(m_aFaceContainer);
			m_aaFaceContainerEDGE.access(*this, m_aFaceContainer);

		//	iterate through the faces and store them at each of their edges.
		//	in order to perform this task at acceptable speed,
		//	the option VRTOPT_STORE_ASSOCIATED_EDGES has to be enabled.
		//	(This is due to the use of GetEdge(...))
		//	if it is disabled here, we'll enable it.
			if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
			{
				LOG("WARNING in edge_store_associated_faces(...): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES." << endl);
				vertex_store_associated_edges(true);
			}

			for(FaceIterator iter = faces_begin(); iter != faces_end(); iter++)
			{
				Face* f = *iter;
			//	iterate through the edges of the face
				int numEdges = f->num_edges();
				for(int i = 0; i < numEdges; ++i)
				{
				//	get the i-th edge
					Edge* e = get_edge(f, i);
					if(e != nullptr)
						m_aaFaceContainerEDGE[e].push_back(f);
				}
			}
		//	store the option
			m_options |= EDGEOPT_STORE_ASSOCIATED_FACES;
		}
	}
	else
	{
		if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
		{
		//	detach edge-face connections
			detach_from_edges(m_aFaceContainer);
			m_options &= (~EDGEOPT_STORE_ASSOCIATED_FACES);
		}
	}
}

void Grid::edge_store_associated_volumes(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
		{
		//	store associated faces
			attach_to_edges(m_aVolumeContainer);
			m_aaVolumeContainerEDGE.access(*this, m_aVolumeContainer);

		//	iterate through the volumes and store them at each of their edges.
		//	in order to perform this task at acceptable speed,
		//	the option VRTOPT_STORE_ASSOCIATED_VOLUMES has to be enabled.
		//	(This is due to the use of GetEdge(...))
		//	if it is disabled here, we'll enable it.
			if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
			{
				LOG("WARNING in Grid::edge_store_associated_volumes(...): auto-enabling VRTOPT_STORE_ASSOCIATED_VOLUMES." << endl);
				vertex_store_associated_volumes(true);
			}

			for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
			{
				Volume* v = *iter;
			//	iterate through the edges of the volume
				int numEdges = v->num_edges();
				for(int i = 0; i < numEdges; ++i)
				{
				//	get the edge-descriptor
					EdgeDescriptor ed;
					v->edge_desc(i, ed);
				//	get the edge that is described by the EdgeDescriptor - if it exists at all.
					Edge* e = get_edge(ed);
					if(e != nullptr)
						m_aaVolumeContainerEDGE[e].push_back(v);
				}
			}

		//	store the option
			m_options |= EDGEOPT_STORE_ASSOCIATED_VOLUMES;
		}
	}
	else
	{
		if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
		{
		//	detach edge-volume connections
			detach_from_edges(m_aVolumeContainer);
			m_options &= (~EDGEOPT_STORE_ASSOCIATED_VOLUMES);
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	FACES
///	creates and removes connectivity data, as specified in optsNew.
void Grid::register_face(Face* f, GridObject* pParent, Volume* createdByVol)
{
	GCM_PROFILE_FUNC();

//	store the element and register it at the pipe.
	m_faceElementStorage.m_attachmentPipe.register_element(f);
	m_faceElementStorage.m_sectionContainer.insert(f, f->container_section());

//	register face at vertices
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
	{
		uint numVrts = f->num_vertices();
		Face::ConstVertexArray vrts = f->vertices();
		for(uint i = 0; i < numVrts; ++i)
			m_aaFaceContainerVERTEX[vrts[i]].push_back(f);
	}

	bool createEdges = option_is_enabled(FACEOPT_AUTOGENERATE_EDGES);;
	const bool edgesStoreFaces = option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES);
	const bool facesStoreEdges = option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES);

//	make sure that the edge container is big enough to hold all edges
	if(facesStoreEdges)
		m_aaEdgeContainerFACE[f].reserve(f->num_edges());

//	create edges if FACEOPT_AUTOGENERATE_EDGES is enabled.
//	register face at associated vertices and edges if the according options are set.
	if(createEdges
		|| edgesStoreFaces
		|| facesStoreEdges)
	{
		GCM_PROFILE(GCM_reg_face_processing_edges);
	//	loop through the edges of the face and check if they already exist in the grid.
		int numEdges = f->num_edges();
		EdgeDescriptor ed;
		for(int i = 0; i < numEdges; ++i)
		{
			f->edge_desc(i, ed);
			Edge* e = find_edge_in_associated_edges(ed.vertex(0), ed);

			if(e == nullptr)
			{
				if(createEdges)
				{
				//	create the edge - regard the parent of f as the parent of the new edge, too.
					e = f->create_edge(i);

					if(facesStoreEdges)
						m_aaEdgeContainerFACE[f].push_back(e);

				//	register the edge
					register_edge(e, pParent, f, createdByVol);
				}
			}
			else
			{
			//	register the face at the edge
				if(edgesStoreFaces)
					m_aaFaceContainerEDGE[e].push_back(f);
			//	register the edge at the face.
				if(facesStoreEdges)
					m_aaEdgeContainerFACE[f].push_back(e);
			}
		}
	}

//	register face at existing volumes. find all volumes that contain face
//	if the option VOLOPT_AUTOGENERATE_FACES is enabled, we can ignore volumes
//	here anyways, since they would be associated with volumes in register_volume
	if(createdByVol != nullptr){
	//	This means that the face was autogenerated by a volume.
	//	In this case the face is already registered in the volumes container.
		if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
			m_aaVolumeContainerFACE[f].push_back(createdByVol);
	}
	else
	{
		if(!option_is_enabled(VOLOPT_AUTOGENERATE_FACES)){
			GCM_PROFILE(GCM_reg_face_processing_volumes);

			int switchVar = 0;
			if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
				switchVar = 1;
			if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
				switchVar += 2;

			if(switchVar > 0)
			{
			//	find associated volumes of f. Don't use get_associated, since it
			//	may lookup m_aaVolumeContainerFACE, which we just want to update...
				vector<Volume*> vVols;
				vVols.reserve(2);

				CollectVolumes(vVols, *this, f, false, true);



				for(auto iter = vVols.begin(); iter != vVols.end(); ++iter)
				{
					Volume* v = *iter;
				//	v contains f
					switch(switchVar)
					{
						case 1: m_aaVolumeContainerFACE[f].push_back(v);
								break;
						case 2: m_aaFaceContainerVOLUME[v].push_back(f);
								break;
						case 3:	m_aaVolumeContainerFACE[f].push_back(v);
								m_aaFaceContainerVOLUME[v].push_back(f);
								break;
					}
				}
			}
		}
	}

	GCM_PROFILE(GCM_notify_face_observers);
//	inform observers about the creation
	NOTIFY_OBSERVERS(m_faceObservers, face_created(this, f, pParent));
	GCM_PROFILE_END();
}

void Grid::register_and_replace_element(Face* f, Face* pReplaceMe)
{
//	check that f and pReplaceMe have the same amount of vertices.
	if(f->num_vertices() != pReplaceMe->num_vertices())
	{
		LOG("WARNING in Grid::register_and_replace_element(Face* f, Face* pReplaceMe): f and pReplaceMe have different numbers of vertices.");
		LOG("  f has not been registered. pReplaceMe has not been replaced.");
		assert(!"ERROR in Grid::register_and_replace_element(Face* f, Face* pReplaceMe): f and pReplaceMe have different numbers of vertices.");
		return;
	}

//	store the element and register it at the pipe.
	m_faceElementStorage.m_attachmentPipe.register_element(f);
	m_faceElementStorage.m_sectionContainer.insert(f, f->container_section());

	pass_on_values(pReplaceMe, f);

	Face::ConstVertexArray vrts = pReplaceMe->vertices();

//	assign vertices
	uint numVrts = f->num_vertices();
	{
		for(uint i = 0; i < numVrts; ++i)
			f->set_vertex(i, vrts[i]);
	}

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_faceObservers, face_created(this, f, pReplaceMe, true));
//	inform observers about the deletion
	NOTIFY_OBSERVERS_REVERSE(m_faceObservers, face_to_be_erased(this, pReplaceMe, f));

//	check if vertices, edges and volumes reference pReplaceMe.
//	if so, correct those references.
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
	{
		for(uint i = 0; i < numVrts; ++i)
			replace(associated_faces_begin(vrts[i]),
					associated_faces_end(vrts[i]),
					pReplaceMe, f);
	}

	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
	//	collect all edges that are associated with pReplaceMe
		vector<Edge*> vEdges;
		CollectEdges(vEdges, *this, pReplaceMe);
		for(uint i = 0; i < vEdges.size(); ++i)
			replace(associated_faces_begin(vEdges[i]),
					associated_faces_end(vEdges[i]),
					pReplaceMe, f);
	}

	if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
	{
	//	collect all volumes that are associated with pReplaceMe
		vector<Volume*> vVolumes;
		CollectVolumes(vVolumes, *this, pReplaceMe);
		for(uint i = 0; i < vVolumes.size(); ++i)
			replace(associated_faces_begin(vVolumes[i]),
					associated_faces_end(vVolumes[i]),
					pReplaceMe, f);
	}

//	now we have to copy all associated elements of pReplaceMe to f
	if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
	{
		m_aaEdgeContainerFACE[f].assign(associated_edges_begin(pReplaceMe),
										associated_edges_end(pReplaceMe));
		m_aaEdgeContainerFACE[pReplaceMe].clear();
	}

	if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		m_aaVolumeContainerFACE[f].assign(associated_volumes_begin(pReplaceMe),
											associated_volumes_end(pReplaceMe));
		m_aaVolumeContainerFACE[pReplaceMe].clear();
	}

//	remove the element from the storage and delete it.
	m_faceElementStorage.m_sectionContainer.erase(get_iterator(pReplaceMe), pReplaceMe->container_section());
	m_faceElementStorage.m_attachmentPipe.unregister_element(pReplaceMe);
	delete pReplaceMe;
}

void Grid::unregister_face(Face* f)
{
//	notify observers that the face is being erased
	NOTIFY_OBSERVERS_REVERSE(m_faceObservers, face_to_be_erased(this, f));

//	remove or disconnect from volumes
	if(num_volumes() > 0)
	{
		if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES) ||
			option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
		{
			vector<Volume*> vVolumes;
			CollectVolumes(vVolumes, *this, f, true);
			auto vIter= vVolumes.begin();

			if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
			{
				if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
					m_aaVolumeContainerFACE[f].clear();

			//	delete the collected volumes
				while(vIter != vVolumes.end())
				{
					Volume* eraseVol = *vIter;
					++vIter;
					erase(eraseVol);
				}
			}
			else
			{
			//	unregister from the collected volumes.
				for(;vIter != vVolumes.end(); ++vIter)
				{
					Volume* vol = *vIter;
					auto iter = find(m_aaFaceContainerVOLUME[vol].begin(),
														m_aaFaceContainerVOLUME[vol].end(), f);
					if(iter != m_aaFaceContainerVOLUME[vol].end())
						m_aaFaceContainerVOLUME[vol].erase(iter);
				}
			}
		}
	}

//	disconnect from edges
	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
	//	iterate through all edges of the face
		uint numEdges = f->num_edges();
		for(uint i = 0; i < numEdges; ++i)
		{
			Edge* e = get_edge(f, i);
			if(e != nullptr)
			{
				auto iter = find(m_aaFaceContainerEDGE[e].begin(),
													m_aaFaceContainerEDGE[e].end(), f);
				if(iter != m_aaFaceContainerEDGE[e].end())
					m_aaFaceContainerEDGE[e].erase(iter);
			}
		}
	}

//	disconnect from vertices
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
	{
		size_t numVrts = f->num_vertices();
		Face::ConstVertexArray vrts = f->vertices();
		for(uint i = 0; i < numVrts; ++i)
		{
			Vertex* vrt = vrts[i];
			auto iter = find(m_aaFaceContainerVERTEX[vrt].begin(),
												m_aaFaceContainerVERTEX[vrt].end(), f);
			if(iter != m_aaFaceContainerVERTEX[vrt].end())
				m_aaFaceContainerVERTEX[vrt].erase(iter);
		}
	}

//	remove the element from the storage
	m_faceElementStorage.m_sectionContainer.erase(get_iterator(f), f->container_section());
	m_faceElementStorage.m_attachmentPipe.unregister_element(f);
}

void Grid::change_face_options(uint optsNew)
{
//	check if associated edge information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, FACEOPT_STORE_ASSOCIATED_EDGES))
	{
		if(!option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
			face_store_associated_edges(true);
	}
	else
	{
		if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
			face_store_associated_edges(false);
	}

//	check if associated volume information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		if(!option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
			face_store_associated_volumes(true);
	}
	else
	{
		if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
			face_store_associated_volumes(false);
	}

//	turn auto-generation of edges on and of
	if(OPTIONS_CONTAIN_OPTION(optsNew, FACEOPT_AUTOGENERATE_EDGES))
	{
		if(!option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
			face_autogenerate_edges(true);
	}
	else
	{
		if(option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
			face_autogenerate_edges(false);
	}
}

void Grid::face_store_associated_edges(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
		{
		//	store associated edges
			attach_to_faces(m_aEdgeContainer);
			m_aaEdgeContainerFACE.access(*this, m_aEdgeContainer);

			bool createEdges = option_is_enabled(FACEOPT_AUTOGENERATE_EDGES);

		//	if EDGEOPT_STORE_ASSOCIATED_FACES is enabled, this is as simple
		//	as to iterate through all edges and store them at their
		//	associated faces.
		//	if not we have to iterate through the faces and store
		//	each of their edges..
		//	If createEdges == true, we want to store the edges sorted.
		//	The else branch thus has to be executed.
			if(!createEdges
				&& option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
			{
			//	iterate through the edges
				for(EdgeIterator iter = edges_begin(); iter != edges_end(); iter++)
				{
					Edge* e = *iter;
				//	iterate through the edges associated faces
					for(auto fIter = associated_faces_begin(e);
						fIter != associated_faces_end(e); fIter++)
					{
						m_aaEdgeContainerFACE[*fIter].push_back(e);
					}
				}
			}
			else
			{
			//	we have to find the edges by hand using GetEdge(...)
				for(FaceIterator iter = faces_begin(); iter != faces_end(); iter++)
				{
					Face* f = *iter;
					uint numEdges = f->num_edges();

					if(createEdges){
						m_aaEdgeContainerFACE[f].resize(numEdges);
						for(uint i = 0; i < numEdges; ++i)
							m_aaEdgeContainerFACE[f][i] = get_edge(f, i);
					}
					else{
						for(uint i = 0; i < numEdges; ++i)
						{
						//	get the i-th edge that is described by the EdgeDescriptor - if it exists at all.
							Edge* e = get_edge(f, i);
							if(e != nullptr)
								m_aaEdgeContainerFACE[f].push_back(e);
						}
					}
				}
			}

		//	store the options
			m_options |= FACEOPT_STORE_ASSOCIATED_EDGES;
		}
	}
	else
	{
		if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
		{
		//	remove face-edge connections
			detach_from_faces(m_aEdgeContainer);
			m_options &= (~FACEOPT_STORE_ASSOCIATED_EDGES);
		}
	}
}

void Grid::face_store_associated_volumes(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
		{
		//	store associated faces
			attach_to_faces(m_aVolumeContainer);
			m_aaVolumeContainerFACE.access(*this, m_aVolumeContainer);

		//	iterate through the volumes and store each of their faces...
			for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
			{
				Volume* v = *iter;

			//	if faces are already stored at their associated volumes, we iterate over them.
				if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
				{
					auto iterEnd = associated_faces_end(v);
					for(auto iter = associated_faces_begin(v); iter != iterEnd; iter++)
						m_aaVolumeContainerFACE[*iter].push_back(v);
				}
				else
				{
					uint numFaces = v->num_faces();
					for(uint i = 0; i < numFaces; ++i)
					{
						Face* f = get_face(v, i);

						if(f)
							m_aaVolumeContainerFACE[f].push_back(v);
					}
				}
			}

		//	store options
			m_options |= FACEOPT_STORE_ASSOCIATED_VOLUMES;
		}
	}
	else
	{
		if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
		{
		//	remove face-edge connections
			detach_from_faces(m_aVolumeContainer);
			m_options &= (~FACEOPT_STORE_ASSOCIATED_VOLUMES);
		}
	}
}

void Grid::face_autogenerate_edges(bool bAutogen)
{
	if(bAutogen)
	{
		if(!option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		{
			bool storeEdges = option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES);
			EdgeDescriptor ed;

		//	generate all missing edges now!
			for(FaceIterator iter = faces_begin(); iter != faces_end(); iter++)
			{
				Face* f = *iter;
			//	check for each edge of the face if it exists. if not create it and link it to the face.
				uint numEdges = f->num_edges();

			//	reserve memory for neighbors
				if(storeEdges){
					m_aaEdgeContainerFACE[f].clear();
					m_aaEdgeContainerFACE[f].reserve(numEdges);
				}

				for(uint i = 0; i < numEdges; ++i)
				{
				//	we can't use get_edge here, since we're manipulating m_aaEdgeContainerFACE on the fly.
					f->edge_desc(i, ed);
					Edge* e = find_edge_in_associated_edges(ed.vertex(0), ed);

					if(e == nullptr)
					{
					//	create the edge
						e = f->create_edge(i);
						if(storeEdges)
							m_aaEdgeContainerFACE[f].push_back(e);
						register_edge(e, f, f, nullptr);
					}
					else if(storeEdges)
						m_aaEdgeContainerFACE[f].push_back(e);
				}
			}

		//	store the option
			m_options |= FACEOPT_AUTOGENERATE_EDGES;

		//	if VOLOPT_AUTOGENERATE_FACES is active, too, then all volume-edges exist.
		//	Now, if VOLOPT_STORE_ASSOCIATED_EDGES is active but
		//	VOLOPT_AUTOGENERATE_EDGES is inactive, then we'll sort the associated-edges
		//	of volumes.
			if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES)
			   && option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES)
			   && (!option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)))
			{
				volume_sort_associated_edge_container();
			}
		}
	}
	else
	{
	//	stop auto-generation
		m_options &= (~FACEOPT_AUTOGENERATE_EDGES);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	VOLUMES
///	creates and removes connectivity data, as specified in optsNew.
void Grid::register_volume(Volume* v, GridObject* pParent)
{
	GCM_PROFILE_FUNC();

//	store the element and register it at the pipe.
	m_volumeElementStorage.m_attachmentPipe.register_element(v);
	m_volumeElementStorage.m_sectionContainer.insert(v, v->container_section());

//	create edges and faces if the according options are enabled.
//	register the volume at the associated vertices, edges and faces, if the according options are enabled.

//	register volume at vertices
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
		uint numVrts = v->num_vertices();
		Volume::ConstVertexArray vrts = v->vertices();
		for(uint i = 0; i < numVrts; ++i)
			m_aaVolumeContainerVERTEX[vrts[i]].push_back(v);
	}

	const bool createEdges = option_is_enabled(VOLOPT_AUTOGENERATE_EDGES);
	const bool createFaces = option_is_enabled(VOLOPT_AUTOGENERATE_FACES);
	const bool createEdgesIndirect = option_is_enabled(VOLOPT_AUTOGENERATE_FACES
													   | FACEOPT_AUTOGENERATE_EDGES);

	const bool edgesStoreVols = option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES);
	const bool volsStoreEdges = option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES);

	const bool facesStoreVols = option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES);
	const bool volsStoreFaces = option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES);


//	if elements are automatically created, then we can directly reserve memory
//	to store neighborhood relations
	if(volsStoreEdges && (createEdges || createEdgesIndirect))
		m_aaEdgeContainerVOLUME[v].reserve(v->num_edges());
	if(volsStoreFaces && createFaces)
		m_aaFaceContainerVOLUME[v].reserve(v->num_faces());

//	create faces first. This makes sense, since faces thus get a chance to
//	generate associated edges. This has the benefit, that the option
//	VOLOPT_AUTOGENERATE_EDGES has no effect in this case.

//	register with faces. autogenerate them if they do not already exist.
	if(createFaces
		|| facesStoreVols
		|| volsStoreFaces)
	{
		GCM_PROFILE(GCM_vol_handle_faces);
	//	iterate through the faces of the volume. create them if demanded.
	//	register faces at the volume and vice versa.
		uint numFaces = v->num_faces();
		FaceDescriptor fd;
		for(uint i = 0; i < numFaces; ++i)
		{
			v->face_desc(i, fd);
			Face* f = find_face_in_associated_faces(fd.vertex(0), fd);

			if(f == nullptr)
			{
				if(createFaces)
				{
					f = v->create_face(i);

					if(volsStoreFaces)
						m_aaFaceContainerVOLUME[v].push_back(f);

				//	register the face. By passing v, we speed up the process.
					register_face(f, pParent, v);
				}
			}
			else
			{
				if(facesStoreVols)
					m_aaVolumeContainerFACE[f].push_back(v);

				if(volsStoreFaces)
					m_aaFaceContainerVOLUME[v].push_back(f);
			}
		}
	}

//	register with edges. autogenerate missing edges if option is set
	if(createEdges
		|| edgesStoreVols
		|| volsStoreEdges)
	{
		GCM_PROFILE(GCM_vol_handle_edges);
	//	loop through all edges of the volume. check if we have to create it.
	//	register the volume with the edges and vice versa after that.
		uint numEdges = v->num_edges();
		EdgeDescriptor ed;
		for(uint i = 0; i < numEdges; ++i)
		{
			v->edge_desc(i, ed);
			Edge* e = find_edge_in_associated_edges(ed.vertex(0), ed);

			if(e == nullptr)
			{
				if(createEdges)
				{
				//	create the edge
					e = v->create_edge(i);

				//	store a reference to the edge, if required
					if(volsStoreEdges)
						m_aaEdgeContainerVOLUME[v].push_back(e);

				//	register the edge.
					register_edge(e, pParent, nullptr, v);
				}
			}
			else
			{
				if(edgesStoreVols){
				//	if the edge was indirectly created during creation of
				//	a side-face above, it already has v in its associated-
				//	volume container. v is the first element in this container
				//	in that case and doesn't have to be added again.
					VolumeContainer& assVols = m_aaVolumeContainerEDGE[e];
					if((assVols.empty()) || (assVols.front() != v))
						assVols.push_back(v);
				}

				if(volsStoreEdges)
					m_aaEdgeContainerVOLUME[v].push_back(e);
			}
		}
	}

	GCM_PROFILE(GCM_notify_volume_observers);
//	inform observers about the creation
	NOTIFY_OBSERVERS(m_volumeObservers, volume_created(this, v, pParent));
	GCM_PROFILE_END();
}

void Grid::register_and_replace_element(Volume* v, Volume* pReplaceMe)
{
//	check that v and pReplaceMe have the same number of vertices.
	if(v->num_vertices() != pReplaceMe->num_vertices())
	{
		LOG("WARNING in Grid::register_and_replace_element(Volume* v, Volume* pReplaceMe): v and pReplaceMe have different numbers of vertices.");
		LOG("  v has not been registered. pReplaceMe has not been replaced.");
		assert(!"ERROR in Grid::register_and_replace_element(Volume* v, Volume* pReplaceMe): v and pReplaceMe have different numbers of vertices.");
		return;
	}

//	store the element and register it at the pipe.
	m_volumeElementStorage.m_attachmentPipe.register_element(v);
	m_volumeElementStorage.m_sectionContainer.insert(v, v->container_section());

	pass_on_values(pReplaceMe, v);

//	assign vertices
	uint numVrts = v->num_vertices();
	Volume::ConstVertexArray vrts = pReplaceMe->vertices();
	for(uint i = 0; i < numVrts; ++i)
		v->set_vertex(i, vrts[i]);

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_volumeObservers, volume_created(this, v, pReplaceMe, true));
//	inform observers about the deletion
	NOTIFY_OBSERVERS_REVERSE(m_volumeObservers, volume_to_be_erased(this, pReplaceMe, v));

//	check if vertices, edges and faces reference pReplaceMe.
//	if so, correct those references.
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
		for(uint i = 0; i < numVrts; ++i)
			replace(associated_volumes_begin(vrts[i]),
					associated_volumes_end(vrts[i]),
					pReplaceMe, v);
	}

	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
	//	collect all edges that are associated with pReplaceMe
		vector<Edge*> vEdges;
		CollectEdges(vEdges, *this, pReplaceMe);
		for(uint i = 0; i < vEdges.size(); ++i)
			replace(associated_volumes_begin(vEdges[i]),
					associated_volumes_end(vEdges[i]),
					pReplaceMe, v);
	}

	if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
	//	collect all faces that are associated with pReplaceMe
		vector<Face*> vFaces;
		CollectFaces(vFaces, *this, pReplaceMe);
		for(uint i = 0; i < vFaces.size(); ++i)
			replace(associated_volumes_begin(vFaces[i]),
					associated_volumes_end(vFaces[i]),
					pReplaceMe, v);
	}

//	now we have to copy all associated elements of pReplaceMe to v
	if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
	{
		m_aaEdgeContainerVOLUME[v].assign(associated_edges_begin(pReplaceMe),
											associated_edges_end(pReplaceMe));
		m_aaEdgeContainerVOLUME[pReplaceMe].clear();
	}

	if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
	{
		m_aaFaceContainerVOLUME[v].assign(associated_faces_begin(pReplaceMe),
											associated_faces_end(pReplaceMe));
		m_aaFaceContainerVOLUME[pReplaceMe].clear();
	}

//	remove the element from the storage and delete it.
	m_volumeElementStorage.m_sectionContainer.erase(get_iterator(pReplaceMe), pReplaceMe->container_section());
	m_volumeElementStorage.m_attachmentPipe.unregister_element(pReplaceMe);
	delete pReplaceMe;
}

void Grid::unregister_volume(Volume* v)
{
//	notify observers that the face is being erased
	NOTIFY_OBSERVERS_REVERSE(m_volumeObservers, volume_to_be_erased(this, v));

//	disconnect from faces
	if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		uint numFaces = v->num_faces();
		vector<Face*> vFaces;
		vFaces.reserve(numFaces);
		CollectFaces(vFaces, *this, v, false);
		numFaces = vFaces.size();
		for(uint i = 0; i < numFaces; ++i)
		{
			Face* f = vFaces[i];

			if(f != nullptr)
			{
				auto iter = find(m_aaVolumeContainerFACE[f].begin(),
														m_aaVolumeContainerFACE[f].end(), v);
				if(iter != m_aaVolumeContainerFACE[f].end())
					m_aaVolumeContainerFACE[f].erase(iter);
			}
		}
	}

//	disconnect from edges
	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		uint numEdges = v->num_edges();
		for(uint i = 0; i < numEdges; ++i)
		{
		//	find the correct entry
			Edge* e = get_edge(v, i);
			if(e != nullptr)
			{
				auto iter = find(m_aaVolumeContainerEDGE[e].begin(),
														m_aaVolumeContainerEDGE[e].end(), v);
				if(iter != m_aaVolumeContainerEDGE[e].end())
					m_aaVolumeContainerEDGE[e].erase(iter);
			}
		}
	}

//	disconnect from vertices
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
	//	iterate through all associated vertices and update their connection-info
		uint numVertices = v->num_vertices();
		Volume::ConstVertexArray vrts = v->vertices();
		for(uint i = 0; i < numVertices; ++i)
		{
		//	find the correct entry
			Vertex* vrt = vrts[i];
			auto iter = find(m_aaVolumeContainerVERTEX[vrt].begin(),
													m_aaVolumeContainerVERTEX[vrt].end(), v);
			if(iter != m_aaVolumeContainerVERTEX[vrt].end())
				m_aaVolumeContainerVERTEX[vrt].erase(iter);
		}
	}

//	remove the element from the storage
	m_volumeElementStorage.m_sectionContainer.erase(get_iterator(v), v->container_section());
	m_volumeElementStorage.m_attachmentPipe.unregister_element(v);
}

void Grid::change_volume_options(uint optsNew)
{
//	check if associated edge information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, VOLOPT_STORE_ASSOCIATED_EDGES))
	{
		if(!option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
			volume_store_associated_edges(true);
	}
	else
	{
		if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
			volume_store_associated_edges(false);
	}

//	check if associated face information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, VOLOPT_STORE_ASSOCIATED_FACES))
	{
		if(!option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
			volume_store_associated_faces(true);
	}
	else
	{
		if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
			volume_store_associated_faces(false);
	}

//	enable or disable edge-auto-generation
	if(OPTIONS_CONTAIN_OPTION(optsNew, VOLOPT_AUTOGENERATE_EDGES))
	{
		if(!option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
			volume_autogenerate_edges(true);
	}
	else
	{
		if(option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
			volume_autogenerate_edges(false);
	}

//	enable or disable face-auto-generation
	if(OPTIONS_CONTAIN_OPTION(optsNew, VOLOPT_AUTOGENERATE_FACES))
	{
		if(!option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
			volume_autogenerate_faces(true);
	}
	else
	{
		if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
			volume_autogenerate_faces(false);
	}
}

void Grid::volume_store_associated_edges(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
		{
		//	store associated edges
			attach_to_volumes(m_aEdgeContainer);
			m_aaEdgeContainerVOLUME.access(*this, m_aEdgeContainer);

			bool createEdges = option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)
								|| option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES);
		//	if EDGEOPT_STORE_ASSOCIATED_VOLUMES is enabled, this is as simple
		//	as to iterate through all edges and store them at their
		//	associated volumes.
		//	if not we have to iterate through the volumes and store
		//	each of their edges..
		//	If createEdges == true, we want to store the edges sorted.
		//	The else branch thus has to be executed.
			if(!createEdges
				&& option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
			{
			//	iterate through the edges
				for(EdgeIterator iter = edges_begin(); iter != edges_end(); iter++)
				{
					Edge* e = *iter;
				//	iterate through the edges associated faces
					for(auto vIter = associated_volumes_begin(e);
						vIter != associated_volumes_end(e); vIter++)
					{
						m_aaEdgeContainerVOLUME[*vIter].push_back(e);
					}
				}
			}
			else
			{
			//	we have to find the edges by hand using GetEdge(...)
				for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
				{
					Volume* v = *iter;
					int numEdges = v->num_edges();

					if(createEdges){
					//	all edges exist. We can insert them sorted.
						m_aaEdgeContainerVOLUME[v].resize(numEdges);
						for(int i = 0; i < numEdges; ++i)
							m_aaEdgeContainerVOLUME[v][i] = get_edge(v, i);
					}
					else{
					//	iterate through the edges of the volume
						for(int i = 0; i < numEdges; ++i)
						{
						//	get the edge-descriptor
							Edge* e = get_edge(v, i);
							if(e != nullptr)
								m_aaEdgeContainerVOLUME[v].push_back(e);
						}
					}
				}
			}

		//	store the option
			m_options |= VOLOPT_STORE_ASSOCIATED_EDGES;
		}
	}
	else
	{
		if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
		{
			//	remove vol-edge connection
			detach_from_volumes(m_aEdgeContainer);
			m_options &= (~VOLOPT_STORE_ASSOCIATED_EDGES);
		}
	}
}

void Grid::volume_store_associated_faces(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
		{
		//	store associated faces
			attach_to_volumes(m_aFaceContainer);
			m_aaFaceContainerVOLUME.access(*this, m_aFaceContainer);

			bool createFaces = option_is_enabled(VOLOPT_AUTOGENERATE_FACES);

		//	if FACEOPT_STORE_ASSOCIATED_VOLUMES is enabled, this is as simple
		//	as to iterate through all faces and store them at their
		//	associated volumes.
		//	if not we have to iterate through the volumes and store
		//	each of their faces..
		//	if createFaces == true then we want to store the faces sorted.
			if(!createFaces
				&& option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
			{
			//	iterate through the edges
				for(FaceIterator iter = faces_begin(); iter != faces_end(); iter++)
				{
					Face* f = *iter;
				//	iterate through the faces associated volumes
					for(auto vIter = associated_volumes_begin(f);
						vIter != associated_volumes_end(f); vIter++)
					{
						m_aaFaceContainerVOLUME[*vIter].push_back(f);
					}
				}
			}
			else
			{
			//	we have to find the edges by hand using GetEdge(...)
				FaceDescriptor fd;

				for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
				{
					Volume* v = *iter;
				//	iterate through the faces of the volume and get the corresponding face from the grid.
					int numFaces = v->num_faces();

					if(createFaces){
						m_aaFaceContainerVOLUME[v].resize(numFaces);
						for(int i = 0; i < numFaces; ++i)
							m_aaFaceContainerVOLUME[v][i] = get_face(v, i);

					}
					else{
						for(int i = 0; i < numFaces; ++i)
						{
							v->face_desc(i, fd);
							Face* f = find_face_in_associated_faces(fd.vertex(0), fd);

							if(f)
								m_aaFaceContainerVOLUME[v].push_back(f);
						}
					}
				}
			}

		//	store options
			m_options |= VOLOPT_STORE_ASSOCIATED_FACES;
		}
	}
	else
	{
		if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
		{
			//	remove vol-edge connection
			detach_from_volumes(m_aFaceContainer);
			m_options &= (~VOLOPT_STORE_ASSOCIATED_FACES);
		}
	}
}

void Grid::volume_autogenerate_edges(bool bAutogen)
{
	if(bAutogen)
	{
		if(!option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
		{
			bool volsStoreEdges = option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES);
			EdgeDescriptor ed;

		//	generate all missing edges now!
			for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
			{
				Volume* v = *iter;

			//	check for each edge of the volume if it exists. if not create it and link it to the edge.
				uint numEdges = v->num_edges();

			//	reserve memory for neighbors
				if(volsStoreEdges){
					m_aaEdgeContainerVOLUME[v].clear();
					m_aaEdgeContainerVOLUME[v].reserve(numEdges);
				}

				for(uint i = 0; i < numEdges; ++i)
				{
				//	we can't use get_edge here, since we modify m_aaEdgeContainerVOLUME on the fly.
					v->edge_desc(i, ed);
					Edge* e = find_edge_in_associated_edges(ed.vertex(0), ed);

					if(e == nullptr)
					{
					//	create the edge
						e = v->create_edge(i);
						if(volsStoreEdges) // has to be performed before register_edge
							m_aaEdgeContainerVOLUME[v].push_back(e);
						register_edge(e, v, nullptr, v);
					}
					else if(volsStoreEdges)
						m_aaEdgeContainerVOLUME[v].push_back(e);
				}
			}

		//	store the option
			m_options |= VOLOPT_AUTOGENERATE_EDGES;
		}
	}
	else
	{
	//	stop auto-generation
		m_options &= (~VOLOPT_AUTOGENERATE_EDGES);
	}
}

void Grid::volume_autogenerate_faces(bool bAutogen)
{
	if(bAutogen)
	{
		if(!option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
		{
			bool volsStoreFaces = option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES);
			FaceDescriptor fd;

		//	generate all missing faces now!
			for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
			{
				Volume* v = *iter;
			//	check for each face of the volume if it exists. if not create it and link it to the face.
				uint numFaces = v->num_faces();

			//	reserve memory for neighbors
				if(volsStoreFaces){
					m_aaFaceContainerVOLUME[v].clear();
					m_aaFaceContainerVOLUME[v].reserve(numFaces);
				}

				for(uint i = 0; i < numFaces; ++i)
				{
				//	we can't use get_face here, since we modify m_aaEdgeContainerVOLUME on the fly.
					v->face_desc(i, fd);
					Face* f = find_face_in_associated_faces(fd.vertex(0), fd);

					if(f == nullptr)
					{
					//	create the face
						f = v->create_face(i);
						if(volsStoreFaces) // has to be performed before register_face
							m_aaFaceContainerVOLUME[v].push_back(f);
						register_face(f, v, v);
					}
					else if(volsStoreFaces)
						m_aaFaceContainerVOLUME[v].push_back(f);
				}
			}

		//	store the option
			m_options |= VOLOPT_AUTOGENERATE_FACES;

		//	if FACEOPT_AUTOGENERATE_EDGES is active, too, then all volume-edges exist.
		//	Now, if VOLOPT_STORE_ASSOCIATED_EDGES is active but
		//	VOLOPT_AUTOGENERATE_EDGES is inactive, then we'll sort the associated-edges
		//	of volumes.
			if(option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)
			   && option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES)
			   && (!option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)))
			{
				volume_sort_associated_edge_container();
			}
		}
	}
	else
	{
	//	stop auto-generation
		m_options &= (~VOLOPT_AUTOGENERATE_FACES);
	}
}


void Grid::volume_sort_associated_edge_container()
{
	if(!option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
		return;

	EdgeContainer tmpCon;
	for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
	{
		Volume* v = *iter;
		uint numEdges = v->num_edges();

	//	since we already have a container of associated edges, we will search it
	//	for matching edges...
	//	to reduce memory allocation, we'll first copy the values to tmpCon
		EdgeContainer& edges = m_aaEdgeContainerVOLUME[v];
		tmpCon.clear();
		for(size_t i = 0; i < edges.size(); ++i)
			tmpCon.push_back(edges[i]);

		edges.clear();

	//	check for each edge-desc, whether a matching edge is contained in edges.
	//	if so push it back to edgesSorted.
		EdgeDescriptor ed;
		for(uint i = 0; i < numEdges; ++i){
			v->edge_desc(i, ed);
			for(size_t j = 0; j < tmpCon.size(); ++j){
				if(tmpCon[j] && CompareVertices(&ed, tmpCon[j])){
					edges.push_back(tmpCon[j]);
					tmpCon[j] = nullptr;
					break;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	replace_vertex
bool Grid::replace_vertex(Vertex* vrtOld, Vertex* vrtNew)
{
//	this bool should be a parameter. However one first would have
//	to add connectivity updates for double-elements in this method,
//	to handle the case when eraseDoubleElements is set to false.
	bool eraseDoubleElements = true;

//	iterate through associated objects of vrtOld.
//	replace in each vrtOld by vrtNew.
//	Delete objects that would be doubled if eraseDoubleElements is set to true.
//	Update connectivity information.
//	erase vrtOld
//	all elements that will be erased, during this method will be
//	unregistered manually on the fly.

//	containers for object collection
	vector<Edge*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVolumes;

	EdgeDescriptor ed;
	FaceDescriptor fd;
	VolumeDescriptor vd;

//	since vrtOld will be erased, we will notify all observers
	NOTIFY_OBSERVERS_REVERSE(m_vertexObservers, vertex_to_be_erased(this, vrtOld));

//	EDGES
	if(num_edges() > 0)
	{
		autoenable_option(VRTOPT_STORE_ASSOCIATED_EDGES,
						"Grid::replace_vertex(...)",
						"VRTOPT_STORE_ASSOCIATED_EDGES");

		auto iter = associated_edges_begin(vrtOld);
		auto iterEnd = associated_edges_end(vrtOld);
		while(iter != iterEnd)
		{
			Edge* e = *iter;
			++iter;

		//	if eraseDoubleElementes is enabled and the new edge would
		//	already exist, we wont replace the vertex in it.
			bool bReplaceVertex = true;

			if(eraseDoubleElements)
			{
			//	create an edge-descriptor of the edge that will be created.
				if(e->vertex(0) == vrtOld)
					ed.set_vertices(e->vertex(1), vrtNew);
				else
					ed.set_vertices(e->vertex(0), vrtNew);

			//	check if this edge already exists.
				Edge* eNew = get_edge(ed);
				if(eNew)
				{
				//	The edge will be removed. Notify observers.
					NOTIFY_OBSERVERS_REVERSE(m_edgeObservers, edge_to_be_erased(this, e));

				//	Since we want to avoid deletion of associated elements
				//	we have to handle the erasure of e by our self.
					bReplaceVertex = false;

				//	before the removal we will replace its entry in associated-
				//	element-lists with the pointer to the existing edge.
					if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
					{
					//	unregister e from the vertex with which e connects vrtOld.
						EdgeContainer& ec = m_aaEdgeContainerVERTEX[ed.vertex(0)];
						auto tmpI = find(ec.begin(), ec.end(), e);
						if(tmpI != ec.end())
							ec.erase(tmpI);
					}

					if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES)
							&& num_faces() > 0)
					{
					//	unregister the edge from all adjacent faces
						CollectFaces(vFaces, *this, e);
						for(uint i = 0; i < vFaces.size(); ++i)
						{
							Face* f = vFaces[i];
							replace(associated_edges_begin(f),
									associated_edges_end(f), e, eNew);
						}
					}

					if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES)
							&& num_volumes() > 0)
					{
					//	unregister the edge from all adjacent volumes
						CollectVolumes(vVolumes, *this, e);
						for(uint i = 0; i < vVolumes.size(); ++i)
						{
							Volume* v = vVolumes[i];
							replace(associated_edges_begin(v),
									associated_edges_end(v), e, eNew);
						}
					}

				//	we can now remove e from the storage.
					m_edgeElementStorage.m_sectionContainer.erase(get_iterator(e), e->container_section());
					m_edgeElementStorage.m_attachmentPipe.unregister_element(e);
					delete e;
				}
			}

		//	perform the replace
			if(bReplaceVertex)
			{
			//	replace the vertex
				if(e->vertex(0) == vrtOld)
					e->set_vertex(0, vrtNew);
				else
					e->set_vertex(1, vrtNew);

			//	register e at vrtNew
				if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
				{
					m_aaEdgeContainerVERTEX[vrtNew].push_back(e);
				}
			}
		}

	//	we will now clear the associated edges array of vrtOld.
		m_aaEdgeContainerVERTEX[vrtOld].clear();
	}

//	FACES
	if(num_faces() > 0)
	{
		autoenable_option(VRTOPT_STORE_ASSOCIATED_FACES,
						"Grid::replace_vertex(...)",
						"VRTOPT_STORE_ASSOCIATED_FACES");

		auto iter = associated_faces_begin(vrtOld);
		auto iterEnd = associated_faces_end(vrtOld);
		while(iter != iterEnd)
		{
			Face* f = *iter;
			uint numVrts = f->num_vertices();
			Face::ConstVertexArray vrts = f->vertices();

			++iter;

		//	if eraseDoubleElementes is enabled and the new face would
		//	already exist, we wont replace the vertex in it.
			bool bReplaceVertex = true;

			if(eraseDoubleElements)
			{
			//	create a face-descriptor of the face that will be created.
				fd.set_num_vertices(numVrts);
				for(uint i = 0; i < numVrts; ++i)
				{
					if(vrts[i] == vrtOld)
						fd.set_vertex(i, vrtNew);
					else
						fd.set_vertex(i, vrts[i]);
				}

			//	check if this face already exists.
				Face* fNew = get_face(fd);
				if(fNew)
				{
				//	The face will be removed. Notify observers.
					NOTIFY_OBSERVERS_REVERSE(m_faceObservers, face_to_be_erased(this, f));

				//	Since we want to avoid deletion of associated elements
				//	we have to handle the erasure of f by our self.
					bReplaceVertex = false;

				//	before the removal we will replace its entry in associated-
				//	element-lists with the pointer to the existing face.
					if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
					{
					//	unregister f from the vertices with which f connects vrtOld.
						for(uint i = 0; i < numVrts; ++i)
						{
							if(vrts[i] != vrtOld)
							{
								FaceContainer& fc = m_aaFaceContainerVERTEX[vrts[i]];
								auto tmpI = find(fc.begin(), fc.end(), f);
								if(tmpI != fc.end())
									fc.erase(tmpI);
							}
						}
					}

					if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES)
						&& num_edges() > 0)
					{
					//	unregister f from adjacent edges.
						CollectEdges(vEdges, *this, f);
						for(uint i = 0; i < vEdges.size(); ++i)
						{
							FaceContainer& fc = m_aaFaceContainerEDGE[vEdges[i]];
							auto tmpI = find(fc.begin(), fc.end(), f);
							if(tmpI != fc.end())
								fc.erase(tmpI);
						}
					}

					if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES)
							&& num_volumes() > 0)
					{
					//	unregister the face from all adjacent volumes
						CollectVolumes(vVolumes, *this, f);
						for(uint i = 0; i < vVolumes.size(); ++i)
						{
							Volume* v = vVolumes[i];
							replace(associated_faces_begin(v),
									associated_faces_end(v), f, fNew);
						}
					}

				//	we can now remove f from the storage.
					m_faceElementStorage.m_sectionContainer.erase(get_iterator(f), f->container_section());
					m_faceElementStorage.m_attachmentPipe.unregister_element(f);
					delete f;
				}
			}

		//	perform the replace
			if(bReplaceVertex)
			{
			//	replace the vertex
				uint numVrts = f->num_vertices();
				for(uint i = 0; i < numVrts; ++i)
				{
					if(vrts[i] == vrtOld)
						f->set_vertex(i, vrtNew);
				}

			//	register f at vrtNew
				if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
				{
					m_aaFaceContainerVERTEX[vrtNew].push_back(f);
				}

			//	register f at existing edges
				if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
				{
				//	only edges which contain vrtNew are relevant.
					uint numEdges = f->num_edges();
					for(uint i = 0; i < numEdges; ++i)
					{
						f->edge_desc(i, ed);
						if(ed.vertex(0) == vrtNew || ed.vertex(1) == vrtNew)
						{
							Edge* tEdge = get_edge(ed);
							if(tEdge)
							{
								FaceContainer& fc = m_aaFaceContainerEDGE[tEdge];
								auto tmpI = find(fc.begin(), fc.end(), f);
								if(tmpI == fc.end())
									fc.push_back(f);
							}
						}
					}
				}
			}
		}
	//	clear associated faces list of vrtOld
		m_aaFaceContainerVERTEX[vrtOld].clear();
	}

//	VOLUMES
	if(num_volumes() > 0)
	{
		autoenable_option(VRTOPT_STORE_ASSOCIATED_VOLUMES,
						"Grid::replace_vertex(...)",
						"VRTOPT_STORE_ASSOCIATED_VOLUMES");

		auto iter = associated_volumes_begin(vrtOld);
		auto iterEnd = associated_volumes_end(vrtOld);
		while(iter != iterEnd)
		{
			Volume* v = *iter;
			uint numVrts = v->num_vertices();
			Volume::ConstVertexArray vrts = v->vertices();
			++iter;

		//	if eraseDoubleElementes is enabled and the new face would
		//	already exist, we wont replace the vertex in it.
			bool bReplaceVertex = true;

			if(eraseDoubleElements)
			{
			//	create a volume-descriptor of the volume that will be created.
				vd.set_num_vertices(numVrts);

				for(uint i = 0; i < numVrts; ++i)
				{
					if(vrts[i] == vrtOld)
						vd.set_vertex(i, vrtNew);
					else
						vd.set_vertex(i, vrts[i]);
				}

			//	check if this volume already exists.
				Volume* vNew = get_volume(vd);
				if(vNew)
				{
				//	The volume will be removed. Notify observers.
					NOTIFY_OBSERVERS_REVERSE(m_volumeObservers, volume_to_be_erased(this, v));

				//	Since we want to avoid deletion of associated elements
				//	we have to handle the erasure of v by our self.
					bReplaceVertex = false;

				//	before the removal we will replace its entry in associated-
				//	element-lists with the pointer to the existing volume.
					if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
					{
					//	unregister v from the vertices with which v connects vrtOld.
						for(uint i = 0; i < numVrts; ++i)
						{
							if(vrts[i] != vrtOld)
							{
								VolumeContainer& vc = m_aaVolumeContainerVERTEX[vrts[i]];
								auto tmpI = find(vc.begin(), vc.end(), v);
								if(tmpI != vc.end())
									vc.erase(tmpI);
							}
						}
					}

					if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES)
						&& num_edges() > 0)
					{
					//	unregister v from adjacent edges.
						CollectEdges(vEdges, *this, v);
						for(uint i = 0; i < vEdges.size(); ++i)
						{
							VolumeContainer& vc = m_aaVolumeContainerEDGE[vEdges[i]];
							auto tmpI = find(vc.begin(), vc.end(), v);
							if(tmpI != vc.end())
								vc.erase(tmpI);
						}
					}

					if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES)
						&& num_faces() > 0)
					{
					//	unregister v from adjacent faces.
						CollectFaces(vFaces, *this, v);
						for(uint i = 0; i < vFaces.size(); ++i)
						{
							VolumeContainer& vc = m_aaVolumeContainerFACE[vFaces[i]];
							auto tmpI = find(vc.begin(), vc.end(), v);
							if(tmpI != vc.end())
								vc.erase(tmpI);
						}
					}


				//	we can now remove v from the storage.
					m_volumeElementStorage.m_sectionContainer.erase(get_iterator(v), v->container_section());
					m_volumeElementStorage.m_attachmentPipe.unregister_element(v);
					delete v;
				}
			}

		//	perform the replace
			if(bReplaceVertex)
			{
			//	replace the vertex
				uint numVrts = v->num_vertices();
				for(uint i = 0; i < numVrts; ++i)
				{
					if(vrts[i] == vrtOld)
						v->set_vertex(i, vrtNew);
				}

			//	register v at vrtNew
				if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
				{
					m_aaVolumeContainerVERTEX[vrtNew].push_back(v);
				}

			//	register v at existing edges
				if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
				{
				//	only edges which contain vrtNew are relevant.
					uint numEdges = v->num_edges();
					for(uint i = 0; i < numEdges; ++i)
					{
						v->edge_desc(i, ed);
						if(ed.vertex(0) == vrtNew || ed.vertex(1) == vrtNew)
						{
							Edge* tEdge = get_edge(ed);
							if(tEdge)
							{
								VolumeContainer& vc = m_aaVolumeContainerEDGE[tEdge];
								auto tmpI = find(vc.begin(), vc.end(), v);
								if(tmpI == vc.end())
									vc.push_back(v);
							}
						}
					}
				}

			//	register v at existing faces
				if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
				{
				//	only faces which contain vrtNew are relevant.
					uint numFaces = v->num_faces();
					for(uint i = 0; i < numFaces; ++i)
					{
						v->face_desc(i, fd);
					//	check whether fd contains vrtNew
						bool bContainsVrtNew = false;
						size_t numFaceVrts = fd.num_vertices();
						Face::ConstVertexArray fvrts = fd.vertices();
						for(uint j = 0; j < numFaceVrts; ++j)
						{
							if(fvrts[j] == vrtNew)
							{
								bContainsVrtNew = true;
								break;
							}
						}

						if(bContainsVrtNew)
						{
							Face* tFace = get_face(fd);
							if(tFace)
							{
								VolumeContainer& vc = m_aaVolumeContainerFACE[tFace];
								auto tmpI = find(vc.begin(), vc.end(), v);
								if(tmpI == vc.end())
									vc.push_back(v);
							}
						}
					}
				}
			}
		}
	}

//	finally erase vrtOld.
	m_vertexElementStorage.m_sectionContainer.erase(get_iterator(vrtOld), vrtOld->container_section());
	m_vertexElementStorage.m_attachmentPipe.unregister_element(vrtOld);
	delete vrtOld;

	return true;
}

////////////////////////////////////////////////////////////////////////
//	replace_vertex_is_valid
bool Grid::replace_vertex_is_valid(Vertex* vrtOld, Vertex* vrtNew)
{
//	iterate through all geometric objects connected to vrtOld and check
//	if they contain vrtNew too. If not return true, else return false.

//	check edges
	if(num_edges() > 0)
	{
	//	we need vertex-edge connectivity
		autoenable_option(VRTOPT_STORE_ASSOCIATED_EDGES,
						"Grid::replace_vertex_is_valid(...)",
						"VRTOPT_STORE_ASSOCIATED_EDGES");

		auto iterEnd = associated_edges_end(vrtOld);
		for(auto iter = associated_edges_begin(vrtOld);
			iter != iterEnd; ++iter)
		{
		//	check whether vrtNew is contained in the edge.
			if((*iter)->vertex(0) == vrtNew)
				return false;
			if((*iter)->vertex(1) == vrtNew)
				return false;
		}
	}

//	check faces
	if(num_faces() > 0)
	{
	//	we need vertex-face connectivity
		autoenable_option(VRTOPT_STORE_ASSOCIATED_FACES,
						"Grid::replace_vertex_is_valid(...)",
						"VRTOPT_STORE_ASSOCIATED_FACES");

		auto iterEnd = associated_faces_end(vrtOld);
		for(auto iter = associated_faces_begin(vrtOld);
			iter != iterEnd; ++iter)
		{
			if(FaceContains(*iter, vrtNew))
				return false;
		}
	}

//	check volumes
	if(num_volumes() > 0)
	{
	//	we need vertex-face connectivity
		autoenable_option(VRTOPT_STORE_ASSOCIATED_VOLUMES,
						"Grid::replace_vertex_is_valid(...)",
						"VRTOPT_STORE_ASSOCIATED_VOLUMES");

		auto iterEnd = associated_volumes_end(vrtOld);
		for(auto iter = associated_volumes_begin(vrtOld);
			iter != iterEnd; ++iter)
		{
			if(VolumeContains(*iter, vrtNew))
				return false;
		}
	}

	return true;
}



////////////////////////////////////////////////////////////////////////////////
//	ASSOCIATED ELEMENTS
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//	ASSOCIATED VERTICES
void Grid::get_associated(SecureVertexContainer& vrts, Edge* e)
{
	vrts.set_external_array(e->vertices(), e->num_vertices());
}

void Grid::get_associated(SecureVertexContainer& vrts, Face* f)
{
	vrts.set_external_array(f->vertices(), f->num_vertices());
}

void Grid::get_associated(SecureVertexContainer& vrts, Volume* v)
{
	vrts.set_external_array(v->vertices(), v->num_vertices());
}


////////////////////////////////////////////////////////////////////////////////
//	ASSOCIATED EDGES
void Grid::get_associated(SecureEdgeContainer& edges, Vertex* v)
{
//	Without the VRTOPT_STORE_ASSOCIATED_... option, this operation would have
//	complexity O(n). This has to be avoided! We thus simply enable the option.
//	This takes some time, however, later queries will greatly benefit.
	if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES)){
	//	only enable the option if edges exist at all
		if(num<Edge>() == 0){
			edges.clear();
			return;
		}

		LOG("WARNING in get_associated(edges, vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES." << endl);
		vertex_store_associated_edges(true);
	}

	EdgeContainer& assEdges = m_aaEdgeContainerVERTEX[v];
	if(assEdges.empty())
		edges.clear();
	else
		edges.set_external_array(&assEdges.front(), assEdges.size());
}

void Grid::get_associated(SecureEdgeContainer& edges, Face* f)
{
//	to improve performance, we first check the grid options.
	if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we can output the associated array directly
		EdgeContainer& assEdges = m_aaEdgeContainerFACE[f];
		if(assEdges.empty())
			edges.clear();
		else
			edges.set_external_array(&assEdges.front(), assEdges.size());
	}
	else{
	//	clear the container
		edges.clear();

	//	if no edges are present, we can leave immediately
		if(num<Edge>() == 0)
			return;

	//	get the edges one by one
		uint numEdges = f->num_edges();
		for(uint i = 0; i < numEdges; ++i){
			Edge* e = get_edge(f, i);
			if(e != nullptr)
				edges.push_back(e);
		}
	}
}

void Grid::get_associated(SecureEdgeContainer& edges, Volume* v)
{
//	to improve performance, we first check the grid options.
	if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we can output the associated array directly
		EdgeContainer& assEdges = m_aaEdgeContainerVOLUME[v];
		if(assEdges.empty())
			edges.clear();
		else
			edges.set_external_array(&assEdges.front(), assEdges.size());
	}
	else{
	//	clear the container
		edges.clear();

	//	if no edges are present, we can leave immediately
		if(num<Edge>() == 0)
			return;

	//	get the edges one by one
		uint numEdges = v->num_edges();
		for(uint i = 0; i < numEdges; ++i){
			Edge* e = get_edge(v, i);
			if(e != nullptr)
				edges.push_back(e);
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	ASSOCIATED FACES
void Grid::get_associated(SecureFaceContainer& faces, Vertex* v)
{
//	Without the VRTOPT_STORE_ASSOCIATED_... option, this operation would have
//	complexity O(n). This has to be avoided! We thus simply enable the option.
//	This takes some time, however, later queries will greatly benefit.
	if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES)){
	//	only enable the option if faces exist at all
		if(num<Face>() == 0){
			faces.clear();
			return;
		}

		LOG("WARNING in get_associated(faces, vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_FACES." << endl);
		vertex_store_associated_faces(true);
	}

	FaceContainer& assFaces = m_aaFaceContainerVERTEX[v];
	if(assFaces.empty())
		faces.clear();
	else
		faces.set_external_array(&assFaces.front(), assFaces.size());
}

void Grid::get_associated(SecureFaceContainer& faces, Edge* e)
{
//	best option: EDGEOPT_STORE_ASSOCIATED_FACES
	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES)){
	//	we can output the associated array directly
		FaceContainer& assFaces = m_aaFaceContainerEDGE[e];
		if(assFaces.empty())
			faces.clear();
		else
			faces.set_external_array(&assFaces.front(), assFaces.size());
	}
	else{
	//	second best: iterate through all faces associated with the first end-point of e
	//	and check for each if it contains e. If so push it into the container.
	//	VRTOPT_STORE_ASSOCIATED_FACES has to be enabled for this. Only continue,
	//	if faces exist at all
		faces.clear();
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES)){
		//	only enable the option if faces exist at all
			if(num<Face>() == 0){
				return;
			}

			LOG("WARNING in get_associated(faces, edge): auto-enabling VRTOPT_STORE_ASSOCIATED_FACES." << endl);
			vertex_store_associated_faces(true);
		}

	//	check as few faces as possible
		Vertex* vrt = e->vertex(0);
/*	//	This check could be beneficial - however, it probably introduces unnecessary overhead.
 		if(m_aaFaceContainerVERTEX[vrt].size() >
			m_aaFaceContainerVERTEX[e->vertex(1)].size())
		{
			vrt = e->vertex(1);
		}
*/

		FaceContainer& assFaces = m_aaFaceContainerVERTEX[vrt];
		for(size_t i = 0; i < assFaces.size(); ++i){
			if(FaceContains(assFaces[i], e))
				faces.push_back(assFaces[i]);
		}
	}
}

void Grid::get_associated(SecureFaceContainer& faces, Volume* v)
{
//	to improve performance, we first check the grid options.
	if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
	{
	//	we can output the associated array directly
		FaceContainer& assFaces = m_aaFaceContainerVOLUME[v];
		if(assFaces.empty())
			faces.clear();
		else
			faces.set_external_array(&assFaces.front(), assFaces.size());
	}
	else{
	//	clear the container
		faces.clear();

	//	if no faces are present, we can leave immediately
		if(num<Face>() == 0)
			return;

	//	get the faces one by one
		uint numFaces = v->num_faces();
		for(uint i = 0; i < numFaces; ++i){
			Face* f = get_face(v, i);
			if(f != nullptr)
				faces.push_back(f);
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	ASSOCIATED VOLUMES
void Grid::get_associated(SecureVolumeContainer& vols, Vertex* v)
{
//	Without the VRTOPT_STORE_ASSOCIATED_... option, this operation would have
//	complexity O(n). This has to be avoided! We thus simply enable the option.
//	This takes some time, however, later queries will greatly benefit.
	if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES)){
	//	only enable the option if volumes exist at all
		if(num<Volume>() == 0){
			vols.clear();
			return;
		}

		LOG("WARNING in get_associated(volumes, vrt): auto-enabling VRTOPT_STORE_ASSOCIATED_VOLUMES." << endl);
		vertex_store_associated_volumes(true);
	}

	VolumeContainer& assVols = m_aaVolumeContainerVERTEX[v];
	if(assVols.empty())
		vols.clear();
	else
		vols.set_external_array(&assVols.front(), assVols.size());
}

void Grid::get_associated(SecureVolumeContainer& vols, Edge* e)
{
//	best option: EDGEOPT_STORE_ASSOCIATED_VOLUMES
	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES)){
	//	we can output the associated array directly
		VolumeContainer& assVols = m_aaVolumeContainerEDGE[e];
		if(assVols.empty())
			vols.clear();
		else
			vols.set_external_array(&assVols.front(), assVols.size());
	}
	else{
	//	second best: iterate through all volumes associated with the first end-point of e
	//	and check for each if it contains e. If so push it into the container.
	//	VRTOPT_STORE_ASSOCIATED_VOLUMES has to be enabled for this. Only continue,
	//	if volumes exist at all
		vols.clear();
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES)){
		//	only enable the option if volumes exist at all
			if(num<Volume>() == 0){
				return;
			}

			LOG("WARNING in get_associated(volumes, edge): auto-enabling VRTOPT_STORE_ASSOCIATED_VOLUMES." << endl);
			vertex_store_associated_volumes(true);
		}

	//	check as few faces as possible
		Vertex* vrt = e->vertex(0);
/*	//	This check could be beneficial - however, it probably introduces unnecessary overhead.
 		if(m_aaVolumeContainerVERTEX[vrt].size() >
			m_aaVolumeContainerVERTEX[e->vertex(1)].size())
		{
			vrt = e->vertex(1);
		}
*/

		VolumeContainer& assVols = m_aaVolumeContainerVERTEX[vrt];
		for(size_t i = 0; i < assVols.size(); ++i){
			if(VolumeContains(assVols[i], e))
				vols.push_back(assVols[i]);
		}
	}
}

void Grid::get_associated(SecureVolumeContainer& vols, Face* f)
{
//	best option: FACEOPT_STORE_ASSOCIATED_VOLUMES
	if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES)){
	//	we can output the associated array directly
		VolumeContainer& assVols = m_aaVolumeContainerFACE[f];
		if(assVols.empty())
			vols.clear();
		else
			vols.set_external_array(&assVols.front(), assVols.size());
	}
	else
		get_associated_vols_raw(vols, f);
}

void Grid::get_associated_vols_raw(SecureVolumeContainer& vols, Face* f)
{

//	iterate through all volumes associated with the first corner of f
//	and check for each if it contains f. If so push it into the container.
//	VRTOPT_STORE_ASSOCIATED_VOLUMES has to be enabled for this. Only continue,
//	if volumes exist at all
	vols.clear();
	if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES)){
	//	only enable the option if volumes exist at all
		if(num<Volume>() == 0){
			return;
		}

		LOG("WARNING in get_associated_vols_raw(volumes, face): auto-enabling VRTOPT_STORE_ASSOCIATED_VOLUMES." << endl);
		vertex_store_associated_volumes(true);
	}

//	check as few faces as possible
	Vertex* vrt = f->vertex(0);

	VolumeContainer& assVols = m_aaVolumeContainerVERTEX[vrt];
	for(size_t i = 0; i < assVols.size(); ++i){
		Volume* v = assVols[i];
		if(VolumeContains(v, f->vertex(1))){
			if(VolumeContains(v, f->vertex(2))){
				if(VolumeContains(v, f))
					vols.push_back(v);
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//	ASSOCIATED SORTED
void Grid::get_associated_sorted(SecureVertexContainer& vrts, Edge* e) const
{
	vrts.set_external_array(e->vertices(), e->num_vertices());
}

void Grid::get_associated_sorted(SecureVertexContainer& vrts, Face* f) const
{
	vrts.set_external_array(f->vertices(), f->num_vertices());
}

void Grid::get_associated_sorted(SecureVertexContainer& vrts, Volume* v) const
{
	vrts.set_external_array(v->vertices(), v->num_vertices());
}


void Grid::get_associated_sorted(SecureEdgeContainer& edges, Vertex*)
{
	edges.set_external_array(nullptr, 0);
}

void Grid::get_associated_sorted(SecureEdgeContainer& edges, Face* f)
{
//	to improve performance, we first check the grid options.
	if(option_is_enabled(FACEOPT_AUTOGENERATE_EDGES
					   | FACEOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we can output the associated array directly
		EdgeContainer& assEdges = m_aaEdgeContainerFACE[f];
		edges.set_external_array(&assEdges.front(), assEdges.size());
	}
	else{
	//	clear the container
		edges.clear();

	//	if no edges are present, we can leave immediately
		if(num<Edge>() == 0)
			return;

	//	get the edges one by one
		uint numEdges = f->num_edges();
		for(uint i = 0; i < numEdges; ++i){
			Edge* e = get_edge(f, i);
			if(e != nullptr)
				edges.push_back(e);
		}
	}
}

void Grid::get_associated_sorted(SecureEdgeContainer& edges, Volume* v)
{
//	to improve performance, we first check the grid options.
	if(option_is_enabled(VOLOPT_AUTOGENERATE_EDGES
							| VOLOPT_STORE_ASSOCIATED_EDGES)
		|| option_is_enabled(VOLOPT_AUTOGENERATE_FACES
							| FACEOPT_AUTOGENERATE_EDGES
							| VOLOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we can output the associated array directly
		EdgeContainer& assEdges = m_aaEdgeContainerVOLUME[v];
		edges.set_external_array(&assEdges.front(), assEdges.size());
	}
	else{
	//	clear the container
		edges.clear();

	//	if no edges are present, we can leave immediately
		if(num<Edge>() == 0)
			return;

	//	get the edges one by one
		uint numEdges = v->num_edges();
		for(uint i = 0; i < numEdges; ++i){
			Edge* e = get_edge(v, i);
			if(e != nullptr)
				edges.push_back(e);
		}
	}
}

void Grid::get_associated_sorted(SecureFaceContainer& faces, Vertex*)
{
	faces.set_external_array(nullptr, 0);
}

void Grid::get_associated_sorted(SecureFaceContainer& faces, Edge*)
{
	faces.set_external_array(nullptr, 0);
}

void Grid::get_associated_sorted(SecureFaceContainer& faces, Volume* v)
{
//	to improve performance, we first check the grid options.
	if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES
					   | VOLOPT_STORE_ASSOCIATED_FACES))
	{
	//	we can output the associated array directly
		FaceContainer& assFaces = m_aaFaceContainerVOLUME[v];
		faces.set_external_array(&assFaces.front(), assFaces.size());
	}
	else{
	//	clear the container
		faces.clear();

	//	if no edges are present, we can leave immediately
		if(num<Face>() == 0)
			return;

	//	get the edges one by one
		uint numFaces = v->num_faces();
		for(uint i = 0; i < numFaces; ++i){
			Face* f = get_face(v, i);
			if(f != nullptr)
				faces.push_back(f);
		}
	}
}

void Grid::get_associated_sorted(SecureVolumeContainer& vols, Vertex*)
{
	vols.set_external_array(nullptr, 0);
}

void Grid::get_associated_sorted(SecureVolumeContainer& vols, Edge*)
{
	vols.set_external_array(nullptr, 0);
}

void Grid::get_associated_sorted(SecureVolumeContainer& vols, Face*)
{
	vols.set_external_array(nullptr, 0);
}

}	//	end of namespace
