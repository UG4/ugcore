//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d16

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


#include "grid.h"
#include "grid_util.h"
#include "common/common.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////
///	this macro helps calling callbacks of different observers.
/**
 * Be sure that callback is a complete function call - including parameters.
 */
#define NOTIFY_OBSERVERS(observerContainer, callback)	{for(Grid::ObserverContainer::iterator iter = observerContainer.begin(); iter != observerContainer.end(); iter++) (*iter)->callback;}

////////////////////////////////////////////////////////////////////////
///	a useful macro that checks if a set of options contains the specified option.
#define OPTIONS_CONTAIN_OPTION(options, option) ((options & option) == option)


namespace ug
{

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	VERTICES
///	creates and removes connectivity data, as specified in optsNew.
void Grid::register_vertex(VertexBase* v, GeometricObject* pParent)
{
//	store the element and register it at the pipe.
	v->m_entryIter = m_elementStorage[VERTEX].m_sectionContainer.insert(static_cast<GeometricObject*>(v), v->shared_pipe_section());
	m_elementStorage[VERTEX].m_attachmentPipe.register_element(v);

//	assign the hash-value
	assign_hash_value(v);

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_vertexObservers, vertex_created(this, v, pParent));
}

void Grid::register_and_replace_element(VertexBase* v, VertexBase* pReplaceMe)
{
	v->m_entryIter = m_elementStorage[VERTEX].m_sectionContainer.insert(static_cast<GeometricObject*>(v), v->shared_pipe_section());
	m_elementStorage[VERTEX].m_attachmentPipe.register_element(v);

//	assign the hash-value
	assign_hash_value(v);

//	pass on values
	pass_on_values(pReplaceMe, v);

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_vertexObservers, vertex_created(this, v, pReplaceMe));
//	inform observers about the replace
	NOTIFY_OBSERVERS(m_vertexObservers, vertex_to_be_replaced(this, pReplaceMe, v));
//	inform observers about the deletion
	NOTIFY_OBSERVERS(m_vertexObservers, vertex_to_be_erased(this, pReplaceMe));

//TODO:	auto-enabling of some options should be optimized (and avoided).
//	all edges, faces and volumes associated with pReplaceMe have to be updated.
//	the option GRIDOPT_VERTEXCENTRIC_INTERCONNECTION has to be enabled.
	if(!option_is_enabled(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION))
	{
		LOG("WARNING in Grid::register_and_replace_element(...) - Vertex: autoenabling grid-option GRIDOPT_VERTEXCENTRIC_INTERCONNECTION.");
		enable_options(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION);
	}

//	update edges
	{
		for(EdgeBaseIterator iter = associated_edges_begin(pReplaceMe);
			iter != associated_edges_end(pReplaceMe); ++iter)
		{
			EdgeBase* e = *iter;
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
	{
		for(FaceIterator iter = associated_faces_begin(pReplaceMe);
			iter != associated_faces_end(pReplaceMe); ++iter)
		{
			Face* f = *iter;
		//	replace the vertex and push f into v's associated faces.
			uint numVrts = f->num_vertices();
			for(uint i = 0; i < numVrts; ++i)
			{
				if(f->vertex(i) == pReplaceMe)
					f->set_vertex(i, v);
			}

			m_aaFaceContainerVERTEX[v].push_back(f);
		}
	}

//	update volumes
	{
		for(VolumeIterator iter = associated_volumes_begin(pReplaceMe);
			iter != associated_volumes_end(pReplaceMe); ++iter)
		{
			Volume* vol = *iter;
		//	replace the vertex and push vol into v's associated volumes.
			uint numVrts = vol->num_vertices();
			for(uint i = 0; i < numVrts; ++i)
			{
				if(vol->vertex(i) == pReplaceMe)
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
	m_elementStorage[VERTEX].m_attachmentPipe.unregister_element(pReplaceMe);
	m_elementStorage[VERTEX].m_sectionContainer.erase(pReplaceMe->m_entryIter, pReplaceMe->shared_pipe_section());
	delete pReplaceMe;
}

void Grid::unregister_vertex(VertexBase* v)
{
//	notify observers that the vertex is being erased
	NOTIFY_OBSERVERS(m_vertexObservers, vertex_to_be_erased(this, v));

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
		//	remove associated volumes
			VolumeIterator iterEnd = associated_volumes_end(v);
			for(VolumeIterator iter = associated_volumes_begin(v); iter != iterEnd;)
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
		//	remove all associated faces.
			FaceIterator iterEnd = associated_faces_end(v);
			for(FaceIterator iter = associated_faces_begin(v); iter != iterEnd;)
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
	//	remove all associated edges.
		EdgeBaseIterator iterEnd = associated_edges_end(v);
		for(EdgeBaseIterator iter = associated_edges_begin(v); iter != iterEnd;)
		{
			EdgeBase* eraseEdge = *iter;
			++iter;
			erase(eraseEdge);
		}
	}

//	remove the element from the storage
	m_elementStorage[VERTEX].m_attachmentPipe.unregister_element(v);
	m_elementStorage[VERTEX].m_sectionContainer.erase(v->m_entryIter, v->shared_pipe_section());
}

void Grid::change_vertex_options(uint optsNew)
{
//	check if associated edge information has to be created or removed.
	if(OPTIONS_CONTAIN_OPTION(optsNew, VRTOPT_STORE_ASSOCIATED_EDGES))
	{
		if(!option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
			vertex_store_associated_edges(true);
	}
	else
	{
		if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
			vertex_store_associated_edges(false);
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

//	store new options
	m_options = m_options | (optsNew & 0x000000FF);
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
			for(EdgeBaseIterator iter = this->edges_begin(); iter != this->edges_end(); iter++)
			{
				EdgeBase* e = *iter;
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
			m_options &= (!VRTOPT_STORE_ASSOCIATED_EDGES);
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
				for(int i = 0; i < numVrts; i++)
					m_aaFaceContainerVERTEX[f->vertex(i)].push_back(f);
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
			m_options &= (!VRTOPT_STORE_ASSOCIATED_FACES);
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
				for(int i = 0; i < numVrts; i++)
					m_aaVolumeContainerVERTEX[v->vertex(i)].push_back(v);
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
			m_options &= (!VRTOPT_STORE_ASSOCIATED_VOLUMES);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	EDGES
///	creates and removes connectivity data, as specified in optsNew.
void Grid::register_edge(EdgeBase* e, GeometricObject* pParent)
{
//	store the element and register it at the pipe.
	e->m_entryIter = m_elementStorage[EDGE].m_sectionContainer.insert(e, e->shared_pipe_section());
	m_elementStorage[EDGE].m_attachmentPipe.register_element(e);

//	register edge at vertices, faces and volumes, if the according options are enabled.
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
		m_aaEdgeContainerVERTEX[e->vertex(0)].push_back(e);
		m_aaEdgeContainerVERTEX[e->vertex(1)].push_back(e);
	}

//	register edges at faces and vice versa
	{
		int switchVar = 0;	// this var will be used to determine what to register where.
		if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
			switchVar = 1;
		if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
			switchVar += 2;

		if(switchVar > 0)
		{
		//	find the faces that contain the edge. search the associated vertices for connected faces.
			FaceIterator iterStart = associated_faces_begin(e->vertex(0));
			FaceIterator iterEnd = associated_faces_end(e->vertex(0));

			EdgeDescriptor ed;
			for(FaceIterator iter = iterStart; iter != iterEnd; iter++)
			{
				Face* f = *iter;
				uint numEdges = f->num_edges();
				for(uint i = 0; i < numEdges; ++i)
				{
					f->edge(i, ed);
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

//	register edges at volumes and vice versa
	{
		int switchVar = 0;	// this var will be used to determine what to register where.
		if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
			switchVar = 1;
		if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
			switchVar += 2;

		if(switchVar > 0)
		{
		//	find the volumes that contain the edge.
			VolumeIterator iterStart = associated_volumes_begin(e->vertex(0));
			VolumeIterator iterEnd = associated_volumes_end(e->vertex(0));

			for(VolumeIterator iter = iterStart; iter != iterEnd; iter++)
			{
				Volume* v = *iter;
				uint numEdges = v->num_edges();
				for(uint i = 0; i < numEdges; ++i)
				{
					EdgeDescriptor ed;
					v->edge(i, ed);
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

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_edgeObservers, edge_created(this, e, pParent));
}

void Grid::register_and_replace_element(EdgeBase* e, EdgeBase* pReplaceMe)
{
//	store the element and register it at the pipe.
	e->m_entryIter = m_elementStorage[EDGE].m_sectionContainer.insert(e, e->shared_pipe_section());
	m_elementStorage[EDGE].m_attachmentPipe.register_element(e);

	pass_on_values(pReplaceMe, e);

//	assign vertices
	e->set_vertex(0, pReplaceMe->vertex(0));
	e->set_vertex(1, pReplaceMe->vertex(1));

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_edgeObservers, edge_created(this, e, pReplaceMe));
//	inform observers about the replace
	NOTIFY_OBSERVERS(m_edgeObservers, edge_to_be_replaced(this, pReplaceMe, e));
//	inform observers about the deletion
	NOTIFY_OBSERVERS(m_edgeObservers, edge_to_be_erased(this, pReplaceMe));

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
	m_elementStorage[EDGE].m_attachmentPipe.unregister_element(pReplaceMe);
	m_elementStorage[EDGE].m_sectionContainer.erase(pReplaceMe->m_entryIter, pReplaceMe->shared_pipe_section());
	delete pReplaceMe;
}

void Grid::unregister_edge(EdgeBase* e)
{
//	notify observers that the edge is being erased
	NOTIFY_OBSERVERS(m_edgeObservers, edge_to_be_erased(this, e));

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
				vector<Volume*>::iterator vIter = vVolumes.begin();

				if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
				{
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
						EdgeContainer::iterator iter = find(m_aaEdgeContainerVOLUME[v].begin(),
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
			vector<Face*>::iterator fIter = vFaces.begin();

			if(option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
			{
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
					EdgeContainer::iterator iter = find(m_aaEdgeContainerFACE[f].begin(),
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
			VertexBase* vrt = e->vertex(i);
			EdgeContainer::iterator iter = find(m_aaEdgeContainerVERTEX[vrt].begin(),
												m_aaEdgeContainerVERTEX[vrt].end(), e);
			if(iter != m_aaEdgeContainerVERTEX[vrt].end())
				m_aaEdgeContainerVERTEX[e->vertex(i)].erase(iter);
		}
	}

//	remove the element from the storage
	m_elementStorage[EDGE].m_attachmentPipe.unregister_element(e);
	m_elementStorage[EDGE].m_sectionContainer.erase(e->m_entryIter, e->shared_pipe_section());
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

//	store new options
	m_options = m_options | (optsNew & 0x0000FF00);
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
					EdgeBase* e = get_edge(f, i);
					if(e != NULL)
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
			m_options &= (!EDGEOPT_STORE_ASSOCIATED_FACES);
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
					v->edge(i, ed);
				//	get the edge that is described by the EdgeDescriptor - if it exists at all.
					EdgeBase* e = get_edge(ed);
					if(e != NULL)
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
			m_options &= (!EDGEOPT_STORE_ASSOCIATED_VOLUMES);
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	FACES
///	creates and removes connectivity data, as specified in optsNew.
void Grid::register_face(Face* f, GeometricObject* pParent)
{
//	store the element and register it at the pipe.
	f->m_entryIter = m_elementStorage[FACE].m_sectionContainer.insert(f, f->shared_pipe_section());
	m_elementStorage[FACE].m_attachmentPipe.register_element(f);

//	register face at vertices
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
	{
		uint numVrts = f->num_vertices();
		for(uint i = 0; i < numVrts; ++i)
			m_aaFaceContainerVERTEX[f->vertex(i)].push_back(f);
	}

//	create edges if FACEOPT_AUTOGENERATE_EDGES is enabled.
//	register face at associated vertices and edges if the according options are set.
	if(option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)
		|| option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES)
		|| option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
	{
	//	loop through the edges of the face and check if they already exist in the grid.
		int numEdges = f->num_edges();
		EdgeDescriptor ed;
		for(int i = 0; i < numEdges; ++i)
		{
			f->edge(i, ed);
			EdgeBase* e = find_edge_in_associated_edges(ed.vertex(0), ed);

			if(e == NULL)
			{
				if(option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
				{
				//	create the edge - regard the parent of f as the parent of the new edge, too.
					e = f->create_edge(i);
					register_edge(e, pParent);
				}
			}
			else
			{
			//	register the face at the edge
				if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
					m_aaFaceContainerEDGE[e].push_back(f);
			//	register the edge at the face.
				if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
					m_aaEdgeContainerFACE[f].push_back(e);
			}
		}
	}

//	register face at existing volumes. find all volumes that contain face
	{
		int switchVar = 0;
		if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
			switchVar = 1;
		if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
			switchVar += 2;

		if(switchVar > 0)
		{
			vector<Volume*> vVols;
			vVols.reserve(2);

			CollectVolumes(vVols, *this, f, false, true);

			for(vector<Volume*>::iterator iter = vVols.begin(); iter != vVols.end(); ++iter)
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

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_faceObservers, face_created(this, f, pParent));
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
	f->m_entryIter = m_elementStorage[FACE].m_sectionContainer.insert(f, f->shared_pipe_section());
	m_elementStorage[FACE].m_attachmentPipe.register_element(f);

	pass_on_values(pReplaceMe, f);

//	assign vertices
	uint numVrts = f->num_vertices();
	{
		for(uint i = 0; i < numVrts; ++i)
			f->set_vertex(i, pReplaceMe->vertex(i));
	}

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_faceObservers, face_created(this, f, pReplaceMe));
//	inform observers about the replace
	NOTIFY_OBSERVERS(m_faceObservers, face_to_be_replaced(this, pReplaceMe, f));
//	inform observers about the deletion
	NOTIFY_OBSERVERS(m_faceObservers, face_to_be_erased(this, pReplaceMe));

//	check if vertices, edges and volumes reference pReplaceMe.
//	if so, correct those references.
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
	{
		for(uint i = 0; i < numVrts; ++i)
			replace(associated_faces_begin(f->vertex(i)),
					associated_faces_end(f->vertex(i)),
					pReplaceMe, f);
	}

	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
	//	collect all edges that are associated with pReplaceMe
		vector<EdgeBase*> vEdges;
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
	m_elementStorage[FACE].m_attachmentPipe.unregister_element(pReplaceMe);
	m_elementStorage[FACE].m_sectionContainer.erase(pReplaceMe->m_entryIter, pReplaceMe->shared_pipe_section());
	delete pReplaceMe;
}

void Grid::unregister_face(Face* f)
{
//	notify observers that the face is being erased
	NOTIFY_OBSERVERS(m_faceObservers, face_to_be_erased(this, f));

//	remove or disconnect from volumes
	if(num_volumes() > 0)
	{
		if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES) ||
			option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
		{
			vector<Volume*> vVolumes;
			CollectVolumes(vVolumes, *this, f, true);
			vector<Volume*>::iterator vIter= vVolumes.begin();

			if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
			{
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
					FaceContainer::iterator iter = find(m_aaFaceContainerVOLUME[vol].begin(),
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
			EdgeBase* e = get_edge(f, i);
			if(e != NULL)
			{
				FaceContainer::iterator iter = find(m_aaFaceContainerEDGE[e].begin(),
													m_aaFaceContainerEDGE[e].end(), f);
				if(iter != m_aaFaceContainerEDGE[e].end())
					m_aaFaceContainerEDGE[e].erase(iter);
			}
		}
	}

//	disconnect from vertices
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
	{
		for(uint i = 0; i < f->num_vertices(); ++i)
		{
			VertexBase* vrt = f->vertex(i);
			FaceContainer::iterator iter = find(m_aaFaceContainerVERTEX[vrt].begin(),
												m_aaFaceContainerVERTEX[vrt].end(), f);
			if(iter != m_aaFaceContainerVERTEX[vrt].end())
				m_aaFaceContainerVERTEX[vrt].erase(iter);
		}
	}

//	remove the element from the storage
	m_elementStorage[FACE].m_attachmentPipe.unregister_element(f);
	m_elementStorage[FACE].m_sectionContainer.erase(f->m_entryIter, f->shared_pipe_section());
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

//	store new options
	m_options = m_options | (optsNew & 0x00FF0000);
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

		//	if EDGEOPT_STORE_ASSOCIATED_FACES is enabled, this is as simple
		//	as to iterate through all edges and store them at their
		//	associated faces.
		//	if not we have to iterate through the faces and store
		//	each of their edges..
			if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
			{
			//	iterate through the edges
				for(EdgeBaseIterator iter = edges_begin(); iter != edges_end(); iter++)
				{
					EdgeBase* e = *iter;
				//	iterate through the edges associated faces
					for(FaceIterator fIter = associated_faces_begin(e);
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
				//	iterate through the edges of the face
					uint numEdges = f->num_edges();
					for(uint i = 0; i < numEdges; ++i)
					{
					//	get the i-th edge that is described by the EdgeDescriptor - if it exists at all.
						EdgeBase* e = get_edge(f, i);
						if(e != NULL)
							m_aaEdgeContainerFACE[f].push_back(e);
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
			m_options &= (!FACEOPT_STORE_ASSOCIATED_EDGES);
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
					FaceIterator iterEnd = associated_faces_end(v);
					for(FaceIterator iter = associated_faces_begin(v); iter != iterEnd; iter++)
						m_aaVolumeContainerFACE[*iter].push_back(v);
				}
				else
				{
					int numFaces = v->num_faces();
					for(int i = 0; i < numFaces; ++i)
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
			m_options &= (!FACEOPT_STORE_ASSOCIATED_VOLUMES);
		}
	}
}

void Grid::face_autogenerate_edges(bool bAutogen)
{
	if(bAutogen)
	{
		if(!option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		{
		//	generate all missing edges now!
			for(FaceIterator iter = faces_begin(); iter != faces_end(); iter++)
			{
				Face* f = *iter;
			//	check for each edge of the face if it exists. if not create it and link it to the face.
				int numEdges = f->num_edges();

				for(int i = 0; i < numEdges; ++i)
				{
					EdgeBase* e = get_edge(f, i);

					if(e == NULL)
					{
					//	create the edge
						e = f->create_edge(i);
						register_edge(e);
					}
				}
			}

		//	store the option
			m_options |= FACEOPT_AUTOGENERATE_EDGES;
		}
	}
	else
	{
	//	stop auto-generation
		m_options &= (!FACEOPT_AUTOGENERATE_EDGES);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
//	VOLUMES
///	creates and removes connectivity data, as specified in optsNew.
void Grid::register_volume(Volume* v, GeometricObject* pParent)
{
//	store the element and register it at the pipe.
	v->m_entryIter = m_elementStorage[VOLUME].m_sectionContainer.insert(v, v->shared_pipe_section());
	m_elementStorage[VOLUME].m_attachmentPipe.register_element(v);

//	create edges and faces if the according options are enabled.
//	register the volume at the associated vertices, edges and faces, if the according options are enabled.

//	register volume at vertices
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
		uint numVrts = v->num_vertices();
		for(uint i = 0; i < numVrts; ++i)
			m_aaVolumeContainerVERTEX[v->vertex(i)].push_back(v);
	}

//	register with edges. autogenerate missing edges if option is set
	if(option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)
		|| option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES)
		|| option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
	{
	//	loop through all edges of the volume. check if we have to create it.
	//	register the volume with the edges and with versa after that.
		uint numEdges = v->num_edges();
		EdgeDescriptor ed;
		for(uint i = 0; i < numEdges; ++i)
		{
			v->edge(i, ed);
			EdgeBase* e = find_edge_in_associated_edges(ed.vertex(0), ed);

			if(e == NULL)
			{
				if(option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
				{
				//	create the edge
					e = v->create_edge(i);
					register_edge(e, pParent);
				}
			}
			else
			{
				if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
					m_aaVolumeContainerEDGE[e].push_back(v);
				if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
					m_aaEdgeContainerVOLUME[v].push_back(e);
			}
		}
	}

//	register with faces. autogenerate them if they do not already exist.
	if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES)
		|| option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES)
		|| option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
	{
	//	iterate through the faces of the volume. create them if demanded.
	//	register faces at the volume and vice versa.
		uint numFaces = v->num_faces();
		FaceDescriptor fd;
		for(uint i = 0; i < numFaces; ++i)
		{
			v->face(i, fd);
			Face* f = find_face_in_associated_faces(fd.vertex(0), fd);

			if(f == NULL)
			{
				if(option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
				{
					f = v->create_face(i);
					register_face(f, pParent);
				}
			}
			else
			{
				if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
					m_aaVolumeContainerFACE[f].push_back(v);
				if(option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
					m_aaFaceContainerVOLUME[v].push_back(f);
			}
		}
	}

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_volumeObservers, volume_created(this, v, pParent));
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
	v->m_entryIter = m_elementStorage[VOLUME].m_sectionContainer.insert(v, v->shared_pipe_section());
	m_elementStorage[VOLUME].m_attachmentPipe.register_element(v);

	pass_on_values(pReplaceMe, v);

//	assign vertices
	uint numVrts = v->num_vertices();
	{
		for(uint i = 0; i < numVrts; ++i)
			v->set_vertex(i, pReplaceMe->vertex(i));
	}

//	inform observers about the creation
	NOTIFY_OBSERVERS(m_volumeObservers, volume_created(this, v, pReplaceMe));
//	inform observers about the replace
	NOTIFY_OBSERVERS(m_volumeObservers, volume_to_be_replaced(this, pReplaceMe, v));
//	inform observers about the deletion
	NOTIFY_OBSERVERS(m_volumeObservers, volume_to_be_erased(this, pReplaceMe));

//	check if vertices, edges and faces reference pReplaceMe.
//	if so, correct those references.
	if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))
	{
		for(uint i = 0; i < numVrts; ++i)
			replace(associated_volumes_begin(v->vertex(i)),
					associated_volumes_end(v->vertex(i)),
					pReplaceMe, v);
	}

	if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
	//	collect all edges that are associated with pReplaceMe
		vector<EdgeBase*> vEdges;
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
	m_elementStorage[VOLUME].m_attachmentPipe.unregister_element(pReplaceMe);
	m_elementStorage[VOLUME].m_sectionContainer.erase(pReplaceMe->m_entryIter, pReplaceMe->shared_pipe_section());
	delete pReplaceMe;
}

void Grid::unregister_volume(Volume* v)
{
//	notify observers that the face is being erased
	NOTIFY_OBSERVERS(m_volumeObservers, volume_to_be_erased(this, v));

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

			if(f != NULL)
			{
				VolumeContainer::iterator iter = find(m_aaVolumeContainerFACE[f].begin(),
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
			EdgeBase* e = get_edge(v, i);
			if(e != NULL)
			{
				VolumeContainer::iterator iter = find(m_aaVolumeContainerEDGE[e].begin(),
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
		for(uint i = 0; i < numVertices; ++i)
		{
		//	find the correct entry
			VertexBase* vrt = v->vertex(i);
			VolumeContainer::iterator iter = find(m_aaVolumeContainerVERTEX[vrt].begin(),
													m_aaVolumeContainerVERTEX[vrt].end(), v);
			if(iter != m_aaVolumeContainerVERTEX[vrt].end())
				m_aaVolumeContainerVERTEX[vrt].erase(iter);
		}
	}

//	remove the element from the storage
	m_elementStorage[VOLUME].m_attachmentPipe.unregister_element(v);
	m_elementStorage[VOLUME].m_sectionContainer.erase(v->m_entryIter, v->shared_pipe_section());
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
//	store new options
	m_options = m_options | (optsNew & 0xFF000000);
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

		//	if EDGEOPT_STORE_ASSOCIATED_VOLUMES is enabled, this is as simple
		//	as to iterate through all edges and store them at their
		//	associated volumes.
		//	if not we have to iterate through the volumes and store
		//	each of their edges..
			if(option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
			{
			//	iterate through the edges
				for(EdgeBaseIterator iter = edges_begin(); iter != edges_end(); iter++)
				{
					EdgeBase* e = *iter;
				//	iterate through the edges associated faces
					for(VolumeIterator vIter = associated_volumes_begin(e);
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
				//	iterate through the edges of the volume
					int numEdges = v->num_edges();
					for(int i = 0; i < numEdges; ++i)
					{
					//	get the edge-descriptor
						EdgeBase* e = get_edge(v, i);
						if(e != NULL)
							m_aaEdgeContainerVOLUME[v].push_back(e);
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
			m_options &= (!VOLOPT_STORE_ASSOCIATED_EDGES);
		}
	}
}

void Grid::volume_store_associated_faces(bool bStoreIt)
{
	if(bStoreIt)
	{
		if(!option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
		{
		//	store associated edges
			attach_to_volumes(m_aFaceContainer);
			m_aaFaceContainerVOLUME.access(*this, m_aFaceContainer);

		//	if FACEOPT_STORE_ASSOCIATED_VOLUMES is enabled, this is as simple
		//	as to iterate through all faces and store them at their
		//	associated volumes.
		//	if not we have to iterate through the volumes and store
		//	each of their faces..
			if(option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
			{
			//	iterate through the edges
				for(FaceIterator iter = faces_begin(); iter != faces_end(); iter++)
				{
					Face* f = *iter;
				//	iterate through the faces associated volumes
					for(VolumeIterator vIter = associated_volumes_begin(f);
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
					for(int i = 0; i < numFaces; ++i)
					{
						v->face(i, fd);
						Face* f = find_face_in_associated_faces(fd.vertex(0), fd);

						if(f)
							m_aaVolumeContainerFACE[f].push_back(v);
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
			m_options &= (!VOLOPT_STORE_ASSOCIATED_FACES);
		}
	}
}

void Grid::volume_autogenerate_edges(bool bAutogen)
{
	if(bAutogen)
	{
		if(!option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
		{
		//	generate all missing edges now!
			for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
			{
				Volume* v = *iter;
			//	check for each edge of the volume if it exists. if not create it and link it to the face.
				uint numEdges = v->num_edges();

				for(uint i = 0; i < numEdges; ++i)
				{
					EdgeBase* e = get_edge(v, i);

					if(e == NULL)
					{
					//	create the edge
						e = v->create_edge(i);
						register_edge(e);
					}
				}
			}

		//	store the option
			m_options |= VOLOPT_STORE_ASSOCIATED_EDGES;
		}
	}
	else
	{
	//	stop auto-generation
		m_options &= (!VOLOPT_STORE_ASSOCIATED_EDGES);
	}
}

void Grid::volume_autogenerate_faces(bool bAutogen)
{
	if(bAutogen)
	{
		if(!option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
		{
		//	generate all missing edges now!
			for(VolumeIterator iter = volumes_begin(); iter != volumes_end(); iter++)
			{
				Volume* v = *iter;
			//	check for each face of the volume if it exists. if not create it and link it to the face.
				uint numFaces = v->num_faces();

				for(uint i = 0; i < numFaces; ++i)
				{
					Face* f = get_face(v, i);

					if(f == NULL)
					{
					//	create the face
						f = v->create_face(i);
						register_face(f);
					}
				}
			}

		//	store the option
			m_options |= VOLOPT_AUTOGENERATE_FACES;
		}
	}
	else
	{
	//	stop auto-generation
		m_options &= (!VOLOPT_AUTOGENERATE_FACES);
	}
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	replace_vertex
bool Grid::replace_vertex(VertexBase* vrtOld, VertexBase* vrtNew)
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
	vector<EdgeBase*> vEdges;
	vector<Face*> vFaces;
	vector<Volume*> vVolumes;

	EdgeDescriptor ed;
	FaceDescriptor fd;
	VolumeDescriptor vd;

//	since vrtOld will be erased, we will notify all observers
	NOTIFY_OBSERVERS(m_vertexObservers, vertex_to_be_erased(this, vrtOld));

//	EDGES
	if(num_edges() > 0)
	{
		autoenable_option(VRTOPT_STORE_ASSOCIATED_EDGES,
						"Grid::replace_vertex(...)",
						"VRTOPT_STORE_ASSOCIATED_EDGES");

		EdgeBaseIterator iter = associated_edges_begin(vrtOld);
		EdgeBaseIterator iterEnd = associated_edges_end(vrtOld);
		while(iter != iterEnd)
		{
			EdgeBase* e = *iter;
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
				EdgeBase* eNew = get_edge(ed);
				if(eNew)
				{
				//	The edge will be removed. Notify observers.
					NOTIFY_OBSERVERS(m_edgeObservers, edge_to_be_erased(this, e));

				//	Since we want to avoid deletion of associated elements
				//	we have to handle the erasure of e by our self.
					bReplaceVertex = false;

				//	before the removal we will replace its entry in associated-
				//	element-lists with the pointer to the existing edge.
					if(option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
					{
					//	unregister e from the vertex with which e connects vrtOld.
						EdgeContainer& ec = m_aaEdgeContainerVERTEX[ed.vertex(0)];
						EdgeContainer::iterator tmpI = find(ec.begin(), ec.end(), e);
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
					m_elementStorage[EDGE].m_attachmentPipe.unregister_element(e);
					m_elementStorage[EDGE].m_sectionContainer.erase(e->m_entryIter, e->shared_pipe_section());
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

		FaceIterator iter = associated_faces_begin(vrtOld);
		FaceIterator iterEnd = associated_faces_end(vrtOld);
		while(iter != iterEnd)
		{
			Face* f = *iter;
			++iter;

		//	if eraseDoubleElementes is enabled and the new face would
		//	already exist, we wont replace the vertex in it.
			bool bReplaceVertex = true;

			if(eraseDoubleElements)
			{
			//	create a face-descriptor of the face that will be created.
				uint numVrts = f->num_vertices();
				fd.set_num_vertices(numVrts);
				for(uint i = 0; i < numVrts; ++i)
				{
					if(f->vertex(i) == vrtOld)
						fd.set_vertex(i, vrtNew);
					else
						fd.set_vertex(i, f->vertex(i));
				}

			//	check if this face already exists.
				Face* fNew = get_face(fd);
				if(fNew)
				{
					LOG("double face" << endl);
				//	The face will be removed. Notify observers.
					NOTIFY_OBSERVERS(m_faceObservers, face_to_be_erased(this, f));

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
							if(f->vertex(i) != vrtOld)
							{
								FaceContainer& fc = m_aaFaceContainerVERTEX[f->vertex(i)];
								FaceContainer::iterator tmpI = find(fc.begin(), fc.end(), f);
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
						LOG("removing face. num associated edges: " << vEdges.size() << endl);
						for(uint i = 0; i < vEdges.size(); ++i)
						{
							FaceContainer& fc = m_aaFaceContainerEDGE[vEdges[i]];
							FaceContainer::iterator tmpI = find(fc.begin(), fc.end(), f);
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
					m_elementStorage[FACE].m_attachmentPipe.unregister_element(f);
					m_elementStorage[FACE].m_sectionContainer.erase(f->m_entryIter, f->shared_pipe_section());
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
					if(f->vertex(i) == vrtOld)
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
						f->edge(i, ed);
						if(ed.vertex(0) == vrtNew || ed.vertex(1) == vrtNew)
						{
							EdgeBase* tEdge = get_edge(ed);
							if(tEdge)
							{
								FaceContainer& fc = m_aaFaceContainerEDGE[tEdge];
								FaceContainer::iterator tmpI = find(fc.begin(), fc.end(), f);
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

		VolumeIterator iter = associated_volumes_begin(vrtOld);
		VolumeIterator iterEnd = associated_volumes_end(vrtOld);
		while(iter != iterEnd)
		{
			Volume* v = *iter;
			++iter;

		//	if eraseDoubleElementes is enabled and the new face would
		//	already exist, we wont replace the vertex in it.
			bool bReplaceVertex = true;

			if(eraseDoubleElements)
			{
			//	create a volume-descriptor of the volume that will be created.
				uint numVrts = v->num_vertices();
				vd.set_num_vertices(numVrts);
				for(uint i = 0; i < numVrts; ++i)
				{
					if(v->vertex(i) == vrtOld)
						vd.set_vertex(i, vrtNew);
					else
						vd.set_vertex(i, v->vertex(i));
				}

			//	check if this volume already exists.
				Volume* vNew = get_volume(vd);
				if(vNew)
				{
				//	The volume will be removed. Notify observers.
					NOTIFY_OBSERVERS(m_volumeObservers, volume_to_be_erased(this, v));

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
							if(v->vertex(i) != vrtOld)
							{
								VolumeContainer& vc = m_aaVolumeContainerVERTEX[v->vertex(i)];
								VolumeContainer::iterator tmpI = find(vc.begin(), vc.end(), v);
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
							VolumeContainer::iterator tmpI = find(vc.begin(), vc.end(), v);
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
							VolumeContainer::iterator tmpI = find(vc.begin(), vc.end(), v);
							if(tmpI != vc.end())
								vc.erase(tmpI);
						}
					}


				//	we can now remove v from the storage.
					m_elementStorage[VOLUME].m_attachmentPipe.unregister_element(v);
					m_elementStorage[VOLUME].m_sectionContainer.erase(v->m_entryIter, v->shared_pipe_section());
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
					if(v->vertex(i) == vrtOld)
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
						v->edge(i, ed);
						if(ed.vertex(0) == vrtNew || ed.vertex(1) == vrtNew)
						{
							EdgeBase* tEdge = get_edge(ed);
							if(tEdge)
							{
								VolumeContainer& vc = m_aaVolumeContainerEDGE[tEdge];
								VolumeContainer::iterator tmpI = find(vc.begin(), vc.end(), v);
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
						v->face(i, fd);
					//	check whether fd contains vrtNew
						bool bContainsVrtNew = false;
						for(uint j = 0; j < fd.num_vertices(); ++j)
						{
							if(fd.vertex(j) == vrtNew)
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
								VolumeContainer::iterator tmpI = find(vc.begin(), vc.end(), v);
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
	m_elementStorage[VERTEX].m_attachmentPipe.unregister_element(vrtOld);
	m_elementStorage[VERTEX].m_sectionContainer.erase(vrtOld->m_entryIter, vrtOld->shared_pipe_section());
	delete vrtOld;

	return true;
}

////////////////////////////////////////////////////////////////////////
//	replace_vertex_is_valid
bool Grid::replace_vertex_is_valid(VertexBase* vrtOld, VertexBase* vrtNew)
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

		EdgeBaseIterator iterEnd = associated_edges_end(vrtOld);
		for(EdgeBaseIterator iter = associated_edges_begin(vrtOld);
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

		FaceIterator iterEnd = associated_faces_end(vrtOld);
		for(FaceIterator iter = associated_faces_begin(vrtOld);
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

		VolumeIterator iterEnd = associated_volumes_end(vrtOld);
		for(VolumeIterator iter = associated_volumes_begin(vrtOld);
			iter != iterEnd; ++iter)
		{
			if(VolumeContains(*iter, vrtNew))
				return false;
		}
	}

	return true;
}
}	//	end of namespace libGrid
