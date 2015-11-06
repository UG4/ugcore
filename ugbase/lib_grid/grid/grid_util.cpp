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

#include "grid_util.h"
#include "grid.h"
#include <iostream>

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CompareVertices
bool CompareVertices(const EdgeVertices* ev1,
					const EdgeVertices* ev2)
{
	if((ev1->vertex(0) == ev2->vertex(0) && ev1->vertex(1) == ev2->vertex(1)) ||
		(ev1->vertex(0) == ev2->vertex(1) && ev1->vertex(1) == ev2->vertex(0)))
		return true;
	return false;
}

bool CompareVertices(const FaceVertices* fv1,
					const FaceVertices* fv2)
{
	uint numVrts = fv1->num_vertices();

	if(numVrts != fv2->num_vertices())
		return false;

	FaceVertices::ConstVertexArray vrts1 = fv1->vertices();
	FaceVertices::ConstVertexArray vrts2 = fv2->vertices();

	for(uint i = 0; i < numVrts; ++i)
	{
		uint j;
		for(j = 0; j < numVrts; ++j)
		{
			if(vrts1[i] == vrts2[j])
				break;
		}

	//	check whether we found a matching vertex
		if(j == numVrts)
			return false;
	}

	return true;
}

bool CompareVertices(const VolumeVertices* vv1,
					const VolumeVertices* vv2)
{
	uint numVrts = vv1->num_vertices();

	if(numVrts != vv2->num_vertices())
		return false;

	VolumeVertices::ConstVertexArray vrts1 = vv1->vertices();
	VolumeVertices::ConstVertexArray vrts2 = vv2->vertices();

	for(uint i = 0; i < numVrts; ++i)
	{
		uint j;
		for(j = 0; j < numVrts; ++j)
		{
			if(vrts1[i] == vrts2[j])
				break;
		}

	//	check whether we found a matching vertex
		if(j == numVrts)
			return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
//	CollectVertices
///	Collects all vertices.
void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid, Vertex* v, bool clearContainer)
{
	// clear container if wanted
	if(clearContainer)
		vVertexOut.clear();

// add vertex pointers to container
	vVertexOut.push_back(v);
}

///	Collects all vertices.
void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid, Edge* e, bool clearContainer)
{
	// clear container if wanted
	if(clearContainer)
		vVertexOut.clear();

	// resize container
	const size_t numVertex = e->num_vertices();

	// add vertex pointers to container
	for(size_t i = 0; i < numVertex; ++i)
		vVertexOut.push_back(e->vertex(i));
}

///	Collects all vertices.
void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid, Face* f, bool clearContainer)
{
	// clear container if wanted
	if(clearContainer)
		vVertexOut.clear();

	// resize container
	const size_t numVertex = f->num_vertices();
	FaceVertices::ConstVertexArray vrts= f->vertices();

	// add vertex pointers to container
	for(size_t i = 0; i < numVertex; ++i)
		vVertexOut.push_back(vrts[i]);
}

///	Collects all vertices.
void CollectVertices(std::vector<Vertex*>& vVertexOut, Grid& grid, Volume* v, bool clearContainer)
{
	// clear container if wanted
	if(clearContainer)
		vVertexOut.clear();

	// resize container
	const size_t numVertex = v->num_vertices();
	VolumeVertices::ConstVertexArray vrts= v->vertices();

	// add vertex pointers to container
	for(size_t i = 0; i < numVertex; ++i)
		vVertexOut.push_back(vrts[i]);
}



///////////////////////////////////////////////////////////////////////////////
//	CollectEdgesSorted
///////////////////////////////////////////////////////////////////////////////

///	Collects all edges of a vertex, thus, none.
void CollectEdgesSorted(vector<Edge*>& vEdgesOut, Grid& grid, Vertex* v, bool clearContainer)
{
	vEdgesOut.clear();
}

///	Collects all edges. (Returns the edge itself)
void CollectEdgesSorted(vector<Edge*>& vEdgesOut, Grid& grid, Edge* e, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

	vEdgesOut.push_back(e);
}

///	Collects all edges which exist in the given grid and which are part of the given face in the order defined by the reference elements.
void CollectEdgesSorted(vector<Edge*>& vEdgesOut, Grid& grid, Face* f, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	to improve performance, we first check the grid options.
	if(grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES
							| FACEOPT_STORE_ASSOCIATED_EDGES))
	{
	//	the edges can be accessed in a sorted way through iterators
		Grid::AssociatedEdgeIterator end = grid.associated_edges_end(f);
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(f);
			iter != end; ++iter)
		{
			vEdgesOut.push_back(*iter);
		}
	}
	else{
	//	if no edges are present, we can leave immediately
		if(grid.num<Edge>() == 0)
			return;
	//	second-best: use GetEdge in order to find the queried edges
		uint numEdges = f->num_edges();
		for(uint i = 0; i < numEdges; ++i)
		{
			Edge* e = grid.get_edge(f, i);
			if(e != NULL)
				vEdgesOut.push_back(e);
		}
	}
}

///	Collects all edges that exist in the given grid are part of the given volume in the order defined by the reference elements.
void CollectEdgesSorted(vector<Edge*>& vEdgesOut, Grid& grid, Volume* v, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	to improve performance, we first check the grid options.
	if(grid.option_is_enabled(VOLOPT_AUTOGENERATE_EDGES
							| VOLOPT_STORE_ASSOCIATED_EDGES)
		|| grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES
							| FACEOPT_AUTOGENERATE_EDGES
							| VOLOPT_STORE_ASSOCIATED_EDGES))
	{
	//	the edges can be accessed in a sorted way through iterators
		Grid::AssociatedEdgeIterator end = grid.associated_edges_end(v);
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
			iter != end; ++iter)
		{
			vEdgesOut.push_back(*iter);
		}
	}
	else{
	//	if no edges are present, we can leave immediately
		if(grid.num<Edge>() == 0)
			return;

	//	second best: use GetEdge in order to find the queried edges.
		uint numEdges = v->num_edges();
		for(uint i = 0; i < numEdges; ++i)
		{
			Edge* e = grid.get_edge(v, i);
			if(e != NULL)
				vEdgesOut.push_back(e);
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
//	CollectEdges
///////////////////////////////////////////////////////////////////////////////

///	Collects all edges which exist in the given grid and which are part of the given vertex.
void CollectEdges(std::vector<Edge*>& vEdgesOut, Grid& grid, Vertex* vrt, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	if no edges are present, we can leave immediately
	if(grid.num<Edge>() == 0)
		return;

	Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(vrt);
	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(vrt);
		iter != iterEnd; ++iter)
	{
		vEdgesOut.push_back(*iter);
	}
}

///	Collects all edges. (Returns the edge itself)
void CollectEdges(vector<Edge*>& vEdgesOut, Grid& grid, Edge* e, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	if no edges are present, we can leave immediately
	if(grid.num<Edge>() == 0)
		return;

	vEdgesOut.push_back(e);
}

///	Collects all edges which exist in the given grid and which are part of the given face.
void CollectEdges(vector<Edge*>& vEdgesOut, Grid& grid, Face* f, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	if no edges are present, we can leave immediately
	if(grid.num<Edge>() == 0)
		return;

//	best-option: FACEOPT_STORE_ASSOCIATED_EDGES
	if(grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
	{
	//	ok. simply push them into the container.
		Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(f);
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(f);
			iter != iterEnd; ++iter)
		{
			vEdgesOut.push_back(*iter);
		}

	//	everything is done. exit.
		return;
	}

//	second-best: use GetEdge in order to find the queried edges
	uint numEdges = f->num_edges();
	for(uint i = 0; i < numEdges; ++i)
	{
		Edge* e = grid.get_edge(f, i);
		if(e != NULL)
			vEdgesOut.push_back(e);
	}
}

///	Collects all edges that exist in the given grid are part of the given volume.
void CollectEdges(vector<Edge*>& vEdgesOut, Grid& grid, Volume* v, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	if no edges are present, we can leave immediately
	if(grid.num<Edge>() == 0)
		return;

//	best option: VOLOPT_STORE_ASSOCIATED_EDGES
	if(grid.option_is_enabled(VOLOPT_STORE_ASSOCIATED_EDGES))
	{
	//	iterate through the associated edges and push them into the container.
		Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(v);
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
			iter != iterEnd; ++iter)
		{
			vEdgesOut.push_back(*iter);
		}

	//	everything is done. exit.
		return;
	}

//	second best: use GetEdge in order to find the queried edges.
	uint numEdges = v->num_edges();
	for(uint i = 0; i < numEdges; ++i)
	{
		Edge* e = grid.get_edge(v, i);
		if(e != NULL)
			vEdgesOut.push_back(e);
	}

}

////////////////////////////////////////////////////////////////////////////////
//	CollectFacesSorted
////////////////////////////////////////////////////////////////////////////////

///	Collects all Faces of a Vertex, thus, none.
void CollectFacesSorted(vector<Face*>& vFacesOut, Grid& grid, Vertex* v, bool clearContainer)
{
	vFacesOut.clear();
}

///	Collects all faces and returns them in the order prescribed by the reference element.
void CollectFacesSorted(vector<Face*>& vFacesOut, Grid& grid, Edge* e, bool clearContainer)
{
	vFacesOut.clear();
}

///	Collects all faces and returns them in the order prescribed by the reference element.
void CollectFacesSorted(vector<Face*>& vFacesOut, Grid& grid, Face* f, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

	if(f != NULL)
		vFacesOut.push_back(f);
}

///	Collects all faces and returns them in the order prescribed by the reference element.
void CollectFacesSorted(vector<Face*>& vFacesOut, Grid& grid, Volume* v, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

//	to improve performance, we first check the grid options.
	if(grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES
							| VOLOPT_STORE_ASSOCIATED_FACES))
	{
	//	the faces can be accessed in a sorted way through iterators
		Grid::AssociatedFaceIterator end = grid.associated_faces_end(v);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v);
			iter != end; ++iter)
		{
			vFacesOut.push_back(*iter);
		}
	}
	else{
	//	if no faces are present, we can leave immediately
		if(grid.num<Face>() == 0)
			return;

		uint numFaces = v->num_faces();
		for(uint i = 0; i < numFaces; ++i)
		{
			Face* f = grid.get_face(v, i);
			if(f != NULL)
				vFacesOut.push_back(f);
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	CollectFaces
///	Collects all faces that exist in the given grid which contain the given vertex.
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, Vertex* vrt, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

//	if no faces are present, we can leave immediately
	if(grid.num<Face>() == 0)
		return;

	Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrt);
	for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt);
		iter != iterEnd; ++iter)
	{
		vFacesOut.push_back(*iter);
	}
}

///	Collects all faces that exist in the given grid which contain the given edge.
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, Edge* e, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

//	if no faces are present, we can leave immediately
	if(grid.num<Face>() == 0)
		return;

//	best option: EDGEOPT_STORE_ASSOCIATED_FACES
	if(grid.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(e);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(e);
			iter != iterEnd; ++iter)
		{
			vFacesOut.push_back(*iter);
		}
		return;
	}

//	second best: iterate through all faces associated with the first end-point of e
//	and check for each if it contains e. If so push it into the container.
	{
		Vertex* v1 = e->vertex(0);
		Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(v1);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v1);
			iter != iterEnd; ++iter)
		{
			if(FaceContains(*iter, e))
				vFacesOut.push_back(*iter);
		}
	}
}

///	Collects all faces. (Returns the face itself)
void CollectFaces(vector<Face*>& vFacesOut, Grid& grid, Face* e, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

//	if no faces are present, we can leave immediately
	if(grid.num<Face>() == 0)
		return;

	vFacesOut.push_back(e);
}

///	Collects all faces that exist in the given grid are part of the given volume.
void CollectFaces(vector<Face*>& vFacesOut, Grid& grid, Volume* v, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

//	if no faces are present, we can leave immediately
	if(grid.num<Face>() == 0)
		return;

//	best option: VOLOPT_STORE_ASSOCIATED_FACES
	if(grid.option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
	{
	//	iterate through the faces and add them to the container
		Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(v);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v);
			iter != iterEnd; ++iter)
		{
			vFacesOut.push_back(*iter);
		}

		return;
	}

//	second best: use FindFace in order to find the queried faces
	uint numFaces = v->num_faces();
	for(uint i = 0; i < numFaces; ++i)
	{
		Face* f = grid.get_face(v, i);
		if(f != NULL)
			vFacesOut.push_back(f);
	}
}

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the two given vertices
bool FaceContains(Face* f, EdgeVertices* ev)
{
	uint numEdges = f->num_edges();
	EdgeDescriptor ed;
	for(uint i = 0; i < numEdges; ++i)
	{
		f->edge_desc(i, ed);
		if(CompareVertices(ev, &ed))
			return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the given vertex
bool FaceContains(FaceVertices* f, Vertex* v)
{
	uint numVrts = f->num_vertices();
	Face::ConstVertexArray vrts = f->vertices();

	for(uint i = 0; i < numVrts; ++i)
	{
		if(vrts[i] == v)
			return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Vertex* vrt, bool clearContainer)
{
	if(clearContainer)
		vVolumesOut.clear();

//	if no volumes are present, we can leave immediately
	if(grid.num<Volume>() == 0)
		return;

	Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(vrt);
	for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(vrt);
		iter != iterEnd; ++iter)
	{
		vVolumesOut.push_back(*iter);
	}
}

////////////////////////////////////////////////////////////////////////
//	CollectVolumes
///	Collects all volumes that exist in the given grid which contain the given edge.
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Edge* e, bool clearContainer)
{
	if(clearContainer)
		vVolumesOut.clear();

//	if no volumes are present, we can leave immediately
	if(grid.num<Volume>() == 0)
		return;

//	best option: EDGEOPT_STORE_ASSOCIATED_VOLUMES
	if(grid.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(e);
		for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(e);
			iter != iterEnd; ++iter)
		{
			vVolumesOut.push_back(*iter);
		}
		return;
	}

//	second best: iterate through all volumes that are associated with the first vertex of e.
//	check for each if it contains e. if so then store it in the container.
	Vertex* v1 = e->vertex(0);
	Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(v1);
	for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(v1);
		iter != iterEnd; ++iter)
	{
		if(VolumeContains(*iter, e))
			vVolumesOut.push_back(*iter);
	}
}

///	Collects all volumes that exist in the given grid which contain the given face.
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, Face* f, bool clearContainer, bool ignoreAssociatedVolumes)
{
	if(clearContainer)
		vVolumesOut.clear();

//	if no volumes are present, we can leave immediately
	if(grid.num<Volume>() == 0)
		return;

	if(!ignoreAssociatedVolumes)
	{
	//	best option: FACEOPT_STORE_ASSOCIATED_VOLUMES
		if(grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
		{
			Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(f);
			for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(f);
				iter != iterEnd; ++iter)
			{
				vVolumesOut.push_back(*iter);
			}
			return;
		}
	}
/*
//	second best: use associated volumes of edges.
	if(grid.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_VOLUMES))
	{
	//	get an edge of the volume. check if associated volumes of that edge
	//	contain the face.
	}
*/

//	worst: iterate through all volumes which are connected to the first vertex of f.
//	check for each if it contains f. If so, store it in the container.
//	to make things a little faster we'll check if the associated volumes contain a second vertex of f.
	Vertex* v1 = f->vertex(0);

	Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(v1);
	for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(v1);
		iter != iterEnd; ++iter)
	{
		Volume* v = *iter;
	//	check if v contains v2
		if(VolumeContains(v, f->vertex(1)))
		{
			if(VolumeContains(v, f->vertex(2)))
			{
				if(VolumeContains(v, f))
					vVolumesOut.push_back(*iter);
			}
		}
	}
}

///	Collects all volumes. (Returns the volume itself)
void CollectVolumes(vector<Volume*>& vVolumesOut, Grid& grid, Volume* v, bool clearContainer)
{
	if(clearContainer)
		vVolumesOut.clear();

//	if no volumes are present, we can leave immediately
	if(grid.num<Volume>() == 0)
		return;

	vVolumesOut.push_back(v);
}

///	Collects all volumes that exist in the given grid which contain the given face.
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, FaceDescriptor& fd, bool clearContainer)
{
	if(clearContainer)
		vVolumesOut.clear();

//	if no volumes are present, we can leave immediately
	if(grid.num<Volume>() == 0)
		return;

//	iterate through all volumes which are connected to the first vertex of fd.
//	check for each if it contains f. If so, store it in the container.
//	to make things a little faster we'll check if the associated volumes contain a second vertex of f.
	Vertex* v1 = fd.vertex(0);

	Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(v1);
	for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(v1);
		iter != iterEnd; ++iter)
	{
		Volume* v = *iter;
	//	check if v contains v2
		if(VolumeContains(v, fd.vertex(1)))
		{
			if(VolumeContains(v, fd.vertex(2)))
			{
				if(VolumeContains(v, &fd))
					vVolumesOut.push_back(*iter);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	VolumeContains
///	returns true if the given volume contains the given vertex
bool VolumeContains(VolumeVertices* v, Vertex* vrt)
{
	uint numVrts = v->num_vertices();
	Volume::ConstVertexArray vrts = v->vertices();

	for(uint i = 0; i < numVrts; ++i)
	{
		if(vrts[i] == vrt)
			return true;
	}
	return false;
}

///	returns true if the given volume contains the given edge
bool VolumeContains(Volume* v, EdgeVertices* ev)
{
	uint numEdges = v->num_edges();
	for(uint i = 0; i < numEdges; ++i)
	{
		EdgeDescriptor ed;
		v->edge_desc(i, ed);
		if(CompareVertices(ev, &ed))
			return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
//	VolumeContains
///	returns true if the given volume contains the given face
bool VolumeContains(Volume* v, FaceVertices* fv)
{
	FaceDescriptor fd;
	uint numFaces = v->num_faces();
	unsigned long hash = hash_key(fv);
	for(uint i = 0; i < numFaces; ++i)
	{
		v->face_desc(i, fd);
		if(hash == hash_key(&fd))
		{
			if(CompareVertices(fv, &fd))
				return true;
		}
	}
	return false;
}

}//	end of namespace libGrid
