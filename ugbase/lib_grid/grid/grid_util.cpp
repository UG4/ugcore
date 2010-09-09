//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d16

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

	for(uint i = 0; i < numVrts; ++i)
	{
		uint j;
		for(j = 0; j < numVrts; ++j)
		{
			if(fv1->vertex(i) == fv2->vertex(j))
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

	for(uint i = 0; i < numVrts; ++i)
	{
		uint j;
		for(j = 0; j < numVrts; ++j)
		{
			if(vv1->vertex(i) == vv2->vertex(j))
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
void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid, VertexBase* v, bool clearContainer)
{
	// clear container if wanted
	if(clearContainer)
		vVertexOut.clear();

	// add vertex pointers to container
		vVertexOut.push_back(v);
}

///	Collects all vertices.
void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid, EdgeBase* e, bool clearContainer)
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
void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid, Face* f, bool clearContainer)
{
	// clear container if wanted
	if(clearContainer)
		vVertexOut.clear();

	// resize container
	const size_t numVertex = f->num_vertices();

	// add vertex pointers to container
	for(size_t i = 0; i < numVertex; ++i)
		vVertexOut.push_back(f->vertex(i));
}

///	Collects all vertices.
void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid, Volume* v, bool clearContainer)
{
	// clear container if wanted
	if(clearContainer)
		vVertexOut.clear();

	// resize container
	const size_t numVertex = v->num_vertices();

	// add vertex pointers to container
	for(size_t i = 0; i < numVertex; ++i)
		vVertexOut.push_back(v->vertex(i));
}



////////////////////////////////////////////////////////////////////////
//	CollectEdgesSorted
///	Collects all edges. (Returns the edge itself)
void CollectEdgesSorted(vector<EdgeBase*>& vEdgesOut, Grid& grid, EdgeBase* e, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

	vEdgesOut.push_back(e);
}

///	Collects all edges which exist in the given grid and which are part of the given face in the order defined by the reference elements.
void CollectEdgesSorted(vector<EdgeBase*>& vEdgesOut, Grid& grid, Face* f, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	second-best: use GetEdge in order to find the queried edges
	uint numEdges = f->num_edges();
	for(uint i = 0; i < numEdges; ++i)
	{
		EdgeBase* e = grid.get_edge(f, i);
		if(e != NULL)
			vEdgesOut.push_back(e);
	}
}

///	Collects all edges that exist in the given grid are part of the given volume in the order defined by the reference elements.
void CollectEdgesSorted(vector<EdgeBase*>& vEdgesOut, Grid& grid, Volume* v, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	second best: use GetEdge in order to find the queried edges.
	uint numEdges = v->num_edges();
	for(uint i = 0; i < numEdges; ++i)
	{
		EdgeBase* e = grid.get_edge(v, i);
		if(e != NULL)
			vEdgesOut.push_back(e);
	}

}

////////////////////////////////////////////////////////////////////////
//	CollectEdges
///	Collects all edges which exist in the given grid and which are part of the given vertex.
void CollectEdges(std::vector<EdgeBase*>& vEdgesOut, Grid& grid, VertexBase* vrt, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

	Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(vrt);
	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(vrt);
		iter != iterEnd; ++iter)
	{
		vEdgesOut.push_back(*iter);
	}
}

///	Collects all edges. (Returns the edge itself)
void CollectEdges(vector<EdgeBase*>& vEdgesOut, Grid& grid, EdgeBase* e, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

	vEdgesOut.push_back(e);
}

///	Collects all edges which exist in the given grid and which are part of the given face.
void CollectEdges(vector<EdgeBase*>& vEdgesOut, Grid& grid, Face* f, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

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
		EdgeBase* e = grid.get_edge(f, i);
		if(e != NULL)
			vEdgesOut.push_back(e);
	}
}

///	Collects all edges that exist in the given grid are part of the given volume.
void CollectEdges(vector<EdgeBase*>& vEdgesOut, Grid& grid, Volume* v, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

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
		EdgeBase* e = grid.get_edge(v, i);
		if(e != NULL)
			vEdgesOut.push_back(e);
	}

}

////////////////////////////////////////////////////////////////////////
//	CollectFaces
///	Collects all faces that exist in the given grid which contain the given vertex.
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, VertexBase* vrt, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

	Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrt);
	for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt);
		iter != iterEnd; ++iter)
	{
		vFacesOut.push_back(*iter);
	}
}

///	Collects all faces that exist in the given grid which contain the given edge.
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, EdgeBase* e, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

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
		VertexBase* v1 = e->vertex(0);
		Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(v1);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v1);
			iter != iterEnd; ++iter)
		{
			if(FaceContains(*iter, e))
				vFacesOut.push_back(*iter);
		}
	}
}

///	Collects all faces that exist in the given grid are part of the given volume.
void CollectFaces(vector<Face*>& vFacesOut, Grid& grid, Volume* v, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

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
		f->edge(i, ed);
		if(CompareVertices(ev, &ed))
			return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the given vertex
bool FaceContains(Face* f, VertexBase* v)
{
	uint numVrts = f->num_vertices();
	for(uint i = 0; i < numVrts; ++i)
	{
		if(f->vertex(i) == v)
			return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, VertexBase* vrt, bool clearContainer)
{
	if(clearContainer)
		vVolumesOut.clear();

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
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, EdgeBase* e, bool clearContainer)
{
	if(clearContainer)
		vVolumesOut.clear();

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
	VertexBase* v1 = e->vertex(0);
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
	VertexBase* v1 = f->vertex(0);

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

///	Collects all volumes that exist in the given grid which contain the given face.
void CollectVolumes(std::vector<Volume*>& vVolumesOut, Grid& grid, FaceDescriptor& fd, bool clearContainer)
{
	if(clearContainer)
		vVolumesOut.clear();

//	iterate through all volumes which are connected to the first vertex of fd.
//	check for each if it contains f. If so, store it in the container.
//	to make things a little faster we'll check if the associated volumes contain a second vertex of f.
	VertexBase* v1 = fd.vertex(0);

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
bool VolumeContains(Volume* v, VertexBase* vrt)
{
	uint numVrts = v->num_vertices();
	for(uint i = 0; i < numVrts; ++i)
	{
		if(v->vertex(i) == vrt)
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
		v->edge(i, ed);
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
	uint hash = hash_key(fv);
	for(uint i = 0; i < numFaces; ++i)
	{
		v->face(i, fd);
		if(hash == hash_key(&fd))
		{
			if(CompareVertices(fv, &fd));
				return true;
		}
	}
	return false;
}

}//	end of namespace libGrid
