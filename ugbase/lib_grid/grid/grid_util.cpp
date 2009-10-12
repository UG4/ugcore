//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d16

#include "grid_util.h"
using namespace std;

namespace ug
{

bool FaceMatches(Face* f, FaceDescriptor& fd, uint faceHash, uint faceDescriptorHash);
bool FaceDescMatches(FaceDescriptor& fd1, FaceDescriptor& fd2,
					uint fd1Hash, uint fd2Hash);
bool VolumeMatches(Volume* v, VolumeDescriptor& vd, uint volHash, uint volDescriptorHash);

inline uint FaceHash(FaceDescriptor& fd)
{
	uint numVrts = fd.num_vertices();
	uint hash = 0;
	for(uint i = 0; i < numVrts; ++i)
		hash += attachment_traits<GeometricObject*, Grid>::get_data_index(fd.vertex(i));
	return hash;
}

inline uint FaceHash(Face* f)
{
	uint numVrts = f->num_vertices();
	uint hash = 0;
	for(uint i = 0; i < numVrts; ++i)
		hash += attachment_traits<GeometricObject*, Grid>::get_data_index(f->vertex(i));
	return hash;
}

inline uint VolumeHash(VolumeDescriptor& vd)
{
	uint numVrts = vd.num_vertices();
	uint hash = 0;
	for(uint i = 0; i < numVrts; ++i)
		hash += attachment_traits<GeometricObject*, Grid>::get_data_index(vd.vertex(i));
	return hash;
}

inline uint VolumeHash(Volume* v)
{
	uint numVrts = v->num_vertices();
	uint hash = 0;
	for(uint i = 0; i < numVrts; ++i)
		hash += attachment_traits<GeometricObject*, Grid>::get_data_index(v->vertex(i));
	return hash;
}

////////////////////////////////////////////////////////////////////////
//	IsVolumeBoundaryFace
///	returns true if the given face is a boundary face.
bool IsVolumeBoundaryFace(Grid& grid, Face* f)
{
//	check if FACEOPT_STORE_ASSOCIATED_VOLUMES is enabled.
//	if so, use it to count the number of adjacent volumes.
	int counter = 0;
	if(grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		for(VolumeIterator iter = grid.associated_volumes_begin(f);
			iter != grid.associated_volumes_end(f); iter++)
		{
			counter++;
		}
	}
	else
	{
	//	iterate over all volumes which are connected to the first vertex
	//	and check if they contain the face...
		VolumeIterator iterEnd = grid.associated_volumes_end(f->vertex(0));
		for(VolumeIterator iter = grid.associated_volumes_begin(f->vertex(0));
			iter != iterEnd; iter++)
		{
			uint numFaces = (*iter)->num_faces();
			for(uint i = 0; i < numFaces; ++i)
			{
				if(VolumeContains(*iter, f))
					counter++;
			}
		}
	}

//	if there are less than two adjacent volumes, the triangle is regarded as a boundary triangle
	if(counter < 2)
		return true;
	return false;
}

bool IsBoundaryEdge2D(Grid& grid, EdgeBase* e)
{
//	get the number of connected faces. if only one face is connected then
//	the edge is considered to be a boundary edge.
	int counter = 0;
	if(grid.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		for(FaceIterator iter = grid.associated_faces_begin(e);
			iter != grid.associated_faces_end(e); ++iter)
			++counter;

		if(counter == 1)
			return true;
	}
	else
	{
	//	fill a vector using a helper function
		vector<Face*> vFaces;
		CollectFaces(vFaces, grid, e, false);
		if(vFaces.size() == 1)
			return true;
	}
	return false;
}

bool IsBoundaryVertex2D(Grid& grid, VertexBase* v)
{
//	check whether one of the associated edges is a boundary edge.
//	if so return true.
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
	{
	//	we have to enable this option, since we need edges in order to detect boundary vertices.
		LOG("WARNING in IsBoundaryVertex2D(...): auto-enabling FACEOPT_AUTOGENERATE_EDGES.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}
	if(!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we have to enable this option, since nothing works without it in reasonable time.
		LOG("WARNING in IsBoundaryVertex2D(...): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES.\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
	}

	for(EdgeBaseIterator iter = grid.associated_edges_begin(v);
		iter != grid.associated_edges_end(v); ++iter)
	{
		if(IsBoundaryEdge2D(grid, *iter))
			return true;
	}

	return false;
}


////////////////////////////////////////////////////////////////////////
//	GetEdge
///	returns the edge of the given grid, that connects the given points.
EdgeBase* FindEdge(Grid& grid, VertexBase* vrt1, VertexBase* vrt2)
{
	EdgeBaseIterator iterEnd = grid.associated_edges_end(vrt1);
	for(EdgeBaseIterator iter = grid.associated_edges_begin(vrt1);
		iter != iterEnd; iter++)
	{
		if(EdgeContains(*iter, vrt1, vrt2))
			return *iter;
	}

	return NULL;
}

////////////////////////////////////////////////////////////////////////
//	GetEdge
///	returns the i-th edge of the given face
EdgeBase* FindEdge(Grid& grid, Face* f, uint ind)
{
	EdgeDescriptor ed;
	f->edge(ind, ed);
	return FindEdge(grid, ed.vertex(0), ed.vertex(1));
}

////////////////////////////////////////////////////////////////////////
//	GetEdge
///	returns the i-th edge of the given volume
EdgeBase* FindEdge(Grid& grid, Volume* v, uint ind)
{
	EdgeDescriptor ed;
	v->edge(ind, ed);
	return FindEdge(grid, ed.vertex(0), ed.vertex(1));
}

////////////////////////////////////////////////////////////////////////
//	CollectEdges
///	Collects all edges which exist in the given grid and which are part of the given face.
void CollectEdges(vector<EdgeBase*>& vEdgesOut, Grid& grid, Face* f, bool clearContainer)
{
	if(clearContainer)
		vEdgesOut.clear();

//	best-option: FACEOPT_STORE_ASSOCIATED_EDGES
	if(grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES))
	{
	//	ok. simply push them into the container.
		EdgeBaseIterator iterEnd = grid.associated_edges_end(f);
		for(EdgeBaseIterator iter = grid.associated_edges_begin(f);
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
		EdgeBase* e = FindEdge(grid, f, i);
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
		EdgeBaseIterator iterEnd = grid.associated_edges_end(v);
		for(EdgeBaseIterator iter = grid.associated_edges_begin(v);
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
		EdgeBase* e = FindEdge(grid, v, i);
		if(e != NULL)
			vEdgesOut.push_back(e);
	}

}

////////////////////////////////////////////////////////////////////////
//	EdgeContains
///	returns true if the given edge contains the given vertices
bool EdgeContains(EdgeBase* e, VertexBase* vrt1, VertexBase* vrt2)
{
	return ((e->vertex(0) == vrt1 && e->vertex(1) == vrt2)
			|| (e->vertex(1) == vrt1 && e->vertex(0) == vrt2));
}

////////////////////////////////////////////////////////////////////////
//	EdgeMatches
///	returns true if the given edge matches the given EdgeDescriptor
bool EdgeMatches(EdgeBase* e, EdgeDescriptor& ed)
{
	return ((e->vertex(0) == ed.vertex(0) && e->vertex(1) == ed.vertex(1))
			|| (e->vertex(1) == ed.vertex(0) && e->vertex(0) == ed.vertex(1)));
}

////////////////////////////////////////////////////////////////////////
//	FindFace
///	returns the face that matches the face descriptor.
Face* FindFace(Grid& grid, FaceDescriptor& fd)
{
//build hash of fd
	uint fdHash = FaceHash(fd);

//	iterate through the faces and search for a match
	VertexBase* vrt = fd.vertex(0);
	FaceIterator iterEnd = grid.associated_faces_end(vrt);
	for(FaceIterator iter = grid.associated_faces_begin(vrt);
		iter != iterEnd; iter++)
	{
		Face* f = *iter;
		if(FaceMatches(f, fd, FaceHash(f), fdHash))
			return f;
	}
	return NULL;
}

////////////////////////////////////////////////////////////////////////
//	GetFace
///	returns the i-th face from the given volume
Face* FindFace(Grid& grid, Volume* v, uint ind, bool ignoreAssociatedFaces)
{
	FaceDescriptor fd;
	v->face(ind, fd);

	if((!ignoreAssociatedFaces)
		&& grid.option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES))
	{
	//	calculate the hash of the face-descriptor
		uint hashFD = FaceHash(fd);

	//	search for a match in the associated faces.
		FaceIterator iterEnd = grid.associated_faces_end(v);
		for(FaceIterator iter = grid.associated_faces_begin(v);
			iter != iterEnd; iter++)
		{
			if(FaceMatches(*iter, fd, FaceHash(*iter), hashFD))
				return *iter;
		}

	//	no match
		return NULL;
	}

	return FindFace(grid, fd);
}

////////////////////////////////////////////////////////////////////////
//	CollectFaces
///	Collects all faces that exist in the given grid which contain the given edge.
void CollectFaces(std::vector<Face*>& vFacesOut, Grid& grid, EdgeBase* e, bool clearContainer)
{
	if(clearContainer)
		vFacesOut.clear();

//	best option: EDGEOPT_STORE_ASSOCIATED_FACES
	if(grid.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES))
	{
		FaceIterator iterEnd = grid.associated_faces_end(e);
		for(FaceIterator iter = grid.associated_faces_begin(e);
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
		FaceIterator iterEnd = grid.associated_faces_end(v1);
		for(FaceIterator iter = grid.associated_faces_begin(v1);
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
		FaceIterator iterEnd = grid.associated_faces_end(v);
		for(FaceIterator iter = grid.associated_faces_begin(v);
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
		Face* f = FindFace(grid, v, i);
		if(f != NULL)
			vFacesOut.push_back(f);
	}
}

////////////////////////////////////////////////////////////////////////
//	FaceMatches
///	returns true if the given face contains exactly the same points as the given descriptor.
bool FaceMatches(Face* f, FaceDescriptor& fd)
{
	return FaceMatches(f, fd, FaceHash(f), FaceHash(fd));
}

bool FaceMatches(Face* f, FaceDescriptor& fd, uint faceHash, uint faceDescriptorHash)
{
	if(faceDescriptorHash != faceHash)
		return false;

	uint numVrts = f->num_vertices();

	if(numVrts != fd.num_vertices())
		return false;

	for(uint i = 0; i < numVrts; ++i)
	{
		bool gotOne = false;
		for(uint j = 0; j < numVrts; ++j)
		{
			if(f->vertex(i) == fd.vertex(j))
			{
				gotOne = true;
				break;
			}
		}
		if(!gotOne)
			return false;
	}

	return true;
}

bool FaceDescMatches(FaceDescriptor& fd1, FaceDescriptor& fd2,
					uint fd1Hash, uint fd2Hash)
{
	if(fd1Hash != fd2Hash)
		return false;

	uint numVrts = fd1.num_vertices();

	if(numVrts != fd2.num_vertices())
		return false;

	for(uint i = 0; i < numVrts; ++i)
	{
		bool gotOne = false;
		for(uint j = 0; j < numVrts; ++j)
		{
			if(fd1.vertex(i) == fd2.vertex(j))
			{
				gotOne = true;
				break;
			}
		}
		if(!gotOne)
			return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
//	FaceContains
///	returns true if the given face contains the two given vertices
bool FaceContains(Face* f, EdgeBase* e)
{
	uint numEdges = f->num_edges();
	EdgeDescriptor ed;
	for(uint i = 0; i < numEdges; ++i)
	{
		f->edge(i, ed);
		if(EdgeMatches(e, ed))
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
//	FindVolume
///	returns the volume that matches the volume descriptor.
Volume* FindVolume(Grid& grid, VolumeDescriptor& vd)
{
//build hash of fd
	uint vdHash = VolumeHash(vd);

//	iterate through the volumes and search for a match
	VertexBase* vrt = vd.vertex(0);
	VolumeIterator iterEnd = grid.associated_volumes_end(vrt);
	for(VolumeIterator iter = grid.associated_volumes_begin(vrt);
		iter != iterEnd; iter++)
	{
		Volume* v = *iter;
		if(VolumeMatches(v, vd, VolumeHash(v), vdHash))
			return v;
	}
	return NULL;
}

bool VolumeMatches(Volume* v, VolumeDescriptor& vd)
{
	return VolumeMatches(v, vd, VolumeHash(v), VolumeHash(vd));
}

bool VolumeMatches(Volume* v, VolumeDescriptor& vd, uint volHash, uint volDescriptorHash)
{
	if(volDescriptorHash != volHash)
		return false;

	uint numVrts = v->num_vertices();

	if(numVrts != vd.num_vertices())
		return false;

	for(uint i = 0; i < numVrts; ++i)
	{
		bool gotOne = false;
		for(uint j = 0; j < numVrts; ++j)
		{
			if(v->vertex(i) == vd.vertex(j))
			{
				gotOne = true;
				break;
			}
		}
		if(!gotOne)
			return false;
	}

	return true;
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
		VolumeIterator iterEnd = grid.associated_volumes_end(e);
		for(VolumeIterator iter = grid.associated_volumes_begin(e);
			iter != iterEnd; ++iter)
		{
			vVolumesOut.push_back(*iter);
		}
		return;
	}

//	second best: iterate through all volumes that are associated with the first vertex of e.
//	check for each if it contains e. if so then store it in the container.
	VertexBase* v1 = e->vertex(0);
	VolumeIterator iterEnd = grid.associated_volumes_end(v1);
	for(VolumeIterator iter = grid.associated_volumes_begin(v1);
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
			VolumeIterator iterEnd = grid.associated_volumes_end(f);
			for(VolumeIterator iter = grid.associated_volumes_begin(f);
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

	VolumeIterator iterEnd = grid.associated_volumes_end(v1);
	for(VolumeIterator iter = grid.associated_volumes_begin(v1);
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

	VolumeIterator iterEnd = grid.associated_volumes_end(v1);
	for(VolumeIterator iter = grid.associated_volumes_begin(v1);
		iter != iterEnd; ++iter)
	{
		Volume* v = *iter;
	//	check if v contains v2
		if(VolumeContains(v, fd.vertex(1)))
		{
			if(VolumeContains(v, fd.vertex(2)))
			{
				if(VolumeContains(v, fd))
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
bool VolumeContains(Volume* v, EdgeBase* e)
{
	uint numEdges = v->num_edges();
	for(uint i = 0; i < numEdges; ++i)
	{
		EdgeDescriptor ed;
		v->edge(i, ed);
		if(EdgeMatches(e, ed))
			return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
//	VolumeContains
///	returns true if the given volume contains the given face
bool VolumeContains(Volume* v, Face* f)
{
	FaceDescriptor fd;
	uint numFaces = v->num_faces();
	uint hash = FaceHash(f);
	for(uint i = 0; i < numFaces; ++i)
	{
		v->face(i, fd);
		if(FaceMatches(f, fd, hash, FaceHash(fd)))
			return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
//	VolumeContains
///	returns true if the given volume contains the given face
bool VolumeContains(Volume* v, FaceDescriptor& fd)
{
	FaceDescriptor vfd;
	uint numFaces = v->num_faces();
	uint hash = FaceHash(fd);
	for(uint i = 0; i < numFaces; ++i)
	{
		v->face(i, vfd);
		if(FaceDescMatches(fd, vfd, hash, FaceHash(vfd)))
			return true;
	}
	return false;
}

}//	end of namespace libGrid
