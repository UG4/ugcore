//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#include "vertex_util.h"
#include "edge_util.h"
#include "../trees/kd_tree_static.h"
#include "misc_util.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
int GetVertexIndex(EdgeBase* e, VertexBase* v)
{
	if(e->vertex(0) == v)
		return 0;
	else if(e->vertex(1) == v)
		return 1;
	return -1;
}

////////////////////////////////////////////////////////////////////////
int GetVertexIndex(Face* f, VertexBase* v)
{
	uint numVrts = f->num_vertices();
	for(uint i = 0; i < numVrts; ++i)
	{
		if(f->vertex(i) == v)
			return i;
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
int GetVertexIndex(Volume* vol, VertexBase* v)
{
	uint numVrts = vol->num_vertices();
	for(uint i = 0; i < numVrts; ++i)
	{
		if(vol->vertex(i) == v)
			return i;
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
VertexBase* GetConnectedVertex(EdgeBase* e, VertexBase* v)
{
	if(e->vertex(0) == v)
		return e->vertex(1);
	else if(e->vertex(1) == v)
		return e->vertex(0);
	return NULL;
}

////////////////////////////////////////////////////////////////////////
int GetConnectedVertexIndex(Face* f, const EdgeDescriptor& ed)
{
	uint numVrts = f->num_vertices();
	for(uint i = 0; i < numVrts; ++i)
	{
		if((f->vertex(i) != ed.vertex(0)) &&
			(f->vertex(i) != ed.vertex(1)))
		{
			return i;
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
void CollectNeighbours(std::vector<VertexBase*>& vNeighborsOut, Grid& grid, VertexBase* v)
{
	vNeighborsOut.clear();

//	check which object-types have to be examined.
	bool searchFaces = true;
	bool searchVolumes = true;

	if(grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		searchFaces = false;

	if(grid.option_is_enabled(VOLOPT_AUTOGENERATE_EDGES) ||
		grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
		searchVolumes = false;

//	search edges first.
	{
		EdgeBaseIterator iterEnd = grid.associated_edges_end(v);
		for(EdgeBaseIterator iter = grid.associated_edges_begin(v);
			iter != iterEnd; ++iter)
		{
			VertexBase* nv = GetConnectedVertex(*iter, v);
		//	check if v is contained in vNeighbors. If not push it in.
			if(find(vNeighborsOut.begin(), vNeighborsOut.end(), nv) == vNeighborsOut.end())
				vNeighborsOut.push_back(nv);
		}
	}

//	search faces if required.
	if(searchFaces && grid.num_faces() > 0)
	{
		FaceIterator iterEnd = grid.associated_faces_end(v);
		for(FaceIterator iter = grid.associated_faces_begin(v);
			iter != iterEnd; ++iter)
		{
			Face* f = *iter;
			uint numVrts = f->num_vertices();

		//	check for each vertex if it is already a member of vNeighbors or not.
			for(uint j = 0; j < numVrts; ++j)
			{
				VertexBase* nv = f->vertex(j);
				if(nv != v)
				{
				//	check if v is contained in vNeighbors. If not push it in.
					if(find(vNeighborsOut.begin(), vNeighborsOut.end(), nv) == vNeighborsOut.end())
						vNeighborsOut.push_back(nv);
				}
			}
		}
	}

//	search volumes if required.
	if(searchVolumes && grid.num_volumes() > 0)
	{
		VolumeIterator iterEnd = grid.associated_volumes_end(v);
		for(VolumeIterator iter = grid.associated_volumes_begin(v);
			iter != iterEnd; ++iter)
		{
			Volume* vol = *iter;
			uint numVrts = vol->num_vertices();

		//	check for each vertex if it is already a member of vNeighbors or not.
			for(uint j = 0; j < numVrts; ++j)
			{
				VertexBase* nv = vol->vertex(j);
				if(nv != v)
				{
				//	check if v is contained in vNeighbors. If not push it in.
					if(find(vNeighborsOut.begin(), vNeighborsOut.end(), nv) == vNeighborsOut.end())
						vNeighborsOut.push_back(nv);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
VertexBase* FindVertexByCoordiante(vector3& coord, VertexBaseIterator iterBegin, VertexIterator iterEnd,
									Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	if(iterBegin == iterEnd)
		return NULL;

	VertexBase* bestVrt = *iterBegin;
	number bestDistSq = VecDistanceSq(coord, aaPos[bestVrt]);

	VertexBaseIterator iter = iterBegin;
	iter++;
	while(iter != iterEnd)
	{
		number distSq = VecDistance(coord, aaPos[*iter]);
		if(distSq < bestDistSq)
		{
			bestDistSq = distSq;
			bestVrt = *iter;
		}

		++iter;
	}

	return bestVrt;
}

////////////////////////////////////////////////////////////////////////
void CalculateBoundingBox(vector3& vMinOut, vector3& vMaxOut, VertexBaseIterator vrtsBegin,
						  VertexBaseIterator vrtsEnd, Grid::VertexAttachmentAccessor<AVector3>& aaPos)
{
    if(vrtsBegin != vrtsEnd)
    {
		vMinOut.x = aaPos[*vrtsBegin].x;
		vMinOut.y = aaPos[*vrtsBegin].y;
		vMinOut.z = aaPos[*vrtsBegin].z;

		vMaxOut = vMinOut;

    	for(VertexBaseIterator iter = vrtsBegin; iter != vrtsEnd; ++iter)
    	{
    		vMinOut.x = std::min(vMinOut.x, aaPos[*iter].x);
    		vMinOut.y = std::min(vMinOut.y, aaPos[*iter].y);
    		vMinOut.z = std::min(vMinOut.z, aaPos[*iter].z);

    		vMaxOut.x = std::max(vMaxOut.x, aaPos[*iter].x);
    		vMaxOut.y = std::max(vMaxOut.y, aaPos[*iter].y);
    		vMaxOut.z = std::max(vMaxOut.z, aaPos[*iter].z);
    	}
    }
}

////////////////////////////////////////////////////////////////////////
//	CalculateBarycenter - mstepnie
/// calculates the barycenter of a set of vertices
vector3 CalculateBarycenter(VertexBaseIterator vrtsBegin, VertexBaseIterator vrtsEnd,
							Grid::VertexAttachmentAccessor<AVector3>& aaPos)
{
	vector3 v = vector3(0,0,0);
	int counter = 0;
	for(VertexBaseIterator iter = vrtsBegin; iter != vrtsEnd; ++iter)
	{
		VecAdd(v,v,aaPos[*iter]);
		counter++;
	}
	
	if(counter>0)
		VecScale(v,v,1.f/counter);
	return v;
}

////////////////////////////////////////////////////////////////////////
//	MergeVertices
///	merges two vertices and restructures the adjacent elements.
void MergeVertices(Grid& grid, VertexBase* v1, VertexBase* v2)
{
//	first we have to check if there are elements that connect the vertices.
//	We have to delete those.
	EraseConnectingElements(grid, v1, v2);

//	create new edges for each edge that is connected with v2.
//	avoid double edges
	if(grid.num_edges() > 0)
	{
		EdgeDescriptor ed;
		EdgeBaseIterator iterEnd = grid.associated_edges_end(v2);
		for(EdgeBaseIterator iter = grid.associated_edges_begin(v2); iter != iterEnd; ++iter)
		{
			EdgeBase* e = *iter;
			if(e->vertex(0) == v2)
				ed.set_vertices(v1, e->vertex(1));
			else
				ed.set_vertices(e->vertex(0), v1);

			if(!grid.get_edge(ed))
				grid.create_by_cloning(e, ed, e);
		}
	}

//	create new faces for each face that is connected to v2
//	avoid double faces.
	if(grid.num_faces() > 0)
	{
		FaceDescriptor fd;
		FaceIterator iterEnd = grid.associated_faces_end(v2);
		for(FaceIterator iter = grid.associated_faces_begin(v2); iter != iterEnd; ++iter)
		{
			Face* f = *iter;
			uint numVrts = f->num_vertices();
			fd.set_num_vertices(numVrts);
			for(uint i = 0; i < numVrts; ++i)
			{
				if(f->vertex(i) == v2)
					fd.set_vertex(i, v1);
				else
					fd.set_vertex(i, f->vertex(i));
			}

			if(!grid.get_face(fd))
				grid.create_by_cloning(f, fd, f);
		}
	}

//	create new volumes for each volume that is connected to v2
	if(grid.num_volumes() > 0)
	{
		VolumeDescriptor vd;
		VolumeIterator iterEnd = grid.associated_volumes_end(v2);
		for(VolumeIterator iter = grid.associated_volumes_begin(v2); iter != iterEnd; ++iter)
		{
			Volume* v = *iter;
			uint numVrts = v->num_vertices();
			vd.set_num_vertices(numVrts);
			for(uint i = 0; i < numVrts; ++i)
			{
				if(v->vertex(i) == v2)
					vd.set_vertex(i, v1);
				else
					vd.set_vertex(i, v->vertex(i));
			}

			assert(!"avoid double volumes! implement FindVolume and use it here.");
			grid.create_by_cloning(v, vd, v);
		}
	}

//	new elements have been created. remove the old ones.
//	it is sufficient to simply erase v2.
	grid.erase(v2);
}

////////////////////////////////////////////////////////////////////////
//TODO:	replace KDTreeStatic by a dynamic kd-tree.
void RemoveDoubles(Grid& grid, const VertexBaseIterator& iterBegin, const VertexBaseIterator& iterEnd, AVector3& aPos, number threshold)
{
	if(!grid.has_vertex_attachment(aPos))
		return;

	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);

	KDTreeStatic<AVector3> kdTree;
	kdTree.create_from_grid(grid, iterBegin, iterEnd, aPos, 20, 10, KDSD_LARGEST);

//	we need temporary attachments:
//	a vector<VertexBase*> attachment, that stores for each vertex all other vertices
//	closer than threshold, which have higher attachment data index.
	typedef Attachment<list<VertexBase*> >	AVertexList;
	AVertexList aVertexList;
	grid.attach_to_vertices(aVertexList);
	Grid::VertexAttachmentAccessor<AVertexList> aaVL(grid, aVertexList);

//	we'll store in this attachment whether a vertex will be merged or not.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);
	{
		for(VertexBaseIterator iter = iterBegin; iter != iterEnd; ++iter)
			aaInt[*iter] = 0;
	}

//	compare squares.
	threshold *= threshold;
//	iterate over all vertices and collect all that have aInt == 0 and are within range.
	for(VertexBaseIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		VertexBase* v = *iter;
		if(aaInt[v] == 0)
		{//	the vertex will not be removed during merge
		//	find all vertices closer than threshold
			list<VertexBase*> neighbours;
			uint numClosest = 3;
			while(numClosest < grid.num_vertices())
			{
				neighbours.clear();
				kdTree.get_neighbourhood(neighbours, aaPos[v], numClosest);

				if(VecDistanceSq(aaPos[neighbours.back()], aaPos[v]) < threshold)
					numClosest *= 2;
				else
					break;
			}

		//	store them in the vertexVec attachment
			if(!neighbours.empty())
			{
				for(list<VertexBase*>::iterator nIter = neighbours.begin();
					nIter != neighbours.end(); ++nIter)
				{
					VertexBase* nv = *nIter;
					if(aaInt[nv] == 0)
					{
						if(nv != v)
						{
							if(VecDistanceSq(aaPos[v], aaPos[nv]) < threshold)
							{
								aaVL[v].push_back(nv);
								aaInt[nv] = 1;
							}
							else
								break;
						}
					}
				}
			}
		}
	}

//	iterate over all vertices again and merge collected ones
	{
		VertexBaseIterator iter = iterBegin;
		while(iter != iterEnd)
		{
			VertexBase* v = *iter;
			if(!aaVL[v].empty())
			{
				list<VertexBase*>::iterator nIter = aaVL[v].begin();
				while(nIter != aaVL[v].end())
				{
					VertexBase* delVrt = *nIter;
					nIter++;
					MergeVertices(grid, v, delVrt);
				}
			}

			++iter;
		}
	}

	grid.detach_from_vertices(aVertexList);
	grid.detach_from_vertices(aInt);
}

////////////////////////////////////////////////////////////////////////
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

}//	end of namespace
