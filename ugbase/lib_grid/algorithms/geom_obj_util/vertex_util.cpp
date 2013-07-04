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
int GetVertexIndex(EdgeVertices* e, VertexBase* v)
{
	if(e->vertex(0) == v)
		return 0;
	else if(e->vertex(1) == v)
		return 1;
	return -1;
}

////////////////////////////////////////////////////////////////////////
int GetVertexIndex(FaceVertices* f, VertexBase* v)
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
int GetVertexIndex(VolumeVertices* vol, VertexBase* v)
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
//	GetConnectedVertex
VertexBase* GetConnectedVertex(EdgeVertices* e, Face* f)
{
	uint numVrts = f->num_vertices();
	for(uint i = 0; i < numVrts; ++i){
		if((f->vertex(i) != e->vertex(0)) &&
			(f->vertex(i) != e->vertex(1)))
			return f->vertex(i);
	}
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
EdgeBase* GetConnectedEdge(Grid& g, VertexBase* vrt, Face* tri)
{
	size_t numEdges = tri->num_edges();
	EdgeDescriptor ed;
	for(size_t i = 0; i < numEdges; ++i){
		tri->edge_desc(i, ed);
		if(!EdgeContains(&ed, vrt))
			return g.get_edge(ed);
	}
	return NULL;
}

////////////////////////////////////////////////////////////////////////
int NumAssociatedEdges(Grid& grid, VertexBase* v)
{
	Grid::edge_traits::secure_container edges;
	grid.associated_elements(edges, v);
	return (int)edges.size();
}

////////////////////////////////////////////////////////////////////////
int NumAssociatedFaces(Grid& grid, VertexBase* v)
{
	Grid::face_traits::secure_container faces;
	grid.associated_elements(faces, v);
	return (int)faces.size();
}

////////////////////////////////////////////////////////////////////////
//	CollectSurfaceNeighborsSorted
bool CollectSurfaceNeighborsSorted(std::vector<VertexBase*>& vNeighborsOut,
								   Grid& grid, VertexBase* v)
{
	vNeighborsOut.clear();
	
//	the algorithm won't work if volumes are connected
	if(grid.num_volumes() > 0){
		if(grid.associated_volumes_begin(v) != grid.associated_volumes_end(v))
			return false;
	}
		
//	the algorithm won't work if no faces are connected
	if(grid.associated_faces_begin(v) == grid.associated_faces_end(v))
		return false;
	
	grid.begin_marking();
	
	if(grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)
	   && grid.option_is_enabled(EDGEOPT_STORE_ASSOCIATED_FACES)){
	//	collect edges in this vector
		vector<EdgeBase*> edges;
	//	start with an arbitrary edge
		EdgeBase* curEdge = *grid.associated_edges_begin(v);
		
		while(curEdge){
			vNeighborsOut.push_back(GetConnectedVertex(curEdge, v));
			grid.mark(curEdge);
			
		//	get associated faces
			Face* f[2];
			if(GetAssociatedFaces(f, grid, curEdge, 2) != 2)
				return false;
			
			curEdge = NULL;
			for(int i = 0; i < 2; ++i){
				if(!grid.is_marked(f[i])){
					CollectEdges(edges, grid, f[i]);
					for(size_t j = 0; j < edges.size(); ++j){
						if(!grid.is_marked(edges[j])){
							if(EdgeContains(edges[j], v)){
								curEdge = edges[j];
								break;
							}
						}
					}
				}
			}
		}
	}
	else{
	//	we can't use GetAssociatedFaces here, since it uses Grid::mark itself if
	//	EDGEOPT_STORE_ASSOCIATED_FACES is not enabled.
	//	Start with an arbitrary face
		Face* f = *grid.associated_faces_begin(v);
		grid.mark(v);

		while(f){
			grid.mark(f);
			
		//	mark one of the edges that is connected to v by marking
		//	the edges endpoints. Make sure that it was not already marked.
			size_t numVrts = f->num_vertices();
			int vind = GetVertexIndex(f, v);
			VertexBase* tvrt = f->vertex((vind + 1)%numVrts);
			if(grid.is_marked(tvrt))
				tvrt = f->vertex((vind + numVrts - 1)%numVrts);
			if(grid.is_marked(tvrt))
				throw(UGError("CollectSurfaceNeighborsSorted: unexpected exit."));
				
			vNeighborsOut.push_back(tvrt);
			grid.mark(tvrt);

		//	iterate through the faces associated with v and find an unmarked one that
		//	contains two marked vertices
			f = NULL;
			Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(v);
			for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v);
				iter != iterEnd; ++iter)
			{
				if(!grid.is_marked(*iter)){
					f = *iter;
					size_t numMarked = 0;
					for(size_t i = 0; i < f->num_vertices(); ++i){
						if(grid.is_marked(f->vertex(i)))
							++numMarked;
					}
					if(numMarked == 2)
						break;
					else
						f = NULL;
				}
			}
		}
	}

	grid.end_marking();
	return true;
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
//	CalculateVertexNormals
bool CalculateVertexNormals(Grid& grid,
							Grid::AttachmentAccessor<VertexBase, APosition>& aaPos,
							Grid::AttachmentAccessor<VertexBase, ANormal>& aaNorm)
{
//	set all normals to zero
	{
		for(VertexBaseIterator iter = grid.begin<VertexBase>();
			iter != grid.end<VertexBase>(); iter++)
			aaNorm[*iter] = vector3(0, 0, 0);
	}
//	loop through all the faces, calculate their normal and add them to their connected points
	{
		for(FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>(); iter++)
		{
			Face* f = *iter;
			vector3 vN;

			CalculateNormal(vN, f, aaPos);

			for(size_t i = 0; i < f->num_vertices(); ++i)
				VecAdd(aaNorm[f->vertex(i)], aaNorm[f->vertex(i)], vN);
		}
	}
//	loop through all the points and normalize their normals
	{
		for(VertexBaseIterator iter = grid.begin<VertexBase>();
			iter != grid.end<VertexBase>(); iter++)
			VecNormalize(aaNorm[*iter], aaNorm[*iter]);
	}
//	done
	return true;
}

bool CalculateVertexNormals(Grid& grid, APosition& aPos, ANormal& aNorm)
{
	if(!grid.has_attachment<VertexBase>(aPos))
		return false;
	if(!grid.has_attachment<VertexBase>(aNorm))
		grid.attach_to<VertexBase>(aNorm);

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<ANormal> aaNorm(grid, aNorm);

	return CalculateVertexNormals(grid, aaPos, aaNorm);
}


////////////////////////////////////////////////////////////////////////
//	MergeVertices
///	merges two vertices and restructures the adjacent elements.
void MergeVertices(Grid& grid, VertexBase* v1, VertexBase* v2)
{
//	make sure that GRIDOPT_VERTEXCENTRIC_INTERCONNECTION is enabled
	if(grid.num_edges() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))){
		LOG("  WARNING in MergeVertices: autoenabling VRTOPT_STORE_ASSOCIATED_EDGES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
	}
	if(grid.num_faces() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))){
		LOG("  WARNING in MergeVertices: autoenabling VRTOPT_STORE_ASSOCIATED_FACES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_FACES);
	}
	if(grid.num_volumes() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))){
		LOG("  WARNING in MergeVertices: autoenabling VRTOPT_STORE_ASSOCIATED_VOLUMES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_VOLUMES);
	}

//	notify the grid, that the two vertices will be merged
	grid.objects_will_be_merged(v1, v1, v2);
	
//	first we have to check if there are elements that connect the vertices.
//	We have to delete those.
	EraseConnectingElements(grid, v1, v2);

//	create new edges for each edge that is connected with v2.
//	avoid double edges
	if(grid.num_edges() > 0)
	{
		EdgeDescriptor ed;
		Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(v2);
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v2); iter != iterEnd; ++iter)
		{
			EdgeBase* e = *iter;
			if(e->vertex(0) == v2)
				ed.set_vertices(v1, e->vertex(1));
			else
				ed.set_vertices(e->vertex(0), v1);

			EdgeBase* existingEdge = grid.get_edge(ed);
			if(!existingEdge)
				grid.create_by_cloning(e, ed, e);
			else
				grid.objects_will_be_merged(existingEdge, existingEdge, e);
		}
	}

//	create new faces for each face that is connected to v2
//	avoid double faces.
	if(grid.num_faces() > 0)
	{
		FaceDescriptor fd;
		Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(v2);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v2); iter != iterEnd; ++iter)
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

			Face* existingFace = grid.get_face(fd);
			if(!existingFace)
				grid.create_by_cloning(f, fd, f);
			else
				grid.objects_will_be_merged(existingFace, existingFace, f);
		}
	}

//	create new volumes for each volume that is connected to v2
	if(grid.num_volumes() > 0)
	{
		VolumeDescriptor vd;
		Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(v2);
		for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(v2); iter != iterEnd; ++iter)
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

			//assert(!"avoid double volumes! implement FindVolume and use it here.");
			grid.create_by_cloning(v, vd, v);
		}
	}

//	new elements have been created. remove the old ones.
//	it is sufficient to simply erase v2.
	grid.erase(v2);
}

////////////////////////////////////////////////////////////////////////
bool IsBoundaryVertex1D(Grid& grid, VertexBase* v,
						Grid::edge_traits::callback cbConsiderEdge)
{
	if(!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we have to enable this option, since nothing works without it in reasonable time.
		LOG("WARNING in IsBoundaryVertex1D(...): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES.\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
	}

//	iterate over associated edges and return true if only one of them
//	should be considered for the polygonal chain
	size_t counter = 0;
	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
		iter != grid.associated_edges_end(v); ++iter)
	{
		if(cbConsiderEdge(*iter)){
			++counter;
			if(counter > 1)
				return false;
		}
	}
	
	return true;
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

	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
		iter != grid.associated_edges_end(v); ++iter)
	{
		if(IsBoundaryEdge2D(grid, *iter))
			return true;
	}

	return false;
}

bool IsBoundaryVertex3D(Grid& grid, VertexBase* v)
{
//	check whether one of the associated edges is a boundary edge.
//	if so return true.
	if(!grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
	{
	//	we have to enable this option, since we need edges in order to detect boundary vertices.
		LOG("  WARNING in IsBoundaryVertex2D(...): auto-enabling VOLOPT_AUTOGENERATE_FACES.\n");
		grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}
	if(!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))
	{
	//	we have to enable this option, since nothing works without it in reasonable time.
		LOG("  WARNING in IsBoundaryVertex2D(...): auto-enabling VRTOPT_STORE_ASSOCIATED_FACES.\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_FACES);
	}

	for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v);
		iter != grid.associated_faces_end(v); ++iter)
	{
		if(IsVolumeBoundaryFace(grid, *iter))
			return true;
	}

	return false;
}

bool LiesOnBoundary(Grid& grid, VertexBase* v)
{
	if(IsBoundaryVertex1D(grid, v))
		return true;

	if((grid.num<Face>() > 0) && IsBoundaryVertex2D(grid, v))
		return true;

	if((grid.num<Volume>() > 0) && IsBoundaryVertex3D(grid, v))
		return true;

	return false;
}


////////////////////////////////////////////////////////////////////////
bool IsRegularSurfaceVertex(Grid& grid, VertexBase* v)
{
//	check how many faces each associated edge has
	Grid::AssociatedEdgeIterator edgesEnd = grid.associated_edges_end(v);
	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
		iter != edgesEnd; ++iter)
	{
		if(NumAssociatedFaces(grid, *iter) != 2)
			return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////
void MarkFixedCreaseVertices(Grid& grid, SubsetHandler& sh,
							int creaseSI, int fixedSI)
{
//	if there are no crease-edges then there is nothing to do.
	if((int)sh.num_subsets() <= creaseSI)
		return;
	if(sh.num<EdgeBase>(creaseSI) == 0)
		return;

//	begin marking
	grid.begin_marking();
//	iterate over all crease-edges
	for(EdgeBaseIterator iter = sh.begin<EdgeBase>(creaseSI);
		iter != sh.end<EdgeBase>(creaseSI); ++iter)
	{
	//	check for both vertices whether they are fixed-vertices
		for(int i = 0; i < 2; ++i)
		{
			VertexBase* v = (*iter)->vertex(i);
		//	if the vertex is not marked (has not been checked yet)
			if(!grid.is_marked(v))
			{
			//	mark it
				grid.mark(v);
			//	count associated crease edges
				int counter = 0;
				Grid::AssociatedEdgeIterator aeIterEnd = grid.associated_edges_end(v);
				for(Grid::AssociatedEdgeIterator aeIter = grid.associated_edges_begin(v);
					aeIter != aeIterEnd; ++aeIter)
				{
					if(sh.get_subset_index(*aeIter) == creaseSI)
					{
					//	the edge is a crease-edge. Increase the counter.
						++counter;
					//	if the counter is higher than 2, the vertex is a fixed vertex.
						if(counter > 2)
						{
							sh.assign_subset(v, fixedSI);
							break;
						}
					}
				}
			}
		}
	}

//	end marking
	grid.end_marking();
}

}//	end of namespace
