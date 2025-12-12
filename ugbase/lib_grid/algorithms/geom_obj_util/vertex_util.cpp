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

#include "vertex_util.h"
#include "edge_util.h"
//ø #include "../trees/kd_tree_static.h"
#include "misc_util.h"

using namespace std;

namespace ug {

////////////////////////////////////////////////////////////////////////
int GetVertexIndex(EdgeVertices* e, Vertex* v)
{
	if(e->vertex(0) == v)
		return 0;
	else if(e->vertex(1) == v)
		return 1;
	return -1;
}

////////////////////////////////////////////////////////////////////////
int GetVertexIndex(FaceVertices* f, Vertex* v)
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
int GetVertexIndex(VolumeVertices* vol, Vertex* v)
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
Vertex* GetConnectedVertex(Edge* e, Vertex* v)
{
	if(e->vertex(0) == v)
		return e->vertex(1);
	else if(e->vertex(1) == v)
		return e->vertex(0);
	return nullptr;
}

////////////////////////////////////////////////////////////////////////
//	GetConnectedVertex
Vertex* GetConnectedVertex(EdgeVertices* e, Face* f)
{
	uint numVrts = f->num_vertices();
	for(uint i = 0; i < numVrts; ++i){
		if((f->vertex(i) != e->vertex(0)) &&
			(f->vertex(i) != e->vertex(1)))
			return f->vertex(i);
	}
	return nullptr;
}

////////////////////////////////////////////////////////////////////////
int GetConnectedVertexIndex(Face* f, const EdgeDescriptor& ed)
{
	return GetConnectedVertexIndex(f, &ed);
}

int GetConnectedVertexIndex(Face* f, const EdgeVertices* e)
{
	uint numVrts = f->num_vertices();
	for(uint i = 0; i < numVrts; ++i)
	{
		if((f->vertex(i) != e->vertex(0)) &&
			(f->vertex(i) != e->vertex(1)))
		{
			return static_cast<int>(i);
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
Edge* GetConnectedEdge(Grid& g, Vertex* vrt, Face* tri)
{
	size_t numEdges = tri->num_edges();
	EdgeDescriptor ed;
	for(size_t i = 0; i < numEdges; ++i){
		tri->edge_desc(i, ed);
		if(!EdgeContains(&ed, vrt))
			return g.get_edge(ed);
	}
	return nullptr;
}

////////////////////////////////////////////////////////////////////////
Vertex* GetSharedVertex(IVertexGroup* vrts0, IVertexGroup* vrts1)
{
	const size_t num0 = vrts0->size();
	const size_t num1 = vrts1->size();

	IVertexGroup::ConstVertexArray v0 = vrts0->vertices();
	IVertexGroup::ConstVertexArray v1 = vrts1->vertices();

	for(size_t i0 = 0; i0 < num0; ++i0){
		Vertex* v = v0[i0];
		for(size_t i1 = 0; i1 < num1; ++i1){
			if(v == v1[i1])
				return v;
		}
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////
size_t GetSharedVertices(
			std::vector<Vertex*>& vrtsOut,
			IVertexGroup* vrts0,
			IVertexGroup* vrts1)
{
	vrtsOut.clear();

	const size_t num0 = vrts0->size();
	const size_t num1 = vrts1->size();

	IVertexGroup::ConstVertexArray v0 = vrts0->vertices();
	IVertexGroup::ConstVertexArray v1 = vrts1->vertices();

	for(size_t i0 = 0; i0 < num0; ++i0){
		Vertex* v = v0[i0];
		for(size_t i1 = 0; i1 < num1; ++i1){
			if(v == v1[i1])
				vrtsOut.push_back(v);
		}
	}

	return vrtsOut.size();
}

////////////////////////////////////////////////////////////////////////
size_t NumSharedVertices(IVertexGroup* vrts0, IVertexGroup* vrts1)
{
	const size_t num0 = vrts0->size();
	const size_t num1 = vrts1->size();

	IVertexGroup::ConstVertexArray v0 = vrts0->vertices();
	IVertexGroup::ConstVertexArray v1 = vrts1->vertices();

	size_t numShared = 0;
	for(size_t i0 = 0; i0 < num0; ++i0){
		Vertex* v = v0[i0];
		for(size_t i1 = 0; i1 < num1; ++i1){
			if(v == v1[i1])
				++numShared;
		}
	}

	return numShared;
}

////////////////////////////////////////////////////////////////////////
int NumAssociatedEdges(Grid& grid, Vertex* v)
{
	Grid::edge_traits::secure_container edges;
	grid.associated_elements(edges, v);
	return static_cast<int>(edges.size());
}

////////////////////////////////////////////////////////////////////////
int NumAssociatedFaces(Grid& grid, Vertex* v)
{
	Grid::face_traits::secure_container faces;
	grid.associated_elements(faces, v);
	return static_cast<int>(faces.size());
}

////////////////////////////////////////////////////////////////////////
//	CollectSurfaceNeighborsSorted
bool CollectSurfaceNeighborsSorted(std::vector<Vertex*>& vNeighborsOut,
								   Grid& grid, Vertex* v)
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
	
	if(grid.option_is_enabled(FaceOptions::FACEOPT_AUTOGENERATE_EDGES)
	   && grid.option_is_enabled(EdgeOptions::EDGEOPT_STORE_ASSOCIATED_FACES)){
	//	collect edges in this vector
		vector<Edge*> edges;
	//	start with an arbitrary edge
		Edge* curEdge = *grid.associated_edges_begin(v);
		
		while(curEdge){
			vNeighborsOut.push_back(GetConnectedVertex(curEdge, v));
			grid.mark(curEdge);
			
		//	get associated faces
			Face* f[2];
			if(GetAssociatedFaces(f, grid, curEdge, 2) != 2)
				return false;
			
			curEdge = nullptr;
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
			Vertex* tvrt = f->vertex((vind + 1)%numVrts);
			if(grid.is_marked(tvrt))
				tvrt = f->vertex((vind + numVrts - 1)%numVrts);
			if(grid.is_marked(tvrt))
				throw(UGError("CollectSurfaceNeighborsSorted: unexpected exit."));
				
			vNeighborsOut.push_back(tvrt);
			grid.mark(tvrt);

		//	iterate through the faces associated with v and find an unmarked one that
		//	contains two marked vertices
			f = nullptr;
			auto iterEnd = grid.associated_faces_end(v);
			for(auto iter = grid.associated_faces_begin(v);
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
						f = nullptr;
				}
			}
		}
	}

	grid.end_marking();
	return true;
}

////////////////////////////////////////////////////////////////////////
Vertex* FindVertexByCoordiante(vector3& coord, VertexIterator iterBegin, VertexIterator iterEnd,
									Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	if(iterBegin == iterEnd)
		return nullptr;

	Vertex* bestVrt = *iterBegin;
	number bestDistSq = VecDistanceSq(coord, aaPos[bestVrt]);

	VertexIterator iter = iterBegin;
	++iter;
	while(iter != iterEnd)
	{
		number distSq = VecDistanceSq(coord, aaPos[*iter]);
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
							Grid::AttachmentAccessor<Vertex, APosition>& aaPos,
							Grid::AttachmentAccessor<Vertex, ANormal>& aaNorm)
{
//	set all normals to zero
	{
		for(VertexIterator iter = grid.begin<Vertex>();
			iter != grid.end<Vertex>(); ++iter)
			aaNorm[*iter] = vector3(0, 0, 0);
	}
//	loop through all the faces, calculate their normal and add them to their connected points
	{
		for(FaceIterator iter = grid.begin<Face>(); iter != grid.end<Face>(); ++iter)
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
		for(VertexIterator iter = grid.begin<Vertex>();
			iter != grid.end<Vertex>(); ++iter)
			VecNormalize(aaNorm[*iter], aaNorm[*iter]);
	}
//	done
	return true;
}

bool CalculateVertexNormals(Grid& grid, APosition& aPos, ANormal& aNorm)
{
	if(!grid.has_attachment<Vertex>(aPos))
		return false;
	if(!grid.has_attachment<Vertex>(aNorm))
		grid.attach_to<Vertex>(aNorm);

	Grid::VertexAttachmentAccessor aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor aaNorm(grid, aNorm);

	return CalculateVertexNormals(grid, aaPos, aaNorm);
}


////////////////////////////////////////////////////////////////////////
//	MergeVertices
///	merges two vertices and restructures the adjacent elements.
void MergeVertices(Grid& grid, Vertex* v1, Vertex* v2)
{
//	make sure that GRIDOPT_VERTEXCENTRIC_INTERCONNECTION is enabled
	if(grid.num_edges() && (!grid.option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_EDGES))){
		UG_LOG("  WARNING in MergeVertices: autoenabling VRTOPT_STORE_ASSOCIATED_EDGES\n");
		grid.enable_options(VertexOptions::VRTOPT_STORE_ASSOCIATED_EDGES);
	}
	if(grid.num_faces() && (!grid.option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_FACES))){
		UG_LOG("  WARNING in MergeVertices: autoenabling VRTOPT_STORE_ASSOCIATED_FACES\n");
		grid.enable_options(VertexOptions::VRTOPT_STORE_ASSOCIATED_FACES);
	}
	if(grid.num_volumes() && (!grid.option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_VOLUMES))){
		UG_LOG("  WARNING in MergeVertices: autoenabling VRTOPT_STORE_ASSOCIATED_VOLUMES\n");
		grid.enable_options(VertexOptions::VRTOPT_STORE_ASSOCIATED_VOLUMES);
	}


	Edge* conEdge = grid.get_edge(v1, v2);
	if(conEdge){
	//	perform an edge-collapse on conEdge
		CollapseEdge(grid, conEdge, v1);
	}
	else{
	//	notify the grid, that the two vertices will be merged
		grid.objects_will_be_merged(v1, v1, v2);
		
	//	we have to check if there are elements that connect the vertices.
	//	We have to delete those.
		EraseConnectingElements(grid, v1, v2);

	//	create new edges for each edge that is connected with v2.
	//	avoid double edges
		if(grid.num_edges() > 0)
		{
			EdgeDescriptor ed;
			auto iterEnd = grid.associated_edges_end(v2);
			for(auto iter = grid.associated_edges_begin(v2); iter != iterEnd; ++iter)
			{
				Edge* e = *iter;
				if(e->vertex(0) == v2)
					ed.set_vertices(v1, e->vertex(1));
				else
					ed.set_vertices(e->vertex(0), v1);

				Edge* existingEdge = grid.get_edge(ed);
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
			auto iterEnd = grid.associated_faces_end(v2);
			for(auto iter = grid.associated_faces_begin(v2); iter != iterEnd; ++iter)
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
			auto iterEnd = grid.associated_volumes_end(v2);
			for(auto iter = grid.associated_volumes_begin(v2); iter != iterEnd; ++iter)
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
}

////////////////////////////////////////////////////////////////////////
bool IsBoundaryVertex1D(Grid& grid, Vertex* v,
						Grid::edge_traits::callback cbConsiderEdge)
{
	if(!grid.option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we have to enable this option, since nothing works without it in reasonable time.
		UG_LOG("WARNING in IsBoundaryVertex1D(...): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES.\n");
		grid.enable_options(VertexOptions::VRTOPT_STORE_ASSOCIATED_EDGES);
	}

//	iterate over associated edges and return true if only one of them
//	should be considered for the polygonal chain
	size_t counter = 0;
	for(auto iter = grid.associated_edges_begin(v);
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
bool IsBoundaryVertex2D(Grid& grid, Vertex* v)
{
//	check whether one of the associated edges is a boundary edge.
//	if so return true.
	if(!grid.option_is_enabled(FaceOptions::FACEOPT_AUTOGENERATE_EDGES))
	{
	//	we have to enable this option, since we need edges in order to detect boundary vertices.
		UG_LOG("WARNING in IsBoundaryVertex2D(...): auto-enabling FACEOPT_AUTOGENERATE_EDGES.\n");
		grid.enable_options(FaceOptions::FACEOPT_AUTOGENERATE_EDGES);
	}
	if(!grid.option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_EDGES))
	{
	//	we have to enable this option, since nothing works without it in reasonable time.
		UG_LOG("WARNING in IsBoundaryVertex2D(...): auto-enabling VRTOPT_STORE_ASSOCIATED_EDGES.\n");
		grid.enable_options(VertexOptions::VRTOPT_STORE_ASSOCIATED_EDGES);
	}

	for(auto iter = grid.associated_edges_begin(v); iter != grid.associated_edges_end(v); ++iter)
	{
		if(IsBoundaryEdge2D(grid, *iter))
			return true;
	}

	return false;
}

bool IsBoundaryVertex3D(Grid& grid, Vertex* v)
{
//	check whether one of the associated edges is a boundary edge.
//	if so return true.
	if(!grid.option_is_enabled(VolumeOptions::VOLOPT_AUTOGENERATE_FACES))
	{
	//	we have to enable this option, since we need edges in order to detect boundary vertices.
		UG_LOG("  WARNING in IsBoundaryVertex2D(...): auto-enabling VOLOPT_AUTOGENERATE_FACES.\n");
		grid.enable_options(VolumeOptions::VOLOPT_AUTOGENERATE_FACES);
	}
	if(!grid.option_is_enabled(VertexOptions::VRTOPT_STORE_ASSOCIATED_FACES))
	{
	//	we have to enable this option, since nothing works without it in reasonable time.
		UG_LOG("  WARNING in IsBoundaryVertex2D(...): auto-enabling VRTOPT_STORE_ASSOCIATED_FACES.\n");
		grid.enable_options(VertexOptions::VRTOPT_STORE_ASSOCIATED_FACES);
	}

	for(auto iter = grid.associated_faces_begin(v); iter != grid.associated_faces_end(v); ++iter)
	{
		if(IsVolumeBoundaryFace(grid, *iter))
			return true;
	}

	return false;
}

bool LiesOnBoundary(Grid& grid, Vertex* v)
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
bool IsRegularSurfaceVertex(Grid& grid, Vertex* v)
{
//	check how many faces each associated edge has
	auto edgesEnd = grid.associated_edges_end(v);
	for(auto iter = grid.associated_edges_begin(v); iter != edgesEnd; ++iter)
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
	if(sh.num_subsets() <= creaseSI)
		return;
	if(sh.num<Edge>(creaseSI) == 0)
		return;

//	begin marking
	grid.begin_marking();
//	iterate over all crease-edges
	for(EdgeIterator iter = sh.begin<Edge>(creaseSI);
		iter != sh.end<Edge>(creaseSI); ++iter)
	{
	//	check for both vertices whether they are fixed-vertices
		for(int i = 0; i < 2; ++i)
		{
			Vertex* v = (*iter)->vertex(i);
		//	if the vertex is not marked (has not been checked yet)
			if(!grid.is_marked(v))
			{
			//	mark it
				grid.mark(v);
			//	count associated crease edges
				int counter = 0;
				auto aeIterEnd = grid.associated_edges_end(v);
				for(auto aeIter = grid.associated_edges_begin(v); aeIter != aeIterEnd; ++aeIter)
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
