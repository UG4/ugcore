/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include "neighborhood.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"

using namespace std;

namespace ug
{
							
////////////////////////////////////////////////////////////////////////
//
void CollectNeighbors(std::vector<Vertex*>& vNeighborsOut,
						Grid& grid, Vertex* vrt, uint nbhType,
						Grid::edge_traits::callback considerEdge,
						Grid::face_traits::callback considerFace,
						Grid::volume_traits::callback considerVol)
{
//	clear the container
	vNeighborsOut.clear();
	
//	begin marking
	grid.begin_marking();
	
//	mark vrt - this makes things easier
	grid.mark(vrt);
	
//	iterate through associated edges
	if(nbhType & NHT_EDGE_NEIGHBORS){
		Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(vrt);
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(vrt);
			iter != iterEnd; ++iter)
		{
			if(considerEdge(*iter)){
				Vertex* neighbour = GetConnectedVertex(*iter, vrt);
				if(!grid.is_marked(neighbour)){
					grid.mark(neighbour);
					vNeighborsOut.push_back(neighbour);
				}
			}
		}
	}

//	iterate through associated faces
	if(nbhType & NHT_FACE_NEIGHBORS){
		Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrt);
		for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt);
			iter != iterEnd; ++iter)
		{
			if(considerFace(*iter)){
				Face* f = *iter;
				size_t numVrts = f->num_vertices();
				Face::ConstVertexArray vrts = f->vertices();
				for(size_t i = 0; i < numVrts; ++i){
					Vertex* neighbour = vrts[i];
					if(!grid.is_marked(neighbour)){
						grid.mark(neighbour);
						vNeighborsOut.push_back(neighbour);
					}
				}
			}
		}
	}

//	iterate through associated volumes
	if(nbhType & NHT_VOLUME_NEIGHBORS){
		Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(vrt);
		for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(vrt);
			iter != iterEnd; ++iter)
		{
			if(considerVol(*iter)){
				Volume* v = *iter;
				size_t numVrts = v->num_vertices();
				Volume::ConstVertexArray vrts = v->vertices();
				for(size_t i = 0; i < numVrts; ++i){
					Vertex* neighbour = vrts[i];
					if(!grid.is_marked(neighbour)){
						grid.mark(neighbour);
						vNeighborsOut.push_back(neighbour);
					}
				}
			}
		}
	}
		
//	end marking
	grid.end_marking();
}
							
////////////////////////////////////////////////////////////////////////
//	CollectNeighbors
void CollectNeighbors(std::vector<Edge*>& vNeighborsOut, Edge* e,
					   Grid& grid, NeighborhoodType nbhType)
{
//	clear the container
	vNeighborsOut.clear();
	
//	default neighbourhood:
	if(nbhType == NHT_DEFAULT)
		nbhType = NHT_VERTEX_NEIGHBORS;

//	edges can only be vertex-neighbours
	if(nbhType != NHT_VERTEX_NEIGHBORS)
		return;

//	begin marking
	grid.begin_marking();
	
//	mark the edge
	grid.mark(e);
	
//	mark the vertices of the edge
	grid.mark(e->vertex(0));
	grid.mark(e->vertex(1));	
	
//	iterate over all edges that are connected to the vertices.
//	if the edge is not yet marked, we have to push it to vNeighboursOut.
	for(uint i = 0; i < 2; ++i)
	{
		Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(e->vertex(i));
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(e->vertex(i));
			iter != iterEnd; ++iter)
		{
			if(!grid.is_marked(*iter))
			{
				vNeighborsOut.push_back(*iter);
				grid.mark(*iter);
			}
		}
	}

//	end marking
	grid.end_marking();
}


////////////////////////////////////////////////////////////////////////
//	CollectNeighbors
void CollectNeighbors(std::vector<Face*>& vNeighborsOut, Face* f,
					   Grid& grid, NeighborhoodType nbhType)
{
//	clear the container
	vNeighborsOut.clear();
	
//	default neighbourhood:
	if(nbhType == NHT_DEFAULT)
		nbhType = NHT_EDGE_NEIGHBORS;

//	faces can't be face-neighbours
	if(nbhType == NHT_FACE_NEIGHBORS)
		return;

//	begin marking
	grid.begin_marking();
	
//	mark the face
	grid.mark(f);
	
//	mark the vertices of the face
	uint numVrts = f->num_vertices();
	Face::ConstVertexArray vrts = f->vertices();
	for(uint i = 0; i < numVrts; ++i)
		grid.mark(vrts[i]);
	
//	in order to get the maximum speed-up, we'll try to use
//	associated elements in grid.
	if((nbhType == NHT_EDGE_NEIGHBORS)
		&& grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_EDGES
								  | EDGEOPT_STORE_ASSOCIATED_FACES
								  | FACEOPT_AUTOGENERATE_EDGES))
	{
	//	iterate through associated edges
		Grid::AssociatedEdgeIterator eEnd = grid.associated_edges_end(f);
		for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(f);
			eIter != eEnd; ++eIter)
		{
		//	iterate through associated folumes of the eace
			Grid::AssociatedFaceIterator fEnd = grid.associated_faces_end(*eIter);
			for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(*eIter);
				iter != fEnd; ++iter)
			{
			//	if the face is not yet marked, then add it to the neighbours
				if(!grid.is_marked(*iter))
				{
					grid.mark(*iter);
					vNeighborsOut.push_back(*iter);
				}
			}
		}
	//	we're done in here. end-marking and return.
		grid.end_marking();
		return;
	}


//	iterate over all faces that are connected to the vertices.
//	if the face shares the elements as required by nbhType and
//	it is not yet marked, we have to push it to vNeighboursOut.
//	to optimize speed we'll check both valid nbhTypes separately.
//	the first case indeed is a subcase of the second
//	(compare counted marked vertices against nbhType)
	switch(nbhType)
	{
	case NHT_VERTEX_NEIGHBORS:
		for(uint i = 0; i < numVrts; ++i)
		{
			Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrts[i]);
			for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrts[i]);
				iter != iterEnd; ++iter)
			{
				if(!grid.is_marked(*iter))
				{
					vNeighborsOut.push_back(*iter);
					grid.mark(*iter);
				}
			}
		}
		break;

	case NHT_EDGE_NEIGHBORS:
		for(uint i = 0; i < numVrts; ++i)
		{
			Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrts[i]);
			for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrts[i]);
				iter != iterEnd; ++iter)
			{
				Face* nf = *iter;
				if(!grid.is_marked(nf))
				{
				//	count the marked vertices that are contained by *iter
				//	if there are more than 1, the faces share an edge
				//	(at least in a regular grid)
					int counter = 0;
					
					size_t numNVrts = nf->num_vertices();
					Face::ConstVertexArray nvrts = nf->vertices();

					for(uint j = 0; j < numNVrts; ++j)
					{
						if(grid.is_marked(nvrts[j]))
						{
							++counter;
							if(counter > 1)
							{
								vNeighborsOut.push_back(nf);
								grid.mark(nf);
								break;
							}
						}
					}
				}
			}
		}
		break;
	default:
		break;
	}

//	end marking
	grid.end_marking();
}

////////////////////////////////////////////////////////////////////////
//	CollectNeighbors
void CollectNeighbors(std::vector<Volume*>& vNeighborsOut, Volume* v,
					   Grid& grid, NeighborhoodType nbhType)
{
//	clear the container
	vNeighborsOut.clear();
	
//	default neighbourhood:
	if(nbhType == NHT_DEFAULT)
		nbhType = NHT_FACE_NEIGHBORS;

//	begin marking
	grid.begin_marking();
	
//	mark the volume
	grid.mark(v);
	
//	mark the vertices of the volume
	uint numVrts = v->num_vertices();
	Volume::ConstVertexArray vrts = v->vertices();
	for(uint i = 0; i < numVrts; ++i)
		grid.mark(vrts[i]);

//	in order to get the maximum speed-up, we'll try to use
//	associated elements in grid.
	if((nbhType == NHT_FACE_NEIGHBORS)
		&& grid.option_is_enabled(VOLOPT_STORE_ASSOCIATED_FACES
								  | FACEOPT_STORE_ASSOCIATED_VOLUMES
								  | VOLOPT_AUTOGENERATE_FACES))
	{
	//	iterate through associated faces
		Grid::AssociatedFaceIterator fEnd = grid.associated_faces_end(v);
		for(Grid::AssociatedFaceIterator fIter = grid.associated_faces_begin(v);
			fIter != fEnd; ++fIter)
		{
		//	iterate through associated volumes of the face
			Grid::AssociatedVolumeIterator vEnd = grid.associated_volumes_end(*fIter);
			for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(*fIter);
				iter != vEnd; ++iter)
			{
			//	if the volume is not yet marked, then add it to the neighbours
				if(!grid.is_marked(*iter))
				{
					grid.mark(*iter);
					vNeighborsOut.push_back(*iter);
				}
			}
		}
	//	we're done in here. end-marking and return.
		grid.end_marking();
		return;
	}

//	iterate over all volumes that are connected to the vertices.
//	if the volume shares the elements as required by nbhType and
//	it is not yet marked, we have to push it to vNeighboursOut.
//	to optimize speed we'll check both valid nbhTypes separately.
//	the first case indeed is a subcase of the second
	if(nbhType & NHT_VERTEX_NEIGHBORS)
		for(uint i = 0; i < numVrts; ++i)
		{
			Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(vrts[i]);
			for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(vrts[i]);
				iter != iterEnd; ++iter)
			{
				if(!grid.is_marked(*iter))
				{
					vNeighborsOut.push_back(*iter);
					grid.mark(*iter);
				}
			}
		}
	else
	{
	//	count the number of shared vertices with volumes, which are connected by
	//	at least one vertex. This case handles both face and edge neighborhoods.
		int minNumSharedVrts = -1;
		if(nbhType & NHT_FACE_NEIGHBORS)
			minNumSharedVrts = 3;
		if(nbhType & NHT_EDGE_NEIGHBORS)
			minNumSharedVrts = 2;

		if(minNumSharedVrts > 0){
			for(uint i = 0; i < numVrts; ++i)
			{
				Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(vrts[i]);
				for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(vrts[i]);
					iter != iterEnd; ++iter)
				{
					Volume* nv = *iter;
					if(!grid.is_marked(nv))
					{
					//	count the marked vertices that are contained by *iter
					//	if there as many as in nbhTypes, the volume is a neighbour.
					//	(at least in a regular grid)
						int counter = 0;
						size_t numNVrts = nv->num_vertices();
						Volume::ConstVertexArray nvrts = nv->vertices();
						for(uint j = 0; j < numNVrts; ++j)
						{
							if(grid.is_marked(nvrts[j]))
							{
								++counter;
								if(counter >= minNumSharedVrts)
								{
									vNeighborsOut.push_back(nv);
									grid.mark(nv);
									break;
								}
							}
						}
					}
				}
			}
		}
	}

//	end marking
	grid.end_marking();
}

void CollectNeighborhood(std::vector<Face*>& facesOut, Grid& grid,
						  Vertex* vrt, size_t range,
						  bool clearContainer)
{
	if(clearContainer)
		facesOut.clear();
	
	vector<Vertex*> candidates;
	size_t rangeBegin = 0;
	size_t rangeEnd = 1;
	
	candidates.push_back(vrt);
	
	grid.begin_marking();
	grid.mark(vrt);
	
//	we iterate over the range
	for(size_t i_range = 0; i_range < range; ++i_range){
	//	iterate from candidatesStart to candidatesEnd
	//	this is important, since we can make sure that only triangles
	//	in the correct range are considered
		for(size_t i_vrt = rangeBegin; i_vrt < rangeEnd; ++i_vrt)
		{
			Vertex* v = candidates[i_vrt];
		//	iterate over associated faces
			for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(v);
				iter != grid.associated_faces_end(v); ++iter)
			{
				Face* f = *iter;
				if(!grid.is_marked(f)){
					grid.mark(f);
					facesOut.push_back(f);
					size_t numVrts = f->num_vertices();
					Face::ConstVertexArray vrts = f->vertices();
					for(size_t i = 0; i < numVrts; ++i){
						if(!grid.is_marked(vrts[i])){
							grid.mark(vrts[i]);
							candidates.push_back(vrts[i]);
						}
					}
				}
			}
		}
		
	//	prepare next iteration
		rangeBegin = rangeEnd;
		rangeEnd = candidates.size();
	}
	
	grid.end_marking();
}

}//	end of namespace
