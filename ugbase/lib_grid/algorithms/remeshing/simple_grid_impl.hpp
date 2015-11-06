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

#ifndef __H__REMESHING__SIMPLE_GRID_IMPL__
#define __H__REMESHING__SIMPLE_GRID_IMPL__

#include <queue>
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/graph/graph.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	ObtainSimpleGrid
template <class TPosAcc, class TIntAcc, class TNormAcc>
bool ObtainSimpleGrid(SimpleGrid& sgOut, Grid& grid,
						Vertex* vrt1, Vertex* vrt2, size_t size,
						TPosAcc& aaPos, TNormAcc& aaNorm,
						TIntAcc& aaInt)
{
//	vVrts will be reused in each call. To avoid unnecessary allocations,
//	we'll reuse this vector.
	static std::vector<Vertex*> vVrts;
	vVrts.clear();
	
//	clear the simple-grid
	sgOut.clear();
	
//	collect vertices and triangles
	grid.begin_marking();
	
//	mark the first two vertices and add them to simple-grid
	grid.mark(vrt1);
	grid.mark(vrt2);
	vVrts.push_back(vrt1);
	vVrts.push_back(vrt2);
	aaInt[vrt1] = 0;
	aaInt[vrt2] = 1;

//	this counter holds the next vertex for which we have to search for
//	associated triangles
	size_t nextVrt = 0;
//	this number holds the first index that is not checked for neighbours.
	size_t vrtsEnd = 2;

//	find the triangles that are adjacent to the edge between vrt1 and vrt2
//	at this point we assume that all associated faces are triangles.
//	If they are not they are simply treated as if they were some.
	Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrt1);
	for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt1);
		iter != iterEnd; ++iter)
	{
		Vertex* vUnmarked = NULL;
		Face* f = *iter;
		int counter = 0;
		
		for(uint j = 0; j < 3; ++j){
			if(grid.is_marked(f->vertex(j)))
				++counter;
			else
				vUnmarked = f->vertex(j);
		}
		
		if(counter > 1){
		//	we found an adjacent triangle. vUnmarked contains the connected vertex
			if(!vUnmarked) goto bail_out;
		//	push the connected vertex to vVrts and assign the index
			aaInt[vUnmarked] = vVrts.size();
			vVrts.push_back(vUnmarked);
		//	add the triangle
			sgOut.triangles.push_back(aaInt[f->vertex(0)]);
			sgOut.triangles.push_back(aaInt[f->vertex(1)]);
			sgOut.triangles.push_back(aaInt[f->vertex(2)]);
		//	mark the face
			grid.mark(f);
		}
	}
	
//	mark the vertices in vVrts that are not yet marked
	for(size_t i = nextVrt; i < vVrts.size(); ++i)
		grid.mark(vVrts[i]);

//	collect all faces in the neighbourhood
	for(size_t i = 0; i < size; ++i)
	{		
		for(; nextVrt < vrtsEnd; ++nextVrt)
		{
			Vertex* vrt = vVrts[nextVrt];
		//	colelct neighbour faces
			Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrt);
			for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt);
				iter != iterEnd; ++iter)
			{
				Face* f = *iter;
			//	if f is unmarked
				if(!grid.is_marked(f)){
				//	add unmarked vertices to vVrts
					for(uint j = 0; j < 3; ++j){
						if(!grid.is_marked(f->vertex(j))){
							aaInt[f->vertex(j)] = vVrts.size();
							grid.mark(f->vertex(j));
							vVrts.push_back(f->vertex(j));
						}
					}

				//	add the triangle
					grid.mark(f);
					sgOut.triangles.push_back(aaInt[f->vertex(0)]);
					sgOut.triangles.push_back(aaInt[f->vertex(1)]);
					sgOut.triangles.push_back(aaInt[f->vertex(2)]);
				}
			}
		}
	//	in the next iteration we'll check all vertices up to this point
		vrtsEnd = vVrts.size();
	}
	
//	copy the vertex-positions and the normals to the grid
	for(size_t i = 0; i < vVrts.size(); ++i)
	{
		sgOut.vertices.push_back(aaPos[vVrts[i]]);
		sgOut.vertexNormals.push_back(aaNorm[vVrts[i]]);
	}
	
//	calculate triangle normals
	CalculateTriangleNormals(sgOut);
	
	grid.end_marking();
	return true;
	
bail_out:
	grid.end_marking();
	return false;
}

////////////////////////////////////////////////////////////////////////
//	ObtainSimpleGrid
template <class TPosAcc, class TIntAcc, class TNormAcc>
bool ObtainSimpleGrid_CollapseEdge(SimpleGrid& sgOut, Grid& grid,
						Edge* e, size_t size,
						TPosAcc& aaPos, TNormAcc& aaNorm,
						TIntAcc& aaInt)
{
//	clear the simple-grid
	sgOut.clear();

//	collect triangles in the neighbourhood of e.
//	Note that faces may (and most likely will) contain some faces twice			
	std::vector<Face*> faces;
	CollectNeighborhood(faces, grid, e->vertex(0), size, false);
	CollectNeighborhood(faces, grid, e->vertex(1), size, false);

//	the first vertex resembles the collapsed edge
	typename TPosAcc::ValueType n;
	VecAdd(n, aaNorm[e->vertex(0)], aaNorm[e->vertex(1)]);
	VecNormalize(n, n);
	sgOut.vertices.push_back(CalculateCenter(e, aaPos));
	sgOut.vertexNormals.push_back(n);

//	now iterate over all associated triangles in the neighbourhood
//	of e and create associated triangles in sgOut.
	grid.begin_marking();
	for(size_t i_face = 0; i_face < faces.size(); ++i_face){
		Face* f = faces[i_face];
	//	make sure that the face is a triangle
		if(f->num_vertices() != 3){
			grid.end_marking();
			return false;
		}
		
	//	avoid multiple insertion of the same face
		if(!grid.is_marked(f)){
			grid.mark(f);

		//	get vertex indices
		//	make sure that the triangles adjacent to e will
		//	not be added to the grid
			int ind[3];
			int edgeVrts = 0;// increase for each vertex that lies on e
			
			for(size_t i = 0; i < 3; ++i){
				Vertex* v = f->vertex(i);
				if((v == e->vertex(0)) || (v == e->vertex(1))){
					ind[i] = 0;
					edgeVrts++;
				}
				else{
				//	get the index of v in sgOut.vertices.
				//	If it hasn't got one, create one.
					if(!grid.is_marked(v)){
					//	NOTE that we add the position even though it is not clear
					//	whether the triangle will be created at all.
						sgOut.vertices.push_back(aaPos[v]);
						sgOut.vertexNormals.push_back(aaNorm[v]);
						aaInt[v] = (int)sgOut.vertices.size() - 1;
						grid.mark(v);
					}
					ind[i] = aaInt[v];
				}
			}
		
		//	add the triangle
			if(edgeVrts < 2){
				sgOut.triangles.push_back(ind[0]);
				sgOut.triangles.push_back(ind[1]);
				sgOut.triangles.push_back(ind[2]);
			}
		}
	}

	grid.end_marking();
	
//	calculate triangle normals
	CalculateTriangleNormals(sgOut);
	
	return true;	
}

}//	end of namespace

#endif
