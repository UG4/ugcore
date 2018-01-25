/*
 * Copyright (c) 2012-2018:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

#ifndef __H__UG_ELEMENT_ANGLES
#define __H__UG_ELEMENT_ANGLES

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#endif

/* system includes */
#include <stddef.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "lib_grid/lib_grid.h"


using namespace std;


namespace ug {



////////////////////////////////////////////////////////////////////////////////////////////
//	CollectAssociatedSides
////////////////////////////////////////////////////////////////////////////////////////////

///	Collects all edges (= 2) which exist in the given face and which share the given vertex.
/**	This method uses Grid::mark **/
UG_API
inline void CollectAssociatedSides(Edge* sidesOut[2], Grid& grid, Face* f, Vertex* vrt)
{
	//PROFILE_BEGIN(CollectAssociatedSides_VERTEX);
	sidesOut[0] = NULL;
	sidesOut[1] = NULL;

	grid.begin_marking();
	for(size_t i = 0; i < f->num_vertices(); ++i){
		grid.mark(f->vertex(i));
	}

	//vector<Edge*> vNeighbourEdgesToVertex;
	//CollectAssociated(vNeighbourEdgesToVertex, grid, vrt, true);

	Grid::AssociatedEdgeIterator iterEnd = grid.associated_edges_end(vrt);
	//Grid::AssociatedEdgeIterator iterEnd = vNeighbourEdgesToVertex.end();
	for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(vrt); iter != iterEnd; ++iter)
	//for(Grid::AssociatedEdgeIterator iter = vNeighbourEdgesToVertex.begin(); iter != iterEnd; ++iter)
	{
		Edge* e = *iter;
		if(grid.is_marked(e->vertex(0)) && grid.is_marked(e->vertex(1))){
			UG_ASSERT(	sidesOut[1] == NULL,
						"Only two edges may be adjacent to a vertex in a face element.");

			if(sidesOut[0] == NULL)
				sidesOut[0] = e;
			else
				sidesOut[1] = e;
		}
	}

	grid.end_marking();
	UG_ASSERT(	sidesOut[1] != NULL,
				"Exactly two edges should be adjacent to a vertex in a face element.")
}

///	Collects all faces (= 2) which exist in the given volume and which share the given edge.
/**	This method uses Grid::mark **/
UG_API
inline void CollectAssociatedSides(Face* sidesOut[2], Grid& grid, Volume* v, Edge* e)
{
	//PROFILE_BEGIN(CollectAssociatedSides_EDGE);
	sidesOut[0] = NULL;
	sidesOut[1] = NULL;

	grid.begin_marking();

	for(size_t i = 0; i < v->num_vertices(); ++i)
		grid.mark(v->vertex(i));

	vector<Face*> vNeighbourFacesToEdge;
	CollectAssociated(vNeighbourFacesToEdge, grid, e, true);

	//Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(e);
	Grid::AssociatedFaceIterator iterEnd = vNeighbourFacesToEdge.end();
	//for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(e); iter != iterEnd; ++iter)
	for(Grid::AssociatedFaceIterator iter = vNeighbourFacesToEdge.begin(); iter != iterEnd; ++iter)
	{
		Face* f = *iter;

	//	check whether all vertices of f are marked
		bool allMarked = true;
		for(size_t i = 0; i < f->num_vertices(); ++i){
			if(!grid.is_marked(f->vertex(i))){
				allMarked = false;
				break;
			}
		}

		if(allMarked){
			if(FaceContains(f, e)){
				UG_ASSERT(	sidesOut[1] == NULL,
							"Only two faces may be adjacent to an edge in a volume element.")

				if(sidesOut[0] == NULL)
					sidesOut[0] = f;
				else
					sidesOut[1] = f;
			}
		}
	}

	grid.end_marking();

	UG_ASSERT(	sidesOut[1] != NULL,
				"Exactly two faces should be adjacent to an edge in a volume element.")
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinAngle
////////////////////////////////////////////////////////////////////////////////////////////

///	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMinAngle(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

///	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Face* f, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();

	/**
	 * in the current implementation this method requires, that all edges
	 * are created for all faces.
	 * TODO: improve this!
	 */

	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING: autoenabling GRIDOPT_AUTOGENERATE_SIDES in GetNeighbours(Face).\n");
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	Get type of vertex attachment in aaPos and define it as ValueType
	typedef typename TAAPosVRT::ValueType ValueType;

//	Initialization
	uint numFaceVrts = f->num_vertices();
	ValueType vNorm1, vNorm2;
	ValueType vDir1, vDir2;
	number minAngle = 180.0;
	number tmpAngle;
	Edge* vNeighbourEdgesToVertex[2];
	Vertex* adjacentVrt1;
	Vertex* adjacentVrt2;

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numFaceVrts; ++vrtIter)
	{
		Vertex* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate vDir vectors of the current two adjacent edges
	//	!!!	Beware of the correct order of the vertices	to get correct angle value !!!
		if(vrt != vNeighbourEdgesToVertex[0]->vertex(0))
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(0);
		else
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(1);

		if(vrt != vNeighbourEdgesToVertex[1]->vertex(0))
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(0);
		else
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(1);

		VecSubtract(vDir1, aaPos[adjacentVrt1], aaPos[vrt]);
		VecSubtract(vDir2, aaPos[adjacentVrt2], aaPos[vrt]);

	//	Normalize
		VecNormalize(vDir1, vDir1);
		VecNormalize(vDir2, vDir2);

	//	Calculate current angle
		tmpAngle = acos(VecDot(vDir1, vDir2));

	//	Check for minimality
		if(tmpAngle < minAngle)
		{
			minAngle = tmpAngle;
		}
	}

//	Transform minAngle from RAD to DEG
	minAngle = 180/PI * minAngle;

	return minAngle;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	INFO: For volume elements the minimal angle corresponds to the smallest dihedral
////////////////////////////////////////////////////////////////////////////////////////////

///	Tetrahedron
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(tet), aaPos);
}

///	Prism
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(prism), aaPos);
}

///	Pyramid
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(pyr), aaPos);
}

///	Hexahedron
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Hexahedron* hex, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(hex), aaPos);
}

///	Volume
template <class TAAPosVRT>
number CalculateMinAngle(Grid& grid, Volume* vol, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, vol, aaPos);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinDihedral
////////////////////////////////////////////////////////////////////////////////////////////

///	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

///	Tetrahedron
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(tet), aaPos);
}

///	Prism
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(prism), aaPos);
}

///	Pyramid
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(pyr), aaPos);
}

///	Hexahedron
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Hexahedron* hex, TAAPosVRT& aaPos)
{
	return CalculateMinDihedral(grid, static_cast<Volume*>(hex), aaPos);
}

///	Volume
template <class TAAPosVRT>
number CalculateMinDihedral(Grid& grid, Volume* v, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();

	/**
	 * in the current implementation this method requires, that all edges
	 * are created for all faces.
	 * TODO: improve this!
	 */

	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING: autoenabling GRIDOPT_AUTOGENERATE_SIDES in GetNeighbours(Face).\n");
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	Initialization
	uint numElementEdges = v->num_edges();
	vector3 vNorm1, vNorm2;
	number minDihedral = 360.0;
	number tmpAngle;
	Face* vNeighbourFacesToEdge[2];

//	Iterate over all element edges
	for(uint eIter = 0; eIter < numElementEdges; ++eIter)
	{
		Edge* e = grid.get_edge(v, eIter);

	//	Get adjacent faces at the current edge and calculate the angle between their normals
	//	!!!	Beware of the correct vExtrDir normals to get correct angle value !!!
		CollectAssociatedSides(vNeighbourFacesToEdge, grid, v, e);

		Vertex* adjacentVrt1;
		Vertex* adjacentVrt2;
		vector3 vDir1;
		vector3 vDir2;
		vector3 vDir3;

	//	Get vertex of each associated face, which is not part of the given edge
		for(size_t i = 0; i < vNeighbourFacesToEdge[0]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(1))
				adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(i);
		}

		for(size_t i = 0; i < vNeighbourFacesToEdge[1]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(1))
				adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(i);
		}

		VecSubtract(vDir1, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		VecSubtract(vDir2, aaPos[adjacentVrt1], aaPos[e->vertex(0)]);
		VecSubtract(vDir3, aaPos[e->vertex(0)], aaPos[adjacentVrt2]);

		VecCross(vNorm1, vDir1, vDir2);
		VecCross(vNorm2, vDir1, vDir3);

		VecNormalize(vNorm1, vNorm1);
		VecNormalize(vNorm2, vNorm2);

		tmpAngle = acos(VecDot(vNorm1, vNorm2));
		tmpAngle = PI - tmpAngle;

	//	Check for minimality
		if(tmpAngle < minDihedral)
		{
			minDihedral = tmpAngle;
		}
	}

//	Transform minAngle from RAD to DEG
	minDihedral = 180/PI * minDihedral;

	return minDihedral;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMaxAngle
////////////////////////////////////////////////////////////////////////////////////////////

///	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

///	Face (Triangles and Quadrilaterals)
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Face* f, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();

	/**
	 * in the current implementation this method requires, that all edges
	 * are created for all faces.
	 * TODO: improve this!
	 */

	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING: autoenabling GRIDOPT_AUTOGENERATE_SIDES in GetNeighbours(Face).\n");
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	Get type of vertex attachment in aaPos and define it as ValueType
	typedef typename TAAPosVRT::ValueType ValueType;

//	Initialization
	uint numFaceVrts = f->num_vertices();
	ValueType vNorm1, vNorm2;
	ValueType vDir1, vDir2;
	number maxAngle = 0;
	number tmpAngle;
	Edge* vNeighbourEdgesToVertex[2];
	Vertex* adjacentVrt1;
	Vertex* adjacentVrt2;

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numFaceVrts; ++vrtIter)
	{
		Vertex* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate vExtrDir vectors of the current two adjacent edges
	//	!!!	Beware of the correct order of the vertices	to get correct angle value !!!
		if(vrt != vNeighbourEdgesToVertex[0]->vertex(0))
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(0);
		else
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(1);

		if(vrt != vNeighbourEdgesToVertex[1]->vertex(0))
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(0);
		else
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(1);

		VecSubtract(vDir1, aaPos[adjacentVrt1], aaPos[vrt]);
		VecSubtract(vDir2, aaPos[adjacentVrt2], aaPos[vrt]);

	//	Normalize
		VecNormalize(vDir1, vDir1);
		VecNormalize(vDir2, vDir2);

	//	Calculate current angle
		tmpAngle = acos(VecDot(vDir1, vDir2));

	//	Check for maximality
		if(tmpAngle > maxAngle)
		{
			maxAngle = tmpAngle;
		}
	}

//	Transform maxAngle from RAD to DEG
	maxAngle = 180/PI * maxAngle;

	return maxAngle;
}

////////////////////////////////////////////////////////////////////////////////////////////
//	INFO: For volume elements the maxAngle corresponds to the largest dihedral
////////////////////////////////////////////////////////////////////////////////////////////

///	Tetrahedron
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(tet), aaPos);
}

///	Prism
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(prism), aaPos);
}

///	Pyramid
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(pyr), aaPos);
}

///	Volume
template <class TAAPosVRT>
number CalculateMaxAngle(Grid& grid, Volume* vol, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, vol, aaPos);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMaxDihedral
////////////////////////////////////////////////////////////////////////////////////////////

///	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

///	Tetrahedron
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(tet), aaPos);
}

///	Prism
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Prism* prism, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(prism), aaPos);
}

///	Pyramid
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos)
{
	return CalculateMaxDihedral(grid, static_cast<Volume*>(pyr), aaPos);
}

///	Volume
template <class TAAPosVRT>
number CalculateMaxDihedral(Grid& grid, Volume* v, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();

	/**
	 * in the current implementation this method requires, that all edges
	 * are created for all faces.
	 * TODO: improve this!
	 */

	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING: autoenabling GRIDOPT_AUTOGENERATE_SIDES in GetNeighbours(Face).\n");
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	Initialization
	uint numElementEdges = v->num_edges();
	vector3 vNorm1, vNorm2;
	number maxDihedral = 0;
	number tmpAngle;
	Face* vNeighbourFacesToEdge[2];

//	Iterate over all element edges
	for(uint eIter = 0; eIter < numElementEdges; ++eIter)
	{
		Edge* e = grid.get_edge(v, eIter);

	//	Get adjacent faces at the current edge and calculate the angle between their normals
	//	!!!	Beware of the correct vExtrDir normals to get correct angle value !!!
		CollectAssociatedSides(vNeighbourFacesToEdge, grid, v, e);

		Vertex* adjacentVrt1;
		Vertex* adjacentVrt2;
		vector3 vDir1;
		vector3 vDir2;
		vector3 vDir3;

		for(size_t i = 0; i < vNeighbourFacesToEdge[0]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(1))
				adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(i);
		}

		for(size_t i = 0; i < vNeighbourFacesToEdge[1]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(1))
				adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(i);
		}

		VecSubtract(vDir1, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		VecSubtract(vDir2, aaPos[adjacentVrt1], aaPos[e->vertex(0)]);
		VecSubtract(vDir3, aaPos[e->vertex(0)], aaPos[adjacentVrt2]);

		VecCross(vNorm1, vDir1, vDir2);
		VecCross(vNorm2, vDir1, vDir3);

		VecNormalize(vNorm1, vNorm1);
		VecNormalize(vNorm2, vNorm2);

		tmpAngle = acos(VecDot(vNorm1, vNorm2));
		tmpAngle = PI - tmpAngle;

	//	Check for maximality
		if(tmpAngle > maxDihedral)
		{
			maxDihedral = tmpAngle;
		}
	}

//	Transform maxDihedral from RAD to DEG
	maxDihedral = 180/PI * maxDihedral;

	return maxDihedral;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateAngles
////////////////////////////////////////////////////////////////////////////////////////////

///	An unimplemented version, so that a compile error occurs if no overload exists.
template <class TElem, class TAAPosVRT>
number CalculateAngles(vector<number>& vAnglesOut, Grid& grid, TElem* elem, TAAPosVRT& aaPos);

/// Face
template <class TAAPosVRT>
void CalculateAngles(vector<number>& vAnglesOut, Grid& grid, Face* f, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();

	/**
	 * in the current implementation this method requires, that all edges
	 * are created for all faces.
	 * TODO: improve this!
	 */

	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING: autoenabling GRIDOPT_AUTOGENERATE_SIDES in GetNeighbours(Face).\n");
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	Get type of vertex attachment in aaPos and define it as ValueType
	typedef typename TAAPosVRT::ValueType ValueType;

//	Initialization
	uint numFaceVrts = f->num_vertices();
	ValueType vNorm1, vNorm2;
	ValueType vDir1, vDir2;
	number tmpAngle;
	Edge* vNeighbourEdgesToVertex[2];
	Vertex* adjacentVrt1;
	Vertex* adjacentVrt2;

//	Iterate over all face vertices
	for(uint vrtIter = 0; vrtIter < numFaceVrts; ++vrtIter)
	{
		Vertex* vrt = f->vertex(vrtIter);

	//	Get adjacent edges at the current vertex and calculate the angle between their normals
		CollectAssociatedSides(vNeighbourEdgesToVertex, grid, f, vrt);

	//	Calculate vDir vectors of the current two adjacent edges
	//	!!!	Beware of the correct order of the vertices	to get correct angle value !!!
		if(vrt != vNeighbourEdgesToVertex[0]->vertex(0))
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(0);
		else
			adjacentVrt1 = vNeighbourEdgesToVertex[0]->vertex(1);

		if(vrt != vNeighbourEdgesToVertex[1]->vertex(0))
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(0);
		else
			adjacentVrt2 = vNeighbourEdgesToVertex[1]->vertex(1);

		VecSubtract(vDir1, aaPos[adjacentVrt1], aaPos[vrt]);
		VecSubtract(vDir2, aaPos[adjacentVrt2], aaPos[vrt]);

	//	Normalize
		VecNormalize(vDir1, vDir1);
		VecNormalize(vDir2, vDir2);

	//	Calculate current angle
		tmpAngle = acos(VecDot(vDir1, vDir2));

	//	Transform minAngle from RAD to DEG
		tmpAngle = 180/PI * tmpAngle;

		vAnglesOut.push_back(tmpAngle);
	}
}

/// Volume
template <class TAAPosVRT>
void CalculateAngles(vector<number>& vAnglesOut, Grid& grid, Volume* v, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();

	/**
	 * in the current implementation this method requires, that all edges
	 * are created for all faces.
	 * TODO: improve this!
	 */

	if(!grid.option_is_enabled(GRIDOPT_AUTOGENERATE_SIDES))
	{
		LOG("WARNING: autoenabling GRIDOPT_AUTOGENERATE_SIDES in GetNeighbours(Face).\n");
		grid.enable_options(GRIDOPT_AUTOGENERATE_SIDES);
	}

//	Initialization
	uint numElementEdges = v->num_edges();
	vector3 vNorm1, vNorm2;
	number tmpAngle;
	Face* vNeighbourFacesToEdge[2];

//	Iterate over all element edges
	for(uint eIter = 0; eIter < numElementEdges; ++eIter)
	{
		Edge* e = grid.get_edge(v, eIter);

	//	Get adjacent faces at the current edge and calculate the angle between their normals
	//	!!!	Beware of the correct vExtrDir normals to get correct angle value !!!
		CollectAssociatedSides(vNeighbourFacesToEdge, grid, v, e);

		Vertex* adjacentVrt1;
		Vertex* adjacentVrt2;
		vector3 vDir1;
		vector3 vDir2;
		vector3 vDir3;

		for(size_t i = 0; i < vNeighbourFacesToEdge[0]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[0]->vertex(i) != e->vertex(1))
				adjacentVrt1 = vNeighbourFacesToEdge[0]->vertex(i);
		}

		for(size_t i = 0; i < vNeighbourFacesToEdge[1]->num_vertices(); ++i)
		{
			if(vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(0) && vNeighbourFacesToEdge[1]->vertex(i) != e->vertex(1))
				adjacentVrt2 = vNeighbourFacesToEdge[1]->vertex(i);
		}

		VecSubtract(vDir1, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		VecSubtract(vDir2, aaPos[adjacentVrt1], aaPos[e->vertex(0)]);
		VecSubtract(vDir3, aaPos[e->vertex(0)], aaPos[adjacentVrt2]);

		VecCross(vNorm1, vDir1, vDir2);
		VecCross(vNorm2, vDir1, vDir3);

		VecNormalize(vNorm1, vNorm1);
		VecNormalize(vNorm2, vNorm2);

		tmpAngle = acos(VecDot(vNorm1, vNorm2));
		tmpAngle = PI - tmpAngle;

	//	Transform maxDihedral from RAD to DEG
		tmpAngle = 180.0/PI * tmpAngle;

		vAnglesOut.push_back(tmpAngle);
	}
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithSmallestMinAngle
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestMinAngle(Grid& grid, TIterator elementsBegin, TIterator elementsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
		return NULL;

//	Initializations
	typename TIterator::value_type elementWithSmallestMinAngle = *elementsBegin;
	number smallestMinAngle = CalculateMinAngle(grid, elementWithSmallestMinAngle, aaPos);
	++elementsBegin;

//	compare all volumes and find that one with smallest minAngle
	for(; elementsBegin != elementsEnd; ++elementsBegin)
	{
		typename TIterator::value_type curElement = *elementsBegin;
		number curSmallestMinAngle = CalculateMinAngle(grid, curElement, aaPos);

		if(curSmallestMinAngle < smallestMinAngle)
		{
			elementWithSmallestMinAngle = curElement;
			smallestMinAngle = curSmallestMinAngle;
		}
	}

	return elementWithSmallestMinAngle;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithLargestMaxAngle
template <class TIterator, class TAAPosVRT>
typename TIterator::value_type
FindElementWithLargestMaxAngle(Grid& grid, TIterator elementsBegin, TIterator elementsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return NULL
	if(elementsBegin == elementsEnd)
		return NULL;

//	Initializations
	typename TIterator::value_type elementWithLargestMaxAngle = *elementsBegin;
	number largestMaxAngle = CalculateMaxAngle(grid, elementWithLargestMaxAngle, aaPos);
	++elementsBegin;

//	compare all volumes and find that one with largest maxAngle
	for(; elementsBegin != elementsEnd; ++elementsBegin)
	{
		typename TIterator::value_type curElement = *elementsBegin;
		number curLargestMaxAngle = CalculateMaxAngle(grid, curElement, aaPos);

		if(curLargestMaxAngle > largestMaxAngle)
		{
			elementWithLargestMaxAngle = curElement;
			largestMaxAngle = curLargestMaxAngle;
		}
	}

	return elementWithLargestMaxAngle;
}


}	 
#endif
