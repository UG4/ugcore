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

#ifndef __H__REMESHING__SIMPLE_GRID__
#define __H__REMESHING__SIMPLE_GRID__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	SimpleGrid
class SimpleGrid
{
	public:
		void clear();
		std::vector<vector3>	vertices;
		std::vector<vector3>	vertexNormals;
		std::vector<int>		triangles;
		std::vector<vector3>	triangleNormals;
};

void PrintSimpleGrid(SimpleGrid& sg);

////////////////////////////////////////////////////////////////////////
//	ObtainSimpleGrid
///	returns a neighbourhood of the edge defined by vrt1 and vrt2 in a SimpleGrid.
/**
 * Uses Grid::mark().
 * Make sure that the queried sub-grid only contains triangles.
 * vrt1 and vrt2 will be the first two vertices in the returned vertex-list.
 * The triangles that contain both vrt1 and vrt2 are the first triangles
 * in the triangle list.
 * The vertices that are contained in those triangles are the next vertices in the
 * vertex list (excluding vrt1 and vrt2).
 *
 * You have to pass an accessor to the positions of the grids vertices (aaPos).
 * Additionally you have to pass an accessor to an integer-attachment at
 * the vertices. This attachment is only intended as a helper during the
 * algorithm. You don't have to initialise it with any special values, nor
 * are the values of this attachment cruical in any consecutive calls.
 *
 * size determines the size of the neighbourhood.
 * If size == 0 then only the triangles that contain the edge (vrt1, vrt2)
 * are copied to sgOut.
 * If size > 0 then each triangle that contains a vertex of the neighbourhood
 * N(size - 1) is contained in N(size).
 */
template <class TPosAcc, class TIntAcc, class TNormAcc>
bool ObtainSimpleGrid(SimpleGrid& sgOut, Grid& grid,
						Vertex* vrt1, Vertex* vrt2, size_t size,
						TPosAcc& aaPos, TNormAcc& aaNorm,
						TIntAcc& aaInt);

////////////////////////////////////////////////////////////////////////
//	ObtainSimpleGrid
///	returns a neighbourhood of the edge e after e has been collapsed.
/**	
 * This algorithm uses Grid::mark.
 *
 * The new vertex which corresponds to the collapsed edge can be found
 * at sgOut.vertices[0].
 *
 * Please note that this method is capable of treating non-surface grids
 * (i.e. grids where edges have more than two adjacent triangles).
 *
 * Make sure that the neighbourhood of e only contains triangles.
 *
 * Please note that the resulting grid is not suited for swap, split
 * or collapse operations.
 */
template <class TPosAcc, class TIntAcc, class TNormAcc>
bool ObtainSimpleGrid_CollapseEdge(SimpleGrid& sgOut, Grid& grid,
						Edge* e, size_t size,
						TPosAcc& aaPos, TNormAcc& aaNorm,
						TIntAcc& aaInt);

////////////////////////////////////////////////////////////////////////
///	caculates the normal of the given triangle and stores it in sg.triangleNormals[triIndex]
/**
 * Make sure that sg.triangleNormas and sg.triangles have the right sizes.
 */
void CalculateTriangleNormal(SimpleGrid& sg, int triIndex);

////////////////////////////////////////////////////////////////////////
///	resizes sg.triangleNormals and calculates them
void CalculateTriangleNormals(SimpleGrid& sg);



////////////////////////////////////////////////////////////////////////
///	the returned degree lies between 0 and 1. The closer to 1 the better.
/**
 * returns the minimal dot-product of each normal of the triangle-corners
 * with the triangles normal.
 */
number GeometricApproximationDegree(SimpleGrid& sg, int triIndex);

////////////////////////////////////////////////////////////////////////
///	sums GeometricApproximationDegree for each triangle.
/**	the returned degree lies between 0 and 1. The closer to 1 the better.*/
number GeometricApproximationDegree(SimpleGrid& sg);

////////////////////////////////////////////////////////////////////////
///	compares the area of the triangle to the length of its edges.
/**	the returned degree lies between 0 and 1. The higher the better.*/
number ShapeQualityDegree(SimpleGrid& sg, int triIndex);

////////////////////////////////////////////////////////////////////////
///	returns the worst quality-degree of the triangles in sg.
/**	the returned degree lies between 0 and 1. The higher the better.*/
number ShapeQualityDegree(SimpleGrid& sg);


////////////////////////////////////////////////////////////////////////
/**	swaps the edge between the first two vertices of the SimpleGrid.
 *	sh has to fullfill the requirements stated in ug::ObtainSimpleGrid.
 *
 *	The manipulated grid is no longer suited for swap, split or
 *	collapse operations.*/
bool SwapEdge(SimpleGrid& sg);

////////////////////////////////////////////////////////////////////////
/**	splits the edge between the first two vertices of the SimpleGrid.
 *	sh has to fullfill the requirements stated in ug::ObtainSimpleGrid.
 *	The new vertexs will be added at the end of the vertex list.
 *	the first two triangles will be altered. Two additional triangles
 *	will be added at the end of the triangle list.
 *
 *	The normal of the new vertex is interpolated from the normals
 *	of the edges end-points.
 *
 *	The manipulated grid is no longer suited for swap, split or
 *	collapse operations.*/
bool SplitEdge(SimpleGrid& sg);

////////////////////////////////////////////////////////////////////////
/**	collapses the edge between the first two vertices of the SimpleGrid.
 *	sh has to fullfill the requirements stated in ug::ObtainSimpleGrid.
 *	two triangles will be removed.
 *
 *	A new vertex will be added at the end of the list. Its normal and
 *	position will be interpolated from the end-points of the edge.
 *	The first two vertices in sg (the endpoints of the collapsed edge)
 *	remain unused.
 *
 *	The manipulated grid is no longer suited for swap, split or
 *	collapse operations.*/
bool CollapseEdge(SimpleGrid& sg);

}//	end of namespace


////////////////////////////////////////////////
//	include implementation
#include "simple_grid_impl.hpp"

#endif
