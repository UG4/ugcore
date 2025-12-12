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

#include "marker_points.h"

#include <sstream>

#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/file_io/file_io.h"

using namespace std;

namespace ug {


MarkerPoint::MarkerPoint() :
	localCoord(0, 0, 0),
	associatedObj(nullptr)
{
}

MarkerPoint::MarkerPoint(const char* strName, const vector3& position,
						 const vector3& normal) :
	name(strName),
	pos(position),
	norm(normal),
	localCoord(0, 0, 0),
	associatedObj(nullptr)
{
}

////////////////////////////////////////////////////////////////////////
///	Loads marker points from a file
bool LoadMarkerPointsFromFile(MarkerPointManager& manager,
							  const char* filename)
{
//	load a grid and copy the points
	Grid grid;
	if(LoadGridFromFile(grid, filename)){
	//	copy the points to the marker-file
		Grid::VertexAttachmentAccessor aaPos(grid, aPosition);
		Grid::VertexAttachmentAccessor<ANormal> aaNorm;
		if(grid.has_vertex_attachment(aNormal))
			aaNorm.access(grid, aNormal);

		int i = 0;
		for(VertexIterator iter = grid.begin<Vertex>();
			iter != grid.end<Vertex>(); ++iter, ++i)
		{
			stringstream ss;
			ss << "marker_" << i;
			if(aaNorm.valid())
				manager.add_marker(MarkerPoint(ss.str().c_str(), aaPos[*iter],
									aaNorm[*iter]));
			else
				manager.add_marker(MarkerPoint(ss.str().c_str(), aaPos[*iter],
									vector3(0, 0, 0)));
		}

		return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////
void SnapMarkerPointToGridVertex(MarkerPoint& markerInOut, Grid& grid,
								 number normalOffset,
								 Grid::VertexAttachmentAccessor<APosition>& aaPos,
								 Grid::VertexAttachmentAccessor<ANormal>* paaNorm)
{
	MarkerPoint& marker = markerInOut;
//	find the closest vertex
//	the specified point
	vector3& p = marker.pos;

//	find the closest point in the grid
	Vertex* vrt;
	{
		VertexIterator iter = grid.begin<Vertex>();
		vrt = *iter;
		number distSq = VecDistanceSq(aaPos[vrt], p);
		++iter;
		for(; iter != grid.end<Vertex>(); ++iter){
			number d = VecDistanceSq(aaPos[*iter], p);
			if(d < distSq){
				distSq = d;
				vrt = *iter;
			}
		}
	}

//	vrt now holds the closest vertex
//	get the normal offset
	vector3 nOff(0, 0, 0);
	if(normalOffset != 0){
		if(paaNorm)
			VecScale(nOff, (*paaNorm)[vrt], normalOffset);
		else{
			CalculateVertexNormal(nOff, grid, vrt, aaPos);
			VecScale(nOff, nOff, normalOffset);
		}
	}

//	the new position
	VecAdd(marker.pos, aaPos[vrt], nOff);

//	the new local coordinate
	marker.localCoord = vector3(0, 0, 0);

//	the new associated geometric object
	marker.associatedObj = vrt;
}

}//	end of namespace
