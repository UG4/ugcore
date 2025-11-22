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

#ifndef __H__LIBGRID__MARKER_POINTS__
#define __H__LIBGRID__MARKER_POINTS__

#include <vector>
#include <string>
#include "common/common.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/common_attachments.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
/**
 * Marker points help when it comes to associating something with a
 * free point.
 *
 * Used in Tetrahedralize to define units.
 */
struct MarkerPoint
{
	MarkerPoint();
	MarkerPoint(const char* name, const vector3& position, const vector3& normal);

	std::string name;

	vector3	pos;
	vector3 norm;
	vector3 localCoord;
	GridObject* associatedObj;
};

////////////////////////////////////////////////////////////////////////
/**
 * the marker point manager handles a group of marker points.
 */
class MarkerPointManager
{
	public:
		void clear()									{m_markers.clear();}
		inline void add_marker()						{m_markers.resize(m_markers.size() + 1);}
		inline void add_markers(size_t num = 1)			{m_markers.resize(m_markers.size() + num);}
		void add_marker(const MarkerPoint& marker)		{m_markers.push_back(marker);}

		inline size_t num_markers()	const				{return m_markers.size();}

		inline MarkerPoint& get_marker(size_t index)				{return m_markers[index];}
		inline const MarkerPoint& get_marker(size_t index) const	{return m_markers[index];}

		inline void set_marker(size_t index,
								const MarkerPoint& marker)			{m_markers[index] = marker;}

		inline MarkerPoint const *get_array()						{return &m_markers.front();}

	protected:
		std::vector<MarkerPoint>	m_markers;
};


////////////////////////////////////////////////////////////////////////
///	Loads marker points from a file
bool LoadMarkerPointsFromFile(MarkerPointManager& manager,
							  const char* filename);

////////////////////////////////////////////////////////////////////////
///	Snaps a marker point to a grid vertex
/**
 * Uses the initial position of the marker-point to determine the
 * closest vertex in the grid. The new position is then calculated
 * by an offset of the vertex-normal to the original vertex-position.
 *
 * If a normal-accessor is supplied, the normal of the target vertex
 * will be taken from it. If not, it will be calculated on the fly.
 *
 * Both aaPos and paaNorm (if supplied) have to access valid vertex
 * attachments of the given grid.
 *
 * Please note that marker.norm and marker.name will remain unchanged.
 * However, the other entries will be changed.
 *
 * Note that this method checks each vertex of the grid and is thus a
 * little slow. Alternate implementations should use a kd-tree or similar.
 */
void SnapMarkerPointToGridVertex(MarkerPoint& markerInOut, Grid& grid,
								 number normalOffset,
								 Grid::VertexAttachmentAccessor<APosition>& aaPos,
								 Grid::VertexAttachmentAccessor<ANormal>* paaNorm = nullptr);

}//	end of namespace

#endif
