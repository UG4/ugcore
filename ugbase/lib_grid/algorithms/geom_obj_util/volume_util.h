/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Martin Stepniewski
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

#ifndef __H__LIB_GRID__VOLUME_UTIL__
#define __H__LIB_GRID__VOLUME_UTIL__

#include <vector>
#include "lib_grid/grid/grid.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/common_attachments.h"

namespace ug
{

/**
 * \brief contains methods to manipulate volumes
 * 
 * \defgroup lib_grid_algorithms_volume_util volume util
 * \ingroup lib_grid_algorithms
 * @{
 */

////////////////////////////////////////////////////////////////////////
//	GetNeighbours - sreiter
///	collects neighbours of the given side of a volume.
/**
 * This algorithm uses Grid::mark.
 *
 * if VOLOPT_AUTOGENERATE_FACES and FACEOPT_STORE_ASSOCIATED_VOLUMES are
 * activated in grid, performance is considerably better.
 *
 * collects all volumes that are adjacent to the given side of v.
 */
UG_API 
void GetNeighbours(std::vector<Volume*>& vVolsOut, Grid& grid, Volume* v,
					int side, bool clearContainer = true);


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateVolumeMinHeight - mstepnie
/// calculates the minimal height of a volume element of type tetrahedron
UG_API
number CalculateMinVolumeHeight(Tetrahedron* tet,
								Grid::VertexAttachmentAccessor<APosition>& aaPos);


////////////////////////////////////////////////////////////////////////
//	CalculateMinTetrahedronHeight - mstepnie
/// calculates the minimal height of a tetrahedral element
//double CalculateMinTetrahedronHeight(Grid& grid, Volume& v);
UG_API 
number CalculateMinTetrahedronHeight(const vector3& a, const vector3& b, 
									 const vector3& c, const vector3& d);


////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronAspectRatio - mstepnie
/// calculates the aspect ratio of a tetrahedral element
UG_API 
number CalculateTetrahedronAspectRatio(Grid& grid, Tetrahedron* tet,
							Grid::VertexAttachmentAccessor<APosition>& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronRootMeanSquareFaceArea - mstepnie
////////////////////////////////////////////////////////////////////////////////////////////
UG_API
number CalculateTetrahedronRootMeanSquareFaceArea(Grid& grid,
												  Tetrahedron* tet,
												  Grid::VertexAttachmentAccessor<APosition>& aaPos);


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronVolToRMSFaceAreaRatio - mstepnie
////////////////////////////////////////////////////////////////////////////////////////////
UG_API
number CalculateTetrahedronVolToRMSFaceAreaRatio(Grid& grid,
											  	 Tetrahedron* tet,
												 Grid::VertexAttachmentAccessor<APosition>& aaPos);


////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron - sreiter
///	returns true if the point lies inside the tetrahedron
UG_API 
inline bool
PointIsInsideTetrahedron(const vector3& v, Tetrahedron* tet,
						 Grid::VertexAttachmentAccessor<APosition>& aaPos);


////////////////////////////////////////////////////////////////////////
/// Returns the number of intersections and the intersection points in intsOut
/**
 * This method evaluates the intersection points of a tetrahedron with a
 * plane. There are 3 or 4 points, the number of the points is returned.
 * If there are 4 intersection points, there order in intsOut corresponds to
 * either the clockwise or counterclockwise ordering of the corners of the
 * quadrilateral.
 *
 * Note that the ordering of the corners of the tetrahedron is important and
 * must correspond to the standard ordering!
 *
 * Proper return values: 0, 3 or 4.
 */
size_t IntersectPlaneWithTetrahedron
(
	vector3 intsOut[4], ///< intersection points
	const vector3& planePoint, ///< a point on the plane
	const vector3& planeNormal, ///< a normal to the plane
	const vector3 t[4] ///< coordinates of the corners of the tetrahedron (order is important!)
);

////////////////////////////////////////////////////////////////////////
///	Checks whether a given point lies in the given volume element
/** \note	This method assumes that the given volume element is convex and
 * 			has planar sides.
 */
template <class TAAPos>
UG_API bool
ContainsPoint(Volume* vol, const vector3& p, TAAPos aaPos);


////////////////////////////////////////////////////////////////////////
///	calculates the center of a volume by averaging the positions of its corners
template<class TVertexPositionAttachmentAccessor>
UG_API 
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const VolumeVertices* vol, TVertexPositionAttachmentAccessor& aaPosVRT);

////////////////////////////////////////////////////////////////////////
///	returns the weighted center of the vertices of the given volume
/** TAAWeightVRT has to be an attachment to the vertices of the grid in which
 * v is contained, with ValueType number (or compatible).
 */
template<class TAAPosVRT, class TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const VolumeVertices* vol, TAAPosVRT& aaPos, TAAWeightVRT& aaWeight);


////////////////////////////////////////////////////////////////////////
///	Refines the volume by connecting its sides with the new center.
/**	Make sure that the specified vertex is belongs to the specified grid, and
 * that its position lies inside the specified volume (self-intersections would
 * occur if it would lie outside).
 * The original volume may optionally be deleted.
 */
void InsertCenterVertex(Grid& g, Volume* vol, Vertex* vrt, bool eraseOldVol);


////////////////////////////////////////////////////////////////////////
///	Converts all volumes in the given range to tetrahedra
template <class TVolIter>
void ConvertToTetrahedra (
		Grid& grid,
		TVolIter volsBegin,
		TVolIter volsEnd);

}//	end of namespace

/// @}
////////////////////////////////////////
//	include implementation
#include "volume_util_impl.hpp"

#endif
