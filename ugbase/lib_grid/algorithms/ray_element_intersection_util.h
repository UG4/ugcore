/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_ray_element_intersection_util
#define __H__UG_ray_element_intersection_util

#include "common/math/ugmath_types.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/grid/grid.h"

namespace ug {

///	2d ray-element intersection for edges
/** sminOut and smaxOut return the relative coordinates of the intersection
 *	with the edge regarding the ray-parameter form.
 *	sminOut and smaxOut usually return the same result. They should only be
 *	different if the ray lies in the same hyperplane as the given edge. Then
 *	they would ontain the point where the ray enters (sminOut) and leaves
 *	(smaxOut) the edge. sminOut and smaxOut are relative coordinates regarding
 *	the ray's parameter form. To obtain a point in 2d use
 *	'PointOnRay(from, dir, sminOut)' or 'PointOnRay(from, dir, smaxOut)'.
 *
 *	\return	true if an intersection was found, false if not.
 *
 *	\note	The case in which the ray lies in the plane of an edge is not yet
 *			implemented correctly.
 *
 *	\todo	if a ray lies in the plane of an edge, sminOut and smaxOut are
 *			currently not correctly computed (see note).
 */
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector2& from,
		const vector2& dir,
		Edge* e,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector2> aaPos,
		number sml = SMALL);

///	2d ray-element intersection for faces
/**	If the method returns true, the point where the ray enters and
 *	leaves the face are returned in sminOut and smaxOut respectively.
 *	sminOut and smaxOut are relative coordinates regarding the ray's parameter
 *	form. To obtain a point in 2d use 'PointOnRay(from, dir, sminOut)'
 *	or 'PointOnRay(from, dir, smaxOut)'.
 *
 *	\return	true if an intersection was found, false if not.
 */
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector2& from,
		const vector2& dir,
		Face* f,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector2> aaPos,
		number sml = SMALL);


///	3d ray-element intersection for faces
/** sminOut and smaxOut return the relative coordinates of the intersection
 *	with the face regarding the ray-parameter form.
 *	sminOut and smaxOut usually return the same result. They should only be
 *	different if the ray lies in the plane of the given face. Then they would
 *	contain the point where the ray enters (sminOut) and leaves (smaxOut) the
 *	face. sminOut and smaxOut are relative coordinates regarding the ray's
 *	parameter form. To obtain a point in 3d use 'PointOnRay(from, dir, sminOut)'
 *	or 'PointOnRay(from, dir, smaxOut)'.
 *
 *	\return	true if an intersection was found, false if not.
 *
 *	\note	The case in which the ray lies in the plane of a face is not yet
 *			implemented correctly.
 *
 *	\todo	if a ray lies in the plane of a face, sminOut and smaxOut are
 *			currently not correctly computed (see note).
 */
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector3& from,
		const vector3& dir,
		Face* f,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		number sml = SMALL);

///	3d ray-element intersection for volumes
/**	If the method returns true, the point where the ray enters and
 *	leaves the volume are returned in sminOut and smaxOut respectively.
 *	sminOut and smaxOut are relative coordinates regarding the ray's parameter
 *	form. To obtain a point in 3d use 'PointOnRay(from, dir, sminOut)'
 *	or 'PointOnRay(from, dir, smaxOut)'.
 *
 *	\return	true if an intersection was found, false if not.
 */
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector3& from,
		const vector3& dir,
		Volume* v,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		number sml = SMALL);

///	3d ray-element intersection for edges
/** Not implemented yet, always returns false
 */
bool RayElementIntersection(
		number& sminOut,
		number& smaxOut,
		const vector3& from,
		const vector3& dir,
		Edge* e,
		Grid& g,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		number sml = SMALL);

}//	end of namespace

#endif