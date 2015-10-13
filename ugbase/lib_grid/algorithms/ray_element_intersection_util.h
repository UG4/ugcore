// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_ray_element_intersection_util
#define __H__UG_ray_element_intersection_util

#include "common/math/ugmath_types.h"
#include "lib_grid/common_attachments.h"
#include "lib_grid/grid/grid.h"

namespace ug{

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

#endif	//__H__UG_ray_element_intersection_util
