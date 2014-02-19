//	Sebastian Reiter (sreiter), Martin Stepniewski (mstepnie)
//	s.b.reiter@googlemail.com, mastep@gmx.de
//	y09 m02 d02

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


////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron - sreiter
///	returns true if the point lies inside the tetrahedron
UG_API 
inline bool
PointIsInsideTetrahedron(const vector3& v, Tetrahedron* tet,
						 Grid::VertexAttachmentAccessor<APosition>& aaPos);


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
///	returns true if the volume is oriented so that all sides point to this outside.
/**
 * Please note that special cases may exist in which the current implementation
 * may not return the correct result. This could be true for badly shaped
 * Hexahedrons or for degenerated elements (elements which have a side-face,
 * which has no height and thus resembles a line-segment).
 *
 * The current implementation checks if all face-normals point away from
 * the geometrical center.
 *
 * \todo this method could be improved by adding specialised versions for
 *		 the different volume types.
 */
template<class TAAPosVRT>
UG_API 
bool
CheckOrientation(Volume* vol, TAAPosVRT& aaPosVRT);

////////////////////////////////////////////////////////////////////////
///	Changes orientation of badly oriented volumes
/**
 * changes the orientation of volumes so that CheckOrientation returns
 * true for all volumes between volsBegin and volsEnd.
 *
 * Make sure that all volumes between volsBegin and volsEnd are registerd
 * volumes of the specified grid.
 *
 * \return number of reoriented volumes.
 */
template<class TVolIterator, class TAAPosVRT>
UG_API 
int
FixOrientation(Grid& grid, TVolIterator volsBegin, TVolIterator volsEnd,
			   TAAPosVRT& aaPosVRT);


////////////////////////////////////////////////////////////////////////
///	Refines the volume by connecting its sides with the new center.
/**	Make sure that the specified vertex is belongs to the specified grid, and
 * that its position lies inside the specified volume (self-intersections would
 * occur if it would lie outside).
 * The original volume may optionally be deleted.
 */
void InsertCenterVertex(Grid& g, Volume* vol, Vertex* vrt, bool eraseOldVol);

}//	end of namespace

/// @}
////////////////////////////////////////
//	include implementation
#include "volume_util_impl.hpp"

#endif
