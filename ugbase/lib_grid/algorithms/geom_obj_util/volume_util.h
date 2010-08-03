//	Sebastian Reiter (sreiter), Martin Stepniewski (mstepnie)
//	s.b.reiter@googlemail.com, mastep@gmx.de
//	y09 m02 d02

#ifndef __H__LIB_GRID__VOLUME_UTIL__
#define __H__LIB_GRID__VOLUME_UTIL__

#include <vector>
#include "lib_grid/lg_base.h"

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
void GetNeighbours(std::vector<Volume*>& vVolsOut, Grid& grid, Volume* v,
					int side, bool clearContainer = true);

////////////////////////////////////////////////////////////////////////
//	CalculateMinTetrahedronHeight - mstepnie
/// calculates the minimal height of a tetrahedral element
//double CalculateMinTetrahedronHeight(Grid& grid, Volume& v);
number CalculateMinTetrahedronHeight(const vector3& a, const vector3& b, 
									 const vector3& c, const vector3& d);


////////////////////////////////////////////////////////////////////////
//	CalculateMinTetrahedronEdge - mstepnie
/// calculates the shortest edge of a tetrahedral element
number CalculateMaxTetrahedronEdgelength(Grid& grid, Volume& v);


////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronAspectRatio - mstepnie
/// calculates the aspect ratio of a tetrahedral element
number CalculateTetrahedronAspectRatio(Grid& grid, Volume& v);


////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronVolume - mstepnie
/// calculates the volume of a tetrahedral element
number CalculateTetrahedronVolume(const vector3& a, const vector3& b,
								  const vector3& c, const vector3& d);

////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron - sreiter
///	returns true if the point lies inside the tetrahedron
inline bool
PointIsInsideTetrahedron(const vector3& v, Tetrahedron* tet,
						 Grid::VertexAttachmentAccessor<APosition>& aaPos);
						 
////////////////////////////////////////////////////////////////////////
template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(Volume* vol, TVertexPositionAttachmentAccessor& aaPosVRT);

}//	end of namespace

/// @}
////////////////////////////////////////
//	include implementation
#include "volume_util_impl.hpp"

#endif
