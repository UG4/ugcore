#ifndef __H__LIB_GRID__GRID_ADAPTION__
#define __H__LIB_GRID__GRID_ADAPTION__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_remeshing
///	@{

////////////////////////////////////////////////////////////////////////
///	Performs local remeshing so that the grid is adapted to the given cylinder.
/**
 * This algorithm uses Grid::mark.
 *
 * The resulting grid features an edge-loop that approximates the
 * intersection of the cylinder with the grid-surface.
 * New faces are inserted accordingly.
 *
 * Please note that volumes-geometries are not supported.
 * 
 * The given selector (sel) will contain the faces that lie inside
 * the cylinder when the algorithm is done.
 *
 * aInt has to be attached to the edges of the grid.
 * It will be used to store temporary values. Initial values are ignored.
 *
 * The algorithm requires the option FACEOPT_AUTOGENERATE_EDGES in grid.
 * If it isn't enabled it will be automatically enabled.
 *
 * An overloaded and somewhat slower version of this method exists,
 * which automatically creates the aInt attachment. It is thus
 * a little more comfortable for the caller.
 *
 * \param rimSnapThreshold:	If a vertex lies closer to the rim than rimSnapThreshold,
 * 							then it will be projected to the rim.
 */
bool AdaptSurfaceGridToCylinder(Selector& selOut, Grid& grid,
						   Vertex* vrtCenter, const vector3& normal,
						   number radius, number rimSnapThreshold, AInt& aInt,
						   APosition& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
///	Performs local remeshing so that the grid is adapted to the given cylinder.
/**
 * This is an overloaded version of AdaptSurfaceGridToCylinder.
 * You don't have to pass an edge-integer-attachment to this version,
 * which makes it a little more comfortable, but at the same time
 * slower. You should worry about this slow-down if you intend to
 * call this method multiple times. You should then consider
 * to call the faster version of AdaptSurfaceGridToCylinder,
 * which takes aInt as a parameter.
 *
 * \param rimSnapThreshold:	If a vertex lies closer to the rim than rimSnapThreshold,
 * 							then it will be projected to the rim.
 */
bool AdaptSurfaceGridToCylinder(Selector& selOut, Grid& grid,
						   Vertex* vrtCenter, const vector3& normal,
						   number radius, number rimSnapThreshold,
						   APosition& aPos = aPosition);

/// @}	// end of add_to_group command

}//	end of namespace

#endif
