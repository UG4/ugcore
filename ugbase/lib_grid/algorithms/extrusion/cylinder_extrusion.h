//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__CYLINDER_EXTRUSION__
#define __H__LIB_GRID__CYLINDER_EXTRUSION__

#include "lib_grid/lg_base.h"

namespace ug
{
/// \addtogroup lib_grid_algorithms_extrusion
///	@{

////////////////////////////////////////////////////////////////////////
///	adapts the grid around the given vertex to a cylinder and extrudes it.
/**
 * If you're not interested in the subsets of the extruded geometry, or if
 * you want to handle them yourself, please use a overloaded version of
 * this method.
 * \todo: 	Add a parameter that allows to specify whether the
 *	 		source-faces of the extrusion shall be erased.
 *
 * \param grid: The Grid
 * \param sh: 	The SubsetHandler
 * \param vrt: 	The Vertex
 * \param direction: the direction
 *
 * \param height: 	The actual extrude-amount is determined by scaling direction
 *					with height. Note that if height is negative, the orientation
 *					of the extruded geometrie will be inverted.
 *
 * \param radius: radius
 *
 * \param rimSnapThreshold:	If a vertex lies closer to the rim than rimSnapThreshold,
 * 							then it will be projected to the rim.
 *
 * \param aaPos: position attachment accessor
 *
 * \param bottomSubInd: default value is -1. Defines the subset into which
 *						the bottom-faces of the cylinder will go
 *
 * \param cylSubInd: 	default value is -1. Defines the subset into which
 *						the bottom-faces of the cylinder will go
 *
 * \param pSel:	You may pass a pointer to a selector that may be used by
 *				this method internally. This makes sense if you call
 *				this method repeatedly, since a repeated allocation and
 *				deallocation can be avoided.
 */
bool ExtrudeCylinder(Grid& grid, SubsetHandler& sh, VertexBase* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int bottomSubInd = -1, int cylSubInd = -1,
					Selector* pSel = NULL);

///	adapts the grid around the given vertex to a cylinder and extrudes it.
/**
 * If you want the method to automatically assign subset-indices to the
 * extruded geometry, then please use an overloaded version of this method.
 *
 * \todo: 	Add a parameter that allows to specify whether the
 *	 		source-faces of the extrusion shall be erased.
 *
 * \param grid: The Grid
 * \param vrt: 	The Vertex
 * \param direction: the direction
 * \param radius: radius
 * \param rimSnapThreshold:	If a vertex lies closer to the rim than rimSnapThreshold,
 * 							then it will be projected to the rim.
 * \param aaPos: position attachment accessor
 *
 * \param height: 	The actual extrude-amount is determined by scaling direction
 *					with height. Note that if height is negative, the orientation
 *					of the extruded geometrie will be inverted.
 *
 * \param pSel:	You may pass a pointer to a selector that may be used by
 *				this method internally. This makes sense if you call
 *				this method repeatedly, since a repeated allocation and
 *				deallocation can be avoided.
 *
 * \param minDot:	Faces whose normal has a dot-product lower that minDot
 *					whith the given direction are not regarded as cylinder-
 *					bottom faces and are thus not extruded.
 */					
bool ExtrudeCylinder(Grid& grid, VertexBase* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					Selector* pSel = NULL);
					
/// @}

}//	end of namespace

#endif
