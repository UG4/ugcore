/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_orientation_util
#define __H__UG_orientation_util

#include "../grid/grid.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
///	inverts the orientation of all elements between elemsBegin and elemsEnd
/**	Make sure that iter_t::value_type is of type Edge*, Face*, or Volume*.*/
template <typename iter_t>
void InvertOrientation(Grid& grid, iter_t elemsBegin,
					   iter_t elemsEnd);

////////////////////////////////////////////////////////////////////////
//	EdgeOrientationMatches
///	checks if the edge-orientation of the edge and the face matches.
/**
 * the match is positive if the face contains the vertices of 'ev'
 * in the same order as ev.
 * please note: if the edge is contained by two faces and both
 * faces have the same edge-orientation as ed, then the face-orientation
 * of the faces differ.
 * \{
 */
UG_API 
bool EdgeOrientationMatches(EdgeVertices* ev, Face* f);

UG_API 
bool OrientationMatches(EdgeVertices* ev, Face* f);
/** \} */

UG_API 
bool OrientationMatches(FaceVertices* fv, Volume* v);

////////////////////////////////////////////////////////////////////////
//	FixOrientation
///	creates uniform orientation of neighboured and boundary faces.
/** This algorithm uses Grid::mark
 *
 * swaps orientation of faces so that all neighboured
 * faces share the same orientation and boundary faces are oriented with
 * an outward normal
 *
 * Value type of TFaceIterator has to be compatible with Face*.
 *
 * Note that all faces between faceBegin and facesEnd have to be members
 * of the specified grid.
 */
template <typename TFaceIterator>
void FixFaceOrientation(Grid& grid, TFaceIterator facesBegin,
						TFaceIterator facesEnd);

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
template<typename TAAPosVRT>
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
template<typename TVolIterator, typename TAAPosVRT>
int
FixOrientation(Grid& grid, TVolIterator volsBegin, TVolIterator volsEnd,
			   TAAPosVRT& aaPosVRT);

}//	end of namespace


////////////////////////////////
// include implementation
#include "orientation_util_impl.hpp"

#endif