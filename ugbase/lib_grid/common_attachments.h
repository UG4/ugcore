/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	In this file attachments are defined that are commonly used by libGrid.
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIB_GRID__COMMON_ATTACHMENTS__
#define __H__LIB_GRID__COMMON_ATTACHMENTS__

#include "common/types.h"
#include "common/math/ugmath_types.h"

#include "lib_grid/attachments/attachment_pipe.h"
//#include "lib_grid/attachments/attachment_info_traits.h"

#include "grid/grid_base_objects.h"

namespace ug {

////////////////////////////////////////////////////////////////////////
//	attachment-types
using ABool = Attachment<bool>;
using AChar = Attachment<char>;
using AByte = Attachment<byte_t>;
using AInt = Attachment<int>;
using AUInt = Attachment<uint>;
using ANumber = Attachment<number>;
using AFloat = Attachment<float>;
using ADouble = Attachment<double>;
using AVector1 = Attachment<vector1>;
using AVector2 = Attachment<vector2>;
using AVector3 = Attachment<vector3>;
using AVector4 = Attachment<vector4>;
using AVertex = Attachment<Vertex*>;
using AEdge = Attachment<Edge*>;
using AFace = Attachment<Face*>;
using AVolume = Attachment<Volume*>;

using APosition1 = AVector1;
using APosition2 = AVector2;
using ANormal2 = AVector3;
using APosition3 = AVector3;
using ANormal3 = AVector3;
using ATexCoord = AVector2;

using APosition = APosition3;
using ANormal = ANormal3;


////////////////////////////////////////////////////////////////////////
//	concrete attachments
///	The standard 3d position type.
UG_API
extern APosition	aPosition;
///	The standard 2d position type
UG_API
extern APosition2	aPosition2;
///	The standard 1d position type
UG_API
extern APosition1	aPosition1;

///	The standard 3d normal type
UG_API
extern ANormal aNormal;

////////////////////////////////////////////////////////////////////////
//	default position attachments for different types
///	this method can be used to retrieve the default position attachments for different types.
/**	Please note that only existing default attachments are returned.
 *
 *	Valid types for TAttachment are:
 *		- APosition (AVector3)		the default 3d position type. Returns aPosition.
 *		- APosition2 (AVector2)		the default 2d position type. Returns aPosition2.
 *		- APosition1 (AVector1)		the default 1d position type. Returns aPosition1.
 */
// ø todo implement or delete | currently not implemented
template <typename TAttachment>
UG_API
inline
TAttachment&
GetDefaultPositionAttachment();


////////////////////////////////////////////////////////////////////////
//	dimension of Position Attachment
///	this function returns the dimension of the position attachment at compile time
template <typename TAPos>
UG_API
inline int GetPositionAttachmentDimension();


}//	end of namespace


////////////////////////////////
//	include implementation
#include "common_attachments_impl.hpp"

#endif
