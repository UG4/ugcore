//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d13

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	In this file attachments are defined that are commonly used by libGrid.
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIB_GRID__COMMON_ATTACHMENTS__
#define __H__LIB_GRID__COMMON_ATTACHMENTS__

#include "common/types.h"
#include "lib_grid/attachments/attachment_pipe.h"
#include "common/math/ugmath.h"
#include "grid/grid_base_objects.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	attachment-types
typedef Attachment<bool>		ABool;
typedef Attachment<byte>		AByte;
typedef Attachment<int>			AInt;
typedef Attachment<uint>		AUInt;
typedef Attachment<number>		ANumber;
typedef	 Attachment<float>		AFloat;
typedef Attachment<double>		ADouble;
typedef Attachment<vector1>		AVector1;
typedef Attachment<vector2>		AVector2;
typedef Attachment<vector3>		AVector3;
typedef Attachment<vector4>		AVector4;
typedef Attachment<Vertex*>	AVertex;
typedef Attachment<Edge*>	AEdge;
typedef Attachment<Face*>		AFace;
typedef Attachment<Volume*>		AVolume;

typedef AVector1	APosition1;
typedef AVector2	APosition2;
typedef AVector3	ANormal2;
typedef AVector3	APosition3;
typedef AVector3	ANormal3;
typedef AVector2	ATexCoord;

typedef APosition3	APosition;
typedef ANormal3	ANormal;

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
extern ANormal		aNormal;

////////////////////////////////////////////////////////////////////////
//	default position attachments for different types
///	this method can be used to retrieve the default position attachments for different types.
/**	Please note that only existing default attachments are returned.
 *
 *	Valid types for TAttachment are:
 *		- APosition (AVector3)		the default 3d position type. Returns aPosition.
 *		- APosition2 (AVector2)		the default 2d position type. Returns aPosition2.
 */
template <class TAttachment>
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
