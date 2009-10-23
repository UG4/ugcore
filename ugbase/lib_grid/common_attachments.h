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
#include "common/util/attachment_pipe.h"
#include "common/math/ugmath.h"
#include "grid/geometric_base_objects.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	attachment-types
typedef Attachment<int>			AInt;
typedef Attachment<uint>		AUInt;
typedef Attachment<number>		ANumber;
typedef	Attachment<float>		AFloat;
typedef Attachment<double>		ADouble;
typedef Attachment<vector2>		AVector2;
typedef Attachment<vector3>		AVector3;
typedef Attachment<vector4>		AVector4;
typedef Attachment<VertexBase*>	AVertexBase;

typedef AVector2	APosition2;
typedef AVector3	ANormal2;
typedef AVector3	APosition3;
typedef AVector3	ANormal3;
typedef AVector2	ATexCoord;

typedef APosition3	APosition;
typedef ANormal3	ANormal;

////////////////////////////////////////////////////////////////////////
//	concrete attachments
extern APosition	aPosition;
extern ANormal		aNormal;

}//	end of namespace

#endif
