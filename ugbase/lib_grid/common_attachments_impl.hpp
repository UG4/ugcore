#ifndef __H__LIB_GRID__COMMON_ATTACHMENTS_IMPL__
#define __H__LIB_GRID__COMMON_ATTACHMENTS_IMPL__

namespace ug
{

template <>
inline
APosition&
GetDefaultPositionAttachment<APosition>()
{
	return aPosition;
}

template <>
inline
APosition2&
GetDefaultPositionAttachment<APosition2>()
{
	return aPosition2;
}

template <>
inline
APosition1&
GetDefaultPositionAttachment<APosition1>()
{
	return aPosition1;
}

template <typename TAPos>
inline int GetPositionAttachmentDimension()
{
	return TAPos::ValueType::Size;
}


}//	end of namespace

#endif
