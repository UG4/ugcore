//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__ATTACHMENT_UTIL_IMPL__
#define __H__LIB_GRID__ATTACHMENT_UTIL_IMPL__

#include "attachment_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	ConvertMathVectorAttachmentValues
template<class TElem, class TSrcAttachment, class TDestAttachment>
bool ConvertMathVectorAttachmentValues(Grid& grid,
							TSrcAttachment& srcAttachment,
							TDestAttachment& destAttachment)
{
	typedef TSrcAttachment ASrc;
	typedef TDestAttachment ADest;

//	make sure that the attachments are attached correctly
	if(!grid.has_attachment<TElem>(srcAttachment))
		return false;
		
	if(!grid.has_attachment<TElem>(destAttachment))
		grid.attach_to<TElem>(destAttachment);
		
//	get data containers
	typename ASrc::ContainerType& cSrc
			= *grid.get_attachment_data_container<TElem>(srcAttachment);
	typename ADest::ContainerType& cDest
			= *grid.get_attachment_data_container<TElem>(destAttachment);

	int srcDim = ASrc::ValueType::Size;
	int destDim = ADest::ValueType::Size;
	
//	iterate through the elements and copy them
//	if destDim is smaller than srcDim we only have to copy
	if(destDim <= srcDim)
	{
		for(uint i = 0; i < cSrc.size(); ++i)
		{
			typename ASrc::ValueType& vSrc = cSrc[i];
			typename ADest::ValueType& vDest = cDest[i];
			
			for(int j = 0; j < destDim; ++j)
				vDest[j] = vSrc[j];
		}
	}
	else
	{
	//	we have to fill the higher dimensions with 0
		for(uint i = 0; i < cSrc.size(); ++i)
		{
			typename ASrc::ValueType& vSrc = cSrc[i];
			typename ADest::ValueType& vDest = cDest[i];
			
			for(int j = 0; j < srcDim; ++j)
				vDest[j] = vSrc[j];
			
			for(int j = srcDim; j < destDim; ++j)
				vDest[j] = 0;
		}
	}
	
	return true;
}

}//	end of namespace

#endif
