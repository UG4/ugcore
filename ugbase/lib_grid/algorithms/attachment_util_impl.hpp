//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__ATTACHMENT_UTIL_IMPL__
#define __H__LIB_GRID__ATTACHMENT_UTIL_IMPL__

#include "attachment_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	SetAttachmentValues
template <class TAttachmentAccessor, class TIter, class TVal>
void SetAttachmentValues(TAttachmentAccessor& aaVal,
						TIter elemsBegin, TIter elemsEnd,
						const TVal& val)
{
	while(elemsBegin != elemsEnd)
	{
		aaVal[*elemsBegin] = val;
		++elemsBegin;
	}
}

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

////////////////////////////////////////////////////////////////////////
template <class TElem, class TAttachment>
bool CopyAttachments(Grid& srcGrid, TAttachment& aSrc,
					Grid& destGrid, TAttachment& aDest)
{
//	make sure the attachments are properly attached.
	if(!srcGrid.has_attachment<TElem>(aSrc))
		return false;
	
	if(!destGrid.has_attachment<TElem>(aDest))
		destGrid.attach_to<TElem>(aDest);
		
//	access the attachments
	Grid::AttachmentAccessor<TElem, TAttachment> aaSrc(srcGrid, aSrc);
	Grid::AttachmentAccessor<TElem, TAttachment> aaDest(destGrid, aDest);
	
//	iterate over the elements
	typename geometry_traits<TElem>::iterator iterSrc = srcGrid.begin<TElem>();
	typename geometry_traits<TElem>::iterator iterDest = destGrid.begin<TElem>();

	while((iterSrc != srcGrid.end<TElem>())
		&&(iterDest != destGrid.end<TElem>()))
	{
		aaDest[*iterDest] = aaSrc[*iterSrc];
		++iterDest;
		++iterSrc;
	}
	
//	done
	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TElemIter, class TAttachment>
bool CopyAttachments(Grid& grid, TElemIter elemsBegin, TElemIter elemsEnd,
					 TAttachment& aSrc, TAttachment& aDest)
{
	typedef typename PtrToValueType<
			typename TElemIter::value_type>::base_type TElem;

//	make sure the attachments are properly attached.
	if(!grid.has_attachment<TElem>(aSrc))
		return false;

	if(!grid.has_attachment<TElem>(aDest))
		grid.attach_to<TElem>(aDest);

//	access the attachments
	Grid::AttachmentAccessor<TElem, TAttachment> aaSrc(grid, aSrc);
	Grid::AttachmentAccessor<TElem, TAttachment> aaDest(grid, aDest);

//	iterate over the elements
	TElemIter iter = elemsBegin;
	while(iter != elemsEnd)
	{
		aaDest[*iter] = aaSrc[*iter];
		++iter;
	}

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////
template <class TIterator, class TAAInt>
void AssignIndices(TIterator begin, TIterator end,
					TAAInt& aaInt, int baseIndex)
{
	for(TIterator iter = begin; iter != end; ++iter)
		aaInt[*iter] = baseIndex++;
}

////////////////////////////////////////////////////////////////////////
template <class TIterator, class TAttAcc>
TIterator FindElementByValue(TIterator begin, TIterator end,
							 const typename TAttAcc::ValueType& val,
							 TAttAcc& aa)
{
	TIterator iter = begin;
	while(iter != end){
		if(aa[*iter] == val)
			break;
		++iter;
	}
	return iter;
}

}//	end of namespace

#endif
