// created by Sebastian Reiter
// y09 m11 d05
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SERIALIZATION_IMPL__
#define __H__LIB_GRID__SERIALIZATION_IMPL__

#include "serialization.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	SerializeAttachment
template <class TElem, class TAttachment>
bool SerializeAttachment(Grid& grid, TAttachment& attachment,
						 std::ostream& out)
{
	return SerializeAttachment<TElem, TAttachment>(
								grid, attachment,
								grid.begin<TElem>(),
								grid.end<TElem>(),
								out);
}

////////////////////////////////////////////////////////////////////////
//	SerializeAttachment
template <class TElem, class TAttachment>
bool SerializeAttachment(Grid& grid, TAttachment& attachment,
						 typename geometry_traits<TElem>::iterator iterBegin,
						 typename geometry_traits<TElem>::iterator iterEnd,
						 std::ostream& out)
{
	if(!grid.has_attachment<TElem>(attachment))
		return false;

//	copy data
	Grid::AttachmentAccessor<TElem, TAttachment> aa(grid, attachment);
	typedef typename TAttachment::ValueType ValueType;
	
//	write a magic number at the beginning and at the end.
	int magicNumber = 8304548;
	out.write((char*)&magicNumber, sizeof(int));

	for(; iterBegin != iterEnd; ++iterBegin)
	{
		out.write((char*)&aa[*iterBegin], sizeof(ValueType));
	}

	out.write((char*)&magicNumber, sizeof(int));

	return true;
}


////////////////////////////////////////////////////////////////////////
//	DeserializeAttachment
template <class TElem, class TAttachment>
bool DeserializeAttachment(Grid& grid, TAttachment& attachment,
						 std::istream& in)
{
	return DeserializeAttachment<TElem, TAttachment>(
					grid, attachment, grid.begin<TElem>(),
					grid.end<TElem>(), in);
}

////////////////////////////////////////////////////////////////////////
//	DeserializeAttachment
template <class TElem, class TAttachment>
bool DeserializeAttachment(Grid& grid, TAttachment& attachment,
						 typename geometry_traits<TElem>::iterator iterBegin,
						 typename geometry_traits<TElem>::iterator iterEnd,
						 std::istream& in)
{
	if(!grid.has_attachment<TElem>(attachment))
		grid.attach_to<TElem>(attachment);

//	copy data
	Grid::AttachmentAccessor<TElem, TAttachment> aa(grid, attachment);
	typedef typename TAttachment::ValueType ValueType;
	
//	compare with the magic number
	int magicNumber = 8304548;
	int tInt;
	in.read((char*)&tInt, sizeof(int));

	if(tInt != magicNumber)
		return false;

	for(; iterBegin != iterEnd; ++iterBegin)
	{
		in.read((char*)&aa[*iterBegin], sizeof(ValueType));
	}

	in.read((char*)&tInt, sizeof(int));

	if(tInt != magicNumber)
		return false;

	return true;
}
}

#endif
