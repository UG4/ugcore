// created by Sebastian Reiter
// y09 m11 d05
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SERIALIZATION_IMPL__
#define __H__LIB_GRID__SERIALIZATION_IMPL__

#include "serialization.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	DataSerializer
inline void GridDataSerializer::serialize(std::ostream& out, VertexBase* vrt) const
{serialize(out, vrt, m_vrtSerializers);}

inline void GridDataSerializer::serialize(std::ostream& out, EdgeBase* edge) const
{serialize(out, edge, m_edgeSerializers);}

inline void GridDataSerializer::serialize(std::ostream& out, Face* face) const
{serialize(out, face, m_faceSerializers);}

inline void GridDataSerializer::serialize(std::ostream& out, Volume* vol) const
{serialize(out, vol, m_volSerializers);}

template<class TGeomObj, class TSerializers>
void GridDataSerializer::serialize(std::ostream& out, TGeomObj* o,
									TSerializers& serializers) const
{
	for(size_t i = 0; i < serializers.size(); ++i){
		serializers[i](out, o);
	}
}

template <class TIterator>
void GridDataSerializer::
serialize(std::ostream& out, TIterator begin, TIterator end) const
{
	for(TIterator iter = begin; iter != end; ++iter)
		serialize(out, *iter);
}

////////////////////////////////////////////////////////////////////////
inline void GridDataDeserializer::deserialize(std::istream& in, VertexBase* vrt) const
{deserialize(in, vrt, m_vrtDeserializers);}

inline void GridDataDeserializer::deserialize(std::istream& in, EdgeBase* edge) const
{deserialize(in, edge, m_edgeDeserializers);}

inline void GridDataDeserializer::deserialize(std::istream& in, Face* face) const
{deserialize(in, face, m_faceDeserializers);}

inline void GridDataDeserializer::deserialize(std::istream& in, Volume* vol) const
{deserialize(in, vol, m_volDeserializers);}

template<class TGeomObj, class TDeserializers>
void GridDataDeserializer::deserialize(std::istream& in, TGeomObj* o,
									   TDeserializers& deserializers) const
{
	for(size_t i = 0; i < deserializers.size(); ++i){
		deserializers[i](in, o);
	}
}

template <class TIterator>
void GridDataDeserializer::
deserialize(std::istream& in, TIterator begin, TIterator end) const
{
	for(TIterator iter = begin; iter != end; ++iter)
		deserialize(in, *iter);
}


////////////////////////////////////////////////////////////////////////
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

//TODO: remove the following test code.
//	test: write a number-value to check whether it is send correctly
/*
	number tNum = 1247.001234;
	out.write((char*)&tNum, sizeof(number));
*/	
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

	if(tInt != magicNumber){
		UG_LOG("  WARNING: magic-number mismatch before read in DeserializeAttachment. Data-salad possible!\n");
		return false;
	}

//TODO: remove the following test code.
//	test: write a number-value to check whether it is send correctly
/*
	number tNum;
	in.read((char*)&tNum, sizeof(number));
	if(tNum != 1247.001234){
		UG_LOG("TEST-NUMBER TRANSMIT FAILED in DeserializeAttachment!\n");
		return false;
	}
*/
	for(; iterBegin != iterEnd; ++iterBegin)
	{
		in.read((char*)&aa[*iterBegin], sizeof(ValueType));
	}
	in.read((char*)&tInt, sizeof(int));

	if(tInt != magicNumber){
		UG_LOG("  WARNING: magic-number mismatch after read in DeserializeAttachment. Data-salad possible!\n");
		return false;
	}

	return true;
}
}

#endif
