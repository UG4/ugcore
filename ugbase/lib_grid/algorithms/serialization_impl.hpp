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
inline void GridDataSerializationHandler::
serialize(BinaryBuffer& out, Vertex* vrt) const
{
	serialize(out, vrt, m_vrtSerializers);
	serialize(out, vrt, m_gridSerializers);
}

inline void GridDataSerializationHandler::
serialize(BinaryBuffer& out, EdgeBase* edge) const
{
	serialize(out, edge, m_edgeSerializers);
	serialize(out, edge, m_gridSerializers);
}

inline void GridDataSerializationHandler::
serialize(BinaryBuffer& out, Face* face) const
{
	serialize(out, face, m_faceSerializers);
	serialize(out, face, m_gridSerializers);
}

inline void GridDataSerializationHandler::
serialize(BinaryBuffer& out, Volume* vol) const
{
	serialize(out, vol, m_volSerializers);
	serialize(out, vol, m_gridSerializers);
}

template <class TIterator>
void GridDataSerializationHandler::
serialize(BinaryBuffer& out, TIterator begin, TIterator end) const
{
	for(TIterator iter = begin; iter != end; ++iter)
		serialize(out, *iter);
}

template<class TGeomObj, class TSerializers>
void GridDataSerializationHandler::
serialize(BinaryBuffer& out, TGeomObj* o,
		  TSerializers& serializers) const
{
//	This method performs the serialization
	for(size_t i = 0; i < serializers.size(); ++i){
		serializers[i]->write_data(out, o);
	}
}


////////////////////////////////////////////////////////////////////////
inline void GridDataSerializationHandler::
deserialize(BinaryBuffer& in, Vertex* vrt)
{
	deserialize(in, vrt, m_vrtSerializers);
	deserialize(in, vrt, m_gridSerializers);
}

inline void GridDataSerializationHandler::
deserialize(BinaryBuffer& in, EdgeBase* edge)
{
	deserialize(in, edge, m_edgeSerializers);
	deserialize(in, edge, m_gridSerializers);
}

inline void GridDataSerializationHandler::
deserialize(BinaryBuffer& in, Face* face)
{
	deserialize(in, face, m_faceSerializers);
	deserialize(in, face, m_gridSerializers);
}

inline void GridDataSerializationHandler::
deserialize(BinaryBuffer& in, Volume* vol)
{
	deserialize(in, vol, m_volSerializers);
	deserialize(in, vol, m_gridSerializers);
}

template <class TIterator>
void GridDataSerializationHandler::
deserialize(BinaryBuffer& in, TIterator begin, TIterator end)
{
	for(TIterator iter = begin; iter != end; ++iter)
		deserialize(in, *iter);
}

template<class TGeomObj, class TDeserializers>
void GridDataSerializationHandler::
deserialize(BinaryBuffer& in, TGeomObj* o,
			TDeserializers& deserializers)
{
//	This method performs the deserialization
	for(size_t i = 0; i < deserializers.size(); ++i){
		deserializers[i]->read_data(in, o);
	}
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SerializeAttachment
template <class TElem, class TAttachment>
bool SerializeAttachment(Grid& grid, TAttachment& attachment,
						 BinaryBuffer& out)
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
						 BinaryBuffer& out)
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
						 BinaryBuffer& in)
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
						 BinaryBuffer& in)
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
