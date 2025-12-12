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

#ifndef __H__LIB_GRID__SERIALIZATION_IMPL__
#define __H__LIB_GRID__SERIALIZATION_IMPL__

#include "serialization.h"

namespace ug {

////////////////////////////////////////////////////////////////////////
//	DataSerializer
inline void GridDataSerializationHandler::
serialize(BinaryBuffer& out, Vertex* vrt) const
{
	serialize(out, vrt, m_vrtSerializers);
	serialize(out, vrt, m_gridSerializers);
}

inline void GridDataSerializationHandler::
serialize(BinaryBuffer& out, Edge* edge) const
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

template <typename TIterator>
void GridDataSerializationHandler::
serialize(BinaryBuffer& out, TIterator begin, TIterator end) const
{
	for(TIterator iter = begin; iter != end; ++iter)
		serialize(out, *iter);
}

template <typename TGeomObj, typename TSerializers>
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
deserialize(BinaryBuffer& in, Edge* edge)
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

template <typename TIterator>
void GridDataSerializationHandler::
deserialize(BinaryBuffer& in, TIterator begin, TIterator end)
{
	for(TIterator iter = begin; iter != end; ++iter)
		deserialize(in, *iter);
}

template<typename TGeomObj, typename TDeserializers>
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
template <typename TElem, typename TAttachment>
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
template <typename TElem, typename TAttachment>
bool SerializeAttachment(Grid& grid, TAttachment& attachment,
						 typename geometry_traits<TElem>::iterator iterBegin,
						 typename geometry_traits<TElem>::iterator iterEnd,
						 BinaryBuffer& out)
{
	if(!grid.has_attachment<TElem>(attachment))
		return false;

//	copy data
	Grid::AttachmentAccessor<TElem, TAttachment> aa(grid, attachment);
	using ValueType = typename TAttachment::ValueType;
	
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
template <typename TElem, typename TAttachment>
bool DeserializeAttachment(Grid& grid, TAttachment& attachment,
						 BinaryBuffer& in)
{
	return DeserializeAttachment<TElem, TAttachment>(
					grid, attachment, grid.begin<TElem>(),
					grid.end<TElem>(), in);
}

////////////////////////////////////////////////////////////////////////
//	DeserializeAttachment
template <typename TElem, typename TAttachment>
bool DeserializeAttachment(Grid& grid, TAttachment& attachment,
						 typename geometry_traits<TElem>::iterator iterBegin,
						 typename geometry_traits<TElem>::iterator iterEnd,
						 BinaryBuffer& in)
{
	if(!grid.has_attachment<TElem>(attachment))
		grid.attach_to<TElem>(attachment);

//	copy data
	Grid::AttachmentAccessor<TElem, TAttachment> aa(grid, attachment);
	using ValueType = typename TAttachment::ValueType;
	
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
