/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIBGRID__SUBSET_HANDLER_INTERFACE_IMPL__
#define __H__LIBGRID__SUBSET_HANDLER_INTERFACE_IMPL__

#include "subset_handler_interface.h"

namespace ug {

/*
template <>
inline AttachmentPipe<Vertex*, ISubsetHandler>&
ISubsetHandler::get_attachment_pipe<Vertex>(int subsetIndex)
{
	assert(subset_attachments_are_enabled() && "ERROR - you have to enable subset-attachments for this subset-handler before executing this method.");
	return *m_vertexAttachmentPipes[subsetIndex];
}

template <>
inline AttachmentPipe<Edge*, ISubsetHandler>&
ISubsetHandler::get_attachment_pipe<Edge>(int subsetIndex)
{
	assert(subset_attachments_are_enabled() && "ERROR - you have to enable subset-attachments for this subset-handler before executing this method.");
	return *m_edgeAttachmentPipes[subsetIndex];
}

template <>
inline AttachmentPipe<Face*, ISubsetHandler>&
ISubsetHandler::get_attachment_pipe<Face>(int subsetIndex)
{
	assert(subset_attachments_are_enabled() && "ERROR - you have to enable subset-attachments for this subset-handler before executing this method.");
	return *m_faceAttachmentPipes[subsetIndex];
}

template <>
inline AttachmentPipe<Volume*, ISubsetHandler>&
ISubsetHandler::get_attachment_pipe<Volume>(int subsetIndex)
{
	assert(subset_attachments_are_enabled() && "ERROR - you have to enable subset-attachments for this subset-handler before executing this method.");
	return *m_volumeAttachmentPipes[subsetIndex];
}
*/

inline int ISubsetHandler::
get_subset_index(Vertex* elem) const
{
	if(elements_are_supported(SubsetHandlerElements::SHE_VERTEX))
		return m_aaSubsetIndexVRT[elem];
	return -1;
}

inline int ISubsetHandler::
get_subset_index(Edge* elem) const
{
	if(elements_are_supported(SubsetHandlerElements::SHE_EDGE))
		return m_aaSubsetIndexEDGE[elem];
	return -1;
}

inline int ISubsetHandler::
get_subset_index(Face* elem) const
{
	if(elements_are_supported(SubsetHandlerElements::SHE_FACE))
		return m_aaSubsetIndexFACE[elem];
	return -1;
}

inline int ISubsetHandler::
get_subset_index(Volume* elem) const
{
	if(elements_are_supported(SubsetHandlerElements::SHE_VOLUME))
		return m_aaSubsetIndexVOL[elem];
	return -1;
}

inline void ISubsetHandler::
subset_assigned(Vertex* v, int subsetIndex)
{
	/*if(subset_attachments_are_enabled())
	{
		if(get_subset_index(v) != -1)
			get_attachment_pipe<Vertex>(get_subset_index(v)).unregister_element(v);

		if(subsetIndex != -1)
			get_attachment_pipe<Vertex>(subsetIndex).register_element(v);
	}*/

	m_aaSubsetIndexVRT[v] = subsetIndex;
}

inline void ISubsetHandler::
subset_assigned(Edge* e, int subsetIndex)
{
	/*if(subset_attachments_are_enabled())
	{
		if(get_subset_index(e) != -1)
			get_attachment_pipe<Edge>(get_subset_index(e)).unregister_element(e);

		if(subsetIndex != -1)
			get_attachment_pipe<Edge>(subsetIndex).register_element(e);
	}*/

	m_aaSubsetIndexEDGE[e] = subsetIndex;
}

inline void
ISubsetHandler::
subset_assigned(Face* f, int subsetIndex)
{
	/*if(subset_attachments_are_enabled())
	{
		if(get_subset_index(f) != -1)
			get_attachment_pipe<Face>(get_subset_index(f)).unregister_element(f);

		if(subsetIndex != -1)
			get_attachment_pipe<Face>(subsetIndex).register_element(f);
	}*/

	m_aaSubsetIndexFACE[f] = subsetIndex;
}

inline void
ISubsetHandler::
subset_assigned(Volume* v, int subsetIndex)
{
	/*if(subset_attachments_are_enabled())
	{
		if(get_subset_index(v) != -1)
			get_attachment_pipe<Volume>(get_subset_index(v)).unregister_element(v);

		if(subsetIndex != -1)
			get_attachment_pipe<Volume>(subsetIndex).register_element(v);
	}*/

	m_aaSubsetIndexVOL[v] = subsetIndex;
}

template <typename TIterator>
void ISubsetHandler::
assign_subset(TIterator iterBegin, TIterator iterEnd, int subsetIndex)
{
	typename TIterator::value_type elem;
	while(iterBegin != iterEnd)
	{
		elem = *iterBegin;
		++iterBegin;
		assign_subset(elem, subsetIndex);
	}
}

inline void
ISubsetHandler::
subset_required(int index)
{
	if(index >= (int)m_subsetInfos.size())
		create_required_subsets(index);
}

inline void
ISubsetHandler::
subset_required(int index) const
{
	if(index >= num_subsets()){
		UG_THROW("Can't create new subsets in const ISubsetHandler. "
				 << "num current subsets: " << num_subsets()
				 << " required subset: " << index);
	}
}

////////////////////////////////////////////////////////////////////////
//	attachment handling
/*
template <typename TGeomObjClass>
void ISubsetHandler::
attach_to(IAttachment& attachment, int subsetIndex)
{
	assert(subset_attachments_are_enabled() && "ERROR - you have to enable subset-attachments for this subset-handler before executing this method.");
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_ID != -1,
			invalid_GeomObjClass);

	subset_info_required(subsetIndex);

	int objType = geometry_traits<TGeomObjClass>::BASE_OBJECT_ID;
	switch(objType)
	{
		case VERTEX:
			get_attachment_pipe<Vertex>(subsetIndex).attach(attachment, 0);
			break;
		case EDGE:
			get_attachment_pipe<Edge>(subsetIndex).attach(attachment, 0);
			break;
		case FACE:
			get_attachment_pipe<Face>(subsetIndex).attach(attachment, 0);
			break;
		case VOLUME:
			get_attachment_pipe<Volume>(subsetIndex).attach(attachment, 0);
			break;
	};
}
*/
/*
template <typename TGeomObjClass, typename TAttachment>
void ISubsetHandler::
attach_to_dv(TAttachment& attachment, int subsetIndex,
			const typename TAttachment::ValueType& defaultValue)
{
	assert(subset_attachments_are_enabled() && "ERROR - you have to enable subset-attachments for this subset-handler before executing this method.");
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_ID != -1,
			invalid_GeomObjClass);

	subset_info_required(subsetIndex);

	int objType = geometry_traits<TGeomObjClass>::BASE_OBJECT_ID;
	switch(objType)
	{
		case VERTEX:
			get_attachment_pipe<Vertex>(subsetIndex).attach(attachment, defaultValue, 0);
			break;
		case EDGE:
			get_attachment_pipe<Edge>(subsetIndex).attach(attachment, defaultValue, 0);
			break;
		case FACE:
			get_attachment_pipe<Face>(subsetIndex).attach(attachment, defaultValue, 0);
			break;
		case VOLUME:
			get_attachment_pipe<Volume>(subsetIndex).attach(attachment, defaultValue, 0);
			break;
	};
}
*/
/*
template <typename TGeomObjClass>
void ISubsetHandler::detach_from(IAttachment& attachment, int subsetIndex)
{
	assert(subset_attachments_are_enabled() && "ERROR - you have to enable subset-attachments for this subset-handler before executing this method.");
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_ID != -1,
				invalid_GeomObjClass);

	assert(subsetIndex >= 0 && subsetIndex < (int)num_subsets() && "bad subset index.");

	int objType = geometry_traits<TGeomObjClass>::BASE_OBJECT_ID;
	switch(objType)
	{
		case VERTEX:
			get_attachment_pipe<Vertex>(subsetIndex).detach(attachment);
			break;
		case EDGE:
			get_attachment_pipe<Edge>(subsetIndex).detach(attachment);
			break;
		case FACE:
			get_attachment_pipe<Face>(subsetIndex).detach(attachment);
			break;
		case VOLUME:
			get_attachment_pipe<Volume>(subsetIndex).detach(attachment);
			break;
	};
}
*/
/*
template <typename TGeomObj, typename TAttachment>
inline typename TAttachment::ContainerType*
ISubsetHandler::get_attachment_data_container(TAttachment& attachment, int subsetIndex)
{
	assert(subset_attachments_are_enabled() && "ERROR - you have to enable subset-attachments for this subset-handler before executing this method.");
	return get_attachment_pipe<TGeomObj>(subsetIndex).get_data_container(attachment);
}
*/

////////////////////////////////////////////////////////////////////////
//	attachments_traits
/*
inline uint
attachment_traits<Vertex*, ISubsetHandler>::
get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem)
{
	return pHandler->get_attachment_data_index(elem);
}

inline void
attachment_traits<Vertex*, ISubsetHandler>::
set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index)
{
	pHandler->set_attachment_data_index(elem, index);
}

inline uint
attachment_traits<Edge*, ISubsetHandler>::
get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem)
{
	return pHandler->get_attachment_data_index(elem);
}

inline void
attachment_traits<Edge*, ISubsetHandler>::
set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index)
{
	pHandler->set_attachment_data_index(elem, index);
}

inline uint
attachment_traits<Face*, ISubsetHandler>::
get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem)
{
	return pHandler->get_attachment_data_index(elem);
}

inline void
attachment_traits<Face*, ISubsetHandler>::
set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index)
{
	pHandler->set_attachment_data_index(elem, index);
}

inline uint
attachment_traits<Volume*, ISubsetHandler>::
get_data_index(ElemHandlerPtr pHandler, ConstElemPtr elem)
{
	return pHandler->get_attachment_data_index(elem);
}

inline void
attachment_traits<Volume*, ISubsetHandler>::
set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index)
{
	pHandler->set_attachment_data_index(elem, index);
}
*/
}//	end of namespace

#endif
