/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__element_storage__
#define __H__UG__element_storage__

#include "grid_base_objects.h"
#include "common/util/section_container.h"

namespace ug
{
///	This struct is used to hold GridObjects and their attachment pipes.
template <typename TElem>
class ElementStorage
{
	public:
		typedef AttachmentPipe<TElem*, ElementStorage> AttachmentPipe; // ø using A = A<T>
		typedef AttachedElementList<AttachmentPipe> AttachedElementList; // ø using A = A<T>
		typedef SectionContainer<TElem*, AttachedElementList > SectionContainer; // ø using A = A<T>

		ElementStorage() : m_attachmentPipe(this)
		{
			m_sectionContainer.get_container().set_pipe(&m_attachmentPipe);
		}
	//	the destructor is important, since destruction order is undefined
	//	and since the AttachedElementList in SectionContainer tries to
	//	unregister itself fomt the assigned pipe.
		~ElementStorage(){
			m_sectionContainer.get_container().set_pipe(nullptr);
		}

		SectionContainer m_sectionContainer;///	holds elements
		AttachmentPipe m_attachmentPipe;///	holds the data of the stored elements.
};


////////////////////////////////////////////////////////////////////////////////
template<>
class attachment_traits<Vertex*, ElementStorage<Vertex> >
{
	public:
		using ElemRef = Vertex*&;
		using ElemPtr = Vertex*;
		using ConstElemPtr = const Vertex*;
		using ElemHandlerPtr = ElementStorage<Vertex>*;
		using ConstElemHandlerPtr = const ElementStorage<Vertex>*;
		using element_iterator = ElementStorage<Vertex>::SectionContainer::iterator;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.begin();}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.end();}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->grid_data_index();}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->set_grid_data_index(index);}
};

template<>
class attachment_traits<Edge*, ElementStorage<Edge> >
{
	public:
		using ElemRef = Edge*&;
		using ElemPtr = Edge*;
		using ConstElemPtr = const Edge*;
		using ElemHandlerPtr = ElementStorage<Edge>*;
		using ConstElemHandlerPtr = const ElementStorage<Edge>*;
		using element_iterator = ElementStorage<Edge>::SectionContainer::iterator;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.begin();}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.end();}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->grid_data_index();}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->set_grid_data_index(index);}
};

template<>
class attachment_traits<Face*, ElementStorage<Face> >
{
	public:
		using ElemRef = Face*&;
		using ElemPtr = Face*;
		using ConstElemPtr = const Face*;
		using ElemHandlerPtr = ElementStorage<Face>*;
		using ConstElemHandlerPtr = const ElementStorage<Face>*;
		using element_iterator = ElementStorage<Face>::SectionContainer::iterator;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.begin();}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.end();}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->grid_data_index();}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->set_grid_data_index(index);}
};

template<>
class attachment_traits<Volume*, ElementStorage<Volume> >
{
	public:
		using ElemRef = Volume*&;
		using ElemPtr = Volume*;
		using ConstElemPtr = const Volume*;
		using ElemHandlerPtr = ElementStorage<Volume>*;
		using ConstElemHandlerPtr = const ElementStorage<Volume>*;
		using element_iterator = ElementStorage<Volume>::SectionContainer::iterator;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.begin();}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.end();}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->grid_data_index();}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->set_grid_data_index(index);}
};


using VertexElementStorage = ElementStorage<Vertex>;
using EdgeElementStorage = ElementStorage<Edge>;
using FaceElementStorage = ElementStorage<Face>;
using VolumeElementStorage = ElementStorage<Volume>;



////////////////////////////////////////////////////////////////////////////////
///	Helper class to access the correct element storage from its element type
template <typename TElem>
struct ElementStorageSelector{
		static inline ElementStorage<TElem>& element_storage(VertexElementStorage& vrts,
									EdgeElementStorage& edges, FaceElementStorage& faces,
									VolumeElementStorage& vols);

		static inline const ElementStorage<TElem>& element_storage(
				const VertexElementStorage& vrts, const EdgeElementStorage& edges,
				const FaceElementStorage& faces, const VolumeElementStorage& vols);
};

template <>
struct ElementStorageSelector<Vertex>{
		static inline ElementStorage<Vertex>& element_storage(VertexElementStorage& vrts,
				EdgeElementStorage& edges, FaceElementStorage& faces, VolumeElementStorage& vols)
		{return vrts;}

		static inline const ElementStorage<Vertex>& element_storage(
				const VertexElementStorage& vrts, const EdgeElementStorage& edges,
				const FaceElementStorage& faces, const VolumeElementStorage& vols)
		{return vrts;}
};

template <>
struct ElementStorageSelector<Edge>{
		static inline ElementStorage<Edge>& element_storage(VertexElementStorage& vrts,
				EdgeElementStorage& edges, FaceElementStorage& faces, VolumeElementStorage& vols)
		{return edges;}

		static inline const ElementStorage<Edge>& element_storage(
				const VertexElementStorage& vrts, const EdgeElementStorage& edges,
				const FaceElementStorage& faces, const VolumeElementStorage& vols)
		{return edges;}
};

template <>
struct ElementStorageSelector<Face>{
		static inline ElementStorage<Face>& element_storage(VertexElementStorage& vrts,
				EdgeElementStorage& edges, FaceElementStorage& faces, VolumeElementStorage& vols)
		{return faces;}

		static inline const ElementStorage<Face>& element_storage(
				const VertexElementStorage& vrts, const EdgeElementStorage& edges,
				const FaceElementStorage& faces, const VolumeElementStorage& vols)
		{return faces;}
};

template <>
struct ElementStorageSelector<Volume>{
		static inline ElementStorage<Volume>& element_storage(VertexElementStorage& vrts,
				EdgeElementStorage& edges, FaceElementStorage& faces, VolumeElementStorage& vols)
		{return vols;}

		static inline const ElementStorage<Volume>& element_storage(
				const VertexElementStorage& vrts, const EdgeElementStorage& edges,
				const FaceElementStorage& faces, const VolumeElementStorage& vols)
		{return vols;}
};


////////////////////////////////////////////////////////////////////////////////
///	Helper class to access the correct element storage from its element type
template <typename TElem>
struct SectionContainerSelector{
		static inline typename ElementStorage<TElem>::SectionContainer& section_container(
				VertexElementStorage::SectionContainer& vrts,
				EdgeElementStorage::SectionContainer& edges,
				FaceElementStorage::SectionContainer& faces,
				VolumeElementStorage::SectionContainer& vols);

		static inline const typename ElementStorage<TElem>::SectionContainer& section_container(
				const VertexElementStorage::SectionContainer& vrts,
				const EdgeElementStorage::SectionContainer& edges,
				const FaceElementStorage::SectionContainer& faces,
				const VolumeElementStorage::SectionContainer& vols);

		static inline typename ElementStorage<TElem>::SectionContainer* section_container(
				VertexElementStorage::SectionContainer* vrts,
				EdgeElementStorage::SectionContainer* edges,
				FaceElementStorage::SectionContainer* faces,
				VolumeElementStorage::SectionContainer* vols);

		static inline const typename ElementStorage<TElem>::SectionContainer* section_container(
				const VertexElementStorage::SectionContainer* vrts,
				const EdgeElementStorage::SectionContainer* edges,
				const FaceElementStorage::SectionContainer* faces,
				const VolumeElementStorage::SectionContainer* vols);
};

template <>
struct SectionContainerSelector<Vertex>{
		static inline ElementStorage<Vertex>::SectionContainer& section_container(
				VertexElementStorage::SectionContainer& vrts,
				EdgeElementStorage::SectionContainer& edges,
				FaceElementStorage::SectionContainer& faces,
				VolumeElementStorage::SectionContainer& vols)
		{return vrts;}

		static inline const ElementStorage<Vertex>::SectionContainer& section_container(
				const VertexElementStorage::SectionContainer& vrts,
				const EdgeElementStorage::SectionContainer& edges,
				const FaceElementStorage::SectionContainer& faces,
				const VolumeElementStorage::SectionContainer& vols)
		{return vrts;}

		static inline ElementStorage<Vertex>::SectionContainer* section_container(
				VertexElementStorage::SectionContainer* vrts,
				EdgeElementStorage::SectionContainer* edges,
				FaceElementStorage::SectionContainer* faces,
				VolumeElementStorage::SectionContainer* vols)
		{return vrts;}

		static inline const ElementStorage<Vertex>::SectionContainer* section_container(
				const VertexElementStorage::SectionContainer* vrts,
				const EdgeElementStorage::SectionContainer* edges,
				const FaceElementStorage::SectionContainer* faces,
				const VolumeElementStorage::SectionContainer* vols)
		{return vrts;}
};

template <>
struct SectionContainerSelector<Edge>{
		static inline ElementStorage<Edge>::SectionContainer& section_container(
				VertexElementStorage::SectionContainer& vrts,
				EdgeElementStorage::SectionContainer& edges,
				FaceElementStorage::SectionContainer& faces,
				VolumeElementStorage::SectionContainer& vols)
		{return edges;}

		static inline const ElementStorage<Edge>::SectionContainer& section_container(
				const VertexElementStorage::SectionContainer& vrts,
				const EdgeElementStorage::SectionContainer& edges,
				const FaceElementStorage::SectionContainer& faces,
				const VolumeElementStorage::SectionContainer& vols)
		{return edges;}

		static inline ElementStorage<Edge>::SectionContainer* section_container(
				VertexElementStorage::SectionContainer* vrts,
				EdgeElementStorage::SectionContainer* edges,
				FaceElementStorage::SectionContainer* faces,
				VolumeElementStorage::SectionContainer* vols)
		{return edges;}

		static inline const ElementStorage<Edge>::SectionContainer* section_container(
				const VertexElementStorage::SectionContainer* vrts,
				const EdgeElementStorage::SectionContainer* edges,
				const FaceElementStorage::SectionContainer* faces,
				const VolumeElementStorage::SectionContainer* vols)
		{return edges;}
};

template <>
struct SectionContainerSelector<Face>{
		static inline ElementStorage<Face>::SectionContainer& section_container(
				VertexElementStorage::SectionContainer& vrts,
				EdgeElementStorage::SectionContainer& edges,
				FaceElementStorage::SectionContainer& faces,
				VolumeElementStorage::SectionContainer& vols)
		{return faces;}

		static inline const ElementStorage<Face>::SectionContainer& section_container(
				const VertexElementStorage::SectionContainer& vrts,
				const EdgeElementStorage::SectionContainer& edges,
				const FaceElementStorage::SectionContainer& faces,
				const VolumeElementStorage::SectionContainer& vols)
		{return faces;}

		static inline ElementStorage<Face>::SectionContainer* section_container(
				VertexElementStorage::SectionContainer* vrts,
				EdgeElementStorage::SectionContainer* edges,
				FaceElementStorage::SectionContainer* faces,
				VolumeElementStorage::SectionContainer* vols)
		{return faces;}

		static inline const ElementStorage<Face>::SectionContainer* section_container(
				const VertexElementStorage::SectionContainer* vrts,
				const EdgeElementStorage::SectionContainer* edges,
				const FaceElementStorage::SectionContainer* faces,
				const VolumeElementStorage::SectionContainer* vols)
		{return faces;}
};

template <>
struct SectionContainerSelector<Volume>{
		static inline ElementStorage<Volume>::SectionContainer& section_container(
				VertexElementStorage::SectionContainer& vrts,
				EdgeElementStorage::SectionContainer& edges,
				FaceElementStorage::SectionContainer& faces,
				VolumeElementStorage::SectionContainer& vols)
		{return vols;}

		static inline const ElementStorage<Volume>::SectionContainer& section_container(
				const VertexElementStorage::SectionContainer& vrts,
				const EdgeElementStorage::SectionContainer& edges,
				const FaceElementStorage::SectionContainer& faces,
				const VolumeElementStorage::SectionContainer& vols)
		{return vols;}

		static inline ElementStorage<Volume>::SectionContainer* section_container(
				VertexElementStorage::SectionContainer* vrts,
				EdgeElementStorage::SectionContainer* edges,
				FaceElementStorage::SectionContainer* faces,
				VolumeElementStorage::SectionContainer* vols)
		{return vols;}

		static inline const ElementStorage<Volume>::SectionContainer* section_container(
				const VertexElementStorage::SectionContainer* vrts,
				const EdgeElementStorage::SectionContainer* edges,
				const FaceElementStorage::SectionContainer* faces,
				const VolumeElementStorage::SectionContainer* vols)
		{return vols;}
};

}//	end of namespace

#endif
