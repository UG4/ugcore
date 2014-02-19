// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 09.12.2011 (m,d,y)

#ifndef __H__UG__element_storage__
#define __H__UG__element_storage__

#include "grid_base_objects.h"
#include "common/util/section_container.h"

namespace ug
{
///	This struct is used to hold GridObjects and their attachment pipes.
template <class TElem>
class ElementStorage
{
	public:
		typedef ug::AttachmentPipe<TElem*, ElementStorage<TElem> >	AttachmentPipe;
		typedef ug::AttachedElementList<AttachmentPipe>	AttachedElementList;
		typedef ug::SectionContainer<TElem*, AttachedElementList >
			SectionContainer;

		ElementStorage() : m_attachmentPipe(this)
		{
			m_sectionContainer.get_container().set_pipe(&m_attachmentPipe);
		}
	//	the destructor is important, since destruction order is undefined
	//	and since the AttachedElementList in SectionContainer tries to
	//	unregister itself fomt the assigned pipe.
		~ElementStorage(){
			m_sectionContainer.get_container().set_pipe(NULL);
		}

		SectionContainer	m_sectionContainer;///	holds elements
		AttachmentPipe		m_attachmentPipe;///	holds the data of the stored elements.
};


////////////////////////////////////////////////////////////////////////////////
template<>
class attachment_traits<Vertex*, ElementStorage<Vertex> >
{
	public:
		typedef Vertex*&		ElemRef;
		typedef Vertex*			ElemPtr;
		typedef const Vertex*	ConstElemPtr;
		typedef ElementStorage<Vertex>*			ElemHandlerPtr;
		typedef const ElementStorage<Vertex>*	ConstElemHandlerPtr;
		typedef ElementStorage<Vertex>::SectionContainer::iterator	element_iterator;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.begin();}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.end();}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->grid_data_index();}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->set_grid_data_index(index);}
};

template<>
class attachment_traits<Edge*, ElementStorage<Edge> >
{
	public:
		typedef Edge*&			ElemRef;
		typedef Edge*			ElemPtr;
		typedef const Edge*		ConstElemPtr;
		typedef ElementStorage<Edge>*		ElemHandlerPtr;
		typedef const ElementStorage<Edge>*	ConstElemHandlerPtr;
		typedef ElementStorage<Edge>::SectionContainer::iterator	element_iterator;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.begin();}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.end();}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->grid_data_index();}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->set_grid_data_index(index);}
};

template<>
class attachment_traits<Face*, ElementStorage<Face> >
{
	public:
		typedef Face*&			ElemRef;
		typedef Face*			ElemPtr;
		typedef const Face*		ConstElemPtr;
		typedef ElementStorage<Face>*			ElemHandlerPtr;
		typedef const ElementStorage<Face>*		ConstElemHandlerPtr;
		typedef ElementStorage<Face>::SectionContainer::iterator	element_iterator;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.begin();}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.end();}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->grid_data_index();}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->set_grid_data_index(index);}
};

template<>
class attachment_traits<Volume*, ElementStorage<Volume> >
{
	public:
		typedef Volume*&			ElemRef;
		typedef Volume*				ElemPtr;
		typedef const Volume*		ConstElemPtr;
		typedef ElementStorage<Volume>*			ElemHandlerPtr;
		typedef const ElementStorage<Volume>*	ConstElemHandlerPtr;
		typedef ElementStorage<Volume>::SectionContainer::iterator	element_iterator;

		static inline element_iterator elements_begin(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.begin();}
		static inline element_iterator elements_end(ElemHandlerPtr pHandler)	{return pHandler->m_sectionContainer.end();}
		static inline uint get_data_index(ConstElemHandlerPtr pHandler, ConstElemPtr elem)	{return elem->grid_data_index();}
		static inline void set_data_index(ElemHandlerPtr pHandler, ElemPtr elem, uint index){elem->set_grid_data_index(index);}
};



typedef ElementStorage<Vertex>	VertexElementStorage;
typedef ElementStorage<Edge>	EdgeElementStorage;
typedef ElementStorage<Face>		FaceElementStorage;
typedef ElementStorage<Volume>		VolumeElementStorage;



////////////////////////////////////////////////////////////////////////////////
///	Helper class to access the correct element storage from its element type
template <class TElem>
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
template <class TElem>
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
