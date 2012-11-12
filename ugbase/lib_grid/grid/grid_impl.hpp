//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d10

#ifndef __H__LIB_GRID__GRID_IMPLEMENTATION__
#define __H__LIB_GRID__GRID_IMPLEMENTATION__

//#include <cassert>
#include "common/common.h"
#include "common/static_assert.h"
#include "grid_util.h"
#include "grid.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	parallelism
bool Grid::
is_parallel() const
{
	return m_distGridMgr.get() != NULL;
}

DistributedGridManager* Grid::
distributed_grid_manager()
{
	return m_distGridMgr.get();
}

const DistributedGridManager* Grid::
distributed_grid_manager() const
{
	return m_distGridMgr.get();
}


////////////////////////////////////////////////////////////////////////
//	create functions
template<class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::create(GeometricObject* pParent)
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::CONTAINER_SECTION != -1
		&&	geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_geometry_type);

	TGeomObj* geomObj = new TGeomObj;
//	int baseObjectType = geometry_traits<GeomObjType>::base_object_type();
//	geomObj->m_elemHandle = m_elementStorage[baseObjectType].m_sectionContainer.insert_element(geomObj, geometry_traits<GeomObjType>::container_section());
//	m_elementStorage[baseObjectType].m_attachmentPipe.register_element(geomObj);

	register_element(geomObj, pParent);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>(get_iterator(geomObj));
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::create(const typename geometry_traits<TGeomObj>::Descriptor& descriptor,
			GeometricObject* pParent)
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::CONTAINER_SECTION != -1
			&&	geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
			invalid_geometry_type);

	TGeomObj* geomObj = new TGeomObj(descriptor);

//	int baseObjectType = geometry_traits<TGeomObj>::base_object_type();
//	geomObj->m_elemHandle = m_elementStorage[baseObjectType].m_sectionContainer.insert_element(geomObj, geometry_traits<GeomObjType>::container_section());
//	m_elementStorage[baseObjectType].m_attachmentPipe.register_element(geomObj);

	register_element(geomObj, pParent);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>(get_iterator(geomObj));
}

template<class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::create_and_replace(typename geometry_traits<TGeomObj>::geometric_base_object* pReplaceMe)
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::CONTAINER_SECTION != -1
		&&	geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_geometry_type);

	TGeomObj* geomObj = new TGeomObj;

	if(geomObj->reference_object_id() == pReplaceMe->reference_object_id())
	{
		register_and_replace_element(geomObj, pReplaceMe);
		return iterator_cast<typename geometry_traits<TGeomObj>::iterator>(get_iterator(geomObj));
	}
	else
	{
		LOG("ERROR in Grid::create_and_replace(...): reference objects do not match!");
		assert(!"ERROR in Grid::create_and_replace(...): reference objects do not match!");
		delete geomObj;
		return end<TGeomObj>();
	}
}

////////////////////////////////////////////////////////////////////////
template <class TGeomObj>
void Grid::reserve(size_t num)
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
				invalid_geometry_type);

	element_storage<TGeomObj>().m_attachmentPipe.reserve(num);
}

////////////////////////////////////////////////////////////////////////
//	erase
template <class GeomObjIter>
void Grid::erase(const GeomObjIter& iterBegin, const GeomObjIter& iterEnd)
{
	GeomObjIter iter = iterBegin;
	while(iter != iterEnd)
	{
		GeomObjIter tmpIter = iter;
		++iter;
		erase(*tmpIter);
	}
}

template <class TGeomObj>
void Grid::clear()
{
	while(begin<TGeomObj>() != end<TGeomObj>())
		erase(*begin<TGeomObj>());
}

////////////////////////////////////////////////////////////////////////
//	Iterators
template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::begin()
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_GeomObj);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(element_storage<TGeomObj>().m_sectionContainer.section_begin(geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::end()
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_GeomObj);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(element_storage<TGeomObj>().m_sectionContainer.section_end(geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::const_iterator
Grid::begin() const
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_GeomObj);

	return iterator_cast<typename geometry_traits<TGeomObj>::const_iterator>
		(element_storage<TGeomObj>().m_sectionContainer.section_begin(geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::const_iterator
Grid::end() const
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_GeomObj);

	return iterator_cast<typename geometry_traits<TGeomObj>::const_iterator>
		(element_storage<TGeomObj>().m_sectionContainer.section_end(geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

template <class TGeomObj>
TGeomObj*
Grid::front()
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_GeomObj);

	return static_cast<TGeomObj*>(*element_storage<TGeomObj>().m_sectionContainer.
										front(geometry_traits<TGeomObj>::CONTAINER_SECTION));
}

template <class TGeomObj>
TGeomObj*
Grid::back()
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_GeomObj);

	return static_cast<TGeomObj*>(*element_storage<TGeomObj>().m_sectionContainer.
										back(geometry_traits<TGeomObj>::CONTAINER_SECTION));
}
////////////////////////////////////////////////////////////////////////
//	element numbers
template <class TGeomObj>
size_t Grid::num() const
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
		invalid_GeomObj);

	int secIndex = geometry_traits<TGeomObj>::CONTAINER_SECTION;

	if(secIndex == -1)
		return element_storage<TGeomObj>().m_sectionContainer.num_elements();

	return element_storage<TGeomObj>().m_sectionContainer.num_elements(secIndex);
}

inline void Grid::
objects_will_be_merged(VertexBase* target, VertexBase* elem1,
						VertexBase* elem2)
{
	for(Grid::ObserverContainer::iterator iter = m_vertexObservers.begin();
		iter != m_vertexObservers.end(); iter++)
	{
		(*iter)->vertices_to_be_merged(this, target, elem1, elem2);
	}
}

inline void Grid::
objects_will_be_merged(EdgeBase* target, EdgeBase* elem1,
						EdgeBase* elem2)
{
	for(Grid::ObserverContainer::iterator iter = m_edgeObservers.begin();
		iter != m_edgeObservers.end(); iter++)
	{
		(*iter)->edges_to_be_merged(this, target, elem1, elem2);
	}
}

inline void Grid::
objects_will_be_merged(Face* target, Face* elem1,
						Face* elem2)
{
	for(Grid::ObserverContainer::iterator iter = m_faceObservers.begin();
		iter != m_faceObservers.end(); iter++)
	{
		(*iter)->faces_to_be_merged(this, target, elem1, elem2);
	}
}

inline void Grid::
objects_will_be_merged(Volume* target, Volume* elem1,
						Volume* elem2)
{
	for(Grid::ObserverContainer::iterator iter = m_volumeObservers.begin();
		iter != m_volumeObservers.end(); iter++)
	{
		(*iter)->volumes_to_be_merged(this, target, elem1, elem2);
	}
}

template <class TGeomObj>
size_t Grid::attachment_container_size() const
{
	return element_storage<TGeomObj>().m_attachmentPipe.num_data_entries();
}

////////////////////////////////////////////////////////////////////////
//	attachment handling
template <class TGeomObjClass>
void Grid::attach_to(IAttachment& attachment, bool passOnValues)
{
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_ID != -1,
			invalid_GeomObjClass);

//	setup the options for this attachment.
	int options = 0;
	if(passOnValues)
		options = 1;

	element_storage<TGeomObjClass>().m_attachmentPipe.attach(attachment, options);
}

inline void Grid::attach_to_all(IAttachment& attachment, bool passOnValues)
{
	attach_to<VertexBase>(attachment, passOnValues);
	attach_to<EdgeBase>(attachment, passOnValues);
	attach_to<Face>(attachment, passOnValues);
	attach_to<Volume>(attachment, passOnValues);
}

inline void Grid::attach_to_all(IAttachment& attachment)
{
	attach_to<VertexBase>(attachment);
	attach_to<EdgeBase>(attachment);
	attach_to<Face>(attachment);
	attach_to<Volume>(attachment);
}

template <class TGeomObjClass, class TAttachment>
void Grid::attach_to_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)
{
	attach_to_dv<TGeomObjClass, TAttachment>(attachment, defaultValue, attachment.default_pass_on_behaviour());
}

template <class TAttachment>
inline void Grid::
attach_to_all_dv(TAttachment& attachment,
				 const typename TAttachment::ValueType& defaultValue)
{
	attach_to_dv<VertexBase>(attachment, defaultValue);
	attach_to_dv<EdgeBase>(attachment, defaultValue);
	attach_to_dv<Face>(attachment, defaultValue);
	attach_to_dv<Volume>(attachment, defaultValue);
}

template <class TGeomObjClass, class TAttachment>
void Grid::attach_to_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)
{
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_ID != -1,
			invalid_GeomObjClass);

//	setup the options for this attachment.
	int options = 0;
	if(passOnValues)
		options = 1;

	element_storage<TGeomObjClass>().m_attachmentPipe.attach(attachment, defaultValue, options);
}

template <class TAttachment>
inline void Grid::
attach_to_all_dv(TAttachment& attachment,
				 const typename TAttachment::ValueType& defaultValue,
				 bool passOnValues)
{
	attach_to_dv<VertexBase>(attachment, defaultValue, passOnValues);
	attach_to_dv<EdgeBase>(attachment, defaultValue, passOnValues);
	attach_to_dv<Face>(attachment, defaultValue, passOnValues);
	attach_to_dv<Volume>(attachment, defaultValue, passOnValues);
}

template <class TGeomObjClass>
void Grid::detach_from(IAttachment& attachment)
{
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_ID != -1,
				invalid_GeomObjClass);

	element_storage<TGeomObjClass>().m_attachmentPipe.detach(attachment);
}

inline void Grid::detach_from_all(IAttachment& attachment)
{
	detach_from<VertexBase>(attachment);
	detach_from<EdgeBase>(attachment);
	detach_from<Face>(attachment);
	detach_from<Volume>(attachment);
}

/*
template <class TGeomObjClass>
util::IAttachmentDataContainer* Grid::get_data_container(util::IAttachment& attachment)
{
	assert(geometry_traits<TGeomObjClass>::base_object_type() != -1
			&& "ERROR in Grid::get_data_container(...). Invalid base_object_type of GeomObjClass!");

	return m_elementStorage[geometry_traits<TGeomObjClass>::base_object_type()].
		m_attachmentPipe.get_data_container(attachment);
}
*/

template <class TGeomObj, class TAttachment>
typename TAttachment::ContainerType*
Grid::get_attachment_data_container(TAttachment& attachment)
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
			invalid_GeomObj);

	return element_storage<TGeomObj>().m_attachmentPipe.get_data_container(attachment);
}

template <class TGeomObj>
typename Grid::traits<TGeomObj>::AttachmentPipe&
Grid::get_attachment_pipe()
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_ID != -1,
			invalid_GeomObj);

	return element_storage<TGeomObj>().m_attachmentPipe;
}

template <class TGeomObj>
uint
Grid::get_attachment_data_index(TGeomObj* pObj) const
{
	typedef typename geometry_traits<TGeomObj>::geometric_base_object BaseObj;
	return attachment_traits<BaseObj*, ElementStorage<BaseObj> >::
			get_data_index(&element_storage<TGeomObj>(), pObj);
}

inline void
Grid::autoenable_option(uint option, const char* caller, const char* optionName)
{
	if(!option_is_enabled(option))
	{
		LOG("WARNING in " << caller << ": auto-enabling " << optionName << "." << std::endl);
		enable_options(option);
	}
}

////////////////////////////////////////////////////////////////////////
//	marks
template <class TIterator>
void Grid::mark(TIterator begin, TIterator end)
{
	for(TIterator iter = begin; iter != end; ++iter)
		mark(*iter);
}

template <class TIterator>
void Grid::unmark(TIterator begin, TIterator end)
{
	for(TIterator iter = begin; iter != end; ++iter)
		unmark(*iter);
}


////////////////////////////////////////////////////////////////////////
template <class TElem>
void Grid::associated_elements(traits<VertexBase>::secure_container& elemsOut, TElem* e)
{
	get_associated(elemsOut, e);
}

template <class TElem>
void Grid::associated_elements(traits<EdgeBase>::secure_container& elemsOut, TElem* e)
{
	get_associated(elemsOut, e);
}

template <class TElem>
void Grid::associated_elements(traits<Face>::secure_container& elemsOut, TElem* e)
{
	get_associated(elemsOut, e);
}

template <class TElem>
void Grid::associated_elements(traits<Volume>::secure_container& elemsOut, TElem* e)
{
	get_associated(elemsOut, e);
}


template <class TElem>
void Grid::associated_elements_sorted(traits<VertexBase>::secure_container& elemsOut, TElem* e)
{
	get_associated_sorted(elemsOut, e);
}

template <class TElem>
void Grid::associated_elements_sorted(traits<EdgeBase>::secure_container& elemsOut, TElem* e)
{
	get_associated_sorted(elemsOut, e);
}

template <class TElem>
void Grid::associated_elements_sorted(traits<Face>::secure_container& elemsOut, TElem* e)
{
	get_associated_sorted(elemsOut, e);
}

template <class TElem>
void Grid::associated_elements_sorted(traits<Volume>::secure_container& elemsOut, TElem* e)
{
	get_associated_sorted(elemsOut, e);
}


template <class TElem>
void Grid::get_associated_sorted(typename traits<TElem>::secure_container& elems, TElem* e)
{
//	we have to retrieve a valid pointer on the element pointer. The only way
//	to receive it, is to return the pointer to the entry in which e is stored in
//	the element storage.
	elems.set_external_array(
		element_storage<typename TElem::geometric_base_object>().m_sectionContainer.
			get_container().get_pointer_to_element(e),
		1);
}


////////////////////////////////////////////////////////////////////////
//	neighbourhood access
template <class TGeomObj>
EdgeBase* Grid::find_edge_in_associated_edges(TGeomObj* obj,
												EdgeVertices& ev)
{
	GRID_PROFILE_FUNC();

	AssociatedEdgeIterator iterEnd = associated_edges_end(obj);
	for(AssociatedEdgeIterator iter = associated_edges_begin(obj);
		iter != iterEnd; ++iter)
	{
		EdgeBase* e = *iter;
		if(CompareVertices(e, &ev))
			return e;
	}
	return NULL;
}

										
template <class TGeomObj>
Face* Grid::find_face_in_associated_faces(TGeomObj* obj,
											FaceVertices& fv)
{
	GRID_PROFILE_FUNC();

	unsigned long key = hash_key(&fv);
	AssociatedFaceIterator iterEnd = associated_faces_end(obj);
	for(AssociatedFaceIterator iter = associated_faces_begin(obj);
		iter != iterEnd; ++iter)
	{
		Face* f = *iter;
		if(key == hash_key(f))
		{
			if(CompareVertices(f, &fv))
				return f;
		}
	}
	
	return NULL;
}
										
template <class TGeomObj>
Volume* Grid::find_volume_in_associated_volumes(TGeomObj* obj,
												VolumeVertices& vv)
{
	GRID_PROFILE_FUNC();

	unsigned long key = hash_key(&vv);
	AssociatedVolumeIterator iterEnd = associated_volumes_end(obj);
	for(AssociatedVolumeIterator iter = associated_volumes_begin(obj);
		iter != iterEnd; ++iter)
	{
		Volume* v = *iter;
		if(key == hash_key(v))
		{
			if(CompareVertices(v, &vv))
				return v;
		}
	}
	
	return NULL;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of Grids AttachmentAccessors

////////////////////////////////////////////////////////////////////////
//	AttachmentAccessor
template <class TElem, class TAttachment>
Grid::AttachmentAccessor<TElem, TAttachment>::
AttachmentAccessor() :
ug::AttachmentAccessor<TElem*, TAttachment, ElementStorage<TElem> >()
{
}

template <class TElem, class TAttachment>
Grid::AttachmentAccessor<TElem, TAttachment>::
AttachmentAccessor(const AttachmentAccessor& aa) :
ug::AttachmentAccessor<TElem*, TAttachment, ElementStorage<TElem> >(aa)
{
}

template <class TElem, class TAttachment>
Grid::AttachmentAccessor<TElem, TAttachment>::
AttachmentAccessor(Grid& grid, TAttachment& a) :
ug::AttachmentAccessor<typename TElem::geometric_base_object*, TAttachment,
					   typename traits<TElem>::ElementStorage>
	(grid.get_attachment_pipe<TElem>(), a)
{
}

template <class TElem, class TAttachment>
Grid::AttachmentAccessor<TElem, TAttachment>::
AttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach) :
ug::AttachmentAccessor<typename TElem::geometric_base_object*, TAttachment,
					   typename traits<TElem>::ElementStorage>()
{
	if(autoAttach){
		if(!grid.has_attachment<TElem>(a))
			grid.attach_to<TElem>(a);
	}

	access(grid, a);
}

////////////////////////////////////////////////////////////////////////
//	VertexAttachmentAccessor
template <class TAttachment>
Grid::VertexAttachmentAccessor<TAttachment>::
VertexAttachmentAccessor() :
	Grid::AttachmentAccessor<VertexBase, TAttachment>()
{
}

template <class TAttachment>
Grid::VertexAttachmentAccessor<TAttachment>::
VertexAttachmentAccessor(const VertexAttachmentAccessor& aa) :
	Grid::AttachmentAccessor<VertexBase, TAttachment>(aa)
{
}

template <class TAttachment>
Grid::VertexAttachmentAccessor<TAttachment>::
VertexAttachmentAccessor(Grid& grid, TAttachment& a) :
	Grid::AttachmentAccessor<VertexBase, TAttachment>(grid, a)
{
}

template <class TAttachment>
Grid::VertexAttachmentAccessor<TAttachment>::
VertexAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach) :
	Grid::AttachmentAccessor<VertexBase, TAttachment>(grid, a, autoAttach)
{
}


////////////////////////////////////////////////////////////////////////
//	EdgeAttachmentAccessor
template <class TAttachment>
Grid::EdgeAttachmentAccessor<TAttachment>::
EdgeAttachmentAccessor() :
	Grid::AttachmentAccessor<EdgeBase, TAttachment>()
{
}

template <class TAttachment>
Grid::EdgeAttachmentAccessor<TAttachment>::
EdgeAttachmentAccessor(const EdgeAttachmentAccessor& aa) :
	Grid::AttachmentAccessor<EdgeBase, TAttachment>(aa)
{
}

template <class TAttachment>
Grid::EdgeAttachmentAccessor<TAttachment>::
EdgeAttachmentAccessor(Grid& grid, TAttachment& a) :
	Grid::AttachmentAccessor<EdgeBase, TAttachment>(grid, a)
{
}

template <class TAttachment>
Grid::EdgeAttachmentAccessor<TAttachment>::
EdgeAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach) :
	Grid::AttachmentAccessor<EdgeBase, TAttachment>(grid, a, autoAttach)
{
}


////////////////////////////////////////////////////////////////////////
//	FaceAttachmentAccessor
template <class TAttachment>
Grid::FaceAttachmentAccessor<TAttachment>::
FaceAttachmentAccessor() :
	Grid::AttachmentAccessor<Face, TAttachment>()
{
}

template <class TAttachment>
Grid::FaceAttachmentAccessor<TAttachment>::
FaceAttachmentAccessor(const FaceAttachmentAccessor& aa) :
	Grid::AttachmentAccessor<Face, TAttachment>(aa)
{
}

template <class TAttachment>
Grid::FaceAttachmentAccessor<TAttachment>::
FaceAttachmentAccessor(Grid& grid, TAttachment& a) :
	Grid::AttachmentAccessor<Face, TAttachment>(grid, a)
{
}

template <class TAttachment>
Grid::FaceAttachmentAccessor<TAttachment>::
FaceAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach) :
	Grid::AttachmentAccessor<Face, TAttachment>(grid, a, autoAttach)
{
}


////////////////////////////////////////////////////////////////////////
//	FaceAttachmentAccessor
template <class TAttachment>
Grid::VolumeAttachmentAccessor<TAttachment>::
VolumeAttachmentAccessor() :
	Grid::AttachmentAccessor<Volume, TAttachment>()
{
}

template <class TAttachment>
Grid::VolumeAttachmentAccessor<TAttachment>::
VolumeAttachmentAccessor(const VolumeAttachmentAccessor& aa) :
	Grid::AttachmentAccessor<Volume, TAttachment>(aa)
{
}

template <class TAttachment>
Grid::VolumeAttachmentAccessor<TAttachment>::
VolumeAttachmentAccessor(Grid& grid, TAttachment& a) :
	Grid::AttachmentAccessor<Volume, TAttachment>(grid, a)
{
}

template <class TAttachment>
Grid::VolumeAttachmentAccessor<TAttachment>::
VolumeAttachmentAccessor(Grid& grid, TAttachment& a, bool autoAttach) :
	Grid::AttachmentAccessor<Volume, TAttachment>(grid, a, autoAttach)
{
}


////////////////////////////////////////////////////////////////////////
//	marks
inline void Grid::mark(GeometricObject* obj)
{
	const int typeID = obj->base_object_id();
	switch(typeID){
		case VERTEX: mark(static_cast<VertexBase*>(obj)); break;
		case EDGE: mark(static_cast<EdgeBase*>(obj)); break;
		case FACE: mark(static_cast<Face*>(obj)); break;
		case VOLUME: mark(static_cast<Volume*>(obj)); break;
	}
}

inline void Grid::mark(VertexBase* obj)
{
	assert(m_bMarking && "ERROR: Grid::mark may only be called between calls to Grid::begin_marking and Grid::end_marking.");
	m_aaMarkVRT[obj] = m_currentMark;
}

inline void Grid::mark(EdgeBase* obj)
{
	assert(m_bMarking && "ERROR: Grid::mark may only be called between calls to Grid::begin_marking and Grid::end_marking.");
	m_aaMarkEDGE[obj] = m_currentMark;
}

inline void Grid::mark(Face* obj)
{
	assert(m_bMarking && "ERROR: Grid::mark may only be called between calls to Grid::begin_marking and Grid::end_marking.");
	m_aaMarkFACE[obj] = m_currentMark;
}

inline void Grid::mark(Volume* obj)
{
	assert(m_bMarking && "ERROR: Grid::mark may only be called between calls to Grid::begin_marking and Grid::end_marking.");
	m_aaMarkVOL[obj] = m_currentMark;
}

inline void Grid::unmark(GeometricObject* obj)
{
	const int typeID = obj->base_object_id();
	switch(typeID){
		case VERTEX: unmark(static_cast<VertexBase*>(obj)); break;
		case EDGE: unmark(static_cast<EdgeBase*>(obj)); break;
		case FACE: unmark(static_cast<Face*>(obj)); break;
		case VOLUME: unmark(static_cast<Volume*>(obj)); break;
	}
}

inline void Grid::unmark(VertexBase* obj)
{
	assert(m_bMarking && "ERROR: Grid::unmark may only be called between calls to Grid::begin_marking and Grid::end_marking.");
	m_aaMarkVRT[obj] = 0;
}

inline void Grid::unmark(EdgeBase* obj)
{
	assert(m_bMarking && "ERROR: Grid::unmark may only be called between calls to Grid::begin_marking and Grid::end_marking.");
	m_aaMarkEDGE[obj] = 0;
}

inline void Grid::unmark(Face* obj)
{
	assert(m_bMarking && "ERROR: Grid::unmark may only be called between calls to Grid::begin_marking and Grid::end_marking.");
	m_aaMarkFACE[obj] = 0;
}

inline void Grid::unmark(Volume* obj)
{
	assert(m_bMarking && "ERROR: Grid::unmark may only be called between calls to Grid::begin_marking and Grid::end_marking.");
	m_aaMarkVOL[obj] = 0;
}

inline bool Grid::is_marked(GeometricObject* obj)
{
	const int typeID = obj->base_object_id();
	switch(typeID){
		case VERTEX: return is_marked(static_cast<VertexBase*>(obj));
		case EDGE: return is_marked(static_cast<EdgeBase*>(obj));
		case FACE: return is_marked(static_cast<Face*>(obj));
		case VOLUME: return is_marked(static_cast<Volume*>(obj));
		default: return false;
	}
}

inline bool Grid::is_marked(VertexBase* obj)
{
	if(m_currentMark == 0)
		return false;
	return m_aaMarkVRT[obj] == m_currentMark;
}

inline bool Grid::is_marked(EdgeBase* obj)
{
	if(m_currentMark == 0)
		return false;
	return m_aaMarkEDGE[obj] == m_currentMark;
}

inline bool Grid::is_marked(Face* obj)
{
	if(m_currentMark == 0)
		return false;
	return m_aaMarkFACE[obj] == m_currentMark;
}

inline bool Grid::is_marked(Volume* obj)
{
	if(m_currentMark == 0)
		return false;
	return m_aaMarkVOL[obj] == m_currentMark;
}

}//	end of namespace libGrid
#endif
