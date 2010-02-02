//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d10

#ifndef __H__LIB_GRID__GRID_IMPLEMENTATION__
#define __H__LIB_GRID__GRID_IMPLEMENTATION__

//#include <cassert>
#include "common/common.h"
#include "common/static_assert.h"
#include "grid_util.h"

namespace ug
{
//	definition of the template-grid functions

////////////////////////////////////////////////////////////////////////
//	create functions
template<class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::create(GeometricObject* pParent)
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION != -1
		&&	geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
		invalid_geometry_type);

	TGeomObj* geomObj = new TGeomObj;
//	int baseObjectType = geometry_traits<GeomObjType>::base_object_type();
//	geomObj->m_elemHandle = m_elementStorage[baseObjectType].m_sectionContainer.insert_element(geomObj, geometry_traits<GeomObjType>::shared_pipe_section());
//	m_elementStorage[baseObjectType].m_attachmentPipe.register_element(geomObj);

	register_element(geomObj, pParent);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>(geomObj->m_entryIter);
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::create(const typename geometry_traits<TGeomObj>::Descriptor& descriptor,
			GeometricObject* pParent)
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION != -1
			&&	geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
			invalid_geometry_type);

	TGeomObj* geomObj = new TGeomObj(descriptor);

//	int baseObjectType = geometry_traits<TGeomObj>::base_object_type();
//	geomObj->m_elemHandle = m_elementStorage[baseObjectType].m_sectionContainer.insert_element(geomObj, geometry_traits<GeomObjType>::shared_pipe_section());
//	m_elementStorage[baseObjectType].m_attachmentPipe.register_element(geomObj);

	register_element(geomObj, pParent);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>(geomObj->m_entryIter);
}

template<class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::create_and_replace(typename geometry_traits<TGeomObj>::geometric_base_object* pReplaceMe)
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION != -1
		&&	geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
		invalid_geometry_type);

	TGeomObj* geomObj = new TGeomObj;

	if(geomObj->reference_object_id() == pReplaceMe->reference_object_id())
	{
		register_and_replace_element(geomObj, pReplaceMe);
		return iterator_cast<typename geometry_traits<TGeomObj>::iterator>(geomObj->m_entryIter);
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

////////////////////////////////////////////////////////////////////////
//	Iterators
template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::begin()
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
		invalid_GeomObj);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(m_elementStorage[geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID].m_sectionContainer.section_begin(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION));
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::end()
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
		invalid_GeomObj);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(m_elementStorage[geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID].m_sectionContainer.section_end(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION));
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::begin() const
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
		invalid_GeomObj);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(m_elementStorage[geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID].m_sectionContainer.section_begin(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION));
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
Grid::end() const
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
		invalid_GeomObj);

	return iterator_cast<typename geometry_traits<TGeomObj>::iterator>
		(m_elementStorage[geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID].m_sectionContainer.section_end(geometry_traits<TGeomObj>::SHARED_PIPE_SECTION));
}


////////////////////////////////////////////////////////////////////////
//	element numbers
template <class TGeomObj>
uint Grid::num() const
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
		invalid_GeomObj);

	int objType = geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID;
	int secIndex = geometry_traits<TGeomObj>::SHARED_PIPE_SECTION;

	if(secIndex == -1)
		return m_elementStorage[objType].m_sectionContainer.num_elements();

	return m_elementStorage[objType].m_sectionContainer.num_elements(secIndex);
}

template <class TGeomObj>
uint Grid::attachment_container_size() const
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
		invalid_GeomObj);

	int objType = geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID;

	return m_elementStorage[objType].m_attachmentPipe.num_data_entries();
}

////////////////////////////////////////////////////////////////////////
//	attachment handling
template <class TGeomObjClass>
void Grid::attach_to(IAttachment& attachment, bool passOnValues)
{
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_TYPE_ID != -1,
			invalid_GeomObjClass);

	int objType = geometry_traits<TGeomObjClass>::BASE_OBJECT_TYPE_ID;
//	setup the options for this attachment.
	int options = 0;
	if(passOnValues)
		options = 1;

	m_elementStorage[objType].m_attachmentPipe.attach(attachment, options);
}

template <class TGeomObjClass, class TAttachment>
void Grid::attach_to_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue)
{
	attach_to_dv<TGeomObjClass, TAttachment>(attachment, defaultValue, attachment.default_pass_on_behaviour());
}

template <class TGeomObjClass, class TAttachment>
void Grid::attach_to_dv(TAttachment& attachment, const typename TAttachment::ValueType& defaultValue, bool passOnValues)
{
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_TYPE_ID != -1,
			invalid_GeomObjClass);

	int objType = geometry_traits<TGeomObjClass>::BASE_OBJECT_TYPE_ID;
//	setup the options for this attachment.
	int options = 0;
	if(passOnValues)
		options = 1;

	m_elementStorage[objType].m_attachmentPipe.attach(attachment, defaultValue, options);
}

template <class TGeomObjClass>
void Grid::detach_from(IAttachment& attachment)
{
	STATIC_ASSERT(geometry_traits<TGeomObjClass>::BASE_OBJECT_TYPE_ID != -1,
				invalid_GeomObjClass);

	int objType = geometry_traits<TGeomObjClass>::BASE_OBJECT_TYPE_ID;
	m_elementStorage[objType].m_attachmentPipe.detach(attachment);
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
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
			invalid_GeomObj);

	return m_elementStorage[geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID].
			m_attachmentPipe.get_data_container(attachment); 
}

template <class TGeomObj>
AttachmentPipe<GeometricObject*, Grid>&
Grid::get_attachment_pipe()
{
	STATIC_ASSERT(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID != -1,
			invalid_GeomObj);

	return m_elementStorage[geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID].m_attachmentPipe;
}

template <class TGeomObj>
uint
Grid::get_attachment_data_index(TGeomObj* pObj) const
{
	return attachment_traits<GeometricObject*, Grid>::get_data_index(this, pObj);
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
//	neighbourhood access
template <class TGeomObj>
EdgeBase* Grid::find_edge_in_associated_edges(TGeomObj* obj,
												EdgeVertices& ev)
{
	EdgeBaseIterator iterEnd = associated_edges_end(obj);
	for(EdgeBaseIterator iter = associated_edges_begin(obj);
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
	unsigned long key = hash_key(&fv);
	FaceIterator iterEnd = associated_faces_end(obj);
	for(FaceIterator iter = associated_faces_begin(obj);
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
	unsigned long key = hash_key(&vv);
	VolumeIterator iterEnd = associated_volumes_end(obj);
	for(VolumeIterator iter = associated_volumes_begin(obj);
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
ug::AttachmentAccessor<GeometricObject*, TAttachment, Grid>()
{
}

template <class TElem, class TAttachment>
Grid::AttachmentAccessor<TElem, TAttachment>::
AttachmentAccessor(const AttachmentAccessor& aa) :
ug::AttachmentAccessor<GeometricObject*, TAttachment, Grid>(aa)
{
}

template <class TElem, class TAttachment>
Grid::AttachmentAccessor<TElem, TAttachment>::
AttachmentAccessor(Grid& grid, TAttachment& a) :
ug::AttachmentAccessor<GeometricObject*, TAttachment, Grid>(grid.get_attachment_pipe<TElem>(), a)
{
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


////////////////////////////////////////////////////////////////////////
//	marks
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
