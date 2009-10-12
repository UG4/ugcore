// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y08 m11 d19

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	In this header the template class GenericSelector and its specializations
//	VertexSelector, EdgeSelector, FaceSelector and VolumeSelector are defined.
//	Additionally the class Selector is defined, which combines the four
//	specializations.
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIBGRID__SELECTOR__
#define __H__LIBGRID__SELECTOR__

#include <list>
#include "grid/grid.h"
#include "common_attachments.h"
#include "common/util/section_container.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class Selector;

////////////////////////////////////////////////////////////////////////
//	GenericSelector
///	template class for element-selectors.
/**
 * a selector is useful if you want to mark an element as selected
 * or if you want to keep a list of selected elements. The Selector
 * automatically avoids multiple selection.
 *
 * Instead of using this class directly, you will most likely use the typedefs
 * VertexSelector, EdgeSelector, FaceSelector and VolumeSelector instead.
 * By specifying special geometric objects through template arguments it is
 * possible to iterate over a subset of selected elements corresponding to the type.
 * Given a VolumeSelector you could for example iterate only over selected
 * Tetrahedrons by calling yourVolumeSelectorInstance.begin<Tetrahedron>().
 * GenericSelector should only be specialized for the geometric base objects
 * VertexBase, EdgeBase, Face and Volume. Behavior is undefined for other
 * specializations! Some of its methods are thus implemented in selector.cpp
 * and instantiated using explicit instantiation.
 */
template <class TElem>
class GenericSelector : public GridObserver
{
	friend class Selector;
	typedef typename geometry_traits<TElem>::iterator	TElemIterator;

	public:
		GenericSelector();
		GenericSelector(Grid& grid);///	calls assign_grid with the given grid.
		GenericSelector(const GenericSelector<TElem>& gSel);

		virtual ~GenericSelector();

		void assign_grid(Grid& grid);///	registers the observer at the grid.
		Grid* get_assigned_grid();

	//	if enabled, all new elements will be automatically enabled. Disabled by default.
		void enable_autoselection(bool bEnable);
		inline bool autoselection_enabled()		{return m_bAutoselectionEnabled;}

	//	if enabled, newly created elements derive their selection status from their parents. Enabled by default.
		void enable_selection_inheritance(bool bEnable);
		inline bool selection_inheritance_enabled()	{return m_bSelectionInheritanceEnabled;}

		void select(TElem* elem);

		template <class TIterator>
		void select(TIterator iterBegin, TIterator iterEnd);

		template <class TSelElem>
		void select_all();
		inline void select_all()	{return select_all<TElem>();}

		void deselect(TElem* elem);

		template <class TIterator>
		void deselect(TIterator iterBegin, TIterator iterEnd);

		template <class TSelElem>
		void clear_selection();
		inline void clear_selection()	{clear_selection<TElem>();}

		template <class TSelElem>
		inline void clear()			{clear_selection<TSelElem>();}
		inline void clear()			{clear_selection<TElem>();}

		bool is_selected(TElem* elem);
		bool is_selected(GeometricObject* elem);

	//	num_selected
		template <class TSelElem>
		uint num_selected();

		inline uint num_selected()	{return num_selected<TElem>();}

		template <class TSelElem>
		inline uint num()	{return num_selected<TSelElem>();}
		inline uint num()	{return num_selected<TElem>();}

	//	empty
		template <class TSelElem>
		inline bool empty()	{return num_selected<TSelElem>() == 0;}

		inline bool empty()	{return empty<TElem>();}

	//	begin
		template <class TSelElem>
		typename geometry_traits<TSelElem>::iterator begin();

		inline typename geometry_traits<TElem>::iterator
		begin()	{return begin<TElem>();}

	//	end
		template <class TSelElem>
		typename geometry_traits<TSelElem>::iterator end();

		inline typename geometry_traits<TElem>::iterator
		end()	{return end<TElem>();}

	//	grid callbacks
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
		virtual void elements_to_be_cleared(Grid* grid);

	//	vertex callbacks
		virtual void vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent = NULL);
		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt);

	//	edge callbacks
		virtual void edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent = NULL);
		virtual void edge_to_be_erased(Grid* grid, EdgeBase* edge);

	//	face callbacks
		virtual void face_created(Grid* grid, Face* face, GeometricObject* pParent = NULL);
		virtual void face_to_be_erased(Grid* grid, Face* face);

	//	volume callbacks
		virtual void volume_created(Grid* grid, Volume* vol, GeometricObject* pParent = NULL);
		virtual void volume_to_be_erased(Grid* grid, Volume* vol);

	protected:
		void elem_created(Grid* grid, TElem* elem, GeometricObject* pParent);
		void elem_to_be_erased(Grid* grid, TElem* elem);

	protected:
		typedef SectionContainer<TElem*, std::list<TElem*> >			ElemSectionContainer;
		typedef Attachment<typename std::list<TElem*>::iterator>		AElemIterator;

	protected:
		const int		m_baseObjectType;
		bool			m_bAutoselectionEnabled;
		bool			m_bSelectionInheritanceEnabled;
		Grid*			m_pGrid;
		AElemIterator 	m_aElemIterator;	/// this attachment will be used to store an iterator into m_selectedElements
		ElemSectionContainer 	m_selectedElements;	/// holds pointers to selected elements.
		std::list<TElem*> 		m_invalidContainer;	/// used to retrieve an iterator which is used to invalidate other iterators.
		Grid::AttachmentAccessor<TElem, AElemIterator>	m_aaElemIterator;
};


////////////////////////////////////////////////////////////////////////
//	some common selectors
typedef GenericSelector<VertexBase>	VertexSelector;
typedef GenericSelector<EdgeBase>	EdgeSelector;
typedef GenericSelector<Face>		FaceSelector;
typedef GenericSelector<Volume>		VolumeSelector;

////////////////////////////////////////////////////////////////////////
//	Selector
///	allows you to select elements of different types with one selector
/**
 * a selector is useful if you want to mark an element as selected
 * or if you want to keep a list of selected elements. The Selector
 * automatically avoids multiple selection.
 */
class Selector : public GridObserver
{
	public:
		Selector();
		Selector(Grid& grid);///	calls assign_grid with the given grid.
		Selector(const Selector& sel);///<	Copy Constructor not yet implemented!
		virtual ~Selector();

		void assign_grid(Grid& grid);///	registers the observer at the grid.
		Grid* get_assigned_grid();

	//	if enabled, all new elements will be automatically enabled. Disabled by default.
		void enable_autoselection(bool bEnable);
		inline bool autoselection_enabled()		{return m_bAutoselectionEnabled;}

	//	if enabled, newly created elements derive their selection status from their parents. Enabled by default.
		void enable_selection_inheritance(bool bEnable);
		inline bool selection_inheritance_enabled()		{return m_bSelectionInheritanceEnabled;}

		inline void select(VertexBase* vrt)		{m_vertexSelector.select(vrt);}

		inline void select(EdgeBase* edge)		{m_edgeSelector.select(edge);}

		inline void select(Face* face)			{m_faceSelector.select(face);}

		inline void select(Volume* vol)			{m_volumeSelector.select(vol);}

		template <class TIterator>
		inline void select(TIterator iterBegin, TIterator iterEnd)
		{
			typename TIterator::value_type v;
			select(v, iterBegin, iterEnd);
		}

		inline void deselect(VertexBase* vrt)	{m_vertexSelector.deselect(vrt);}
		inline void deselect(VertexBaseIterator vrtsBegin, VertexBaseIterator vrtsEnd)
						{m_vertexSelector.deselect(vrtsBegin, vrtsEnd);}

		inline void deselect(EdgeBase* edge)	{m_edgeSelector.deselect(edge);}
		inline void deselect(EdgeBaseIterator edgesBegin, EdgeBaseIterator edgesEnd)
						{m_edgeSelector.deselect(edgesBegin, edgesEnd);}

		inline void deselect(Face* face)		{m_faceSelector.deselect(face);}
		inline void deselect(FaceIterator facesBegin, FaceIterator facesEnd)
						{m_faceSelector.deselect(facesBegin, facesEnd);}

		inline void deselect(Volume* vol)		{m_volumeSelector.deselect(vol);}
		inline void deselect(VolumeIterator volsBegin, VolumeIterator volsEnd)
						{m_volumeSelector.deselect(volsBegin, volsEnd);}

		template <class TSelElem>
		void select_all();

		inline void select_all()	{select_all<GeometricObject>();}

		template <class TSelElem>
		void clear_selection();

		inline void clear_selection()	{clear_selection<GeometricObject>();}

		template <class TSelElem>
		inline void clear()	{clear_selection<TSelElem>();}
		inline void clear()	{clear_selection<GeometricObject>();}

		inline bool is_selected(VertexBase* vrt)	{return m_vertexSelector.is_selected(vrt);}
		inline bool is_selected(EdgeBase* edge)		{return m_edgeSelector.is_selected(edge);}
		inline bool is_selected(Face* face)			{return m_faceSelector.is_selected(face);}
		inline bool is_selected(Volume* vol)		{return m_volumeSelector.is_selected(vol);}

		template <class TSelElem>
		inline uint num_selected()					{TSelElem* t; return num_selected<TSelElem>(t);}

		template <class TSelElem>
		inline uint num()							{return num_selected<TSelElem>();}

		template <class TSelElem>
		inline bool empty()							{TSelElem* t; return empty<TSelElem>(t);}

		inline bool empty()							{return empty<VertexBase>() && empty<EdgeBase>() && empty<Face>() && empty<Volume>();}

		template <class TSelElem>
		typename geometry_traits<TSelElem>::iterator
		begin();

		template <class TSelElem>
		typename geometry_traits<TSelElem>::iterator
		end();

		inline VertexBaseIterator vertices_begin()	{return m_vertexSelector.begin();}
		inline VertexBaseIterator vertices_end()	{return m_vertexSelector.end();}
		inline EdgeBaseIterator edges_begin()		{return m_edgeSelector.begin();}
		inline EdgeBaseIterator edges_end()			{return m_edgeSelector.end();}
		inline FaceIterator faces_begin()			{return m_faceSelector.begin();}
		inline FaceIterator faces_end()				{return m_faceSelector.end();}
		inline VolumeIterator volumes_begin()		{return m_volumeSelector.begin();}
		inline VolumeIterator volumes_end()			{return m_volumeSelector.end();}

	//	geometric-object-collection
		GeometricObjectCollection get_geometric_object_collection();

	//	grid callbacks
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);

	protected:
		template <class TIterator>
		inline void select(VertexBase* dummy, TIterator iterBegin, TIterator iterEnd)
		{
			m_vertexSelector.select(iterBegin, iterEnd);
		}

		template <class TIterator>
		inline void select(EdgeBase* dummy, TIterator iterBegin, TIterator iterEnd)
		{
			m_edgeSelector.select(iterBegin, iterEnd);
		}

		template <class TIterator>
		inline void select(Face* dummy, TIterator iterBegin, TIterator iterEnd)
		{
			m_faceSelector.select(iterBegin, iterEnd);
		}

		template <class TIterator>
		inline void select(Volume* dummy, TIterator iterBegin, TIterator iterEnd)
		{
			m_volumeSelector.select(iterBegin, iterEnd);
		}

	//	argument is only needed for compiler-selection of the right methods.
	//	pType will not be used in this method an can hold any value.
	//	select all
		template <class TSelElem>
		inline void select_all(GeometricObject* pType)	{m_vertexSelector.select_all<VertexBase>();
														 m_edgeSelector.select_all<EdgeBase>();
														 m_faceSelector.select_all<Face>();
														 m_volumeSelector.select_all<Volume>();}

		template <class TSelElem>
		inline void select_all(VertexBase* pType)		{m_vertexSelector.select_all<TSelElem>();}

		template <class TSelElem>
		inline void select_all(EdgeBase* pType)			{m_edgeSelector.select_all<TSelElem>();}

		template <class TSelElem>
		inline void select_all(Face* pType)				{m_faceSelector.select_all<TSelElem>();}

		template <class TSelElem>
		inline void select_all(Volume* pType)			{m_volumeSelector.select_all<TSelElem>();}

	//	clear selection
		template <class TSelElem>
		inline void clear_selection(GeometricObject* pType)	{m_vertexSelector.clear_selection<VertexBase>();
															 m_edgeSelector.clear_selection<EdgeBase>();
															 m_faceSelector.clear_selection<Face>();
															 m_volumeSelector.clear_selection<Volume>();}

		template <class TSelElem>
		inline void clear_selection(VertexBase* pType)		{m_vertexSelector.clear_selection<TSelElem>();}

		template <class TSelElem>
		inline void clear_selection(EdgeBase* pType)		{m_edgeSelector.clear_selection<TSelElem>();}

		template <class TSelElem>
		inline void clear_selection(Face* pType)			{m_faceSelector.clear_selection<TSelElem>();}

		template <class TSelElem>
		inline void clear_selection(Volume* pType)			{m_volumeSelector.clear_selection<TSelElem>();}

	//	num selected
		template <class TSelElem>
		inline uint num_selected(VertexBase* pType)			{return m_vertexSelector.num_selected<TSelElem>();}

		template <class TSelElem>
		inline uint num_selected(EdgeBase* pType)			{return m_edgeSelector.num_selected<TSelElem>();}

		template <class TSelElem>
		inline uint num_selected(Face* pType)				{return m_faceSelector.num_selected<TSelElem>();}

		template <class TSelElem>
		inline uint num_selected(Volume* pType)				{return m_volumeSelector.num_selected<TSelElem>();}

	//	empty
		template <class TSelElem>
		inline bool empty(VertexBase* pType)				{return m_vertexSelector.empty<TSelElem>();}

		template <class TSelElem>
		inline bool empty(EdgeBase* pType)					{return m_edgeSelector.empty<TSelElem>();}

		template <class TSelElem>
		inline bool empty(Face* pType)						{return m_faceSelector.empty<TSelElem>();}

		template <class TSelElem>
		inline bool empty(Volume* pType)					{return m_volumeSelector.empty<TSelElem>();}

	//	begin
		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(VertexBase* pType)	{return m_vertexSelector.begin<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(EdgeBase* pType)		{return m_edgeSelector.begin<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(Face* pType)			{return m_faceSelector.begin<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(Volume* pType)		{return m_volumeSelector.begin<TSelElem>();}

	//	end
		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(VertexBase* pType)	{return m_vertexSelector.end<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(EdgeBase* pType)	{return m_edgeSelector.end<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(Face* pType)		{return m_faceSelector.end<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(Volume* pType)		{return m_volumeSelector.end<TSelElem>();}


	protected:
		Grid*			m_pGrid;
		bool			m_bAutoselectionEnabled;
		bool			m_bSelectionInheritanceEnabled;
		VertexSelector	m_vertexSelector;
		EdgeSelector	m_edgeSelector;
		FaceSelector	m_faceSelector;
		VolumeSelector	m_volumeSelector;
};

}//	end of namespace

//	include the implementation
#include "selector_impl.hpp"

#endif
