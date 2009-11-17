// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y08 m11 d19

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	In this header the template class GenericElementSelector and its specializations
//	VertexSelector, EdgeSelector, FaceSelector and VolumeSelector are defined.
//	Additionally the class Selector is defined, which combines the four
//	specializations.
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIBGRID__SELECTOR__
#define __H__LIBGRID__SELECTOR__

#include "grid/grid.h"
#include "selection_policies.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class Selector;

////////////////////////////////////////////////////////////////////////
//	GenericElementSelector
///	template class for element-selectors.
/**
 * a selector is useful if you want to mark an element as selected
 * or if you want to keep a list of selected elements. The Selector
 * automatically avoids multiple selection.
 *
 * Instead of using this class directly, you will most likely use the typedefs
 * VertexSelector, EdgeSelector, FaceSelector or VolumeSelector instead.
 * By specifying special geometric objects through template arguments it is
 * possible to iterate over a subset of selected elements corresponding to the type.
 * Given a VolumeSelector you could for example iterate only over selected
 * Tetrahedrons by calling yourVolumeSelectorInstance.begin<Tetrahedron>().
 * GenericElementSelector should only be specialized for the geometric base objects
 * VertexBase, EdgeBase, Face and Volume. Behavior is undefined for other
 * specializations! Some of its methods are thus implemented in selector.cpp
 * and instantiated using explicit instantiation.
 *
 * Please note that some of the functionality of the selector is implemented
 * in the SelectionPolicy, that is passed as a template parameter.
 */
template <class TElem, class SelectionPolicy>
class GenericElementSelector : public GridObserver, public SelectionPolicy
{
	friend class Selector;
	public:
		typedef typename geometry_traits<TElem>::iterator	TElemIterator;
		typedef typename SelectionPolicy::GridRefType			TGridRef;
		typedef typename SelectionPolicy::GridPtrType			TGridPtr;
		typedef typename SelectionPolicy::GridRefType			GridRefType;
		typedef typename SelectionPolicy::GridPtrType			GridPtrType;

	
	public:
	//	We're using some methods of SelectionPolicy.
	//	To avoid that a method from this class hides a method
	//	from the SelectionPolicy, we'll explicitly connect them here.
		using SelectionPolicy::select;
		using SelectionPolicy::deselect;
		using SelectionPolicy::clear_selection;
		using SelectionPolicy::is_selected;
		using SelectionPolicy::num_selected;

	//	begin() and end() are not mentioned here, since they have differing
	//	function arguments for different grid types.

	public:
		GenericElementSelector();
		GenericElementSelector(TGridRef grid);///<	calls assign_grid with the given grid.
		GenericElementSelector(const GenericElementSelector<TElem, SelectionPolicy>& gSel);

		virtual ~GenericElementSelector();

	///	registers the observer at the grid.
		void assign_grid(TGridRef grid);
		TGridPtr get_assigned_grid();

	///	if enabled, all new elements will be automatically enabled. Disabled by default.
		void enable_autoselection(bool bEnable);
		inline bool autoselection_enabled()		{return m_bAutoselectionEnabled;}

	///	if enabled, newly created elements derive their selection status from their parents. Enabled by default.
		void enable_selection_inheritance(bool bEnable);
		inline bool selection_inheritance_enabled()	{return m_bSelectionInheritanceEnabled;}
		
	///	calls select(TElem*) on all elements between iterBegin and iterEnd.
		template <class TIterator>
		void select(TIterator iterBegin, TIterator iterEnd);

	///	calls select(TElem*) on all elements of type TSelElem of the grid
		template <class TSelElem>
		void select_all();
		
	//	calls select(TElem*) on all elements
		inline void select_all()	{return select_all<TElem>();}
		
		template <class TIterator>
		void deselect(TIterator iterBegin, TIterator iterEnd);
		
		inline void clear_selection()	{SelectionPolicy::template clear_selection<TElem>();}

		template <class TSelElem>
		inline void clear()			{SelectionPolicy::template clear_selection<TSelElem>();}
		inline void clear()			{SelectionPolicy::template clear_selection<TElem>();}
		
		inline uint num_selected()	{return SelectionPolicy::template num_selected<TElem>();}

		template <class TSelElem>
		inline uint num()	{return SelectionPolicy::template num_selected<TSelElem>();}
		inline uint num()	{return SelectionPolicy::template num_selected<TElem>();}
		
	//	empty
		template <class TSelElem>
		inline bool empty()	{return SelectionPolicy::template num_selected<TSelElem>() == 0;}

		inline bool empty()	{return empty<TElem>();}
		
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
		bool			m_bAutoselectionEnabled;
		bool			m_bSelectionInheritanceEnabled;
		TGridPtr		m_pGrid;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GenericSelector
///	base-implementation of a selector
/**
 * ...
 */
template <class TElementSelectors>
class GenericSelector : public GridObserver
{
	public:
		typedef GenericSelector<TElementSelectors>		ClassType;
		typedef typename TElementSelectors::TGridRef	TGridRef;
		typedef typename TElementSelectors::TGridPtr	TGridPtr;
		typedef typename TElementSelectors::TVertexSelector	TVertexSelector;
		typedef typename TElementSelectors::TEdgeSelector	TEdgeSelector;
		typedef typename TElementSelectors::TFaceSelector	TFaceSelector;
		typedef typename TElementSelectors::TVolumeSelector	TVolumeSelector;

	public:
		GenericSelector();
		GenericSelector(TGridRef grid);///	calls assign_grid with the given grid.
		GenericSelector(const ClassType& sel);///<	Copy Constructor not yet implemented!
		virtual ~GenericSelector();

		void assign_grid(TGridRef grid);///	registers the observer at the grid.
		TGridPtr get_assigned_grid();

	//	if enabled, all new elements will be automatically enabled. Disabled by default.
		void enable_autoselection(bool bEnable);
		inline bool autoselection_enabled()		{return m_bAutoselectionEnabled;}

	//	if enabled, newly created elements derive their selection status from their parents. Enabled by default.
		void enable_selection_inheritance(bool bEnable);
		inline bool selection_inheritance_enabled()		{return m_bSelectionInheritanceEnabled;}

		void select(GeometricObject* obj);
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

		void deselect(GeometricObject* obj);
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

		bool is_selected(GeometricObject* obj);
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

	protected:
		TGridPtr		m_pGrid;
		bool			m_bAutoselectionEnabled;
		bool			m_bSelectionInheritanceEnabled;
		TVertexSelector	m_vertexSelector;
		TEdgeSelector	m_edgeSelector;
		TFaceSelector	m_faceSelector;
		TVolumeSelector	m_volumeSelector;
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	some common element selectors for Grid and MultiGrid
typedef GenericElementSelector<VertexBase, GridSelectionPolicy<VertexBase> >
		VertexSelector;
typedef GenericElementSelector<EdgeBase, GridSelectionPolicy<EdgeBase> >
		EdgeSelector;
typedef GenericElementSelector<Face, GridSelectionPolicy<Face> >
		FaceSelector;
typedef GenericElementSelector<Volume, GridSelectionPolicy<Volume> >
		VolumeSelector;

struct ElementSelectors
{
	typedef Grid&			TGridRef;
	typedef Grid*			TGridPtr;
	typedef VertexSelector	TVertexSelector;
	typedef EdgeSelector	TEdgeSelector;
	typedef FaceSelector	TFaceSelector;
	typedef VolumeSelector	TVolumeSelector;
};

typedef GenericElementSelector<VertexBase, MultiGridSelectionPolicy<VertexBase> >
		MGVertexSelector;
typedef GenericElementSelector<EdgeBase, MultiGridSelectionPolicy<EdgeBase> >
		MGEdgeSelector;
typedef GenericElementSelector<Face, MultiGridSelectionPolicy<Face> >
		MGFaceSelector;
typedef GenericElementSelector<Volume, MultiGridSelectionPolicy<Volume> >
		MGVolumeSelector;

struct MGElementSelectors
{
	typedef MultiGrid&			TGridRef;
	typedef MultiGrid*			TGridPtr;
	typedef MGVertexSelector	TVertexSelector;
	typedef MGEdgeSelector		TEdgeSelector;
	typedef MGFaceSelector		TFaceSelector;
	typedef MGVolumeSelector	TVolumeSelector;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	Selector and MGSelector

////////////////////////////////////////////////////////////////////////
//	Selector
///	allows to select elements of different types with one selector
/**
 * a selector is useful if you want to mark an element as selected
 * or if you want to keep a list of selected elements. The Selector
 * automatically avoids multiple selection.
 */
class Selector : public GenericSelector<ElementSelectors>
{
	public:
		typedef GenericSelector<ElementSelectors>	BaseClass;
		typedef BaseClass::TGridRef					TGridRef;

	public:
		Selector() : BaseClass()	{}
		Selector(TGridRef grid) : BaseClass(grid)		{}
		Selector(const Selector& sel) : BaseClass(sel)	{}

	//	begin
		template <class TSelElem>
		inline
		typename geometry_traits<TSelElem>::iterator
		begin()
		{
			TSelElem* pTmp;
			return begin<TSelElem>(pTmp);
		}

	//	end
		template <class TSelElem>
		inline
		typename geometry_traits<TSelElem>::iterator
		end()
		{
			TSelElem* pTmp;
			return end<TSelElem>(pTmp);
		}

	//	convenience begin and end
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

	protected:
	//	begin
		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(const VertexBase* pType)	{return m_vertexSelector.begin<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(const EdgeBase* pType)		{return m_edgeSelector.begin<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(const Face* pType)			{return m_faceSelector.begin<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(const Volume* pType)		{return m_volumeSelector.begin<TSelElem>();}

	//	end
		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(const VertexBase* pType)	{return m_vertexSelector.end<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(const EdgeBase* pType)	{return m_edgeSelector.end<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(const Face* pType)		{return m_faceSelector.end<TSelElem>();}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(const Volume* pType)		{return m_volumeSelector.end<TSelElem>();}
};


////////////////////////////////////////////////////////////////////////
//	MGSelector
///	allows to select elements of different types with one selector
/**
 * a selector is useful if you want to mark an element as selected
 * or if you want to keep a list of selected elements. The Selector
 * automatically avoids multiple selection.
 */
class MGSelector : public GenericSelector<MGElementSelectors>
{
	public:
		typedef GenericSelector<MGElementSelectors>	BaseClass;
		typedef BaseClass::TGridRef					TGridRef;

	public:
		using BaseClass::num_selected;
		using BaseClass::num;

	public:
		MGSelector() : BaseClass()	{}
		MGSelector(TGridRef grid) : BaseClass(grid)		{}
		MGSelector(const MGSelector& sel) : BaseClass(sel)	{}

	//	begin
		template <class TSelElem>
		inline
		typename geometry_traits<TSelElem>::iterator
		begin(int level)									{TSelElem* pTmp; return begin<TSelElem>(pTmp, level);}

	//	end
		template <class TSelElem>
		inline
		typename geometry_traits<TSelElem>::iterator
		end(int level)										{TSelElem* pTmp; return end<TSelElem>(pTmp, level);}

	//	convenience begin and end
		inline VertexBaseIterator vertices_begin(int level)	{return m_vertexSelector.begin(level);}
		inline VertexBaseIterator vertices_end(int level)	{return m_vertexSelector.end(level);}
		inline EdgeBaseIterator edges_begin(int level)		{return m_edgeSelector.begin(level);}
		inline EdgeBaseIterator edges_end(int level)		{return m_edgeSelector.end(level);}
		inline FaceIterator faces_begin(int level)			{return m_faceSelector.begin(level);}
		inline FaceIterator faces_end(int level)			{return m_faceSelector.end(level);}
		inline VolumeIterator volumes_begin(int level)		{return m_volumeSelector.begin(level);}
		inline VolumeIterator volumes_end(int level)		{return m_volumeSelector.end(level);}

	//	level-support for num-queries
		template <class TSelElem>
		inline uint num_selected(int level)					{TSelElem* t; return num_selected<TSelElem>(level, t);}

		template <class TSelElem>
		inline uint num(int level)							{return num_selected<TSelElem>(level);}

	//	num levels
		inline uint num_levels()	{return std::max(std::max(m_vertexSelector.num_levels(), m_edgeSelector.num_levels()),
													std::max(m_faceSelector.num_levels(), m_volumeSelector.num_levels()));}
		
	//	multi-level-geometric-object-collection
		MultiLevelGeometricObjectCollection
		get_multi_level_geometric_object_collection();

	protected:
	//	begin
		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(int level, const VertexBase* pType)	{return m_vertexSelector.begin<TSelElem>(level);}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(int level, const EdgeBase* pType)		{return m_edgeSelector.begin<TSelElem>(level);}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(int level, const Face* pType)			{return m_faceSelector.begin<TSelElem>(level);}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		begin(int level, const Volume* pType)		{return m_volumeSelector.begin<TSelElem>(level);}

	//	end
		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(int level, const VertexBase* pType)	{return m_vertexSelector.end<TSelElem>(level);}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(int level, const EdgeBase* pType)	{return m_edgeSelector.end<TSelElem>(level);}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(int level, const Face* pType)		{return m_faceSelector.end<TSelElem>(level);}

		template <class TSelElem>
		inline typename geometry_traits<TSelElem>::iterator
		end(int level, const Volume* pType)		{return m_volumeSelector.end<TSelElem>(level);}

	//	num selected
		template <class TSelElem>
		inline uint num_selected(int level, const VertexBase* pType)	{return m_vertexSelector.num_selected<TSelElem>(level);}

		template <class TSelElem>
		inline uint num_selected(int level, const EdgeBase* pType)		{return m_edgeSelector.num_selected<TSelElem>(level);}

		template <class TSelElem>
		inline uint num_selected(int level, const Face* pType)			{return m_faceSelector.num_selected<TSelElem>(level);}

		template <class TSelElem>
		inline uint num_selected(int level, const Volume* pType)		{return m_volumeSelector.num_selected<TSelElem>(level);}
};

}//	end of namespace

////////////////////////////////
//	include implementation
#include "selector_impl.hpp"

#endif
