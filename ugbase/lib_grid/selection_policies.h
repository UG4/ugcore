// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m11 d09

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	The selection-policies help with the creation of element-selectors
//	for different grid types.
//	They handle the storage of selected elements and feature methods
//	to access those elements.
//	The policies normally won't be used by a user of the library.
//	They are used internally in the declaration of specialized
//	selectors.
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIB_GRID__SELECTION_POLICIES__
#define __H__LIB_GRID__SELECTION_POLICIES__

#include <list>
#include "grid/grid.h"
#include "common_attachments.h"
#include "common/util/section_container.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class MultiGrid;
class Selector;
class MGSelector;

////////////////////////////////////////////////////////////////////////
//	GridSelectionPolicy
///	SelectionPolicy for the Grid class.
/**
 * The selection policy for the Grid class implements the select,
 * deselect, clear_selection, is_selected and num_selected methods,
 * that are called by a user of a selector and by the selector itself.
 * Moreover it features begin and end methods, that return iterators
 * to the selected elements.
 *
 * TElem may be either VertexBase, EdgeBase, Face or Volume.
 */
template <class TElem>
class GridSelectionPolicy
{
	friend class Selector;
	
	public:
		typedef typename geometry_traits<TElem>::iterator	TElemIterator;
		typedef Grid* GridPtrType;
		typedef Grid& GridRefType;
		
	public:
		GridSelectionPolicy();
		~GridSelectionPolicy();
		
		void select(TElem* elem);

		void deselect(TElem* elem);

		template <class TSelElem>
		void clear_selection();

		bool is_selected(TElem* elem);
		bool is_selected(GeometricObject* elem);

	//	num_selected
		template <class TSelElem>
		uint num_selected();

	//	begin
		template <class TSelElem>
		typename geometry_traits<TSelElem>::iterator
		begin();

		inline typename geometry_traits<TElem>::iterator
		begin()	{return begin<TElem>();}

	//	end
		template <class TSelElem>
		typename geometry_traits<TSelElem>::iterator
		end();

		inline typename geometry_traits<TElem>::iterator
		end()	{return end<TElem>();}

	protected:
		void new_grid(Grid* gridNew);
		
		inline void init_element(TElem* elem)
		{m_aaElemIterator[elem] = m_invalidContainer.begin();}

	protected:
		typedef SectionContainer<TElem*, std::list<TElem*> >			ElemSectionContainer;
		typedef Attachment<typename std::list<TElem*>::iterator>		AElemIterator;
		
	protected:
		const int				m_baseObjectType;
		Grid*					m_pGrid;
		ElemSectionContainer 	m_selectedElements;	/// holds pointers to selected elements.
		std::list<TElem*> 		m_invalidContainer;	/// used to retrieve an iterator which is used to invalidate other iterators.

		AElemIterator 	m_aElemIterator;	/// this attachment will be used to store an iterator into m_selectedElements
		Grid::AttachmentAccessor<TElem, AElemIterator>	m_aaElemIterator;
};


////////////////////////////////////////////////////////////////////////
//	MultiGridSelectionPolicy
///	SelectionPolicy for the MultiGrid class.
/**
 * The selection policy for the MultiGrid class implements the select,
 * deselect, clear_selection, is_selected and num_selected methods,
 * that are called by a user of a selector and by the selector itself.
 * Moreover it features begin and end methods, that return iterators
 * to the selected elements.
 *
 * TElem may be either VertexBase, EdgeBase, Face or Volume.
 */
template <class TElem>
class MultiGridSelectionPolicy
{
	friend class MGSelector;
	
	public:
		typedef typename geometry_traits<TElem>::iterator	TElemIterator;
		typedef MultiGrid* GridPtrType;
		typedef MultiGrid& GridRefType;
		
	public:
		MultiGridSelectionPolicy();
		~MultiGridSelectionPolicy();
		
	//	selection and deselection
		void select(TElem* elem);

		void deselect(TElem* elem);

	//	clear selection
		template <class TSelElem>
		void clear_selection();

		template <class TSelElem>
		void clear_selection(int level);

		inline void clear_selection(int level)		{clear_selection<TElem>(level);}

	//	is selected
		bool is_selected(TElem* elem);
		bool is_selected(GeometricObject* elem);

	//	num_selected
		template <class TSelElem>
		uint num_selected();

		template <class TSelElem>
		uint num_selected(int level);

		inline uint num_selected(int level)	{return num_selected<TElem>(level);}

	//	begin
		template <class TSelElem>
		typename geometry_traits<TSelElem>::iterator
		begin(int level);

		inline typename geometry_traits<TElem>::iterator
		begin(int level)	{return begin<TElem>(level);}

	//	end
		template <class TSelElem>
		typename geometry_traits<TSelElem>::iterator
		end(int level);

		inline typename geometry_traits<TElem>::iterator
		end(int level)	{return end<TElem>(level);}

	//	num levels
		inline uint num_levels()	{return m_vSections.size();}

	protected:
		typedef SectionContainer<TElem*, std::list<TElem*> >			ElemSectionContainer;
		typedef std::vector<ElemSectionContainer*>						SectionContainerVec;
		typedef Attachment<typename std::list<TElem*>::iterator>		AElemIterator;

	protected:
		void new_grid(MultiGrid* gridNew);
		
		inline void init_element(TElem* elem)
		{m_aaElemIterator[elem] = m_invalidContainer.begin();}

		void grow_sections(int newSize);

		inline ElemSectionContainer& get_section(int level)	{grow_sections(level+1); return *m_vSections[level];}
		
	protected:
		const int				m_baseObjectType;
		MultiGrid*				m_pMultiGrid;
		SectionContainerVec 	m_vSections;	/// holds pointers to section-containers.
		std::list<TElem*> 		m_invalidContainer;	/// used to retrieve an iterator which is used to invalidate other iterators.

		AElemIterator 	m_aElemIterator;	/// this attachment will be used to store an iterator into m_selectedElements
		Grid::AttachmentAccessor<TElem, AElemIterator>	m_aaElemIterator;
};

}//	end of namespace

////////////////////////////////
//	include implementation
#include "selection_policies_impl.hpp"

#endif
