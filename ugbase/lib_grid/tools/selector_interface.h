// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

#ifndef __H__LIBGRID__SELECTOR_INTERFACE__
#define __H__LIBGRID__SELECTOR_INTERFACE__

#include <cassert>
#include "lib_grid/grid/grid.h"

namespace ug
{
///	\addtogroup lib_grid
/// \{

////////////////////////////////////////////////////////////////////////
//	SelectorElements
///	Use these constants to specify which elements shall be supported by a Selector.
/**
 * You may combine the constants using or-operations.
 */
enum SelectorElements
{
	SE_NONE = 0,
	SE_VERTEX = 1,
	SE_EDGE = 1<<1,
	SE_FACE = 1<<2,
	SE_VOLUME = 1 << 3,
	SE_ALL = SE_VERTEX | SE_EDGE | SE_FACE | SE_VOLUME
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	ISelector
///	base-implementation of a selector
/**
 * A selector is a useful class, that allows the user to mark
 * elements of a grid as selected or deselected.
 * The selection status is maintained even if new elements
 * are created or old ones deleted. Features like
 * autoselection and selection_inheritance allow users to
 * follow the creation and removal of elements in all kind of
 * algorithms.
 *
 * Derived classes like Selector and MGSelector or
 * VertexSelector, EdgeSelector, FaceSelector and VolumeSelector
 * keep selected elements in linked lists and allow the user
 * to iterate over them. If a selected element is removed
 * from the underlying grid, it is automatically removed
 * from the element-lists of all associated selectors.
 *
 * Please note that the selector has to be registered at a
 * grid before it may be used. Derived classes may use the
 * protected method set_grid to do so, or pass a grid directly to
 * the constructor.
 *
 * This is a base implementation that can not be instatiated.
 * It features most cruical methods like is_selected,
 * select, deselect etc.
 * 
 * autoselection and selection_inheritance can be enabled here.
 * If autoselection is enabled, then all elements that are created
 * while the selector is registered at a grid will automatically
 * be selected.
 *
 * If selection_inheritance is enabled, then all elements that are
 * created while the selector is registered at a grid, and whose
 * parents (from the same grid) are selected will be automatically
 * selected, too. selection_inheritance is ignored if autoselection
 * is enabled.
 *
 * Support for different element-types can be enabled / disabled
 * by derived classes using set_supported_elements,
 * enable_element_support or disable_element_support.
 *
 * Derived classes have to store the selected elements in linked lists
 * of type ISelector::SectionContainer and should feature begin and
 * end itertors to the sections of those lists.
 *
 * ISelector communicates with derived classes via the methods
 * add_to_list and erase_from_list.
 *
 * Especially for the implementation of erase_from_list the protected
 * member function get_iterator may be useful, which returns the
 * iterator that is associated with the selected element.
 *
 * If a derived class wants to change the selection-status of an
 * element behind the scenes (avoids select, deselect), then it
 * may use the protected methods mark_selected, mark_deselected.
 */
class ISelector : public GridObserver
{
	public:
		ISelector(uint supportedElements = SE_ALL);
		ISelector(Grid& grid, uint supportedElements = SE_ALL);
		virtual ~ISelector();

	//	virtual methods
		virtual void clear() = 0;

	//	selection
		inline void select(GeometricObject* elem);

		template <class TElem>
		inline void select(TElem* elem);

		template <class TIterator>
		inline void select(TIterator iterBegin, TIterator iterEnd);
		
	//	deselection
		inline void deselect(GeometricObject* elem);
		
		template <class TElem>
		inline void deselect(TElem* elem);
		
		template <class TIterator>
		inline void deselect(TIterator iterBegin, TIterator iterEnd);

	//	selection status
		inline bool is_selected(GeometricObject* elem) const;
		inline bool is_selected(VertexBase* vrt) const	{if(!elements_are_supported(SE_VERTEX)) return false; return *m_aaIterVRT[vrt] != NULL;}
		inline bool is_selected(EdgeBase* edge) const	{if(!elements_are_supported(SE_EDGE)) return false; return *m_aaIterEDGE[edge] != NULL;}
		inline bool is_selected(Face* face) const		{if(!elements_are_supported(SE_FACE)) return false; return *m_aaIterFACE[face] != NULL;}
		inline bool is_selected(Volume* vol) const		{if(!elements_are_supported(SE_VOLUME)) return false; return *m_aaIterVOL[vol] != NULL;}

	//	non-virtual methods.
		inline Grid* get_assigned_grid() const		{return m_pGrid;}

	///	returns true if the given element-types are supported.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		inline bool elements_are_supported(uint shElements) const;

	//	if enabled, all new elements will be automatically enabled. Disabled by default.
		void enable_autoselection(bool bEnable);
		inline bool autoselection_enabled()		{return m_bAutoselectionEnabled;}

	//	if enabled, newly created elements derive their selection status from their parents. Enabled by default.
		void enable_selection_inheritance(bool bEnable);
		inline bool selection_inheritance_enabled()		{return m_bSelectionInheritanceEnabled;}

	/**	restricts subset inheritance so that new elements derive their
	 * 	selection status only from parents with the same base-type.
	 * 	Disabled by default.
	 * 	NOTE: strict inheritance only has an effect if
	 * 	selection inheritance is enabled.*/
		void enable_strict_inheritance(bool bEnable);
		inline bool strict_inheritance_enabled()	{return m_bStrictInheritanceEnabled;}
	//	grid callbacks
	/*
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
	*/
		virtual void grid_to_be_destroyed(Grid* grid);
		virtual void elements_to_be_cleared(Grid* grid);

	//	element callbacks
		virtual void vertex_created(Grid* grid, VertexBase* vrt,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

		virtual void edge_created(Grid* grid, EdgeBase* e,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

		virtual void face_created(Grid* grid, Face* f,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

		virtual void volume_created(Grid* grid, Volume* vol,
									GeometricObject* pParent = NULL,
									bool replacesParent = false);

		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt,
										 VertexBase* replacedBy = NULL);

		virtual void edge_to_be_erased(Grid* grid, EdgeBase* e,
										 EdgeBase* replacedBy = NULL);

		virtual void face_to_be_erased(Grid* grid, Face* f,
										 Face* replacedBy = NULL);

		virtual void volume_to_be_erased(Grid* grid, Volume* vol,
										 Volume* replacedBy = NULL);

		virtual void vertices_to_be_merged(Grid* grid, VertexBase* target,
										 VertexBase* elem1, VertexBase* elem2);

		virtual void edges_to_be_merged(Grid* grid, EdgeBase* target,
										 EdgeBase* elem1, EdgeBase* elem2);

		virtual void faces_to_be_merged(Grid* grid, Face* target,
										 Face* elem1, Face* elem2);

		virtual void volumes_to_be_merged(Grid* grid, Volume* target,
										 Volume* elem1, Volume* elem2);

	protected:
		typedef ug::SectionContainer<GeometricObject*, std::list<GeometricObject*> >	SectionContainer;
		typedef SectionContainer::iterator iterator;
		typedef Attachment<iterator> AIterator;

	protected:
		virtual iterator add_to_list(VertexBase* elem) = 0;
		virtual iterator add_to_list(EdgeBase* elem) = 0;
		virtual iterator add_to_list(Face* elem) = 0;
		virtual iterator add_to_list(Volume* elem) = 0;

		virtual void erase_from_list(VertexBase* elem) = 0;
		virtual void erase_from_list(EdgeBase* elem) = 0;
		virtual void erase_from_list(Face* elem) = 0;
		virtual void erase_from_list(Volume* elem) = 0;

	protected:
	///	performs grid registration / deregistration and initialisation of the observer.
	/**	If you call this method with NULL, deregistration and cleanup is performed.
	 *
	 *	Please note: sine set_grid calls virtual methods it shouldn't
	 *	be invoked from any constructors / destructors.*/
		void set_grid(Grid* grid);

	///	set the type of elements that shall be handled by the SubsetHandler.
	/**	Pass an or-combination of constants enumerated in SubsetHandlerElements.
	 *	\sa SubsetHandler::enable_element_support*/
		void set_supported_elements(uint shElements);

	///	enable support for element-types. Does not invalidate previous settings.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		void enable_element_support(uint shElements);

	///	disable support for element-types.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		void disable_element_support(uint shElements);

		inline void mark_selected(VertexBase* elem, iterator iter)	{assert(elements_are_supported(SE_VERTEX)); m_aaIterVRT[elem] = iter;}
		inline void mark_selected(EdgeBase* elem, iterator iter)	{assert(elements_are_supported(SE_EDGE)); m_aaIterEDGE[elem] = iter;}
		inline void mark_selected(Face* elem, iterator iter)		{assert(elements_are_supported(SE_FACE)); m_aaIterFACE[elem] = iter;}
		inline void mark_selected(Volume* elem, iterator iter)		{assert(elements_are_supported(SE_VOLUME)); m_aaIterVOL[elem] = iter;}

		inline void mark_deselected(VertexBase* elem)	{assert(elements_are_supported(SE_VERTEX)); m_aaIterVRT[elem] = m_invalidIterator;}
		inline void mark_deselected(EdgeBase* elem)		{assert(elements_are_supported(SE_EDGE)); m_aaIterEDGE[elem] = m_invalidIterator;}
		inline void mark_deselected(Face* elem)			{assert(elements_are_supported(SE_FACE)); m_aaIterFACE[elem] = m_invalidIterator;}
		inline void mark_deselected(Volume* elem)		{assert(elements_are_supported(SE_VOLUME)); m_aaIterVOL[elem] = m_invalidIterator;}

		inline iterator get_iterator(VertexBase* elem)	{assert(elements_are_supported(SE_VERTEX)); return m_aaIterVRT[elem];}
		inline iterator get_iterator(EdgeBase* elem)	{assert(elements_are_supported(SE_EDGE)); return m_aaIterEDGE[elem];}
		inline iterator get_iterator(Face* elem)		{assert(elements_are_supported(SE_FACE)); return m_aaIterFACE[elem];}
		inline iterator get_iterator(Volume* elem)		{assert(elements_are_supported(SE_VOLUME)); return m_aaIterVOL[elem];}

	///	helper for GridObserver callbacks.
		template <class TElem>
		void elems_to_be_merged(Grid* grid, TElem* target,
								TElem* elem1, TElem* elem2);

	private:
		ISelector(const ISelector& sel){};///<	Copy Constructor not yet implemented!

	protected:
		Grid*	m_pGrid;
		uint	m_supportedElements;
		bool	m_bAutoselectionEnabled;
		bool	m_bSelectionInheritanceEnabled;
		bool	m_bStrictInheritanceEnabled;
		
		AIterator 	m_aIterator;	/// this attachment will be used to store an iterator into m_selectedElements
		Grid::AttachmentAccessor<VertexBase, AIterator>	m_aaIterVRT;
		Grid::AttachmentAccessor<EdgeBase, AIterator>	m_aaIterEDGE;
		Grid::AttachmentAccessor<Face, AIterator>		m_aaIterFACE;
		Grid::AttachmentAccessor<Volume, AIterator>		m_aaIterVOL;
		
		std::list<GeometricObject*>	m_invalidContainer;
		iterator					m_invalidIterator;
};

/// \}

}//	end of namespace

////////////////////////////////
//	include implementation
#include "selector_interface_impl.hpp"

#endif
