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

#ifndef __H__LIBGRID__SELECTOR_INTERFACE__
#define __H__LIBGRID__SELECTOR_INTERFACE__

#include <cassert>
#include "lib_grid/grid/grid.h"
#include "common/types.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/parallel_grid_layout.h"
	#include "pcl/pcl_interface_communicator.h"
#endif


namespace ug
{
/** \ingroup lib_grid_tools
 *  \{ */

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
class UG_API ISelector : public GridObserver
{
	public:
		typedef byte_t	status_t;

		enum{
			DESELECTED = 0,
			SELECTED = 1
		};

	public:
		ISelector(uint supportedElements = SE_ALL);
		ISelector(Grid& grid, uint supportedElements = SE_ALL);
		virtual ~ISelector();

	//	virtual methods
		virtual void clear() = 0;

	//	selection
	///	selects an element
	/**	You may optionally pass a status-flag. Note that 0 is reserved for
	 * non-selected elements. status thus has to be bigger then 0.
	 *
	 * If the element is already selected, then select only alters the status-
	 * flag of the element. It however does not change the elements position
	 * in the list of selected elements.
	 * \{
	 */
		inline void select(GridObject* elem, byte_t status);
		inline void select(GridObject* elem)
		{select(elem, 1);}

		template <class TElem>
		inline void select(TElem* elem, byte_t status);
		template <class TElem>
		inline void select(TElem* elem)
		{select(elem, 1);}

		template <class TIterator>
		inline void select(TIterator iterBegin, TIterator iterEnd, byte_t status = 1);
	/**	\} */

	///	selects an element
	/**	In this context 'mark' is simply a synonym for 'select' and simply forwards
	 * to the corresponding 'select' method.
	 * \{ */
	 	template <class TElem>
		inline void mark(TElem* elem)
		{select(elem);}

		template <class TElem>
		inline void mark(TElem* elem, byte_t status)
		{select(elem, status);}
	/** \} */


	//	deselection
		inline void deselect(GridObject* elem);
		
		template <class TElem>
		inline void deselect(TElem* elem);
		
		template <class TIterator>
		inline void deselect(TIterator iterBegin, TIterator iterEnd);

	///	deselects an element
	/**	In this context 'unmark' is simply a synonym for 'deselect'.*/
	 	template <class TElem>
		inline void unmark(TElem* elem)
		{deselect(elem);}

	//	selection status
	///	returns the selection state of the specified elelent
	/** \{ */
		inline byte_t get_selection_status(GridObject* elem) const;
		inline byte_t get_selection_status(Vertex* vrt) const	{if(!elements_are_supported(SE_VERTEX)) return 0; return m_aaSelVRT[vrt];}
		inline byte_t get_selection_status(Edge* edge) const	{if(!elements_are_supported(SE_EDGE)) return 0; return m_aaSelEDGE[edge];}
		inline byte_t get_selection_status(Face* face) const		{if(!elements_are_supported(SE_FACE)) return 0; return m_aaSelFACE[face];}
		inline byte_t get_selection_status(Volume* vol) const		{if(!elements_are_supported(SE_VOLUME)) return 0; return m_aaSelVOL[vol];}
	/** \} */

	///	returns the selection state of the specified elelent
	/** In this context, 'get_mark' is simply a synonym for 'get_selection_status'
	 * and simply forwards to the corresponding method.*/
		template <class TElem>
		inline byte_t get_mark(TElem* elem) const
		{return get_selection_status(elem);}

	///	returns true if an element is selected
		template <class TElem>
		inline bool is_selected(TElem* elem) const		{return get_selection_status(elem) != 0;}

	//	non-virtual methods.
		inline Grid* grid() const		{return m_pGrid;}

	///	returns a geometric object collection, containing all selected objects
		virtual GridObjectCollection get_grid_objects() const = 0;

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
		virtual void vertex_created(Grid* grid, Vertex* vrt,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void edge_created(Grid* grid, Edge* e,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void face_created(Grid* grid, Face* f,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void volume_created(Grid* grid, Volume* vol,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void vertex_to_be_erased(Grid* grid, Vertex* vrt,
										 Vertex* replacedBy = NULL);

		virtual void edge_to_be_erased(Grid* grid, Edge* e,
										 Edge* replacedBy = NULL);

		virtual void face_to_be_erased(Grid* grid, Face* f,
										 Face* replacedBy = NULL);

		virtual void volume_to_be_erased(Grid* grid, Volume* vol,
										 Volume* replacedBy = NULL);

		virtual void vertices_to_be_merged(Grid* grid, Vertex* target,
										 Vertex* elem1, Vertex* elem2);

		virtual void edges_to_be_merged(Grid* grid, Edge* target,
										 Edge* elem1, Edge* elem2);

		virtual void faces_to_be_merged(Grid* grid, Face* target,
										 Face* elem1, Face* elem2);

		virtual void volumes_to_be_merged(Grid* grid, Volume* target,
										 Volume* elem1, Volume* elem2);

	///	returns true if the selector contains vertices
		virtual bool contains_vertices() const = 0;

	///	returns true if the selector contains edges
		virtual bool contains_edges() const = 0;

	///	returns true if the selector contains faces
		virtual bool contains_faces() const = 0;

	///	returns true if the selector contains volumes
		virtual bool contains_volumes() const = 0;

	///	broadcasts the current selection
	/**	This method is only interesting for parallel algorithms.
	 * By setting 'deselect' to true, all copies of an element will be deselected if
	 * one or more copies already were deselected.
	 *
	 * If 'deselect' is false (this is the default), the selection-states of all
	 * copies of an element will be united and applied to all copies of that element.
	 *
	 * Set 'includeGhosts' to true, to also broadcast the selection-state to
	 * ghost-copies. Default is false.
	 */
		virtual void broadcast_selection_states(bool deselect = false,
												bool includeGhosts = false);

	protected:
		typedef Grid::traits<Vertex>::AttachedElementList	AttachedVertexList;
		typedef Grid::traits<Edge>::AttachedElementList		AttachedEdgeList;
		typedef Grid::traits<Face>::AttachedElementList			AttachedFaceList;
		typedef Grid::traits<Volume>::AttachedElementList		AttachedVolumeList;

		typedef Grid::traits<Vertex>::SectionContainer		VertexSectionContainer;
		typedef Grid::traits<Edge>::SectionContainer		EdgeSectionContainer;
		typedef Grid::traits<Face>::SectionContainer			FaceSectionContainer;
		typedef Grid::traits<Volume>::SectionContainer			VolumeSectionContainer;

	protected:
		virtual void add_to_list(Vertex* elem) = 0;
		virtual void add_to_list(Edge* elem) = 0;
		virtual void add_to_list(Face* elem) = 0;
		virtual void add_to_list(Volume* elem) = 0;

		virtual void erase_from_list(Vertex* elem) = 0;
		virtual void erase_from_list(Edge* elem) = 0;
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
	//	Protected non-virtual to avoid virtual calls during construction
		void set_supported_elements(uint shElements);

	///	enable support for element-types. Does not invalidate previous settings.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
	//	Protected non-virtual to avoid virtual calls during construction
		void enable_element_support(uint shElements);

	///	disable support for element-types.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
	//	Protected non-virtual to avoid virtual calls during construction
		void disable_element_support(uint shElements);

		inline void mark_selected(Vertex* elem, byte_t status)	{assert(elements_are_supported(SE_VERTEX)); m_aaSelVRT[elem] = status;}
		inline void mark_selected(Edge* elem, byte_t status)		{assert(elements_are_supported(SE_EDGE)); m_aaSelEDGE[elem] = status;}
		inline void mark_selected(Face* elem, byte_t status)			{assert(elements_are_supported(SE_FACE)); m_aaSelFACE[elem] = status;}
		inline void mark_selected(Volume* elem, byte_t status)		{assert(elements_are_supported(SE_VOLUME)); m_aaSelVOL[elem] = status;}

		inline void mark_deselected(Vertex* elem)	{assert(elements_are_supported(SE_VERTEX)); m_aaSelVRT[elem] = 0;}
		inline void mark_deselected(Edge* elem)		{assert(elements_are_supported(SE_EDGE)); m_aaSelEDGE[elem] = 0;}
		inline void mark_deselected(Face* elem)			{assert(elements_are_supported(SE_FACE)); m_aaSelFACE[elem] = 0;}
		inline void mark_deselected(Volume* elem)		{assert(elements_are_supported(SE_VOLUME)); m_aaSelVOL[elem] = 0;}

	///	helper for GridObserver callbacks.
		template <class TElem>
		void elems_to_be_merged(Grid* grid, TElem* target,
								TElem* elem1, TElem* elem2);

	private:
		ISelector(const ISelector& sel){};///<	Copy Constructor not yet implemented!

		#ifdef UG_PARALLEL
			template <class TIntfcCom>
			void broadcast_selection_states(bool deselect, bool includeGhosts,
											TIntfcCom& icom);
		#endif

	protected:
		Grid*	m_pGrid;
		uint	m_supportedElements;
		bool	m_bAutoselectionEnabled;
		bool	m_bSelectionInheritanceEnabled;
		bool	m_bStrictInheritanceEnabled;
		
	//	will use a default constructor
		typedef Attachment<unsigned char> 	AUChar;
		AUChar								m_aSelected;

		Grid::AttachmentAccessor<Vertex, AUChar>	m_aaSelVRT;
		Grid::AttachmentAccessor<Edge, AUChar>		m_aaSelEDGE;
		Grid::AttachmentAccessor<Face, AUChar>			m_aaSelFACE;
		Grid::AttachmentAccessor<Volume, AUChar>		m_aaSelVOL;

		#ifdef UG_PARALLEL
			pcl::InterfaceCommunicator<VertexLayout>	m_icomVRT;
			pcl::InterfaceCommunicator<EdgeLayout>		m_icomEDGE;
			pcl::InterfaceCommunicator<FaceLayout>		m_icomFACE;
			pcl::InterfaceCommunicator<VolumeLayout>	m_icomVOL;
		#endif
};

/// \}

}//	end of namespace

////////////////////////////////
//	include implementation
#include "selector_interface_impl.hpp"

#endif
