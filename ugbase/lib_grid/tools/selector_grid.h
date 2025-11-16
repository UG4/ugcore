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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	...
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIBGRID__SELECTOR_GRID__
#define __H__LIBGRID__SELECTOR_GRID__

#include <cassert>
#include "selector_interface.h"

namespace ug
{

/** \ingroup lib_grid_tools
 *  \{ */

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	Selector
///	specialization of ISelector for a grid of class Grid.
/** A selector is a useful class, that allows the user to mark
 * elements of a grid as selected or deselected.
 * The selection status is maintained even if new elements
 * are created or old ones deleted. Features like
 * autoselection and selection_inheritance allow users to
 * follow the creation and removal of elements in all kind of
 * algorithms.
 *
 * Please note that the selector has to be registered at a
 * grid before it may be used. You may register it using the constructor
 * or the method assign_grid.
 *
 * This is a specialization of ISelector for the Grid class.
 *
 * The following methods are the most used:
 *	- select, deselect, is_selected (see ISelector)
 *	- begin, end, num, clear.
 *
 * You may specify the element-type on which begin, end, num and clear
 * operate via a template parameter:
 *
 *
 * \code
 * Grid g;
 * Selector sel(g);
 *
 * // ... create elements and select some
 *
 * // number of selected vertices
 * int nSelVrts = sel.num<Vertex>();
 *
 * // number of selected triangles
 * int nSelTris = sel.num<Triangle>();
 *
 * // iteration over all faces
 * for(FaceIterator iter = sel.begin<Face>();
 *		iter != sel.end<Face>(); ++iter){
 * // ...
 * }
 * \endcode
 */
class UG_API Selector : public ISelector
{
	public:
		using BaseClass = ISelector;
		using grid_type = Grid;

	///	The traits class holds some important types for each element-type
		template <class TElem>
		struct traits{
			using iterator = typename geometry_traits<TElem>::iterator;
			using const_iterator = typename geometry_traits<TElem>::const_iterator;
			using level_iterator = typename geometry_traits<TElem>::iterator;
			using const_level_iterator = typename geometry_traits<TElem>::const_iterator;
		};

	public:
		Selector(uint supportedElements = SE_ALL);
		Selector(Grid& grid, uint supportedElements = SE_ALL);
		virtual ~Selector();

		void assign_grid(Grid& grid);
		void assign_grid(Grid* grid);

	///	set the type of elements that shall be handled by the Selector.
	/**	Pass an or-combination of constants enumerated in SelectorElements.
	 *	\sa Selector::enable_element_support*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		void set_supported_elements(uint shElements);

	///	enable support for element-types. Does not invalidate previous settings.
	/**	pass an or-combination of constants enumerated in SelectorElements.*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		void enable_element_support(uint shElements);

	///	disable support for element-types.
	/**	pass an or-combination of constants enumerated in SelectorElements.*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		void disable_element_support(uint shElements);

		virtual void clear();

		template <class TElem>
		inline void clear();

		template <class TElem>
		inline size_t num() const;
		
		inline size_t num() const;
		
	//	empty
		inline bool empty() const;

		template <class TElem>
		inline bool empty() const;

	//	begin
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		begin();

		template <class TElem>
		inline typename geometry_traits<TElem>::const_iterator
		begin() const;
		
	//	end
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		end();
		
		template <class TElem>
		inline typename geometry_traits<TElem>::const_iterator
		end() const;

	//	convenience begin and end
		inline VertexIterator vertices_begin()	{return begin<Vertex>();}
		inline VertexIterator vertices_end()	{return end<Vertex>();}
		inline EdgeIterator edges_begin()		{return begin<Edge>();}
		inline EdgeIterator edges_end()			{return end<Edge>();}
		inline FaceIterator faces_begin()			{return begin<Face>();}
		inline FaceIterator faces_end()				{return end<Face>();}
		inline VolumeIterator volumes_begin()		{return begin<Volume>();}
		inline VolumeIterator volumes_end()			{return end<Volume>();}

	///	returns the first selected element of the given type.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <class TElem> TElem* front();
		
	///	returns the last selected element of the given type.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <class TElem> TElem* back();

	//	geometric-object-collection
		virtual GridObjectCollection get_grid_objects() const;

	//	callbacks that allows us to clean-up
	//	derived from GridObserver
		//virtual void unregistered_from_grid(Grid* grid);
		virtual void grid_to_be_destroyed(Grid* grid);

	////////////////////////////////////////
	//	for compatibility with MGSelector
	///	always returns 1
		inline size_t num_levels() const;
		
	///	calls num();
		inline uint num(size_t) const;
	///	calls num<TElem>();
		template <class TElem> inline size_t num(size_t) const;
		
	//	empty
	///	calls empty();
		inline bool empty(size_t) const;
	//	calls empty<TElem>();
		template <class TElem>
		inline bool empty(size_t) const;

	//	begin
	///	calls begin<TElem>();
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		begin(size_t);

	//	end
	///	calls end<TElem>();
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		end(size_t);

	///	returns true if the selector contains vertices
		virtual bool contains_vertices() const	{return num<Vertex>() > 0;}

	///	returns true if the selector contains edges
		virtual bool contains_edges() const		{return num<Edge>() > 0;}

	///	returns true if the selector contains faces
		virtual bool contains_faces() const		{return num<Face>() > 0;}

	///	returns true if the selector contains volumes
		virtual bool contains_volumes() const	{return num<Volume>() > 0;}

	protected:
		using ISelector::AttachedVertexList;
		using ISelector::AttachedEdgeList;
		using ISelector::AttachedFaceList;
		using ISelector::AttachedVolumeList;

		using ISelector::VertexSectionContainer;
		using ISelector::EdgeSectionContainer;
		using ISelector::FaceSectionContainer;
		using ISelector::VolumeSectionContainer;

	protected:
		void clear_lists();

		virtual void add_to_list(Vertex* elem);
		virtual void add_to_list(Edge* elem);
		virtual void add_to_list(Face* elem);
		virtual void add_to_list(Volume* elem);

		virtual void erase_from_list(Vertex* elem);
		virtual void erase_from_list(Edge* elem);
		virtual void erase_from_list(Face* elem);
		virtual void erase_from_list(Volume* elem);

	///	returns the iterator at which the given element lies in the section container
	/**	This method may only be called if the element is indeed selected
	 * \{
	 */
		inline VertexSectionContainer::iterator
		get_iterator(Vertex* o)
		{
			assert(is_selected(o) && "object not selected.");
			return section_container<Vertex>().get_container().get_iterator(o);
		}

		inline EdgeSectionContainer::iterator
		get_iterator(Edge* o)
		{
			assert(is_selected(o) && "object not selected");
			return section_container<Edge>().get_container().get_iterator(o);
		}

		inline FaceSectionContainer::iterator
		get_iterator(Face* o)
		{
			assert(is_selected(o) && "object not selected");
			return section_container<Face>().get_container().get_iterator(o);
		}

		inline VolumeSectionContainer::iterator
		get_iterator(Volume* o)
		{
			assert(is_selected(o) && "object not selected");
			return section_container<Volume>().get_container().get_iterator(o);
		}
	/**	\}	*/

	///	returns the section container for the given type, subset and level
		template <class TElem> inline
		typename Grid::traits<TElem>::SectionContainer&
		section_container();

	///	returns the const section container for the given type, subset and level
		template <class TElem> inline
		const typename Grid::traits<TElem>::SectionContainer&
		section_container() const;

		template <class TElem>
		inline int get_section_index() const;
		
	private:
		Selector(const Selector& sel){};///<	Copy Constructor not yet implemented!

	protected:
		VertexSectionContainer 	m_vertices;
		EdgeSectionContainer 	m_edges;
		FaceSectionContainer 	m_faces;
		VolumeSectionContainer 	m_volumes;
};

/** \} */
}//	end of namespace

////////////////////////////////
//	include implementation
#include "selector_grid_impl.hpp"

#endif
