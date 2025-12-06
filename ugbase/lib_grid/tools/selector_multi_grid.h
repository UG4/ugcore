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

#ifndef __H__LIBGRID__SELECTOR_MULTI_GRID__
#define __H__LIBGRID__SELECTOR_MULTI_GRID__

#include <cassert>
#include "selector_interface.h"
#include "../multi_grid.h"
#include "../lib_grid_messages.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class MultiGrid;

/** \ingroup lib_grid_tools
 *  \{ */

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	MGSelector
///	specialization of ISelector for grids of class MultiGrid.
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
 * This is a specialization of ISelector for the MultiGrid class.
 *
 * The following methods are the most used:
 *	- select, deselect, is_selected (see ISelector)
 *	- begin, end, num, clear.
 *
 * Note that the number of levels in the MGSelector always matches the number
 * of levels in the associated multigrid. This is guaranteed through a callback
 * mechanism.
 *
 * You may specify the element-type on which begin, end, num and clear
 * operate via a template parameter, and the level via a
 * normal int-parameter:
 *
 * \code
 * MultiGrid mg;
 * MGSelector sel(mg);
 *
 * // ... create elements and select some
 *
 * // number of selected vertices on level 1
 * int nSelVrts = sel.num<Vertex>(1);
 *
 * // total number of selected triangles
 * int nSelTris = sel.num<Triangle>();
 *
 * // iteration over all faces
 * for(uint i = 0; i < sel.num_levels(); ++i){
 *	for(FaceIterator iter = sel.begin<Face>(i);
 *		iter != sel.end<Face>(i); ++iter){
 * // ...
 *	}
 * }
 * \endcode
 */

class UG_API MGSelector : public ISelector
{
	public:
		using BaseClass = ISelector;
		using grid_type = MultiGrid;
		
	/** This iterator is used by MGSelector to provide iteration across all levels.
	 * The TMGSelector and TLevelIterator template argument allows to use this
	 * iterator for const and non-const use.*/
		template <typename TElem, typename TMGSelector, typename TLevelIterator>
		class MGSelectionIterator
		{
			public:
				using this_type = MGSelectionIterator;
				using iterator_category = std::forward_iterator_tag;
				using difference_type = size_t;
				using pointer = TElem**;
				using value_type = TElem*;
				using reference = value_type&;

				MGSelectionIterator() : m_sel(nullptr), m_lvl(0)	{}

			///	copy constructor that allows creation of const-iterators from non-const iterators
				/*explicit*/ MGSelectionIterator(const MGSelectionIterator<TElem, MGSelector,
							        		typename geometry_traits<TElem>::iterator>& iter)
				{
					m_sel = iter.m_sel;
					m_lvl = iter.m_lvl;
					m_iter = iter.m_iter;
				}

				this_type operator ++ ()	{increment(); return *this;}
				this_type operator ++ (int unused)	{this_type i = *this; increment(); return i;}

				bool operator == (const this_type& iter) const {return equal(iter);}
				bool operator != (const this_type& iter) const {return !equal(iter);}

				value_type operator * () {return dereference();}
				value_type operator * () const {return dereference();}

			private:
				friend class MGSelector;
		//		friend class MGSelectionIterator<TElem, const MGSelector,
		//										 typename MGSelector::traits<TElem>::const_level_iterator>;
				using level_iterator = TLevelIterator;

				MGSelectionIterator(TMGSelector* sel, int lvl,
									level_iterator iter)
				{
					m_sel = sel;
					m_lvl = lvl;
					m_iter = iter;
				}

				inline bool equal(const this_type& other) const
				{
					return m_iter == other.m_iter;
				}

			///	returns next valid iterator
				void increment()
				{
					++m_iter;
					while((m_iter == m_sel->template end<TElem>(m_lvl))
						  && (m_lvl + 1 < m_sel->num_levels()))
					{
						++m_lvl;
						m_iter = m_sel->template begin<TElem>(m_lvl);
					}
				}

			///	dereference
				inline value_type dereference() const
				{
					return *m_iter;
				}

			private:
				TMGSelector*	m_sel;
				size_t			m_lvl;
				level_iterator	m_iter;
		};

	///	The traits class holds some important types for each element-type
		template <typename TElem>
		struct traits{
			using value_t = TElem*;
			using level_iterator = typename geometry_traits<TElem>::iterator;
			using const_level_iterator = typename geometry_traits<TElem>::const_iterator;
			using iterator = MGSelectionIterator<TElem, MGSelector, level_iterator>;
			using const_iterator = MGSelectionIterator<TElem, const MGSelector, const_level_iterator>;
		};


		explicit MGSelector(byte_t supportedElements = SE_ALL);
		explicit MGSelector(MultiGrid& grid, byte_t supportedElements = SE_ALL);

		~MGSelector() override;

		void assign_grid(MultiGrid& grid);
		void assign_grid(MultiGrid* grid);
		[[nodiscard]] inline MultiGrid* multi_grid() const {return m_pMultiGrid;}

	///	set the type of elements that shall be handled by the Selector.
	/**	Pass an or-combination of constants enumerated in SelectorElements.
	 *	\sa Selector::enable_element_support*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		inline void set_supported_elements(byte_t shElements);

	///	enable support for element-types. Does not invalidate previous settings.
	/**	pass an or-combination of constants enumerated in SelectorElements.*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		inline void enable_element_support(byte_t shElements);

	///	disable support for element-types.
	/**	pass an or-combination of constants enumerated in SelectorElements.*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		void disable_element_support(byte_t shElements);

		[[nodiscard]] inline size_t num_levels() const	{return m_levels.size();}
		[[nodiscard]] inline size_t top_level() const
		{
			size_t l = m_levels.size();
			if(l == 0)
				return 0;
			else
				return l - 1;
		}

		void clear() override;

		template <typename TElem>
		inline void clear();

		void clear(int level);

		template <typename TElem>
		inline void clear(int level);

		template <typename TElem>
		[[nodiscard]] inline size_t num(int level) const;
		
		[[nodiscard]] inline size_t num(int level) const;

		template <typename TElem>
		[[nodiscard]] inline size_t num() const;
		
		[[nodiscard]] inline size_t num() const;

	//	empty
		[[nodiscard]] inline bool empty(int level) const;

		template <typename TElem>
		[[nodiscard]] inline bool empty(int level) const;

		[[nodiscard]] inline bool empty() const;

		template <typename TElem>
		[[nodiscard]] inline bool empty() const;

	//	begin
		template <typename TElem>
		inline typename traits<TElem>::iterator
		begin();

		template <typename TElem>
		[[nodiscard]] inline typename traits<TElem>::const_iterator
		begin() const;

		template <typename TElem>
		inline typename traits<TElem>::level_iterator
		begin(int level);

		template <typename TElem>
		[[nodiscard]] inline typename traits<TElem>::const_level_iterator
		begin(int level) const;
		
	//	end
		template <typename TElem>
		inline typename traits<TElem>::iterator
		end();

		template <typename TElem>
		inline typename traits<TElem>::const_iterator
		end() const;

		template <typename TElem>
		inline typename traits<TElem>::level_iterator
		end(int level);
		
		template <typename TElem>
		[[nodiscard]] inline typename traits<TElem>::const_level_iterator
		end(int level) const;

	//	convenience begin and end
		inline traits<Vertex>::level_iterator vertices_begin(int level) {return begin<Vertex>(level);}
		inline traits<Vertex>::level_iterator vertices_end(int level)	{return end<Vertex>(level);}
		inline traits<Edge>::level_iterator edges_begin(int level)		{return begin<Edge>(level);}
		inline traits<Edge>::level_iterator edges_end(int level)		{return end<Edge>(level);}
		inline traits<Face>::level_iterator faces_begin(int level)		{return begin<Face>(level);}
		inline traits<Face>::level_iterator faces_end(int level)		{return end<Face>(level);}
		inline traits<Volume>::level_iterator volumes_begin(int level)	{return begin<Volume>(level);}
		inline traits<Volume>::level_iterator volumes_end(int level)	{return end<Volume>(level);}

	///	returns the first selected element of the given type on the specified level.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <typename TElem> TElem* front(int level);
		
	///	returns the last selected element of the given type on the specified level.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <typename TElem> TElem* back(int level);
		
	//	geometric-object-collection
		[[nodiscard]] GridObjectCollection get_grid_objects() const override;

	//	callbacks that allows us to clean-up
	//	derived from GridObserver
	/*
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
	*/
		void grid_to_be_destroyed(Grid* grid) override;
		
	///	returns true if the selector contains vertices
		[[nodiscard]] bool contains_vertices() const override {return num<Vertex>() > 0;}

	///	returns true if the selector contains edges
		[[nodiscard]] bool contains_edges() const override {return num<Edge>() > 0;}

	///	returns true if the selector contains faces
		[[nodiscard]] bool contains_faces() const override {return num<Face>() > 0;}

	///	returns true if the selector contains volumes
		[[nodiscard]] bool contains_volumes() const override {return num<Volume>() > 0;}

	protected:
		void clear_lists();

		void add_to_list(Vertex* elem) override;
		void add_to_list(Edge* elem) override;
		void add_to_list(Face* elem) override;
		void add_to_list(Volume* elem) override;

		void erase_from_list(Vertex* elem) override;
		void erase_from_list(Edge* elem) override;
		void erase_from_list(Face* elem) override;
		void erase_from_list(Volume* elem) override;

	protected:
		using ISelector::AttachedVertexList;
		using ISelector::AttachedEdgeList;
		using ISelector::AttachedFaceList;
		using ISelector::AttachedVolumeList;

		using ISelector::VertexSectionContainer;
		using ISelector::EdgeSectionContainer;
		using ISelector::FaceSectionContainer;
		using ISelector::VolumeSectionContainer;

		struct Level{
			VertexSectionContainer	m_vertices;
			EdgeSectionContainer	m_edges;
			FaceSectionContainer	m_faces;
			VolumeSectionContainer	m_volumes;
		};

		using LevelVec = std::vector<Level*>;

	protected:
	///	returns the section container for the given type, subset and level
		template <typename TElem> inline
		typename Grid::traits<TElem>::SectionContainer&
		section_container(int level);

	///	returns the const section container for the given type, subset and level
		template <typename TElem> inline
		const typename Grid::traits<TElem>::SectionContainer&
		section_container(int level) const;
		
		template <typename TElem>
		[[nodiscard]] inline int get_section_index() const;

		inline void level_required(int newSize);
		void add_level();

	///	This method should only be called if a complete cleanup is required.
		void cleanup();

	///	returns the iterator at which the given element lies in the section container
	/**	This method may only be called if the element is indeed selected
	 * \{
	 */
		inline VertexSectionContainer::iterator
		get_level_iterator(Vertex* o)
		{
			assert(is_selected(o) && "object not selected.");
			return section_container<Vertex>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline EdgeSectionContainer::iterator
		get_level_iterator(Edge* o)
		{
			assert(is_selected(o) && "object not selected");
			return section_container<Edge>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline FaceSectionContainer::iterator
		get_level_iterator(Face* o)
		{
			assert(is_selected(o) && "object not selected");
			return section_container<Face>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline VolumeSectionContainer::iterator
		get_level_iterator(Volume* o)
		{
			assert(is_selected(o) && "object not selected");
			return section_container<Volume>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}
	/**	\}	*/

	///	callback for multigrid messages
		void multigrid_changed(const GridMessage_MultiGridChanged& gm);

	private:
		MGSelector(const MGSelector& sel){};///<	Copy Constructor not yet implemented!

	protected:
		MultiGrid*	m_pMultiGrid;
		LevelVec 	m_levels;
		VertexIterator m_tmpVBegin;
		VertexIterator m_tmpVEnd;

	//	callback-id (automatically unregisters callback, when the selector is deleted).
		MessageHub::SPCallbackId	m_callbackId;

	//	we use a shared attachment for the entry-lists of all section containers
		AttachedVertexList::AEntry	m_aSharedEntryVRT;
		AttachedEdgeList::AEntry	m_aSharedEntryEDGE;
		AttachedFaceList::AEntry	m_aSharedEntryFACE;
		AttachedVolumeList::AEntry	m_aSharedEntryVOL;
};

/** \} */

}//	end of namespace

////////////////////////////////
//	include implementation
#include "selector_multi_grid_impl.hpp"

#endif
