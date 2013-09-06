// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

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
 * int nSelVrts = sel.num<VertexBase>(1);
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
		typedef ISelector	BaseClass;
		typedef MultiGrid	grid_type;
		
	/** This iterator is used by MGSelector to provide iteration across all levels.
	 * The TMGSelector and TLevelIterator template argument allows to use this
	 * iterator for const and non-const use.*/
		template <class TElem, class TMGSelector, class TLevelIterator>
		class MGSelectionIterator
		{
			public:
				typedef MGSelectionIterator										this_type;
				typedef std::forward_iterator_tag	iterator_category;
				typedef size_t						difference_type;
				typedef TElem**						pointer;
				typedef TElem*						value_type;
				typedef value_type&					reference;

				MGSelectionIterator() : m_sel(NULL), m_lvl(0)	{}

			///	copy constructor that allows creation of const-iterators from non-const iterators
				MGSelectionIterator(const MGSelectionIterator<TElem, MGSelector,
									typename geometry_traits<TElem>::iterator>& iter)
				{
					m_sel = iter.m_sel;
					m_lvl = iter.m_lvl;
					m_iter = iter.m_iter;
				}

				this_type operator ++()	{increment(); return *this;}
				this_type operator ++(int unused)	{this_type i = *this; increment(); return i;}

				bool operator ==(const this_type& iter) const {return equal(iter);}
				bool operator !=(const this_type& iter) const {return !equal(iter);}

				value_type operator *()	{return dereference();}

			private:
				friend class MGSelector;
		//		friend class MGSelectionIterator<TElem, const MGSelector,
		//										 typename MGSelector::traits<TElem>::const_level_iterator>;
				typedef TLevelIterator level_iterator;

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
		template <class TElem>
		struct traits{
			typedef TElem*											value_t;
			typedef typename geometry_traits<TElem>::iterator		level_iterator;
			typedef typename geometry_traits<TElem>::const_iterator	const_level_iterator;
			typedef MGSelectionIterator<TElem, MGSelector,
										level_iterator>				iterator;
			typedef MGSelectionIterator<TElem, const MGSelector,
										const_level_iterator>		const_iterator;
		};


		MGSelector(uint supportedElements = SE_ALL);
		MGSelector(MultiGrid& grid, uint supportedElements = SE_ALL);
		virtual ~MGSelector();

		void assign_grid(MultiGrid& grid);
		void assign_grid(MultiGrid* grid);
		inline MultiGrid* multi_grid()	{return m_pMultiGrid;}

	///	set the type of elements that shall be handled by the Selector.
	/**	Pass an or-combination of constants enumerated in SelectorElements.
	 *	\sa Selector::enable_element_support*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		inline void set_supported_elements(uint shElements);

	///	enable support for element-types. Does not invalidate previous settings.
	/**	pass an or-combination of constants enumerated in SelectorElements.*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		inline void enable_element_support(uint shElements);

	///	disable support for element-types.
	/**	pass an or-combination of constants enumerated in SelectorElements.*/
	//	forwards to protected ISelector method. This rather complicated setup
	//	is required to avoid virtual method calls during construction.
		void disable_element_support(uint shElements);

		inline size_t num_levels() const	{return m_levels.size();}
		inline size_t top_level() const
		{
			size_t l = m_levels.size();
			if(l == 0)
				return 0;
			else
				return l - 1;
		}

		virtual void clear();

		template <class TElem>
		inline void clear();

		void clear(int level);

		template <class TElem>
		inline void clear(int level);

		template <class TElem>
		inline size_t num(int level) const;
		
		inline size_t num(int level) const;

		template <class TElem>
		inline size_t num() const;
		
		inline size_t num() const;

	//	empty
		inline bool empty(int level) const;

		template <class TElem>
		inline bool empty(int level) const;

		inline bool empty() const;

		template <class TElem>
		inline bool empty() const;

	//	begin
		template <class TElem>
		inline typename traits<TElem>::iterator
		begin();

		template <class TElem>
		inline typename traits<TElem>::const_iterator
		begin() const;

		template <class TElem>
		inline typename traits<TElem>::level_iterator
		begin(int level);

		template <class TElem>
		inline typename traits<TElem>::const_level_iterator
		begin(int level) const;
		
	//	end
		template <class TElem>
		inline typename traits<TElem>::iterator
		end();

		template <class TElem>
		inline typename traits<TElem>::const_iterator
		end() const;

		template <class TElem>
		inline typename traits<TElem>::level_iterator
		end(int level);
		
		template <class TElem>
		inline typename traits<TElem>::const_level_iterator
		end(int level) const;

	//	convenience begin and end
		inline traits<VertexBase>::level_iterator	vertices_begin(int level)	{return begin<VertexBase>(level);}
		inline traits<VertexBase>::level_iterator	vertices_end(int level)		{return end<VertexBase>(level);}
		inline traits<EdgeBase>::level_iterator		edges_begin(int level)		{return begin<EdgeBase>(level);}
		inline traits<EdgeBase>::level_iterator		edges_end(int level)		{return end<EdgeBase>(level);}
		inline traits<Face>::level_iterator			faces_begin(int level)		{return begin<Face>(level);}
		inline traits<Face>::level_iterator			faces_end(int level)		{return end<Face>(level);}
		inline traits<Volume>::level_iterator		volumes_begin(int level)	{return begin<Volume>(level);}
		inline traits<Volume>::level_iterator		volumes_end(int level)		{return end<Volume>(level);}

	///	returns the first selected element of the given type on the specified level.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <class TElem> TElem* front(int level);
		
	///	returns the last selected element of the given type on the specified level.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <class TElem> TElem* back(int level);
		
	//	geometric-object-collection
		virtual GeometricObjectCollection get_geometric_objects() const;

	//	callbacks that allows us to clean-up
	//	derived from GridObserver
	/*
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
	*/
		virtual void grid_to_be_destroyed(Grid* grid);
		
	///	returns true if the selector contains vertices
		virtual bool contains_vertices() const	{return num<VertexBase>() > 0;}

	///	returns true if the selector contains edges
		virtual bool contains_edges() const		{return num<EdgeBase>() > 0;}

	///	returns true if the selector contains faces
		virtual bool contains_faces() const		{return num<Face>() > 0;}

	///	returns true if the selector contains volumes
		virtual bool contains_volumes() const	{return num<Volume>() > 0;}

	protected:
		void clear_lists();

		virtual void add_to_list(VertexBase* elem);
		virtual void add_to_list(EdgeBase* elem);
		virtual void add_to_list(Face* elem);
		virtual void add_to_list(Volume* elem);

		virtual void erase_from_list(VertexBase* elem);
		virtual void erase_from_list(EdgeBase* elem);
		virtual void erase_from_list(Face* elem);
		virtual void erase_from_list(Volume* elem);

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
		typedef std::vector<Level*>	LevelVec;

	protected:
	///	returns the section container for the given type, subset and level
		template <class TElem> inline
		typename Grid::traits<TElem>::SectionContainer&
		section_container(int level);

	///	returns the const section container for the given type, subset and level
		template <class TElem> inline
		const typename Grid::traits<TElem>::SectionContainer&
		section_container(int level) const;
		
		template <class TElem>
		inline int get_section_index() const;

		inline void level_required(int newSize);
		void add_level();

	///	This method should only be called if a complete cleanup is required.
		void cleanup();

	///	returns the iterator at which the given element lies in the section container
	/**	This method may only be called if the element is indeed selected
	 * \{
	 */
		inline VertexSectionContainer::iterator
		get_level_iterator(VertexBase* o)
		{
			assert((is_selected(o) >= 0) && "object not selected.");
			return section_container<VertexBase>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline EdgeSectionContainer::iterator
		get_level_iterator(EdgeBase* o)
		{
			assert((is_selected(o) >= 0) && "object not selected");
			return section_container<EdgeBase>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline FaceSectionContainer::iterator
		get_level_iterator(Face* o)
		{
			assert((is_selected(o) >= 0) && "object not selected");
			return section_container<Face>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline VolumeSectionContainer::iterator
		get_level_iterator(Volume* o)
		{
			assert((is_selected(o) >= 0) && "object not selected");
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
		VertexBaseIterator m_tmpVBegin;
		VertexBaseIterator m_tmpVEnd;

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
