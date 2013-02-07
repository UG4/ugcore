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
		
	public:
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

		inline uint num_levels()	{return m_levels.size();}

		virtual void clear();

		template <class TElem>
		inline void clear();

		void clear(int level);

		template <class TElem>
		inline void clear(int level);

		template <class TElem>
		inline uint num(int level);
		
		inline uint num(int level);

		template <class TElem>
		inline uint num();
		
		inline uint num();

	//	empty
		inline bool empty(int level);

		template <class TElem>
		inline bool empty(int level);

		inline bool empty();

		template <class TElem>
		inline bool empty();

	//	begin
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		begin(int level);

		template <class TElem>
		inline typename geometry_traits<TElem>::const_iterator
		begin(int level) const;
		
	//	end
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		end(int level);
		
		template <class TElem>
		inline typename geometry_traits<TElem>::const_iterator
		end(int level) const;

	//	convenience begin and end
		inline VertexBaseIterator vertices_begin(int level)	{return begin<VertexBase>(level);}
		inline VertexBaseIterator vertices_end(int level)	{return end<VertexBase>(level);}
		inline EdgeBaseIterator edges_begin(int level)		{return begin<EdgeBase>(level);}
		inline EdgeBaseIterator edges_end(int level)		{return end<EdgeBase>(level);}
		inline FaceIterator faces_begin(int level)			{return begin<Face>(level);}
		inline FaceIterator faces_end(int level)			{return end<Face>(level);}
		inline VolumeIterator volumes_begin(int level)		{return begin<Volume>(level);}
		inline VolumeIterator volumes_end(int level)		{return end<Volume>(level);}

	///	returns the first selected element of the given type on the specified level.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <class TElem> TElem* front(int level);
		
	///	returns the last selected element of the given type on the specified level.
	/**	Make sure that elements of the given type exist!
	 *	Behaviour is undefined, if not.*/
		template <class TElem> TElem* back(int level);
		
	//	geometric-object-collection
		virtual GeometricObjectCollection get_geometric_objects();

	//	callbacks that allows us to clean-up
	//	derived from GridObserver
	/*
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
	*/
		virtual void grid_to_be_destroyed(Grid* grid);
		
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
		get_iterator(VertexBase* o)
		{
			assert((is_selected(o) >= 0) && "object not selected.");
			return section_container<VertexBase>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline EdgeSectionContainer::iterator
		get_iterator(EdgeBase* o)
		{
			assert((is_selected(o) >= 0) && "object not selected");
			return section_container<EdgeBase>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline FaceSectionContainer::iterator
		get_iterator(Face* o)
		{
			assert((is_selected(o) >= 0) && "object not selected");
			return section_container<Face>(m_pMultiGrid->get_level(o)).
				get_container().get_iterator(o);
		}

		inline VolumeSectionContainer::iterator
		get_iterator(Volume* o)
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
