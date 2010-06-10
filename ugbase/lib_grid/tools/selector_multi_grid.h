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

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class MultiGrid;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	MGSelector
///	specialization of ISelector for grids of class MultiGrid.
/**
 * A selector is a useful class, that allows the user to mark
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

class MGSelector : public ISelector
{
	public:
		typedef ISelector	BaseClass;

	public:
		MGSelector(uint supportedElements = SE_ALL);
		MGSelector(MultiGrid& grid, uint supportedElements = SE_ALL);
		virtual ~MGSelector();

		void assign_grid(MultiGrid& grid);
		inline MultiGrid* get_assigned_multi_grid()	{return m_pMultiGrid;}

	///	set the type of elements that shall be handled by the Selector.
	/**	Pass an or-combination of constants enumerated in SelectorElements.
	 *	\sa Selector::enable_element_support*/
	//	forwards to protected ISelector method
		inline void set_supported_elements(uint shElements)	{ISelector::set_supported_elements(shElements);}

	///	enable support for element-types. Does not invalidate previous settings.
	/**	pass an or-combination of constants enumerated in SelectorElements.*/
	//	forwards to protected ISelector method
		inline void enable_element_support(uint shElements)	{ISelector::enable_element_support(shElements);}

	///	disable support for element-types.
	/**	pass an or-combination of constants enumerated in SelectorElements.*/
	//	forwards to protected ISelector method
		void disable_element_support(uint shElements)		{ISelector::disable_element_support(shElements);}

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

	//	end
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		end(int level);

	//	convenience begin and end
		inline VertexBaseIterator vertices_begin(int level)	{return begin<VertexBase>(level);}
		inline VertexBaseIterator vertices_end(int level)	{return end<VertexBase>(level);}
		inline EdgeBaseIterator edges_begin(int level)		{return begin<EdgeBase>(level);}
		inline EdgeBaseIterator edges_end(int level)		{return end<EdgeBase>(level);}
		inline FaceIterator faces_begin(int level)			{return begin<Face>(level);}
		inline FaceIterator faces_end(int level)			{return end<Face>(level);}
		inline VolumeIterator volumes_begin(int level)		{return begin<Volume>(level);}
		inline VolumeIterator volumes_end(int level)		{return end<Volume>(level);}

	//	multi-level-geometric-object-collection
		MultiLevelGeometricObjectCollection
		get_multi_level_geometric_object_collection();

	//	callbacks that allows us to clean-up
	//	derived from GridObserver
	/*
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
	*/
		virtual void grid_to_be_destroyed(Grid* grid);
		
	protected:
		iterator begin(int objID, int level, int section);
		iterator end(int objID, int level, int section);

		void clear_lists();

		virtual iterator add_to_list(VertexBase* elem);
		virtual iterator add_to_list(EdgeBase* elem);
		virtual iterator add_to_list(Face* elem);
		virtual iterator add_to_list(Volume* elem);

		virtual void erase_from_list(VertexBase* elem);
		virtual void erase_from_list(EdgeBase* elem);
		virtual void erase_from_list(Face* elem);
		virtual void erase_from_list(Volume* elem);

	protected:
		struct Level{
			SectionContainer m_elements[NUM_GEOMETRIC_BASE_OBJECTS];
		};
		typedef std::vector<Level*>	LevelVec;

	protected:
		template <class TElem>
		inline SectionContainer& get_section_container(int level);

		template <class TElem>
		inline int get_section_index();

		void level_required(int newSize);

	private:
		MGSelector(const MGSelector& sel){};///<	Copy Constructor not yet implemented!

	protected:
		MultiGrid*	m_pMultiGrid;
		LevelVec 	m_levels;
		VertexBaseIterator m_tmpVBegin;
		VertexBaseIterator m_tmpVEnd;
};

}//	end of namespace

////////////////////////////////
//	include implementation
#include "selector_multi_grid_impl.hpp"

#endif
