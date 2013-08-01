// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

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
 * int nSelVrts = sel.num<VertexBase>();
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
		typedef ISelector	BaseClass;
		typedef Grid		grid_type;

	///	The traits class holds some important types for each element-type
		template <class TElem>
		struct traits{
			typedef typename geometry_traits<TElem>::iterator		iterator;
			typedef typename geometry_traits<TElem>::const_iterator	const_iterator;
			typedef typename geometry_traits<TElem>::iterator		level_iterator;
			typedef typename geometry_traits<TElem>::const_iterator	const_level_iterator;
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
		inline VertexBaseIterator vertices_begin()	{return begin<VertexBase>();}
		inline VertexBaseIterator vertices_end()	{return end<VertexBase>();}
		inline EdgeBaseIterator edges_begin()		{return begin<EdgeBase>();}
		inline EdgeBaseIterator edges_end()			{return end<EdgeBase>();}
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
		virtual GeometricObjectCollection get_geometric_objects();

	//	callbacks that allows us to clean-up
	//	derived from GridObserver
		//virtual void unregistered_from_grid(Grid* grid);
		virtual void grid_to_be_destroyed(Grid* grid);

	////////////////////////////////////////
	//	for compatibility with MGSelector
	///	always returns 1
		inline size_t num_levels();
		
	///	calls num();
		inline uint num(size_t);
	///	calls num<TElem>();
		template <class TElem> inline size_t num(size_t);
		
	//	empty
	///	calls empty();
		inline bool empty(size_t);
	//	calls empty<TElem>();
		template <class TElem>
		inline bool empty(size_t);

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

		virtual void add_to_list(VertexBase* elem);
		virtual void add_to_list(EdgeBase* elem);
		virtual void add_to_list(Face* elem);
		virtual void add_to_list(Volume* elem);

		virtual void erase_from_list(VertexBase* elem);
		virtual void erase_from_list(EdgeBase* elem);
		virtual void erase_from_list(Face* elem);
		virtual void erase_from_list(Volume* elem);

	///	returns the iterator at which the given element lies in the section container
	/**	This method may only be called if the element is indeed selected
	 * \{
	 */
		inline VertexSectionContainer::iterator
		get_iterator(VertexBase* o)
		{
			assert((is_selected(o) >= 0) && "object not selected.");
			return section_container<VertexBase>().get_container().get_iterator(o);
		}

		inline EdgeSectionContainer::iterator
		get_iterator(EdgeBase* o)
		{
			assert((is_selected(o) >= 0) && "object not selected");
			return section_container<EdgeBase>().get_container().get_iterator(o);
		}

		inline FaceSectionContainer::iterator
		get_iterator(Face* o)
		{
			assert((is_selected(o) >= 0) && "object not selected");
			return section_container<Face>().get_container().get_iterator(o);
		}

		inline VolumeSectionContainer::iterator
		get_iterator(Volume* o)
		{
			assert((is_selected(o) >= 0) && "object not selected");
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
