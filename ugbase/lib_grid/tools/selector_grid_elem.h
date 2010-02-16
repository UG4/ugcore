// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d16

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	...
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIBGRID__SELECTOR_GRID_ELEM__
#define __H__LIBGRID__SELECTOR_GRID_ELEM__

#include <cassert>
#include "selector_interface.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	TElemSelector
///	specialization of ISelector for a subset of the elements in a grid of class Grid.
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
 * This is a specialization of ISelector for the class Grid that only
 * works on one element type (either VertexBase, EdgeBase, Face or
 * Volume). Normally you will use the typedefs VertexSelector,
 * EdgeSelector, FaceSelector or VolumeSelector instead of this class.
 *
 * The following methods are the most used:
 *	- select, deselect, is_selected (see ISelector)
 *	- begin, end, num, clear.
 *
 * You may specify the element-type on which begin, end, num and clear
 * operate via a template parameter.
 *
 * \code
 * Grid g;
 * FaceSelector fsel(g);
 *
 * // ... create elements and select some
 *
 * // number of selected triangles
 * int nSelTris = fsel.num<Triangle>();
 *
 * // iteration over all faces
 *	for(FaceIterator iter = fsel.begin(i);
 *		iter != fsel.end(i); ++iter){
 * // ...
 *	}
 * \endcode
 */
template <class TBaseElem>
class TElemSelector : public ISelector
{
	public:
		typedef ISelector	BaseClass;
		typedef typename geometry_traits<TBaseElem>::iterator BaseElemIterator;
	public:
		TElemSelector();
		TElemSelector(Grid& grid);
		virtual ~TElemSelector()	{}

		void assign_grid(Grid& grid);

		virtual void clear();

		template <class TElem>
		inline void clear();

		template <class TElem>
		inline uint num();
		
		inline uint num();
		
	//	empty
		inline bool empty();

		template <class TElem>
		inline bool empty();

	//	begin
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		begin();

		inline BaseElemIterator begin();

	//	end
		template <class TElem>
		inline typename geometry_traits<TElem>::iterator
		end();

		inline BaseElemIterator end();


	//	geometric-object-collection
		GeometricObjectCollection get_geometric_object_collection();

	//	callbacks that allows us to clean-up
	//	derived from GridObserver
		virtual void unregistered_from_grid(Grid* grid);

	protected:
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
		template <class TElem>
		inline SectionContainer& get_section_container();

		template <class TElem>
		inline int get_section_index();

	private:
		TElemSelector(const TElemSelector& sel){};///<	Copy Constructor not yet implemented!

	protected:
		SectionContainer 	m_elements;
};


////////////////////////////////////////////////////////////////////////
//	typedefs of the four element-selectors
typedef TElemSelector<VertexBase>	VertexSelector;
typedef TElemSelector<EdgeBase>		EdgeSelector;
typedef TElemSelector<Face>			FaceSelector;
typedef TElemSelector<Volume>		VolumeSelector;

}//	end of namespace

////////////////////////////////
//	include implementation
#include "selector_grid_elem_impl.hpp"

#endif
