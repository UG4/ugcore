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
#include "selector_grid.h"

namespace ug
{

/** \ingroup lib_grid_tools
 *  \{ */

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
 * works on one element type (either Vertex, EdgeBase, Face or
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
class UG_API TElemSelector : public Selector
{
	public:
		typedef typename geometry_traits<TBaseElem>::iterator iterator;
		typedef typename geometry_traits<TBaseElem>::const_iterator const_iterator;


	public:
		TElemSelector()	: Selector(SE_NONE)
		{
			int typId = geometry_traits<TBaseElem>::BASE_OBJECT_ID;
			switch(typId){
				case VERTEX: this->enable_element_support(SE_VERTEX); break;
				case EDGE: this->enable_element_support(SE_EDGE); break;
				case FACE: this->enable_element_support(SE_FACE); break;
				case VOLUME: this->enable_element_support(SE_VOLUME); break;
				default: break;
			}
		}

		TElemSelector(Grid& grid) : Selector(grid, SE_NONE)
		{
			int typId = geometry_traits<TBaseElem>::BASE_OBJECT_ID;
			switch(typId){
				case VERTEX: this->enable_element_support(SE_VERTEX); break;
				case EDGE: this->enable_element_support(SE_EDGE); break;
				case FACE: this->enable_element_support(SE_FACE); break;
				case VOLUME: this->enable_element_support(SE_VOLUME); break;
				default: break;
			}
		}

		using Selector::begin;// support for iteration over Edge, ConstrainedEdge, ...
		inline iterator begin()						{return Selector::begin<TBaseElem>();}
		inline const_iterator begin() const			{return Selector::begin<TBaseElem>();}

		using Selector::end;// support for iteration over Edge, ConstrainedEdge, ...
		inline iterator end()						{return Selector::end<TBaseElem>();}
		inline const_iterator end() const			{return Selector::end<TBaseElem>();}

};

////////////////////////////////////////////////////////////////////////
//	typedefs of the four element-selectors

typedef TElemSelector<Vertex>	VertexSelector;
typedef TElemSelector<EdgeBase>		EdgeSelector;
typedef TElemSelector<Face>			FaceSelector;
typedef TElemSelector<Volume>		VolumeSelector;

/** \} */
}//	end of namespace

#endif
