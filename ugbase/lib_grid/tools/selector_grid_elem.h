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

#ifndef __H__LIBGRID__SELECTOR_GRID_ELEM__
#define __H__LIBGRID__SELECTOR_GRID_ELEM__

// #include <cassert>
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
 * works on one element type (either Vertex, Edge, Face or
 * Volume). Normally you will use the type definitions VertexSelector,
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

template <typename TBaseElem>
class UG_API TElemSelector : public Selector
{
	public:
		using iterator = typename geometry_traits<TBaseElem>::iterator;
		using const_iterator = typename geometry_traits<TBaseElem>::const_iterator;


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

		explicit TElemSelector(Grid& grid) : Selector(grid, SE_NONE)
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
		[[nodiscard]] inline const_iterator begin() const			{return Selector::begin<TBaseElem>();}

		using Selector::end;// support for iteration over Edge, ConstrainedEdge, ...
		inline iterator end()						{return Selector::end<TBaseElem>();}
		[[nodiscard]] inline const_iterator end() const			{return Selector::end<TBaseElem>();}

};

////////////////////////////////////////////////////////////////////////
//	type definitions of the four element-selectors

using VertexSelector = TElemSelector<Vertex>;
using EdgeSelector = TElemSelector<Edge>;
using FaceSelector = TElemSelector<Face>;
using VolumeSelector = TElemSelector<Volume>;

/** \} */
}//	end of namespace

#endif
