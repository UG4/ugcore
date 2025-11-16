/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__element_visitor_util__
#define __H__UG__element_visitor_util__

#include <vector>
#include "grid/grid.h"
#include "common/util/metaprogramming_util.h"

namespace ug
{

///	Visits all elements between begin and end and executes the visitorCallback on them
template <class TIter>
void VisitAll(const TIter begin, const TIter end,
			  boost::function<void (typename TIter::value_type)> visitorCallback)
{
	TIter iter = begin;
	while(iter != end){
		typename TIter::value_type val = *iter;
		++iter;
		visitorCallback(val);
	}
}

///	Visits all boundary elements of the area specified through the iterators.
/** WARNING: This method uses Grid::mark
 * TIter::value_type has to be compatible to Edge*, Face* or Volume*.
 * Lets say that TIter::value_type equals TElem*.
 *
 * You have to specify a callback cbBelongsToArea with the signature
 * 'void (TElem*)'. which tells for each element of TIter::value_type, whether
 * it lies in the area or not.
 * Note that this callback should optimally return true for all elements between
 * begin and end, and false for all other elements.
 *
 * You furthermore have to specify a callback cbVisitSide, with the signature
 * 'void (TElem::side*)'. This callback is executed for all sides of the
 * given area.
 */
template <class TIter>
void VisitAreaBoundary(Grid& g, const TIter begin, const TIter end,
		      boost::function<bool (typename TIter::value_type)> cbBelongsToArea,
			  boost::function<void (typename Pointer2Value<typename TIter::value_type>::type::side)> cbVisitSide)
{
	using TElem =  typename Pointer2Value<typename TIter::value_type>::type;
	using TSide = typename TElem::side;

	g.begin_marking();

	std::vector<TSide*> sides;
	std::vector<TElem*> elems;
	TIter iter = begin;
	while(iter != facesEnd){
		TElem* elem = *iter;
		++iter;
		CollectAssociated(sides, g, elem);
		for(size_t i = 0; i < sides.size(); ++i){
			TSide* side = sides[i];
			if(!g.is_marked(side)){
				g.mark(side);
			//	collect associated elems and check whether all belong to the area
				CollectAssociated(elems, g, side);
				int numInArea = 0;
				for(size_t i_elem = 0; i_elem < elems.size(); ++i_elem){
					if(cbBelongsToArea(elems[i_elem]))
						++numInArea;
				}

				if(numInArea == 1){
					cbVisitSide(side);
				}
			}
		}
	}

	grid.end_marking();
}

}//	end of namespace

#endif
