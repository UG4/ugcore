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

#ifndef __H__UG__MULTI_GRID_UTIL__
#define __H__UG__MULTI_GRID_UTIL__

//#include "lib_grid/lg_base.h"

namespace ug {

/**	iterates over the multi-grid and assigns all surface-elements to surfaceViewOut.*/
template <typename TElem>
void CollectSurfaceViewElements(ISubsetHandler& surfaceViewOut,
								MultiGrid& mg,
								MGSubsetHandler& mgsh,
								bool clearContainer = true);

/**	calls CollectSurfaceViewElements and then assigns all associated elements
 *	of lower dimension to the surface-view, too.
 *
 *	TSurfaceView has to be a SubsetHandler or MGSubsetHandler compatible type.*/
template <typename TElem, typename TSurfaceView>
void CreateSurfaceView(TSurfaceView& surfaceViewOut,
						MultiGrid& mg,
						MGSubsetHandler& mgsh);

///	returns true, if the element lies one level below the surface
/**	If checkSides == false, then only children of the same base type as TElem
 * will be regarded. If checkSides == true, also children of the elements sides
 * will be regarded, too.
 * If all regarded children lie on the surface (i.e. do not have children them selfs),
 * then the element is regarded as a surface element.
 */
template <typename TElem>
bool IsSubSurfaceElement(MultiGrid& mg, TElem* e, bool checkSides = false);


}//	end of namespace

////////////////////////////////
//	include implementation
#include "multi_grid_util_impl.hpp"

#endif
