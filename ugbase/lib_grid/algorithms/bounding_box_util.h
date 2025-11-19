/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_bounding_box_util
#define __H__UG_bounding_box_util

#include "common/math/misc/shapes.h"

namespace ug{

///	calculates the smallest axis aligned box that contains the given vertex
template <class TAAPos>
AABox<typename TAAPos::ValueType>
CalculateBoundingBox (Vertex* e, TAAPos aaPos)
{
	return AABox<typename TAAPos::ValueType>(aaPos[e], aaPos[e]);
}

///	calculates the smallest axis aligned box that contains the given element
template <class TElem, class TAAPos>
AABox<typename TAAPos::ValueType>
CalculateBoundingBox (TElem* e, TAAPos aaPos)
{
	using box_t = AABox<typename TAAPos::ValueType>;
	typename TElem::ConstVertexArray vrts = e->vertices();
	box_t box(aaPos[vrts[0]], aaPos[vrts[0]]);
	for(size_t i = 1; i < e->num_vertices(); ++i){
		box = box_t(box, aaPos[vrts[i]]);
	}
	return box;
}

}//	end of namespace

#endif