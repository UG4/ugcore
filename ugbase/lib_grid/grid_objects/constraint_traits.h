/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_LIB_GRID__CONSTRAINT_TRAITS__
#define __H__UG_LIB_GRID__CONSTRAINT_TRAITS__

#include "grid_objects_0d.h"
#include "grid_objects_1d.h"
#include "grid_objects_2d.h"

namespace ug{
/**	constraint traits provide the associated constrained and constraining grid
 * grid object types to a given grid object type.*/
template <class TElem>
struct constraint_traits{
	using elem_t = TElem;
	using constrained_t = void;
	using constraining_t = void;
};

template <>
struct constraint_traits<Vertex>{
	using elem_t = Vertex;
	using constrained_t = ConstrainedVertex;
	using constraining_t = void;
};

template <>
struct constraint_traits<Edge>{
	using elem_t = Edge;
	using constrained_t = ConstrainedEdge;
	using constraining_t = ConstrainingEdge;
};

template <>
struct constraint_traits<Face>{
	using elem_t = Face;
	using constrained_t = ConstrainedFace;
	using constraining_t = ConstrainingFace;
};

template <>
struct constraint_traits<Triangle>{
	using elem_t = Triangle;
	using constrained_t = ConstrainedTriangle;
	using constraining_t = ConstrainingTriangle;
};

template <>
struct constraint_traits<Quadrilateral>{
	using elem_t = Quadrilateral;
	using constrained_t = ConstrainedQuadrilateral;
	using constraining_t = ConstrainingQuadrilateral;
};

}//	end of namespace

#endif
