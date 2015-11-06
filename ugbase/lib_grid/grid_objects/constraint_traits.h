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
	typedef TElem	elem_t;
	typedef void	constrained_t;
	typedef void	constraining_t;
};

template <>
struct constraint_traits<Vertex>{
	typedef Vertex				elem_t;
	typedef ConstrainedVertex	constrained_t;
	typedef void				constraining_t;
};

template <>
struct constraint_traits<Edge>{
	typedef Edge				elem_t;
	typedef ConstrainedEdge		constrained_t;
	typedef ConstrainingEdge	constraining_t;
};

template <>
struct constraint_traits<Face>{
	typedef Face				elem_t;
	typedef ConstrainedFace		constrained_t;
	typedef ConstrainingFace	constraining_t;
};

template <>
struct constraint_traits<Triangle>{
	typedef Triangle				elem_t;
	typedef ConstrainedTriangle		constrained_t;
	typedef ConstrainingTriangle	constraining_t;
};

template <>
struct constraint_traits<Quadrilateral>{
	typedef Quadrilateral				elem_t;
	typedef ConstrainedQuadrilateral	constrained_t;
	typedef ConstrainingQuadrilateral	constraining_t;
};

}//	end of namespace

#endif
