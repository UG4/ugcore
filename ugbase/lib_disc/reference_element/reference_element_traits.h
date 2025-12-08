/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_TRAITS__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__REFERENCE_ELEMENT_TRAITS__

#include "reference_element.h"
#include "lib_grid/grid_objects/grid_objects.h"

namespace ug{

/// traits for reference elements
/**
 * The traits class provides for a Grid-element type the corresponding
 * reference-element type. It is used to determine the reference element type
 * at compile time.
 */
template <typename TElem>
struct reference_element_traits;


///////////////////////////////////////////////////////////////////////////////
// RegularVertex
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Vertex>
{
	using reference_element_type = ReferenceVertex;
	static constexpr int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<RegularVertex>
	: reference_element_traits<Vertex> {};

template <>
struct reference_element_traits<ConstrainedVertex>
	: reference_element_traits<Vertex> {};

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<RegularEdge>
{
	using reference_element_type = ReferenceEdge;
	static constexpr int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<Edge>
{
	using reference_element_type = ReferenceEdge;
	static constexpr int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<ConstrainedEdge>
	: reference_element_traits<RegularEdge>{};

template <>
struct reference_element_traits<ConstrainingEdge>
	: reference_element_traits<RegularEdge>{};

///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Triangle>
{
	using reference_element_type = ReferenceTriangle;
	static constexpr int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<ConstrainedTriangle>
	: reference_element_traits<Triangle> {};

template <>
struct reference_element_traits<ConstrainingTriangle>
	: reference_element_traits<Triangle> {};

///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Quadrilateral>
{
	using reference_element_type = ReferenceQuadrilateral;
	static constexpr int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<ConstrainedQuadrilateral>
	: reference_element_traits<Quadrilateral>{};

template <>
struct reference_element_traits<ConstrainingQuadrilateral>
	: reference_element_traits<Quadrilateral>{};

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Tetrahedron>
{
	using reference_element_type = ReferenceTetrahedron;
	static constexpr int dim = reference_element_type::dim;
};

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Pyramid>
{
	using reference_element_type = ReferencePyramid;
	static constexpr int dim = reference_element_type::dim;
};

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Prism>
{
	using reference_element_type = ReferencePrism;
	static constexpr int dim = reference_element_type::dim;
};

///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Hexahedron>
{
	using reference_element_type = ReferenceHexahedron;
	static constexpr int dim = reference_element_type::dim;
};

///////////////////////////////////////////////////////////////////////////////
// Octahedron
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Octahedron>
{
	using reference_element_type = ReferenceOctahedron;
	static constexpr int dim = reference_element_type::dim;
};

}

#endif