/*
 * reference_element_traits.h
 *
 *  Created on: 30.06.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_TRAITS__
#define __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_TRAITS__

#include "reference_vertex.h"
#include "reference_edge.h"
#include "reference_triangle.h"
#include "reference_quadrilateral.h"
#include "reference_tetrahedron.h"
#include "reference_pyramid.h"
#include "reference_prism.h"
#include "reference_hexahedron.h"

namespace ug{

/// returns the reference element dimension at run-time
inline int ReferenceElementDimension(ReferenceObjectID roid)
{
	switch(roid)
	{
		case ROID_VERTEX: return 0;
		case ROID_EDGE: return 1;
		case ROID_TRIANGLE: return 2;
		case ROID_QUADRILATERAL: return 2;
		case ROID_TETRAHEDRON: return 3;
		case ROID_PYRAMID: return 3;
		case ROID_PRISM: return 3;
		case ROID_HEXAHEDRON: return 3;
		default: throw(UGFatalError("ReferenceObjectId not found."));
	}
}

/// traits for reference elements
/**
 * The traits class provides for a Grid-element type the corresponding
 * reference-element type. It is used to determine the reference element type
 * at compile time.
 */
template <class TElem>
struct reference_element_traits;


///////////////////////////////////////////////////////////////////////////////
// Vertex
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<VertexBase>
{
	typedef ReferenceVertex reference_element_type;
	static const int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<Vertex>
	: public reference_element_traits<VertexBase> {};

template <>
struct reference_element_traits<HangingVertex>
	: public reference_element_traits<VertexBase> {};

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Edge>
{
	typedef ReferenceEdge reference_element_type;
	static const int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<ConstrainedEdge>
	: public reference_element_traits<Edge>{};

template <>
struct reference_element_traits<ConstrainingEdge>
	: public reference_element_traits<Edge>{};

///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Triangle>
{
	typedef ReferenceTriangle reference_element_type;
	static const int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<ConstrainedTriangle>
	: public reference_element_traits<Triangle> {};

template <>
struct reference_element_traits<ConstrainingTriangle>
	: public reference_element_traits<Triangle> {};

///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Quadrilateral>
{
	typedef ReferenceQuadrilateral reference_element_type;
	static const int dim = reference_element_type::dim;
};

template <>
struct reference_element_traits<ConstrainedQuadrilateral>
	: public reference_element_traits<Quadrilateral>{};

template <>
struct reference_element_traits<ConstrainingQuadrilateral>
	: public reference_element_traits<Quadrilateral>{};

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Tetrahedron>
{
	typedef ReferenceTetrahedron reference_element_type;
	static const int dim = reference_element_type::dim;
};

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Pyramid>
{
	typedef ReferencePyramid reference_element_type;
	static const int dim = reference_element_type::dim;
};

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Prism>
{
	typedef ReferencePrism reference_element_type;
	static const int dim = reference_element_type::dim;
};

///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <>
struct reference_element_traits<Hexahedron>
{
	typedef ReferenceHexahedron reference_element_type;
	static const int dim = reference_element_type::dim;
};

}

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_TRAITS__ */
