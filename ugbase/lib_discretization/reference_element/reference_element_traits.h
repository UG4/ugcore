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
class reference_element_traits;


///////////////////////////////////////////////////////////////////////////////
// Vertex
///////////////////////////////////////////////////////////////////////////////

template <>
class reference_element_traits<VertexBase>
{
	public:
		typedef ReferenceVertex reference_element_type;
};

template <>
class reference_element_traits<Vertex>
{
	public:
		typedef ReferenceVertex reference_element_type;
};

template <>
class reference_element_traits<HangingVertex>
{
	public:
		typedef ReferenceVertex reference_element_type;
};

///////////////////////////////////////////////////////////////////////////////
// Edge
///////////////////////////////////////////////////////////////////////////////

template <>
class reference_element_traits<Edge>
{
	public:
		typedef ReferenceEdge reference_element_type;
};

template <>
class reference_element_traits<ConstrainedEdge>
{
	public:
		typedef ReferenceEdge reference_element_type;
};

template <>
class reference_element_traits<ConstrainingEdge>
{
	public:
		typedef ReferenceEdge reference_element_type;
};

///////////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////////

template <>
class reference_element_traits<Triangle>
{
	public:
		typedef ReferenceTriangle reference_element_type;
};

template <>
class reference_element_traits<ConstrainedTriangle>
{
	public:
		typedef ReferenceTriangle reference_element_type;
};

template <>
class reference_element_traits<ConstrainingTriangle>
{
	public:
		typedef ReferenceTriangle reference_element_type;
};

///////////////////////////////////////////////////////////////////////////////
// Quadrilateral
///////////////////////////////////////////////////////////////////////////////

template <>
class reference_element_traits<Quadrilateral>
{
	public:
		typedef ReferenceQuadrilateral reference_element_type;
};

template <>
class reference_element_traits<ConstrainedQuadrilateral>
{
	public:
		typedef ReferenceQuadrilateral reference_element_type;
};

template <>
class reference_element_traits<ConstrainingQuadrilateral>
{
	public:
		typedef ReferenceQuadrilateral reference_element_type;
};

///////////////////////////////////////////////////////////////////////////////
// Tetrahedron
///////////////////////////////////////////////////////////////////////////////

template <>
class reference_element_traits<Tetrahedron>
{
	public:
		typedef ReferenceTetrahedron reference_element_type;
};

///////////////////////////////////////////////////////////////////////////////
// Pyramid
///////////////////////////////////////////////////////////////////////////////

template <>
class reference_element_traits<Pyramid>
{
	public:
		typedef ReferencePyramid reference_element_type;
};

///////////////////////////////////////////////////////////////////////////////
// Prism
///////////////////////////////////////////////////////////////////////////////

template <>
class reference_element_traits<Prism>
{
	public:
		typedef ReferencePrism reference_element_type;
};

///////////////////////////////////////////////////////////////////////////////
// Hexahedron
///////////////////////////////////////////////////////////////////////////////

template <>
class reference_element_traits<Hexahedron>
{
	public:
		typedef ReferenceHexahedron reference_element_type;
};

}

#endif /* __H__LIBDISCRETIZATION__REFERENCE_ELEMENT__REFERENCE_ELEMENT_TRAITS__ */
