/*
 * reference_element.cpp
 *
 *  Created on: 22.07.2010
 *      Author: andreasvogel
 */

#include "reference_element.h"


// register elements at factory
namespace {
using namespace ug;

std::vector<const ReferenceElement* > ReferenceElementFactory::m_vElem =
	std::vector<const ReferenceElement* >();

///////////////
// Vertex
///////////////
/*
ReferenceElementWrapper<ReferenceVertex> refVertex;
DimReferenceElementWrapper<ReferenceVertex, 2> dimRefVertex;

static const bool registered_1 = ReferenceElementFactory::register_reference_element(refVertex);
static const bool registered_2 = DimReferenceElementFactory<2>::register_reference_element(dimRefVertex);
*/
///////////////
// Edge
///////////////

ReferenceElementWrapper<ReferenceEdge> refEdge;
DimReferenceElementWrapper<ReferenceEdge, 1> dimRefEdge;

static const bool registered_3 = ReferenceElementFactory::register_reference_element(refEdge);
static const bool registered_4 = DimReferenceElementFactory<1>::register_reference_element(dimRefEdge);

///////////////
// Triangle
///////////////
ReferenceElementWrapper<ReferenceTriangle> refTriangle;
DimReferenceElementWrapper<ReferenceTriangle, 2> dimRefTriangle;

static const bool registered_5 = ReferenceElementFactory::register_reference_element(refTriangle);
static const bool registered_6 = DimReferenceElementFactory<2>::register_reference_element(dimRefTriangle);

///////////////
// Quadrilateral
///////////////
ReferenceElementWrapper<ReferenceQuadrilateral> refQuadrilateral;
DimReferenceElementWrapper<ReferenceQuadrilateral, 2> dimRefQuadrilateral;

static const bool registered_7 = ReferenceElementFactory::register_reference_element(refQuadrilateral);
static const bool registered_8 = DimReferenceElementFactory<2>::register_reference_element(dimRefQuadrilateral);

///////////////
// Tetrahedron
///////////////
ReferenceElementWrapper<ReferenceTetrahedron> refTetrahedron;
DimReferenceElementWrapper<ReferenceTetrahedron, 3> dimRefTetrahedron;

static const bool registered_9 = ReferenceElementFactory::register_reference_element(refTetrahedron);
static const bool registered_10 = DimReferenceElementFactory<3>::register_reference_element(dimRefTetrahedron);

};
