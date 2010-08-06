/*
 * reference_element.cpp
 *
 *  Created on: 22.07.2010
 *      Author: andreasvogel
 */

#include "reference_element.h"


// register elements at factory
namespace ug{
std::vector<const ReferenceElement* > ReferenceElementFactory::m_vElem =
	std::vector<const ReferenceElement* >();
}

namespace {
using namespace ug;

///////////////
// Vertex
///////////////

ReferenceElementWrapper<ReferenceVertex> refVertex;

static const bool registered_1 = ReferenceElementFactory::register_reference_element(refVertex);

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

///////////////
// Prism
///////////////
ReferenceElementWrapper<ReferencePrism> refPrism;
DimReferenceElementWrapper<ReferencePrism, 3> dimRefPrism;

static const bool registered_15 = ReferenceElementFactory::register_reference_element(refPrism);
static const bool registered_16 = DimReferenceElementFactory<3>::register_reference_element(dimRefPrism);

///////////////
// Pyramid
///////////////
ReferenceElementWrapper<ReferencePyramid> refPyramid;
DimReferenceElementWrapper<ReferencePyramid, 3> dimRefPyramid;

static const bool registered_11 = ReferenceElementFactory::register_reference_element(refPyramid);
static const bool registered_12 = DimReferenceElementFactory<3>::register_reference_element(dimRefPyramid);

///////////////
// Hexahedron
///////////////
ReferenceElementWrapper<ReferenceHexahedron> refHexahedron;
DimReferenceElementWrapper<ReferenceHexahedron, 3> dimRefHexahedron;

static const bool registered_13 = ReferenceElementFactory::register_reference_element(refHexahedron);
static const bool registered_14 = DimReferenceElementFactory<3>::register_reference_element(dimRefHexahedron);

};
