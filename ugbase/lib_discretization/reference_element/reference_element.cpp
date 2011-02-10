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

bool RegisterStandardDimReferenceElements()
{
///////////////
// Vertex
///////////////

static ReferenceElementWrapper<ReferenceVertex> refVertex;

static const bool registered_1 = ReferenceElementFactory::register_reference_element(refVertex);

///////////////
// Edge
///////////////

static ReferenceElementWrapper<ReferenceEdge> refEdge;
static DimReferenceElementWrapper<ReferenceEdge, 1> dimRefEdge;

static const bool registered_3 = ReferenceElementFactory::register_reference_element(refEdge);
static const bool registered_4 = DimReferenceElementFactory<1>::register_reference_element(dimRefEdge);

///////////////
// Triangle
///////////////
static ReferenceElementWrapper<ReferenceTriangle> refTriangle;
static DimReferenceElementWrapper<ReferenceTriangle, 2> dimRefTriangle;

static const bool registered_5 = ReferenceElementFactory::register_reference_element(refTriangle);
static const bool registered_6 = DimReferenceElementFactory<2>::register_reference_element(dimRefTriangle);

///////////////
// Quadrilateral
///////////////
static ReferenceElementWrapper<ReferenceQuadrilateral> refQuadrilateral;
static DimReferenceElementWrapper<ReferenceQuadrilateral, 2> dimRefQuadrilateral;

static const bool registered_7 = ReferenceElementFactory::register_reference_element(refQuadrilateral);
static const bool registered_8 = DimReferenceElementFactory<2>::register_reference_element(dimRefQuadrilateral);

///////////////
// Tetrahedron
///////////////
static ReferenceElementWrapper<ReferenceTetrahedron> refTetrahedron;
static DimReferenceElementWrapper<ReferenceTetrahedron, 3> dimRefTetrahedron;

static const bool registered_9 = ReferenceElementFactory::register_reference_element(refTetrahedron);
static const bool registered_10 = DimReferenceElementFactory<3>::register_reference_element(dimRefTetrahedron);

///////////////
// Prism
///////////////
static ReferenceElementWrapper<ReferencePrism> refPrism;
static DimReferenceElementWrapper<ReferencePrism, 3> dimRefPrism;

static const bool registered_15 = ReferenceElementFactory::register_reference_element(refPrism);
static const bool registered_16 = DimReferenceElementFactory<3>::register_reference_element(dimRefPrism);

///////////////
// Pyramid
///////////////
static ReferenceElementWrapper<ReferencePyramid> refPyramid;
static DimReferenceElementWrapper<ReferencePyramid, 3> dimRefPyramid;

static const bool registered_11 = ReferenceElementFactory::register_reference_element(refPyramid);
static const bool registered_12 = DimReferenceElementFactory<3>::register_reference_element(dimRefPyramid);

///////////////
// Hexahedron
///////////////
static ReferenceElementWrapper<ReferenceHexahedron> refHexahedron;
static DimReferenceElementWrapper<ReferenceHexahedron, 3> dimRefHexahedron;

static const bool registered_13 = ReferenceElementFactory::register_reference_element(refHexahedron);
static const bool registered_14 = DimReferenceElementFactory<3>::register_reference_element(dimRefHexahedron);


return (registered_1 &&
		registered_3 &&
		registered_4 &&
		registered_5 &&
		registered_6 &&
		registered_7 &&
		registered_8 &&
		registered_9 &&
		registered_10 &&
		registered_11 &&
		registered_12 &&
		registered_13 &&
		registered_14 &&
		registered_15 &&
		registered_16);
}

};
