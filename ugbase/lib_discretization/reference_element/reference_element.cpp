/*
 * reference_element.cpp
 *
 *  Created on: 22.07.2010
 *      Author: andreasvogel
 */

#include "reference_element.h"


namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Reference Element
///////////////////////////////////////////////////////////////////////////////
void ReferenceElement::print_info() const
{
	using namespace std;

	string GeomObjects[4] ={"Corner", "Edge", "Face", "Volume"};

	cout << "Reference Element Info: " << endl;
	cout << "----------------------- " << endl;

	cout << "Size: " << size() << endl;
	cout << "Dimension where Reference Element lives: " << dimension() << endl;

	for(int i = dimension(); i>=0 ;i--)
		cout << "Number of " << GeomObjects[i] << "s: " << num(i) << endl;

	for(int dim_i = dimension(); dim_i>=0 ;dim_i--)
	{
		for(size_t i=0; i < num(dim_i); i++)
		{
			cout << GeomObjects[dim_i] << " with id '" << i << "' contains the following GeomObjects:" << endl;
			for(int dim_j=dim_i; dim_j>= 0; dim_j--)
			{
				cout << num(dim_i,i,dim_j) << " " << GeomObjects[dim_j] << "s with id: ";
				for(size_t j=0; j< num(dim_i,i,dim_j); j++)
				{
					cout << id(dim_i,i,dim_j,j) << " ";
				}
				cout << endl;
			}
		}
	}
}

template <int d>
void DimReferenceElement<d>::print_info() const
{
	using namespace std;

	ReferenceElement::print_info();

	cout << "corners:\n";
	for(size_t i = 0; i< num(0); i++)
	{
		cout << i << ":" << corner(i) << "\n";
	}
}


///////////////////////////////////////////////////////////////////////////////
// Reference Element Factory
///////////////////////////////////////////////////////////////////////////////

// register elements at factory
std::vector<const ReferenceElement* > ReferenceElementProvider::m_vElem =
	std::vector<const ReferenceElement* >();

bool RegisterStandardDimReferenceElements()
{
///////////////
// Vertex
///////////////

static ReferenceElementWrapper<ReferenceVertex> refVertex;

static const bool registered_1 = ReferenceElementProvider::add(refVertex);

///////////////
// Edge
///////////////

static ReferenceElementWrapper<ReferenceEdge> refEdge;
static DimReferenceElementWrapper<ReferenceEdge, 1> dimRefEdge;

static const bool registered_3 = ReferenceElementProvider::add(refEdge);
static const bool registered_4 = DimReferenceElementProvider<1>::add(dimRefEdge);

///////////////
// Triangle
///////////////
static ReferenceElementWrapper<ReferenceTriangle> refTriangle;
static DimReferenceElementWrapper<ReferenceTriangle, 2> dimRefTriangle;

static const bool registered_5 = ReferenceElementProvider::add(refTriangle);
static const bool registered_6 = DimReferenceElementProvider<2>::add(dimRefTriangle);

///////////////
// Quadrilateral
///////////////
static ReferenceElementWrapper<ReferenceQuadrilateral> refQuadrilateral;
static DimReferenceElementWrapper<ReferenceQuadrilateral, 2> dimRefQuadrilateral;

static const bool registered_7 = ReferenceElementProvider::add(refQuadrilateral);
static const bool registered_8 = DimReferenceElementProvider<2>::add(dimRefQuadrilateral);

///////////////
// Tetrahedron
///////////////
static ReferenceElementWrapper<ReferenceTetrahedron> refTetrahedron;
static DimReferenceElementWrapper<ReferenceTetrahedron, 3> dimRefTetrahedron;

static const bool registered_9 = ReferenceElementProvider::add(refTetrahedron);
static const bool registered_10 = DimReferenceElementProvider<3>::add(dimRefTetrahedron);

///////////////
// Prism
///////////////
static ReferenceElementWrapper<ReferencePrism> refPrism;
static DimReferenceElementWrapper<ReferencePrism, 3> dimRefPrism;

static const bool registered_15 = ReferenceElementProvider::add(refPrism);
static const bool registered_16 = DimReferenceElementProvider<3>::add(dimRefPrism);

///////////////
// Pyramid
///////////////
static ReferenceElementWrapper<ReferencePyramid> refPyramid;
static DimReferenceElementWrapper<ReferencePyramid, 3> dimRefPyramid;

static const bool registered_11 = ReferenceElementProvider::add(refPyramid);
static const bool registered_12 = DimReferenceElementProvider<3>::add(dimRefPyramid);

///////////////
// Hexahedron
///////////////
static ReferenceElementWrapper<ReferenceHexahedron> refHexahedron;
static DimReferenceElementWrapper<ReferenceHexahedron, 3> dimRefHexahedron;

static const bool registered_13 = ReferenceElementProvider::add(refHexahedron);
static const bool registered_14 = DimReferenceElementProvider<3>::add(dimRefHexahedron);


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

///////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
///////////////////////////////////////////////////////////////////////////////
template class DimReferenceElement<1>;
template class DimReferenceElement<2>;
template class DimReferenceElement<3>;

}; // end namespace ug

