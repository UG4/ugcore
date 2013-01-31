/*
 * reference_element.cpp
 *
 *  Created on: 22.07.2010
 *      Author: andreasvogel
 */

#include "reference_element.h"
#include "common/util/provider.h"


namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Reference Element
///////////////////////////////////////////////////////////////////////////////
ReferenceElement::ReferenceElement()
{
	for(size_t d = 0; d < MAXDIM+1; ++d)
	{
		m_vNum[d] = 0;

		for(size_t obj = 0; obj < MAXOBJECTS; ++obj)
		{
			m_vRefElemType[d][obj] = ROID_UNKNOWN;
			for(size_t d2 = 0; d2 < MAXDIM+1; ++d2)
				m_vSubNum[d][obj][d2] = 0;
		}
	}

	for(size_t d = 0; d < MAXDIM+1; ++d)
		for(size_t obj = 0; obj < MAXOBJECTS; ++obj)
			for(size_t d2 = 0; d2 < MAXDIM+1; ++d2)
				for(size_t obj2 = 0; obj2 < MAXOBJECTS; ++obj2)
					m_id[d][obj][d2][obj2] = -1;

	for(size_t roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
		m_vNumRefElem[roid] = 0;
}

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
// Reference Element Provider
///////////////////////////////////////////////////////////////////////////////

// register elements at factory
const ReferenceElement* ReferenceElementProvider::m_vElem[NUM_REFERENCE_OBJECTS];

ReferenceElementProvider::
ReferenceElementProvider()
{
	static bool bInit = false;

	if(!bInit)
	{
		bInit = true;
		bool bRes = true;
		bRes &= add_elem(Provider<ReferenceVertex>::get());
		// not adding for reference vertex

		bRes &= add_elem(Provider<ReferenceEdge>::get());
		bRes &= add_dim_elem<1>(Provider<ReferenceEdge>::get());

		bRes &= add_elem(Provider<ReferenceTriangle>::get());
		bRes &= add_dim_elem<2>(Provider<ReferenceTriangle>::get());

		bRes &= add_elem(Provider<ReferenceQuadrilateral>::get());
		bRes &= add_dim_elem<2>(Provider<ReferenceQuadrilateral>::get());

		bRes &= add_elem(Provider<ReferenceTetrahedron>::get());
		bRes &= add_dim_elem<3>(Provider<ReferenceTetrahedron>::get());

		bRes &= add_elem(Provider<ReferencePrism>::get());
		bRes &= add_dim_elem<3>(Provider<ReferencePrism>::get());

		bRes &= add_elem(Provider<ReferencePyramid>::get());
		bRes &= add_dim_elem<3>(Provider<ReferencePyramid>::get());

		bRes &= add_elem(Provider<ReferenceHexahedron>::get());
		bRes &= add_dim_elem<3>(Provider<ReferenceHexahedron>::get());

		if(!bRes) UG_THROW("Error while registering Reference Elements");
	}
}

bool ReferenceElementProvider::add_elem(const ReferenceElement& elem)
{
	const ReferenceObjectID roid = elem.reference_object_id();
	UG_ASSERT(roid >= 0, "roid ="<<roid<<" wrong")
	UG_ASSERT(roid < NUM_REFERENCE_OBJECTS, "roid ="<<roid<<" wrong")
	m_vElem[roid] = &elem;
	return true;
}

const ReferenceElement& ReferenceElementProvider::get_elem(ReferenceObjectID roid)
{
	UG_ASSERT(roid >= 0, "roid ="<<roid<<" wrong")
	UG_ASSERT(roid < NUM_REFERENCE_OBJECTS, "roid ="<<roid<<" wrong")
	UG_ASSERT(m_vElem[roid] != NULL, "Null pointer for roid ="<<roid);
	return *m_vElem[roid];
}

template <int dim>
bool ReferenceElementProvider::add_dim_elem(const DimReferenceElement<dim>& elem)
{
	const ReferenceObjectID roid = elem.reference_object_id();
	UG_ASSERT(roid >= 0, "roid ="<<roid<<" wrong")
	UG_ASSERT(roid < NUM_REFERENCE_OBJECTS, "roid ="<<roid<<" wrong")
	static const DimReferenceElement<dim>** vDimElem = get_vector<dim>();
	vDimElem[roid] = &elem;
	return true;
}


///////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
///////////////////////////////////////////////////////////////////////////////
template class DimReferenceElement<1>;
template class DimReferenceElement<2>;
template class DimReferenceElement<3>;



///////////////////////////////////////////////////////////////////////////////
//	ReferenceVertex
///////////////////////////////////////////////////////////////////////////////
ReferenceVertex::ReferenceVertex()
{
	// dimension
	m_dim = 0;

	// size
	m_size = 1.0;

	//number of Geometric Objects
	m_vNum[POINT] = 1;

	// number of Geometric Objects
	m_vSubNum[POINT][0][POINT] = 1;

	m_vRefElemType[POINT][0] = ROID_VERTEX;

	//reset m_id to -1
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<MAXOBJECTS; ++j)
			for(int k=0; k<=dim; ++k)
				for(size_t l=0; l<MAXOBJECTS; l++)
				{
					m_id[i][j][k][l] = -1;
				}

	//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<m_vNum[i]; ++j)
		{
			m_id[i][j][i][0] = j;
		}

	// Reference Element Types
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		m_vNumRefElem[i] = 0;
	}
	m_vNumRefElem[ROID_VERTEX] = 1;
}

///////////////////////////////////////////////////////////////////////////////
//	ReferenceEdge
///////////////////////////////////////////////////////////////////////////////

ReferenceEdge::ReferenceEdge()
{
	// dimension
	m_dim = 1;

	// size
	m_size = 1.0;

	//number of Geometric Objects
	m_vNum[POINT] = 2;
	m_vNum[EDGE] = 1;

	// number of Geometric Objects
	m_vSubNum[EDGE][0][POINT] = 2;
	m_vSubNum[EDGE][0][EDGE] = 1;

	m_vRefElemType[EDGE][0] = ROID_EDGE;

	for(size_t i = 0; i < m_vNum[POINT]; ++i)
	{
		m_vSubNum[POINT][i][POINT] = 1;
		m_vSubNum[POINT][i][EDGE] = 1;

		m_vRefElemType[POINT][i] = ROID_VERTEX;
	}

	//reset m_id to -1
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<MAXOBJECTS; ++j)
			for(int k=0; k<=dim; ++k)
				for(size_t l=0; l<MAXOBJECTS; l++)
				{
					m_id[i][j][k][l] = -1;
				}

	//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<m_vNum[i]; ++j)
		{
			m_id[i][j][i][0] = j;
		}

	// Points <-> Face
	for(size_t i=0; i<m_vNum[POINT]; ++i)
	{
		m_id[EDGE][0][POINT][i] = i;
		m_id[POINT][i][EDGE][0] = 0;
	}

	// Reference Corners
	m_vCorner[0][0] = 0.0;
	m_vCorner[1][0] = 1.0;

	m_vCoInt[0][0] = 0;
	m_vCoInt[1][0] = 1;

	// Reference Element Types
	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		m_vNumRefElem[i] = 0;
	}
	m_vNumRefElem[ROID_VERTEX] = 2;
	m_vNumRefElem[ROID_EDGE] = 1;
}

///////////////////////////////////////////////////////////////////////////////
//	ReferenceTriangle
///////////////////////////////////////////////////////////////////////////////

ReferenceTriangle::ReferenceTriangle()
{
	// dimension
	m_dim = 2;

	// size
	m_size = 0.5;

	//number of Geometric Objects
 	m_vNum[POINT] = 3;
 	m_vNum[EDGE] = 3;
 	m_vNum[FACE] = 1;

	// number of Geometric Objects
 	m_vSubNum[FACE][0][POINT] = 3;
 	m_vSubNum[FACE][0][EDGE] = 3;
 	m_vSubNum[FACE][0][FACE] = 1;

 	m_vRefElemType[FACE][0] = ROID_TRIANGLE;

 	for(size_t i = 0; i < m_vNum[EDGE]; ++i)
 	{
	 	m_vSubNum[EDGE][i][POINT] = 2;
	 	m_vSubNum[EDGE][i][EDGE] = 1;
	 	m_vSubNum[EDGE][i][FACE] = 1;

	 	m_vRefElemType[EDGE][i] = ROID_EDGE;
 	}

 	for(size_t i = 0; i < m_vNum[POINT]; ++i)
 	{
	 	m_vSubNum[POINT][i][POINT] = 1;
	 	m_vSubNum[POINT][i][EDGE] = 2;
	 	m_vSubNum[POINT][i][FACE] = 1;

	 	m_vRefElemType[POINT][i] = ROID_VERTEX;
 	}

	//reset m_id to -1
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<MAXOBJECTS; ++j)
			for(int k=0; k<=dim; ++k)
				for(size_t l=0; l<MAXOBJECTS; l++)
				{
				 	m_id[i][j][k][l] = -1;
				}

	//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<m_vNum[i]; ++j)
		{
		 	m_id[i][j][i][0] = j;
		}

	//Edges <-> Face
	for(size_t i=0; i<m_vNum[EDGE]; ++i)
	{
	 	m_id[FACE][0][EDGE][i] = i;
	 	m_id[EDGE][i][FACE][0] = 0;
	}

	// Points <-> Face
	for(size_t i=0; i<m_vNum[POINT]; ++i)
	{
	 	m_id[FACE][0][POINT][i] = i;
	 	m_id[POINT][i][FACE][0] = 0;
	}

	// Points of Edges
	// edge 0 = (0,1)
 	m_id[EDGE][0][POINT][0] = 0;
 	m_id[EDGE][0][POINT][1] = 1;
 	// edge 1 = (1,2)
 	m_id[EDGE][1][POINT][0] = 1;
 	m_id[EDGE][1][POINT][1] = 2;
	// edge 2 = (2,0)
 	m_id[EDGE][2][POINT][0] = 2;
 	m_id[EDGE][2][POINT][1] = 0;

 	// Edges of Point
 	m_id[POINT][0][EDGE][0] = 2;
 	m_id[POINT][0][EDGE][1] = 0;

 	m_id[POINT][1][EDGE][0] = 0;
 	m_id[POINT][1][EDGE][1] = 1;

 	m_id[POINT][2][EDGE][0] = 1;
 	m_id[POINT][2][EDGE][1] = 2;


	// Reference Corners
 	m_vCorner[0] = MathVector<dim>(0.0, 0.0);
 	m_vCorner[1] = MathVector<dim>(1.0, 0.0);
 	m_vCorner[2] = MathVector<dim>(0.0, 1.0);

 	m_vCoInt[0] = MathVector<dim,int>(0, 0);
 	m_vCoInt[1] = MathVector<dim,int>(1, 0);
 	m_vCoInt[2] = MathVector<dim,int>(0, 1);

 	// Reference Element Types
 	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
 	{
		m_vNumRefElem[i] = 0;
 	}
 	m_vNumRefElem[ROID_VERTEX] = 3;
 	m_vNumRefElem[ROID_EDGE] = 3;
 	m_vNumRefElem[ROID_TRIANGLE] = 1;
}

///////////////////////////////////////////////////////////////////////////////
//	ReferenceQuadrilateral
///////////////////////////////////////////////////////////////////////////////

ReferenceQuadrilateral::ReferenceQuadrilateral()
{
	// dimension
	m_dim = 2;

	// size
	m_size = 1.0;

	//number of Geometric Objects
 	m_vNum[POINT] = 4;
 	m_vNum[EDGE] = 4;
 	m_vNum[FACE] = 1;

	// number of Geometric Objects
 	m_vSubNum[FACE][0][POINT] = 4;
 	m_vSubNum[FACE][0][EDGE] = 4;
 	m_vSubNum[FACE][0][FACE] = 1;
 	m_vRefElemType[FACE][0] = ROID_QUADRILATERAL;

 	for(size_t i = 0; i < m_vNum[EDGE]; ++i)
 	{
 		m_vSubNum[EDGE][i][EDGE] = 1;
	 	m_vSubNum[EDGE][i][POINT] = 2;
	 	m_vSubNum[EDGE][i][FACE] = 1;

	 	m_vRefElemType[EDGE][i] = ROID_EDGE;
 	}

 	for(size_t i = 0; i < m_vNum[EDGE]; ++i)
 	{
 		m_vSubNum[POINT][i][POINT] = 1;
 		m_vSubNum[POINT][i][EDGE] = 2;
 		m_vSubNum[POINT][i][FACE] = 1;

	 	m_vRefElemType[POINT][i] = ROID_VERTEX;
 	}

	//reset m_id to -1
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<MAXOBJECTS; ++j)
			for(int k=0; k<=dim; ++k)
				for(size_t l=0; l<MAXOBJECTS; l++)
				{
				 	m_id[i][j][k][l] = -1;
				}

	//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
	for(int d=0; d<=dim; ++d)
		for(size_t j=0; j<m_vNum[d]; ++j)
		{
		 	m_id[d][j][d][0] = j;
		}

	//Edges <-> Face
	for(size_t i=0; i<m_vNum[EDGE]; ++i)
	{
	 	m_id[FACE][0][EDGE][i] = i;
	 	m_id[EDGE][i][FACE][0] = 0;
	}

	// Points <-> Face
	for(size_t i=0; i<m_vNum[POINT]; ++i)
	{
	 	m_id[FACE][0][POINT][i] = i;
	 	m_id[POINT][i][FACE][0] = 0;
	}

	// Points of Edges
	// edge 0 = (0,1)
 	m_id[EDGE][0][POINT][0] = 0;
 	m_id[EDGE][0][POINT][1] = 1;
	// edge 1 = (1,2)
 	m_id[EDGE][1][POINT][0] = 1;
 	m_id[EDGE][1][POINT][1] = 2;
	// edge 2 = (2,3)
 	m_id[EDGE][2][POINT][0] = 2;
 	m_id[EDGE][2][POINT][1] = 3;
	// edge 3 = (3,0)
 	m_id[EDGE][3][POINT][0] = 3;
 	m_id[EDGE][3][POINT][1] = 0;

 	// Edges of Point
 	m_id[POINT][0][EDGE][0] = 3;
 	m_id[POINT][0][EDGE][1] = 0;

 	m_id[POINT][1][EDGE][0] = 0;
 	m_id[POINT][1][EDGE][1] = 1;

 	m_id[POINT][2][EDGE][0] = 1;
 	m_id[POINT][2][EDGE][1] = 2;

 	m_id[POINT][3][EDGE][0] = 2;
 	m_id[POINT][3][EDGE][1] = 3;

	// Reference Corners
 	m_vCorner[0] = MathVector<dim>(0.0, 0.0);
 	m_vCorner[1] = MathVector<dim>(1.0, 0.0);
 	m_vCorner[2] = MathVector<dim>(1.0, 1.0);
 	m_vCorner[3] = MathVector<dim>(0.0, 1.0);

 	m_vCoInt[0] = MathVector<dim,int>(0, 0);
 	m_vCoInt[1] = MathVector<dim,int>(1, 0);
 	m_vCoInt[2] = MathVector<dim,int>(1, 1);
 	m_vCoInt[3] = MathVector<dim,int>(0, 1);

 	// Reference Element Types
 	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
 	{
		m_vNumRefElem[i] = 0;
 	}
 	m_vNumRefElem[ROID_VERTEX] = 4;
 	m_vNumRefElem[ROID_EDGE] = 4;
 	m_vNumRefElem[ROID_QUADRILATERAL] = 1;
}

///////////////////////////////////////////////////////////////////////////////
//	ReferenceTetrahedron
///////////////////////////////////////////////////////////////////////////////

ReferenceTetrahedron::ReferenceTetrahedron()
{
	// dimension
	m_dim = 3;

	// size
	m_size = 1.0/6.0;

	//number of Geometric Objects
 	m_vNum[POINT] = 4;
 	m_vNum[EDGE] = 6;
 	m_vNum[FACE] = 4;
 	m_vNum[VOLUME] = 1;

	// number of Geometric Objects
 	m_vSubNum[VOLUME][0][POINT] = 4;
 	m_vSubNum[VOLUME][0][EDGE] = 6;
 	m_vSubNum[VOLUME][0][FACE] = 4;
 	m_vSubNum[VOLUME][0][VOLUME] = 1;
	m_vRefElemType[VOLUME][0] = ROID_TETRAHEDRON;

 	for(size_t i = 0; i < m_vNum[FACE]; ++i)
 	{
		m_vSubNum[FACE][i][POINT] = 3;
		m_vSubNum[FACE][i][EDGE] = 3;
		m_vSubNum[FACE][i][FACE] = 1;
		m_vSubNum[FACE][i][VOLUME] = 1;

		m_vRefElemType[FACE][i] = ROID_TRIANGLE;
 	}

 	for(size_t i = 0; i < m_vNum[EDGE]; ++i)
 	{
	 	m_vSubNum[EDGE][i][POINT] = 2;
	 	m_vSubNum[EDGE][i][EDGE] = 1;
	 	m_vSubNum[EDGE][i][FACE] = 2;
	 	m_vSubNum[EDGE][i][VOLUME] = 1;

	 	m_vRefElemType[EDGE][i] = ROID_EDGE;
 	}

 	for(size_t i = 0; i < m_vNum[POINT]; ++i)
 	{
	 	m_vSubNum[POINT][i][POINT] = 1;
	 	m_vSubNum[POINT][i][EDGE] = 3;
	 	m_vSubNum[POINT][i][FACE] = 3;
	 	m_vSubNum[POINT][i][VOLUME] = 1;

	 	m_vRefElemType[POINT][i] = ROID_VERTEX;
 	}

	//reset m_id to -1
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<MAXOBJECTS; ++j)
			for(int k=0; k<=dim; ++k)
				for(size_t l=0; l<MAXOBJECTS; l++)
				{
				 	m_id[i][j][k][l] = -1;
				}

	//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<m_vNum[i]; ++j)
		{
		 	m_id[i][j][i][0] = j;
		}

	// Face <-> Volume
	for(size_t i=0; i<m_vNum[FACE]; ++i)
	{
	 	m_id[VOLUME][0][FACE][i] = i;
	 	m_id[FACE][i][VOLUME][0] = 0;
	}

	// Edge <-> Volume
	for(size_t i=0; i<m_vNum[EDGE]; ++i)
	{
	 	m_id[VOLUME][0][EDGE][i] = i;
	 	m_id[EDGE][i][VOLUME][0] = 0;
	}

	// Point <-> Volume
	for(size_t i=0; i<m_vNum[POINT]; ++i)
	{
	 	m_id[VOLUME][0][POINT][i] = i;
	 	m_id[POINT][i][VOLUME][0] = 0;
	}

	// Points <-> Faces
 	m_id[FACE][0][POINT][0] = 0;
 	m_id[FACE][0][POINT][1] = 2;
 	m_id[FACE][0][POINT][2] = 1;

 	m_id[FACE][1][POINT][0] = 1;
 	m_id[FACE][1][POINT][1] = 2;
 	m_id[FACE][1][POINT][2] = 3;

 	m_id[FACE][2][POINT][0] = 0;
 	m_id[FACE][2][POINT][1] = 3;
 	m_id[FACE][2][POINT][2] = 2;

 	m_id[FACE][3][POINT][0] = 0;
 	m_id[FACE][3][POINT][1] = 1;
 	m_id[FACE][3][POINT][2] = 3;


 	m_id[POINT][0][FACE][0] = 0;
 	m_id[POINT][0][FACE][1] = 2;
 	m_id[POINT][0][FACE][2] = 3;

 	m_id[POINT][1][FACE][0] = 0;
 	m_id[POINT][1][FACE][1] = 3;
 	m_id[POINT][1][FACE][2] = 1;

 	m_id[POINT][2][FACE][0] = 0;
 	m_id[POINT][2][FACE][1] = 1;
 	m_id[POINT][2][FACE][2] = 2;

 	m_id[POINT][3][FACE][0] = 3;
 	m_id[POINT][3][FACE][1] = 2;
 	m_id[POINT][3][FACE][2] = 1;

 	// Edges <-> Faces
 	m_id[FACE][0][EDGE][0] = 2;
 	m_id[FACE][0][EDGE][1] = 1;
 	m_id[FACE][0][EDGE][2] = 0;

 	m_id[FACE][1][EDGE][0] = 1;
 	m_id[FACE][1][EDGE][1] = 5;
 	m_id[FACE][1][EDGE][2] = 4;

 	m_id[FACE][2][EDGE][0] = 3;
 	m_id[FACE][2][EDGE][1] = 5;
 	m_id[FACE][2][EDGE][2] = 2;

 	m_id[FACE][3][EDGE][0] = 0;
 	m_id[FACE][3][EDGE][1] = 4;
 	m_id[FACE][3][EDGE][2] = 3;

 	m_id[EDGE][0][FACE][0] = 0;
 	m_id[EDGE][0][FACE][1] = 3;

 	m_id[EDGE][1][FACE][0] = 0;
 	m_id[EDGE][1][FACE][1] = 1;

 	m_id[EDGE][2][FACE][0] = 0;
 	m_id[EDGE][2][FACE][1] = 2;

 	m_id[EDGE][3][FACE][0] = 3;
 	m_id[EDGE][3][FACE][1] = 2;

 	m_id[EDGE][4][FACE][0] = 1;
 	m_id[EDGE][4][FACE][1] = 3;

 	m_id[EDGE][5][FACE][0] = 2;
 	m_id[EDGE][5][FACE][1] = 1;


	// Points of Edges
	// edge 0 = (0,1)
 	m_id[EDGE][0][POINT][0] = 0;
 	m_id[EDGE][0][POINT][1] = 1;
 	// edge 1 = (1,2)
 	m_id[EDGE][1][POINT][0] = 1;
 	m_id[EDGE][1][POINT][1] = 2;
	// edge 2 = (2,0)
 	m_id[EDGE][2][POINT][0] = 2;
 	m_id[EDGE][2][POINT][1] = 0;
	// edge 3 = (0,3)
 	m_id[EDGE][3][POINT][0] = 0;
 	m_id[EDGE][3][POINT][1] = 3;
	// edge 4 = (1,3)
 	m_id[EDGE][4][POINT][0] = 1;
 	m_id[EDGE][4][POINT][1] = 3;
	// edge 5 = (2,3)
 	m_id[EDGE][5][POINT][0] = 2;
 	m_id[EDGE][5][POINT][1] = 3;

 	// Edges of Point
 	m_id[POINT][0][EDGE][0] = 2;
 	m_id[POINT][0][EDGE][1] = 0;
 	m_id[POINT][0][EDGE][2] = 3;

 	m_id[POINT][1][EDGE][0] = 0;
 	m_id[POINT][1][EDGE][1] = 1;
 	m_id[POINT][1][EDGE][2] = 4;

 	m_id[POINT][2][EDGE][0] = 1;
 	m_id[POINT][2][EDGE][1] = 2;
 	m_id[POINT][2][EDGE][2] = 5;

 	m_id[POINT][3][EDGE][0] = 3;
 	m_id[POINT][3][EDGE][1] = 4;
 	m_id[POINT][3][EDGE][2] = 5;

	// Reference Corners
 	m_vCorner[0] = MathVector<dim>(0.0, 0.0, 0.0);
 	m_vCorner[1] = MathVector<dim>(1.0, 0.0, 0.0);
 	m_vCorner[2] = MathVector<dim>(0.0, 1.0, 0.0);
 	m_vCorner[3] = MathVector<dim>(0.0, 0.0, 1.0);

	m_vCoInt[0] = MathVector<dim,int>(0, 0, 0);
 	m_vCoInt[1] = MathVector<dim,int>(1, 0, 0);
 	m_vCoInt[2] = MathVector<dim,int>(0, 1, 0);
 	m_vCoInt[3] = MathVector<dim,int>(0, 0, 1);

 	// Reference Element Types
 	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
 	{
		m_vNumRefElem[i] = 0;
 	}
 	m_vNumRefElem[ROID_VERTEX] = 4;
 	m_vNumRefElem[ROID_EDGE] = 6;
 	m_vNumRefElem[ROID_TRIANGLE] = 4;
 	m_vNumRefElem[ROID_TETRAHEDRON] = 1;
}

///////////////////////////////////////////////////////////////////////////////
//	ReferencePyramid
///////////////////////////////////////////////////////////////////////////////

ReferencePyramid::ReferencePyramid()
{
	// dimension
	m_dim = 3;

	// size
	m_size = 1.0/3.0;

	//number of Geometric Objects
 	m_vNum[POINT] = 5;
 	m_vNum[EDGE] = 8;
 	m_vNum[FACE] = 5;
 	m_vNum[VOLUME] = 1;

	// number of Geometric Objects
 	m_vSubNum[VOLUME][0][POINT] = 5;
 	m_vSubNum[VOLUME][0][EDGE] = 8;
 	m_vSubNum[VOLUME][0][FACE] = 5;
 	m_vSubNum[VOLUME][0][VOLUME] = 1;
	m_vRefElemType[VOLUME][0] = ROID_PYRAMID;

	m_vSubNum[FACE][0][POINT] = 4;
	m_vSubNum[FACE][0][EDGE] = 4;
	m_vSubNum[FACE][0][FACE] = 1;
	m_vSubNum[FACE][0][VOLUME] = 1;
	m_vRefElemType[FACE][0] = ROID_QUADRILATERAL;

	for(size_t i = 1; i < m_vNum[FACE]; ++i)
 	{
		m_vSubNum[FACE][i][POINT] = 3;
		m_vSubNum[FACE][i][EDGE] = 3;
		m_vSubNum[FACE][i][FACE] = 1;
		m_vSubNum[FACE][i][VOLUME] = 1;

		m_vRefElemType[FACE][i] = ROID_TRIANGLE;
 	}

 	for(size_t i = 0; i < m_vNum[EDGE]; ++i)
 	{
	 	m_vSubNum[EDGE][i][POINT] = 2;
	 	m_vSubNum[EDGE][i][EDGE] = 1;
	 	m_vSubNum[EDGE][i][FACE] = 2;
	 	m_vSubNum[EDGE][i][VOLUME] = 1;

	 	m_vRefElemType[EDGE][i] = ROID_EDGE;
 	}

 	for(size_t i = 0; i < m_vNum[POINT] - 1; ++i)
 	{
	 	m_vSubNum[POINT][i][POINT] = 1;
	 	m_vSubNum[POINT][i][EDGE] = 3;
	 	m_vSubNum[POINT][i][FACE] = 3;
	 	m_vSubNum[POINT][i][VOLUME] = 1;

	 	m_vRefElemType[POINT][i] = ROID_VERTEX;
 	}
 	m_vSubNum[POINT][4][POINT] = 1;
 	m_vSubNum[POINT][4][EDGE] = 4;
 	m_vSubNum[POINT][4][FACE] = 4;
 	m_vSubNum[POINT][4][VOLUME] = 1;

 	m_vRefElemType[POINT][4] = ROID_VERTEX;

	//reset m_id to -1
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<MAXOBJECTS; ++j)
			for(int k=0; k<=dim; ++k)
				for(size_t l=0; l<MAXOBJECTS; l++)
				{
				 	m_id[i][j][k][l] = -1;
				}

	//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<m_vNum[i]; ++j)
		{
		 	m_id[i][j][i][0] = j;
		}

	// Face <-> Volume
	for(size_t i=0; i<m_vNum[FACE]; ++i)
	{
	 	m_id[VOLUME][0][FACE][i] = i;
	 	m_id[FACE][i][VOLUME][0] = 0;
	}

	// Edge <-> Volume
	for(size_t i=0; i<m_vNum[EDGE]; ++i)
	{
	 	m_id[VOLUME][0][EDGE][i] = i;
	 	m_id[EDGE][i][VOLUME][0] = 0;
	}

	// Point <-> Volume
	for(size_t i=0; i<m_vNum[POINT]; ++i)
	{
	 	m_id[VOLUME][0][POINT][i] = i;
	 	m_id[POINT][i][VOLUME][0] = 0;
	}

	// Points <-> Faces
 	m_id[FACE][0][POINT][0] = 0;
 	m_id[FACE][0][POINT][1] = 3;
 	m_id[FACE][0][POINT][2] = 2;
 	m_id[FACE][0][POINT][3] = 1;

 	m_id[FACE][1][POINT][0] = 0;
 	m_id[FACE][1][POINT][1] = 1;
 	m_id[FACE][1][POINT][2] = 4;

 	m_id[FACE][2][POINT][0] = 1;
 	m_id[FACE][2][POINT][1] = 2;
 	m_id[FACE][2][POINT][2] = 4;

 	m_id[FACE][3][POINT][0] = 2;
 	m_id[FACE][3][POINT][1] = 3;
 	m_id[FACE][3][POINT][2] = 4;

 	m_id[FACE][4][POINT][0] = 3;
 	m_id[FACE][4][POINT][1] = 0;
 	m_id[FACE][4][POINT][2] = 4;

 	m_id[POINT][0][FACE][0] = 0;
 	m_id[POINT][0][FACE][1] = 4;//1;
 	m_id[POINT][0][FACE][2] = 1;//4;

 	m_id[POINT][1][FACE][0] = 0;
 	m_id[POINT][1][FACE][1] = 1;
 	m_id[POINT][1][FACE][2] = 2;

 	m_id[POINT][2][FACE][0] = 0;
 	m_id[POINT][2][FACE][1] = 2;
 	m_id[POINT][2][FACE][2] = 3;

 	m_id[POINT][3][FACE][0] = 0;
 	m_id[POINT][3][FACE][1] = 3;
 	m_id[POINT][3][FACE][2] = 4;

 	m_id[POINT][4][FACE][0] = 1;
 	m_id[POINT][4][FACE][1] = 2;
 	m_id[POINT][4][FACE][2] = 3;
 	m_id[POINT][4][FACE][3] = 4;

 	// Edges <-> Faces
 	m_id[FACE][0][EDGE][0] = 3;
 	m_id[FACE][0][EDGE][1] = 2;
 	m_id[FACE][0][EDGE][2] = 1;
 	m_id[FACE][0][EDGE][3] = 0;

 	m_id[FACE][1][EDGE][0] = 0;
 	m_id[FACE][1][EDGE][1] = 5;
 	m_id[FACE][1][EDGE][2] = 4;

 	m_id[FACE][2][EDGE][0] = 1;
 	m_id[FACE][2][EDGE][1] = 6;
 	m_id[FACE][2][EDGE][2] = 5;

 	m_id[FACE][3][EDGE][0] = 2;
 	m_id[FACE][3][EDGE][1] = 7;
 	m_id[FACE][3][EDGE][2] = 6;

 	m_id[FACE][4][EDGE][0] = 3;
 	m_id[FACE][4][EDGE][1] = 4;
 	m_id[FACE][4][EDGE][2] = 7;

 	m_id[EDGE][0][FACE][0] = 0;
 	m_id[EDGE][0][FACE][1] = 1;

 	m_id[EDGE][1][FACE][0] = 0;
 	m_id[EDGE][1][FACE][1] = 2;

 	m_id[EDGE][2][FACE][0] = 0;
 	m_id[EDGE][2][FACE][1] = 3;

 	m_id[EDGE][3][FACE][0] = 0;
 	m_id[EDGE][3][FACE][1] = 4;

 	m_id[EDGE][4][FACE][0] = 1;
 	m_id[EDGE][4][FACE][1] = 4;

 	m_id[EDGE][5][FACE][0] = 2;
 	m_id[EDGE][5][FACE][1] = 1;

 	m_id[EDGE][6][FACE][0] = 3;
 	m_id[EDGE][6][FACE][1] = 2;

 	m_id[EDGE][7][FACE][0] = 4;
 	m_id[EDGE][7][FACE][1] = 3;

 	// Points of Edges
	// edge 0 = (0,1)
 	m_id[EDGE][0][POINT][0] = 0;
 	m_id[EDGE][0][POINT][1] = 1;
 	// edge 1 = (1,2)
 	m_id[EDGE][1][POINT][0] = 1;
 	m_id[EDGE][1][POINT][1] = 2;
	// edge 2 = (2,3)
 	m_id[EDGE][2][POINT][0] = 2;
 	m_id[EDGE][2][POINT][1] = 3;
	// edge 3 = (3,0)
 	m_id[EDGE][3][POINT][0] = 3;
 	m_id[EDGE][3][POINT][1] = 0;
	// edge 4 = (0,4)
 	m_id[EDGE][4][POINT][0] = 0;
 	m_id[EDGE][4][POINT][1] = 4;
	// edge 5 = (1,4)
 	m_id[EDGE][5][POINT][0] = 1;
 	m_id[EDGE][5][POINT][1] = 4;
	// edge 6 = (2,4)
 	m_id[EDGE][6][POINT][0] = 2;
 	m_id[EDGE][6][POINT][1] = 4;
	// edge 7 = (3,4)
 	m_id[EDGE][7][POINT][0] = 3;
 	m_id[EDGE][7][POINT][1] = 4;

 	// Edges of Point
 	m_id[POINT][0][EDGE][0] = 3;
 	m_id[POINT][0][EDGE][1] = 0;
 	m_id[POINT][0][EDGE][2] = 4;

 	m_id[POINT][1][EDGE][0] = 0;
 	m_id[POINT][1][EDGE][1] = 1;//5;
 	m_id[POINT][1][EDGE][2] = 5;//1;

 	m_id[POINT][2][EDGE][0] = 1;
 	m_id[POINT][2][EDGE][1] = 2;//6;
 	m_id[POINT][2][EDGE][2] = 6;//2;

 	m_id[POINT][3][EDGE][0] = 2;
 	m_id[POINT][3][EDGE][1] = 3;//7;
 	m_id[POINT][3][EDGE][2] = 7;//3;

 	m_id[POINT][4][EDGE][0] = 4;
 	m_id[POINT][4][EDGE][1] = 5;
 	m_id[POINT][4][EDGE][2] = 6;
 	m_id[POINT][4][EDGE][3] = 7;

 	// Reference Corners
 	m_vCorner[0] = MathVector<dim>(0.0, 0.0, 0.0);
 	m_vCorner[1] = MathVector<dim>(1.0, 0.0, 0.0);
 	m_vCorner[2] = MathVector<dim>(1.0, 1.0, 0.0);
 	m_vCorner[3] = MathVector<dim>(0.0, 1.0, 0.0);
 	m_vCorner[4] = MathVector<dim>(0.0, 0.0, 1.0);

 	m_vCoInt[0] = MathVector<dim,int>(0, 0, 0);
 	m_vCoInt[1] = MathVector<dim,int>(1, 0, 0);
 	m_vCoInt[2] = MathVector<dim,int>(1, 1, 0);
 	m_vCoInt[3] = MathVector<dim,int>(0, 1, 0);
 	m_vCoInt[4] = MathVector<dim,int>(0, 0, 1);

	// Reference Element Types
 	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
 	{
		m_vNumRefElem[i] = 0;
 	}
 	m_vNumRefElem[ROID_VERTEX] = 5;
 	m_vNumRefElem[ROID_EDGE] = 18;
 	m_vNumRefElem[ROID_QUADRILATERAL] = 1;
 	m_vNumRefElem[ROID_TRIANGLE] = 4;
 	m_vNumRefElem[ROID_PYRAMID] = 1;
}

///////////////////////////////////////////////////////////////////////////////
//	ReferencePrism
///////////////////////////////////////////////////////////////////////////////

ReferencePrism::ReferencePrism()
{
	// dimension
	m_dim = 3;

	// size
	m_size = 0.5;

	//number of Geometric Objects
 	m_vNum[POINT] = 6;
 	m_vNum[EDGE] = 9;
 	m_vNum[FACE] = 5;
 	m_vNum[VOLUME] = 1;

	// number of Geometric Objects
 	m_vSubNum[VOLUME][0][POINT] = 6;
 	m_vSubNum[VOLUME][0][EDGE] = 9;
 	m_vSubNum[VOLUME][0][FACE] = 5;
 	m_vSubNum[VOLUME][0][VOLUME] = 1;
	m_vRefElemType[VOLUME][0] = ROID_PRISM;

 	m_vSubNum[FACE][0][POINT] = 3;
 	m_vSubNum[FACE][0][EDGE] = 3;
 	m_vSubNum[FACE][0][FACE] = 1;
 	m_vSubNum[FACE][0][VOLUME] = 1;
	m_vRefElemType[FACE][0] = ROID_TRIANGLE;

 	m_vSubNum[FACE][1][POINT] = 4;
 	m_vSubNum[FACE][1][EDGE] = 4;
 	m_vSubNum[FACE][1][FACE] = 1;
 	m_vSubNum[FACE][1][VOLUME] = 1;
	m_vRefElemType[FACE][1] = ROID_QUADRILATERAL;

 	m_vSubNum[FACE][2][POINT] = 4;
 	m_vSubNum[FACE][2][EDGE] = 4;
 	m_vSubNum[FACE][2][FACE] = 1;
 	m_vSubNum[FACE][2][VOLUME] = 1;
	m_vRefElemType[FACE][2] = ROID_QUADRILATERAL;

 	m_vSubNum[FACE][3][POINT] = 4;
 	m_vSubNum[FACE][3][EDGE] = 4;
 	m_vSubNum[FACE][3][FACE] = 1;
 	m_vSubNum[FACE][3][VOLUME] = 1;
	m_vRefElemType[FACE][3] = ROID_QUADRILATERAL;

 	m_vSubNum[FACE][4][POINT] = 3;
 	m_vSubNum[FACE][4][EDGE] = 3;
 	m_vSubNum[FACE][4][FACE] = 1;
 	m_vSubNum[FACE][4][VOLUME] = 1;
	m_vRefElemType[FACE][4] = ROID_TRIANGLE;

 	for(size_t i = 0; i < m_vNum[EDGE]; ++i)
 	{
	 	m_vSubNum[EDGE][i][POINT] = 2;
	 	m_vSubNum[EDGE][i][EDGE] = 1;
	 	m_vSubNum[EDGE][i][FACE] = 2;
	 	m_vSubNum[EDGE][i][VOLUME] = 1;

	 	m_vRefElemType[EDGE][i] = ROID_EDGE;
 	}

 	for(size_t i = 0; i < m_vNum[POINT]; ++i)
 	{
	 	m_vSubNum[POINT][i][POINT] = 1;
	 	m_vSubNum[POINT][i][EDGE] = 3;
	 	m_vSubNum[POINT][i][FACE] = 3;
	 	m_vSubNum[POINT][i][VOLUME] = 1;

	 	m_vRefElemType[POINT][i] = ROID_VERTEX;
 	}

	//reset m_id to -1
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<MAXOBJECTS; ++j)
			for(int k=0; k<=dim; ++k)
				for(size_t l=0; l<MAXOBJECTS; l++)
				{
				 	m_id[i][j][k][l] = -1;
				}

	//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<m_vNum[i]; ++j)
		{
		 	m_id[i][j][i][0] = j;
		}

	// Face <-> Volume
	for(size_t i=0; i<m_vNum[FACE]; ++i)
	{
	 	m_id[VOLUME][0][FACE][i] = i;
	 	m_id[FACE][i][VOLUME][0] = 0;
	}

	// Edge <-> Volume
	for(size_t i=0; i<m_vNum[EDGE]; ++i)
	{
	 	m_id[VOLUME][0][EDGE][i] = i;
	 	m_id[EDGE][i][VOLUME][0] = 0;
	}

	// Point <-> Volume
	for(size_t i=0; i<m_vNum[POINT]; ++i)
	{
	 	m_id[VOLUME][0][POINT][i] = i;
	 	m_id[POINT][i][VOLUME][0] = 0;
	}

	// Points <-> Faces
 	m_id[FACE][0][POINT][0] = 0;
 	m_id[FACE][0][POINT][1] = 2;
 	m_id[FACE][0][POINT][2] = 1;

 	m_id[FACE][1][POINT][0] = 0;
 	m_id[FACE][1][POINT][1] = 1;
 	m_id[FACE][1][POINT][2] = 4;
 	m_id[FACE][1][POINT][3] = 3;

 	m_id[FACE][2][POINT][0] = 1;
 	m_id[FACE][2][POINT][1] = 2;
 	m_id[FACE][2][POINT][2] = 5;
 	m_id[FACE][2][POINT][3] = 4;

 	m_id[FACE][3][POINT][0] = 2;
 	m_id[FACE][3][POINT][1] = 0;
 	m_id[FACE][3][POINT][2] = 3;
 	m_id[FACE][3][POINT][3] = 5;

 	m_id[FACE][4][POINT][0] = 3;
 	m_id[FACE][4][POINT][1] = 4;
 	m_id[FACE][4][POINT][2] = 5;


 	m_id[POINT][0][FACE][0] = 0;
 	m_id[POINT][0][FACE][1] = 3;
 	m_id[POINT][0][FACE][2] = 1;

 	m_id[POINT][1][FACE][0] = 0;
 	m_id[POINT][1][FACE][1] = 1;
 	m_id[POINT][1][FACE][2] = 2;

 	m_id[POINT][2][FACE][0] = 0;
 	m_id[POINT][2][FACE][1] = 2;
 	m_id[POINT][2][FACE][2] = 3;

 	m_id[POINT][3][FACE][0] = 1;
 	m_id[POINT][3][FACE][1] = 3;
 	m_id[POINT][3][FACE][2] = 4;

 	m_id[POINT][4][FACE][0] = 2;
 	m_id[POINT][4][FACE][1] = 1;
 	m_id[POINT][4][FACE][2] = 4;

 	m_id[POINT][5][FACE][0] = 3;
 	m_id[POINT][5][FACE][1] = 2;
 	m_id[POINT][5][FACE][2] = 4;

 	// Edges <-> Faces
 	m_id[FACE][0][EDGE][0] = 2;
 	m_id[FACE][0][EDGE][1] = 1;
 	m_id[FACE][0][EDGE][2] = 0;

 	m_id[FACE][1][EDGE][0] = 0;
 	m_id[FACE][1][EDGE][1] = 4;
 	m_id[FACE][1][EDGE][2] = 6;
 	m_id[FACE][1][EDGE][3] = 3;

 	m_id[FACE][2][EDGE][0] = 1;
 	m_id[FACE][2][EDGE][1] = 5;
 	m_id[FACE][2][EDGE][2] = 7;
 	m_id[FACE][2][EDGE][3] = 4;

 	m_id[FACE][3][EDGE][0] = 2;
 	m_id[FACE][3][EDGE][1] = 3;
 	m_id[FACE][3][EDGE][2] = 8;
 	m_id[FACE][3][EDGE][3] = 5;

 	m_id[FACE][4][EDGE][0] = 6;
 	m_id[FACE][4][EDGE][1] = 7;
 	m_id[FACE][4][EDGE][2] = 8;

 	m_id[EDGE][0][FACE][0] = 0;
 	m_id[EDGE][0][FACE][1] = 1;

 	m_id[EDGE][1][FACE][0] = 0;
 	m_id[EDGE][1][FACE][1] = 2;

 	m_id[EDGE][2][FACE][0] = 0;
 	m_id[EDGE][2][FACE][1] = 3;

 	m_id[EDGE][3][FACE][0] = 1;
 	m_id[EDGE][3][FACE][1] = 3;

 	m_id[EDGE][4][FACE][0] = 2;
 	m_id[EDGE][4][FACE][1] = 1;

 	m_id[EDGE][5][FACE][0] = 3;
 	m_id[EDGE][5][FACE][1] = 2;

 	m_id[EDGE][6][FACE][0] = 1;
 	m_id[EDGE][6][FACE][1] = 4;

 	m_id[EDGE][7][FACE][0] = 2;
 	m_id[EDGE][7][FACE][1] = 4;

 	m_id[EDGE][8][FACE][0] = 3;
 	m_id[EDGE][8][FACE][1] = 4;

	// Points of Edges
	// edge 0 = (0,1)
 	m_id[EDGE][0][POINT][0] = 0;
 	m_id[EDGE][0][POINT][1] = 1;
 	// edge 1 = (1,2)
 	m_id[EDGE][1][POINT][0] = 1;
 	m_id[EDGE][1][POINT][1] = 2;
	// edge 2 = (2,0)
 	m_id[EDGE][2][POINT][0] = 2;
 	m_id[EDGE][2][POINT][1] = 0;
	// edge 3 = (0,3)
 	m_id[EDGE][3][POINT][0] = 0;
 	m_id[EDGE][3][POINT][1] = 3;
	// edge 4 = (1,3)
 	m_id[EDGE][4][POINT][0] = 1;
 	m_id[EDGE][4][POINT][1] = 4;
	// edge 5 = (2,3)
 	m_id[EDGE][5][POINT][0] = 2;
 	m_id[EDGE][5][POINT][1] = 5;
	// edge 6 = (3,4)
 	m_id[EDGE][6][POINT][0] = 3;
 	m_id[EDGE][6][POINT][1] = 4;
	// edge 7 = (4,5)
 	m_id[EDGE][7][POINT][0] = 4;
 	m_id[EDGE][7][POINT][1] = 5;
	// edge 8 = (5,3)
 	m_id[EDGE][8][POINT][0] = 5;
 	m_id[EDGE][8][POINT][1] = 3;

 	// Edges of Point
 	m_id[POINT][0][EDGE][0] = 2;
 	m_id[POINT][0][EDGE][1] = 0;
 	m_id[POINT][0][EDGE][2] = 3;

 	m_id[POINT][1][EDGE][0] = 0;
 	m_id[POINT][1][EDGE][1] = 1;
 	m_id[POINT][1][EDGE][2] = 4;

 	m_id[POINT][2][EDGE][0] = 1;
 	m_id[POINT][2][EDGE][1] = 2;
 	m_id[POINT][2][EDGE][2] = 5;

 	m_id[POINT][3][EDGE][0] = 3;
 	m_id[POINT][3][EDGE][1] = 6;
 	m_id[POINT][3][EDGE][2] = 8;

 	m_id[POINT][4][EDGE][0] = 4;
 	m_id[POINT][4][EDGE][1] = 7;
 	m_id[POINT][4][EDGE][2] = 6;

 	m_id[POINT][5][EDGE][0] = 5;
 	m_id[POINT][5][EDGE][1] = 8;
 	m_id[POINT][5][EDGE][2] = 7;

 	// Reference Corners
 	m_vCorner[0] = MathVector<dim>(0.0, 0.0, 0.0);
 	m_vCorner[1] = MathVector<dim>(1.0, 0.0, 0.0);
 	m_vCorner[2] = MathVector<dim>(0.0, 1.0, 0.0);
 	m_vCorner[3] = MathVector<dim>(0.0, 0.0, 1.0);
 	m_vCorner[4] = MathVector<dim>(1.0, 0.0, 1.0);
 	m_vCorner[5] = MathVector<dim>(0.0, 1.0, 1.0);

 	m_vCoInt[0] = MathVector<dim,int>(0, 0, 0);
 	m_vCoInt[1] = MathVector<dim,int>(1, 0, 0);
 	m_vCoInt[2] = MathVector<dim,int>(0, 1, 0);
 	m_vCoInt[3] = MathVector<dim,int>(0, 0, 1);
 	m_vCoInt[4] = MathVector<dim,int>(1, 0, 1);
 	m_vCoInt[5] = MathVector<dim,int>(0, 1, 1);

 	// Reference Element Types
 	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
 	{
		m_vNumRefElem[i] = 0;
 	}
 	m_vNumRefElem[ROID_VERTEX] = 6;
 	m_vNumRefElem[ROID_EDGE] = 9;
 	m_vNumRefElem[ROID_TRIANGLE] = 2;
 	m_vNumRefElem[ROID_QUADRILATERAL] = 3;
 	m_vNumRefElem[ROID_PRISM] = 1;
}

///////////////////////////////////////////////////////////////////////////////
//	ReferenceHexahedron
///////////////////////////////////////////////////////////////////////////////

ReferenceHexahedron::ReferenceHexahedron()
{
	// dimension
	m_dim = 3;

	// size
	m_size = 1.0;

	//number of Geometric Objects
 	m_vNum[POINT] = 8;
 	m_vNum[EDGE] = 12;
 	m_vNum[FACE] = 6;
 	m_vNum[VOLUME] = 1;

	// number of Geometric Objects
 	m_vSubNum[VOLUME][0][POINT] = 8;
 	m_vSubNum[VOLUME][0][EDGE] = 12;
 	m_vSubNum[VOLUME][0][FACE] = 6;
 	m_vSubNum[VOLUME][0][VOLUME] = 1;
	m_vRefElemType[VOLUME][0] = ROID_HEXAHEDRON;

 	for(size_t i = 0; i < m_vNum[FACE]; ++i)
 	{
		m_vSubNum[FACE][i][POINT] = 4;
		m_vSubNum[FACE][i][EDGE] = 4;
		m_vSubNum[FACE][i][FACE] = 1;
		m_vSubNum[FACE][i][VOLUME] = 1;

		m_vRefElemType[FACE][i] = ROID_QUADRILATERAL;
 	}

 	for(size_t i = 0; i < m_vNum[EDGE]; ++i)
 	{
	 	m_vSubNum[EDGE][i][POINT] = 2;
	 	m_vSubNum[EDGE][i][EDGE] = 1;
	 	m_vSubNum[EDGE][i][FACE] = 2;
	 	m_vSubNum[EDGE][i][VOLUME] = 1;

	 	m_vRefElemType[EDGE][i] = ROID_EDGE;
 	}

 	for(size_t i = 0; i < m_vNum[POINT]; ++i)
 	{
	 	m_vSubNum[POINT][i][POINT] = 1;
	 	m_vSubNum[POINT][i][EDGE] = 3;
	 	m_vSubNum[POINT][i][FACE] = 3;
	 	m_vSubNum[POINT][i][VOLUME] = 1;

	 	m_vRefElemType[POINT][i] = ROID_VERTEX;
 	}

	//reset m_id to -1
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<MAXOBJECTS; ++j)
			for(int k=0; k<=dim; ++k)
				for(size_t l=0; l<MAXOBJECTS; l++)
				{
				 	m_id[i][j][k][l] = -1;
				}

	//self references: (i.e. Point <-> Point, Edge <-> Edge, etc.)
	for(int i=0; i<=dim; ++i)
		for(size_t j=0; j<m_vNum[i]; ++j)
		{
		 	m_id[i][j][i][0] = j;
		}

	// Face <-> Volume
	for(size_t i=0; i<m_vNum[FACE]; ++i)
	{
	 	m_id[VOLUME][0][FACE][i] = i;
	 	m_id[FACE][i][VOLUME][0] = 0;
	}

	// Edge <-> Volume
	for(size_t i=0; i<m_vNum[EDGE]; ++i)
	{
	 	m_id[VOLUME][0][EDGE][i] = i;
	 	m_id[EDGE][i][VOLUME][0] = 0;
	}

	// Point <-> Volume
	for(size_t i=0; i<m_vNum[POINT]; ++i)
	{
	 	m_id[VOLUME][0][POINT][i] = i;
	 	m_id[POINT][i][VOLUME][0] = 0;
	}

	// Points <-> Faces
 	m_id[FACE][0][POINT][0] = 0;
 	m_id[FACE][0][POINT][1] = 3;
 	m_id[FACE][0][POINT][2] = 2;
 	m_id[FACE][0][POINT][3] = 1;

 	m_id[FACE][1][POINT][0] = 0;
 	m_id[FACE][1][POINT][1] = 1;
 	m_id[FACE][1][POINT][2] = 5;
 	m_id[FACE][1][POINT][3] = 4;

 	m_id[FACE][2][POINT][0] = 1;
 	m_id[FACE][2][POINT][1] = 2;
 	m_id[FACE][2][POINT][2] = 6;
 	m_id[FACE][2][POINT][3] = 5;

 	m_id[FACE][3][POINT][0] = 2;
 	m_id[FACE][3][POINT][1] = 3;
 	m_id[FACE][3][POINT][2] = 7;
 	m_id[FACE][3][POINT][3] = 6;

 	m_id[FACE][4][POINT][0] = 3;
 	m_id[FACE][4][POINT][1] = 0;
 	m_id[FACE][4][POINT][2] = 4;
 	m_id[FACE][4][POINT][3] = 7;

 	m_id[FACE][5][POINT][0] = 4;
 	m_id[FACE][5][POINT][1] = 5;
 	m_id[FACE][5][POINT][2] = 6;
 	m_id[FACE][5][POINT][3] = 7;

 	m_id[POINT][0][FACE][0] = 0;
 	m_id[POINT][0][FACE][1] = 4;
 	m_id[POINT][0][FACE][2] = 1;

 	m_id[POINT][1][FACE][0] = 0;
 	m_id[POINT][1][FACE][1] = 1;
 	m_id[POINT][1][FACE][2] = 2;

 	m_id[POINT][2][FACE][0] = 0;
 	m_id[POINT][2][FACE][1] = 2;
 	m_id[POINT][2][FACE][2] = 3;

 	m_id[POINT][3][FACE][0] = 0;
 	m_id[POINT][3][FACE][1] = 3;
 	m_id[POINT][3][FACE][2] = 4;

 	m_id[POINT][4][FACE][0] = 1;
 	m_id[POINT][4][FACE][1] = 4;
 	m_id[POINT][4][FACE][2] = 5;

 	m_id[POINT][5][FACE][0] = 2;
 	m_id[POINT][5][FACE][1] = 1;
 	m_id[POINT][5][FACE][2] = 5;

 	m_id[POINT][6][FACE][0] = 3;
 	m_id[POINT][6][FACE][1] = 2;
 	m_id[POINT][6][FACE][2] = 5;

 	m_id[POINT][7][FACE][0] = 4;
 	m_id[POINT][7][FACE][1] = 3;
 	m_id[POINT][7][FACE][2] = 5;

 	// Edges <-> Faces
 	m_id[FACE][0][EDGE][0] = 3;
 	m_id[FACE][0][EDGE][1] = 2;
 	m_id[FACE][0][EDGE][2] = 1;
 	m_id[FACE][0][EDGE][3] = 0;

 	m_id[FACE][1][EDGE][0] = 0;
 	m_id[FACE][1][EDGE][1] = 5;
 	m_id[FACE][1][EDGE][2] = 8;
 	m_id[FACE][1][EDGE][3] = 4;

 	m_id[FACE][2][EDGE][0] = 1;
 	m_id[FACE][2][EDGE][1] = 6;
 	m_id[FACE][2][EDGE][2] = 9;
 	m_id[FACE][2][EDGE][3] = 5;

 	m_id[FACE][3][EDGE][0] = 2;
 	m_id[FACE][3][EDGE][1] = 7;
 	m_id[FACE][3][EDGE][2] = 10;
 	m_id[FACE][3][EDGE][3] = 6;

 	m_id[FACE][4][EDGE][0] = 3;
 	m_id[FACE][4][EDGE][1] = 4;
 	m_id[FACE][4][EDGE][2] = 11;
 	m_id[FACE][4][EDGE][3] = 7;

 	m_id[FACE][5][EDGE][0] = 8;
 	m_id[FACE][5][EDGE][1] = 9;
 	m_id[FACE][5][EDGE][2] = 10;
 	m_id[FACE][5][EDGE][3] = 11;

 	m_id[EDGE][0][FACE][0] = 0;
 	m_id[EDGE][0][FACE][1] = 1;

 	m_id[EDGE][1][FACE][0] = 0;
 	m_id[EDGE][1][FACE][1] = 2;

 	m_id[EDGE][2][FACE][0] = 0;
 	m_id[EDGE][2][FACE][1] = 3;

 	m_id[EDGE][3][FACE][0] = 0;
 	m_id[EDGE][3][FACE][1] = 4;

 	m_id[EDGE][4][FACE][0] = 1;
 	m_id[EDGE][4][FACE][1] = 4;

 	m_id[EDGE][5][FACE][0] = 2;
 	m_id[EDGE][5][FACE][1] = 1;

 	m_id[EDGE][6][FACE][0] = 3;
 	m_id[EDGE][6][FACE][1] = 2;

 	m_id[EDGE][7][FACE][0] = 4;
 	m_id[EDGE][7][FACE][1] = 3;

 	m_id[EDGE][8][FACE][0] = 1;
 	m_id[EDGE][8][FACE][1] = 5;

 	m_id[EDGE][9][FACE][0] = 2;
 	m_id[EDGE][9][FACE][1] = 5;

 	m_id[EDGE][10][FACE][0] = 3;
 	m_id[EDGE][10][FACE][1] = 5;

 	m_id[EDGE][11][FACE][0] = 4;
 	m_id[EDGE][11][FACE][1] = 5;

 	// Points of Edges
	// edge 0 = (0,1)
 	m_id[EDGE][0][POINT][0] = 0;
 	m_id[EDGE][0][POINT][1] = 1;
 	// edge 1 = (1,2)
 	m_id[EDGE][1][POINT][0] = 1;
 	m_id[EDGE][1][POINT][1] = 2;
	// edge 2 = (2,3)
 	m_id[EDGE][2][POINT][0] = 2;
 	m_id[EDGE][2][POINT][1] = 3;
	// edge 3 = (3,0)
 	m_id[EDGE][3][POINT][0] = 3;
 	m_id[EDGE][3][POINT][1] = 0;
	// edge 4 = (0,4)
 	m_id[EDGE][4][POINT][0] = 0;
 	m_id[EDGE][4][POINT][1] = 4;
	// edge 5 = (1,5)
 	m_id[EDGE][5][POINT][0] = 1;
 	m_id[EDGE][5][POINT][1] = 5;
	// edge 6 = (2,6)
 	m_id[EDGE][6][POINT][0] = 2;
 	m_id[EDGE][6][POINT][1] = 6;
	// edge 7 = (3,7)
 	m_id[EDGE][7][POINT][0] = 3;
 	m_id[EDGE][7][POINT][1] = 7;
	// edge 8 = (4,5)
 	m_id[EDGE][8][POINT][0] = 4;
 	m_id[EDGE][8][POINT][1] = 5;
	// edge 9 = (5,6)
 	m_id[EDGE][9][POINT][0] = 5;
 	m_id[EDGE][9][POINT][1] = 6;
	// edge 10 = (6,7)
 	m_id[EDGE][10][POINT][0] = 6;
 	m_id[EDGE][10][POINT][1] = 7;
	// edge 11 = (7,4)
 	m_id[EDGE][11][POINT][0] = 7;
 	m_id[EDGE][11][POINT][1] = 4;

 	// Edges of Point
 	m_id[POINT][0][EDGE][0] = 3;
 	m_id[POINT][0][EDGE][1] = 0;
 	m_id[POINT][0][EDGE][2] = 4;

 	m_id[POINT][1][EDGE][0] = 0;
 	m_id[POINT][1][EDGE][1] = 1;
 	m_id[POINT][1][EDGE][2] = 5;

 	m_id[POINT][2][EDGE][0] = 1;
 	m_id[POINT][2][EDGE][1] = 2;
 	m_id[POINT][2][EDGE][2] = 6;

 	m_id[POINT][3][EDGE][0] = 2;
 	m_id[POINT][3][EDGE][1] = 3;
 	m_id[POINT][3][EDGE][2] = 7;

 	m_id[POINT][4][EDGE][0] = 4;
 	m_id[POINT][4][EDGE][1] = 8;
 	m_id[POINT][4][EDGE][2] = 11;

 	m_id[POINT][5][EDGE][0] = 5;
 	m_id[POINT][5][EDGE][1] = 9;
 	m_id[POINT][5][EDGE][2] = 8;

 	m_id[POINT][6][EDGE][0] = 6;
 	m_id[POINT][6][EDGE][1] = 10;
 	m_id[POINT][6][EDGE][2] = 9;

 	m_id[POINT][7][EDGE][0] = 7;
 	m_id[POINT][7][EDGE][1] = 11;
 	m_id[POINT][7][EDGE][2] = 10;

 	// Reference Corners
 	m_vCorner[0] = MathVector<dim>(0.0, 0.0, 0.0);
 	m_vCorner[1] = MathVector<dim>(1.0, 0.0, 0.0);
 	m_vCorner[2] = MathVector<dim>(1.0, 1.0, 0.0);
 	m_vCorner[3] = MathVector<dim>(0.0, 1.0, 0.0);
 	m_vCorner[4] = MathVector<dim>(0.0, 0.0, 1.0);
 	m_vCorner[5] = MathVector<dim>(1.0, 0.0, 1.0);
 	m_vCorner[6] = MathVector<dim>(1.0, 1.0, 1.0);
 	m_vCorner[7] = MathVector<dim>(0.0, 1.0, 1.0);

	m_vCoInt[0] = MathVector<dim,int>(0, 0, 0);
 	m_vCoInt[1] = MathVector<dim,int>(1, 0, 0);
 	m_vCoInt[2] = MathVector<dim,int>(1, 1, 0);
 	m_vCoInt[3] = MathVector<dim,int>(0, 1, 0);
 	m_vCoInt[4] = MathVector<dim,int>(0, 0, 1);
 	m_vCoInt[5] = MathVector<dim,int>(1, 0, 1);
 	m_vCoInt[6] = MathVector<dim,int>(1, 1, 1);
 	m_vCoInt[7] = MathVector<dim,int>(0, 1, 1);

 	// Reference Element Types
 	for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
 	{
		m_vNumRefElem[i] = 0;
 	}
 	m_vNumRefElem[ROID_VERTEX] = 8;
 	m_vNumRefElem[ROID_EDGE] = 12;
 	m_vNumRefElem[ROID_QUADRILATERAL] = 6;
 	m_vNumRefElem[ROID_HEXAHEDRON] = 1;
}


}; // end namespace ug

