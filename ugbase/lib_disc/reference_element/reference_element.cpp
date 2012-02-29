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

		if(!bRes) throw(UGFatalError("Error while registering Reference Elements"));
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

}; // end namespace ug

