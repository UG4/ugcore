//	created by Andreas Vogel

#ifndef __H__UG__LIB_DISC__DOMAIN_UTIL_GENERAL_IMPL__
#define __H__UG__LIB_DISC__DOMAIN_UTIL_GENERAL_IMPL__

#include "domain_util.h"

#include <string>
#include <sstream>

#include "lib_disc/reference_element/reference_element.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// CollectCornerCoordinates
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//	returns the corner coordinates of a geometric object
template <typename TElem, typename TAAPos>
void CollectCornerCoordinates(	std::vector<typename TAAPos::ValueType>& vCornerCoordsOut,
								TElem* elem, const TAAPos& aaPos, bool clearContainer)
{
	if(clearContainer)
		vCornerCoordsOut.clear();

	// number of vertices of element
	const size_t numVertices = NumVertices(elem);

	// loop vertices
	for(size_t i = 0; i < numVertices; ++i)
	{
		// get element
		VertexBase* vert = GetVertex(elem, i);

		// write corner coordinates
		vCornerCoordsOut.push_back(aaPos[vert]);
	}
}

//	returns the corner coordinates of a geometric object
template <typename TElem, typename TAAPos>
void CollectCornerCoordinates(	std::vector<typename TAAPos::ValueType>& vCornerCoordsOut,
								const TElem& elem, const TAAPos& aaPos, bool clearContainer)
{
//	cast constness away
	TElem* pElem = const_cast<TElem*>(&elem);

//	forward
	return CollectCornerCoordinates(vCornerCoordsOut, pElem, aaPos, clearContainer);
}

///	returns the corner coordinates of a geometric object
template <typename TElem, typename TDomain>
void CollectCornerCoordinates(	std::vector<typename TDomain::position_type>& vCornerCoordsOut,
								const TElem& elem, const TDomain& domain, bool clearContainer)
{
	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.position_accessor();

	CollectCornerCoordinates(vCornerCoordsOut, elem, aaPos, clearContainer);
}

////////////////////////////////////////////////////////////////////////
///	returns the size of a geometric object
template <typename TElem, typename TPosition>
number ElementSize(const TElem& elem, const Grid::VertexAttachmentAccessor<Attachment<TPosition> >& aaPos)
{
	// corner coords
	std::vector<TPosition> vCornerCoords;

	// load corner coords
	CollectCornerCoordinates(vCornerCoords, elem, aaPos);

	// get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type TRefElem;

	// dimension of Positions
	static const int dim = TPosition::Size;

	// return Element Size
	return ElementSize<TRefElem, dim>(&vCornerCoords[0]);
}

///	returns the size of a geometric object
template <typename TElem, typename TDomain>
number ElementSize(const TElem& elem, const TDomain& domain)
{
	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.position_accessor();

	return ElementSize(elem, aaPos);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Class ElemGlobCornerCoords
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/** Gets and keeps the corner coordinates of an element
 *
 * This class can be used with the Provider class to get the
 * corner coordinates of an element in an element loop.
 */
template <typename TDomain, typename TElem>
class ElemGlobCornerCoords
{
//	Position type
	typedef typename TDomain::position_type position_type;

// Reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

// The element
	TElem * elem;

// Array of the vertex coordinates
	position_type m_vCornerCoord[ref_elem_type::numCorners];
	
	public:

// Total number of the corners:
	static const size_t numCorners = ref_elem_type::numCorners;
	
// Get and save the positions of the vertices
	void update(TDomain *domain, TElem *the_elem, bool forced = false)
	{
	//	has the element been changed?
		if ((! forced) && the_elem == elem) return; // no
		
		elem = the_elem;
		
	//	get the element accessor
		const typename TDomain::position_accessor_type &aaPos = domain->position_accessor();
		
	//	extract the coordinates
		for(size_t i = 0; i < numCorners; ++i)
			m_vCornerCoord[i] = aaPos[the_elem->vertex(i)];
	};

// Return the element
	const TElem * Element() {return elem;};

//	Return corner coords as a vector
	position_type *vGlobalCorner() {return m_vCornerCoord;};

// Return corner coords of one corner
	position_type &GlobalCorner(size_t i) {return m_vCornerCoord[i];};

// Constructor
	ElemGlobCornerCoords () : elem(0) {};
};

////////////////////////////////////////////////////////////////////////
//	ElementDiameter
////////////////////////////////////////////////////////////////////////

template <typename TElem, typename TDomain>
number ElementDiameterSq(TElem& elem, TDomain& domain)
{
	return ElementDiameterSq(*domain.grid(), domain.position_accessor(), &elem);
}

template <typename TElem, typename TDomain>
number ElementDiameter(TElem& elem, TDomain& domain)
{
	return ElementDiameter(*domain.grid(), domain.position_accessor(), &elem);
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOMAIN_UTIL_GENERAL_IMPL__ */
