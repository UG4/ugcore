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
		Vertex* vert = GetVertex(elem, i);

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

///	returns the corner coordinates of a geometric object
template <typename TElem, typename TDomain>
void FillCornerCoordinates(	typename TDomain::position_type vCornerCoordsOut[],
                           	const TElem& elem, const TDomain& domain)
{
	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.position_accessor();

	const Vertex* const* vVertex = const_cast<TElem*>(&elem)->vertices();

	// write corner coordinates
	for(size_t i = 0; i < TElem::NUM_VERTICES; ++i)
		vCornerCoordsOut[i] = aaPos[vVertex[i]];
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
