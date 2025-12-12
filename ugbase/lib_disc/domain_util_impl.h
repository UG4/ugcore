/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__DOMAIN_UTIL_GENERAL_IMPL__
#define __H__UG__LIB_DISC__DOMAIN_UTIL_GENERAL_IMPL__

#include "domain_util.h"

#include <string>
#include <sstream>

#include "lib_disc/reference_element/reference_element.h" // for reference_element_traits

namespace ug {

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
	auto* pElem = const_cast<TElem*>(&elem);

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
template <typename TDomain>
void CollectCornerCoordinates(int base_object_id,
			std::vector<typename TDomain::position_type>& vCornerCoordsOut,
			GridObject& elem, const TDomain& domain, bool clearContainer)
{
	switch(base_object_id)
	{
		case VERTEX: CollectCornerCoordinates<Vertex, TDomain>(vCornerCoordsOut, *dynamic_cast< Vertex*>(&elem), domain, clearContainer); return;
		case EDGE: CollectCornerCoordinates<Edge, TDomain>(vCornerCoordsOut, *dynamic_cast< Edge*>(&elem), domain, clearContainer); return;
		case FACE: CollectCornerCoordinates<Face, TDomain>(vCornerCoordsOut, *dynamic_cast< Face*>(&elem), domain, clearContainer); return;
		case VOLUME: CollectCornerCoordinates<Volume, TDomain>(vCornerCoordsOut, *dynamic_cast< Volume*>(&elem), domain, clearContainer); return;
		default: UG_THROW("CollectCornerCoordinates, base_object_id " << base_object_id << "unknown.");
	}
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


///	returns the coordinates of a vertex (specialization)
template <typename TDomain>
void FillCornerCoordinates(	typename TDomain::position_type vCornerCoordsOut[],
                           	const RegularVertex& vtx, const TDomain& domain)
{
	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.position_accessor();

	// write vertex coordinates
	vCornerCoordsOut[0] = aaPos[&vtx];
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
	using TRefElem = typename reference_element_traits<TElem>::reference_element_type;

	// dimension of Positions
	static constexpr int dim = TPosition::Size;

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
number ElementDiameterSq(const TElem& elem, TDomain& domain)
{
	return ElementDiameterSq(*domain.grid(), domain.position_accessor(), &elem);
}

template <typename TElem, typename TDomain>
number ElementDiameter(const TElem& elem, TDomain& domain)
{
	return ElementDiameter(*domain.grid(), domain.position_accessor(), &elem);
}

} // end namespace ug

#endif