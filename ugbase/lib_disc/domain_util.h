/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Andreas Vogel, Jan Friebertshäuser
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

#ifndef __H__UG__LIB_DISC__DOMAIN_UTIL__
#define __H__UG__LIB_DISC__DOMAIN_UTIL__

// other lib_discretization headers
#include "domain.h"

namespace ug {

/// \addtogroup lib_disc_domain
/// \{

///	Loads a domain from a grid-file.
/**	By optionally specifying a procId, you can make sure that the domain is only
 * loaded on one process. Pass -1, if you want to load it on all processes.
 * Note that the procId is only important in parallel environments.
 * \{
 */
template <typename TDomain>
void LoadDomain(TDomain& domain, const char* filename);

template <typename TDomain>
void LoadDomain(TDomain& domain, const char* filename, int procId);
/**	\} */

///	Saves the domain to a grid-file.
template <typename TDomain>
void SaveDomain(TDomain& domain, const char* filename);


////////////////////////////////////////////////////////////////////////
///	returns the corner coordinates of a geometric object
/** Returns the corner coordinated of a geometric object in a vector
 *
 * This function collects the corner coordinates for a given geometric object
 * in the order prescribed by the reference elements
 *
 * \param[out]	vCornerCoordsOut	vector of corner coordinates
 * \param[in] 	elem				Geometric Object
 * \param[in] 	aaPos				AttachmentAccessor for Positions
 * \param[in]	clearContainer		empty container before filling
 */
template <typename TElem, typename TAAPos>
void CollectCornerCoordinates(	std::vector<typename TAAPos::ValueType>& vCornerCoordsOut,
								const TElem& elem, const TAAPos& aaPos,
								bool clearContainer = true);


////////////////////////////////////////////////////////////////////////
///	returns the corner coordinates of a geometric object
/** Returns the corner coordinated of a geometric object in a vector
 *
 * This function collects the corner coordinates for a given geometric object
 * in the order prescribed by the reference elements
 *
 * \param[out]	vCornerCoordsOut	vector of corner coordinates
 * \param[in] 	elem				Geometric Object
 * \param[in] 	domain				Domain
 * \param[in]	clearContainer		empty container before filling
 */
template <typename TElem, typename TDomain>
void CollectCornerCoordinates(	std::vector<typename TDomain::position_type>& vCornerCoordsOut,
								const TElem& elem, const TDomain& domain,
								bool clearContainer = true);

////////////////////////////////////////////////////////////////////////
///	returns the corner coordinates of a geometric object
///////////////////////////////////////////////////////////////////////////
/** Returns the corner coordinates for a given geometric object in a vector
 *
 * This function calls CollectCornerCoordinates based on the base object id
 *
 * @param[out] vCornerCoordsOut
 * @param[in] base_object_id
 * @param[in] elem
 * @param[in] domain
 * @param[in] clearContainer
 */
template <typename TDomain>
void CollectCornerCoordinates(int base_object_id,
			std::vector<typename TDomain::position_type>& vCornerCoordsOut,
			GridObject& elem, const TDomain& domain, bool clearContainer);

////////////////////////////////////////////////////////////////////////
///	returns the size of a geometric object
/** Returns the size of a geometric object
 *
 * This function returns the size of a geometric object.
 *
 * \param[in] 	elem				Geometric Object
 * \param[in] 	aaPos				AttachmentAccessor for Positions
 * \return 		number				Size of Element
 */
template <typename TElem, typename TPosition>
number ElementSize(const TElem& elem,
                   const Grid::VertexAttachmentAccessor<Attachment<TPosition> >& aaPos);

////////////////////////////////////////////////////////////////////////
///	returns the size of a geometric object
/** Returns the size of a geometric object
 *
 * This function returns the size of a geometric object.
 *
 * \param[in] 	elem				Geometric Object
 * \param[in] 	domain				Domain
 * \return 		number				Size of Element
 */
template <typename TElem, typename TDomain>
number ElementSize(const TElem& elem, const TDomain& domain);

////////////////////////////////////////////////////////////////////////
//	ElementDiameter
///	returns the maximal squared distance between to element vertices
template <typename TElem, typename TDomain>
number ElementDiameterSq(const TElem& elem, TDomain& domain);

///	returns the maximal distance between to element vertices
template <typename TElem, typename TDomain>
number ElementDiameter(const TElem& elem, TDomain& domain);

///	returns the maximal diameter of all elements between iterBegin and iterEnd.
template <typename TDomain>
number MaxElementDiameter(TDomain& domain, int level);

///	returns the minimal diameter of all elements between iterBegin and iterEnd.
template <typename TDomain>
number MinElementDiameter(TDomain& domain, int level);

// end group lib_disc_domain
/// \}

} // end namespace ug

////////////////////////////////
//	include implementation
#include "domain_util_impl.h"

#endif