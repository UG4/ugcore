//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m06 d30

#ifndef __H__UG__LIB_DISC__DOMAIN_UTIL__
#define __H__UG__LIB_DISC__DOMAIN_UTIL__

// other lib_discretization headers
#include "domain.h"

namespace ug{

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

// end group lib_disc_domain
/// \}

} // end namespace ug

////////////////////////////////
//	include implementation
#include "domain_util_impl.h"

#endif /* __H__UG__LIB_DISC__DOMAIN_UTIL__ */
