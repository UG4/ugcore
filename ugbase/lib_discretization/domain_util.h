//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m06 d30

#ifndef __H__LIBDISCRETIZATION__DOMAIN_UTIL__
#define __H__LIBDISCRETIZATION__DOMAIN_UTIL__

// extern libraries
#include <vector>

// other ug4 libraries
#include "lib_grid/lg_base.h"

// other lib_discretization headers
#include "./domain.h"

namespace ug{

/// \ingroup lib_disc_domain
/// @{

////////////////////////////////////////////////////////////////////////
///	creates a domain from a grid-file and distributes it over multiple processes.
/** After the grid is loaded it will be distributed to the specified processes.
 *
 *	\param domainOut	prepared Domain
 *	\param shTopViewOut	prepared Surface View
 *	\param filename		File to reas Domain from
 *	\param numProcs		number of processes
 *
 *	\param keepSrcGrid: If set to true a copy of the whole grid
 *		(including pre-refinement) will be kept on process 0.
 *		Vertical interfaces will be created.
 *
 *	\param numPreRefinements: you may specify the amount of refinement
 *		performed before the grid is distributed.
 *
 *	\param numPostRefinements: you may specify the amount of refinement
 *		performed after the grid was distributed.
 *
 *	\param writeProcessGrids: specify if distributed grids should be written
 *		to file
 *
 *	\param autoAssignInnerObjectsToSubset, autoAssignBoundaryObjectsToSubset:
 *		if both are > -2, inner and boundary elements will automatically be
 *		assigned to the given subsets.
 *
 *	\todo: add a parameter by which the grid-partitioning algorithm can be
 *		selected.
 */
template <class TDomain>
bool PrepareDomain(TDomain& domainOut, SubsetHandler& shTopViewOut,
					const char* filename,
					int numProcs,
					bool keepSrcGrid,
					size_t numPreRefinements = 0,
					size_t numPostRefinements = 0,
					bool writeProcessGrids = true,
					int autoAssignInnerObjectsToSubset = -2,
					int autoAssignBoundaryObjectsToSubset = -2);

////////////////////////////////////////////////////////////////////////
///	returns the current dimension of the subset
/** Returns the dimension of geometric objects, that are contained in the subset
 *
 * This function returns the dimension of the subset. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained in the subset
 *
 * \param[in]	sh			SubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		dimension	Dimension of Subset
 * 				-1 			if no Dimension accessible
 */
inline int DimensionOfSubset(const SubsetHandler& sh, int si);

////////////////////////////////////////////////////////////////////////
///	returns the current dimension of the subset
/** Returns the dimension of geometric objects, that are contained in the subset
 *
 * This function returns the dimension of the subset. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained in the subset
 *
 * \param[in]	sh			MultiGridSubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		dimension	Dimension of Subset
 * 				-1 			if no Dimension accessible
 */
inline int DimensionOfSubset(const MGSubsetHandler& sh, int si);

////////////////////////////////////////////////////////////////////////
///	returns the current dimension of the subset
/** Returns the dimension of geometric objects, that are contained in the subset
 *
 * This function returns the dimension of the subset. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained in the subset
 *
 * \param[in]	sh			ISubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		dimension	Dimension of Subset
 * 				-1 			if no Dimension accessible
 */
inline int DimensionOfSubset(const ISubsetHandler& sh, int si);

////////////////////////////////////////////////////////////////////////
///	returns the current dimension for all subset
/** Returns the dimension of geometric objects, that are contained in the subset handler
 *
 * This function returns the dimension of the subsets. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained the union of all subset
 *
 * \param[in]	sh			ISubsetHandler
 *
 * \return		dimension	Dimension of Subset
 * 				-1 			if no Dimension accessible
 */
inline int DimensionOfSubsets(const ISubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
///	returns the current dimension of the subset
/** Returns the dimension of geometric objects, that are contained in the subset
 *
 * This function returns the dimension of the subset. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained in the subset
 *
 * \param[in]	domain		Domain
 * \param[in]	si			Subset Index
 *
 * \return		dimension	Dimension of Subset
 * 				-1 			if no Dimension accessible
 */
template <typename TDomain>
inline int DimensionOfSubset(const TDomain& domain, int si);

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
/// writes domain to *.ugx file
/** Writes a domain to *.ugx format
 *
 * This function writes a domain to a ugx-file.
 *
 * \param[in] 	filename		Filename
 * \param[in]	domain			Domain that is written to file
 * \return 		true			if successful
 * 				false			if error occurred
 */
template <typename TDomain>
bool WriteDomainToUGX(const char* filename, const TDomain& domain);

} // end namespace ug

/// @}

////////////////////////////////
//	include implementation
#include "domain_util_general_impl.h"
#ifdef UG_PARALLEL
	#include "domain_util_parallel_impl.hpp"
#else
	#include "domain_util_impl.hpp"
#endif

#endif /* __H__LIBDISCRETIZATION__DOMAIN_UTIL__ */
