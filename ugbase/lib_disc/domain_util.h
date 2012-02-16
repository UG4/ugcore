//	created by Sebastian Reiter, Andreas Vogel
//	s.b.reiter@googlemail.com
//	y10 m06 d30

#ifndef __H__UG__LIB_DISC__DOMAIN_UTIL__
#define __H__UG__LIB_DISC__DOMAIN_UTIL__

// extern libraries
#include <vector>

// other ug4 libraries
#include "lib_grid/lg_base.h"

// other lib_discretization headers
#include "domain.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
#endif

namespace ug{

///	Loads a domain from a grid-file.
/**	By optionally specifying a procId, you can make sure that the domain is only
 * loaded on one process. Pass -1, if you want to load it on all processes.
 * Note that the procId is only important in parallel environments.
 * \{
 */
template <typename TDomain>
bool LoadDomain(TDomain& domain, const char* filename);

template <typename TDomain>
bool LoadDomain(TDomain& domain, const char* filename, int procId);
/**	\} */

///	Saves the domain to a grid-file.
template <typename TDomain>
bool SaveDomain(TDomain& domain, const char* filename);


////////////////////////////////////////////////////////////////////////
/// returns if a subset is a regular grid
/**
 * This function returns if a subset contains constrained/constraining elements
 * such as hanging vertices, contrained edges/faces. In this case, the subset
 * does not form a regular grid.
 *
 *
 * \param[in]	sh			SubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		true		if subset is regular grid
 * 				false 		if subset is non-regular grid
 */
inline bool SubsetIsRegularGrid(const SubsetHandler& sh, int si);

////////////////////////////////////////////////////////////////////////
/// returns if a subset is a regular grid
/**
 * This function returns if a subset contains constrained/constraining elements
 * such as hanging vertices, contrained edges/faces. In this case, the subset
 * does not form a regular grid.
 *
 *
 * \param[in]	sh			SubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		true		if subset is regular grid
 * 				false 		if subset is non-regular grid
 */
inline bool SubsetIsRegularGrid(const MGSubsetHandler& sh, int si);

////////////////////////////////////////////////////////////////////////
/// returns if a subset is a regular grid
/**
 * This function returns if a subset contains constrained/constraining elements
 * such as hanging vertices, contrained edges/faces. In this case, the subset
 * does not form a regular grid.
 *
 *
 * \param[in]	sh			SubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		true		if subset is regular grid
 * 				false 		if subset is non-regular grid
 */
inline bool SubsetIsRegularGrid(const ISubsetHandler& sh, int si);

/// abbreviations for return types
enum {DIM_SUBSET_EMPTY_GRID = -1};

////////////////////////////////////////////////////////////////////////
///	Returns the dimension of geometric objects, that are contained in the subset
/**
 * This function returns the dimension of the subset. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained in the subset
 * If a ParallelCommunicator is passed, the highest dimension within all
 * procs in the ProcessCommunicator is returned.
 *
 * \param[in]	sh			ISubsetHandler
 * \param[in]	si			Subset Index
 * \param[in]	pProcCom	ParallelCommunicator (optinal)
 *
 * \return		dimension					Dimension of Subset
 * 				DIM_SUBSET_EMPTY_GRID		if empty Grid given
 */
inline int DimensionOfSubset(const ISubsetHandler& sh, int si
#ifdef UG_PARALLEL
                             , const pcl::ProcessCommunicator* pProcCom = NULL
#endif
							);

////////////////////////////////////////////////////////////////////////
///	Returns the dimension of geometric objects, that are contained in the subset handler
/**
 * This function returns the dimension of the subsets. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained the union of all subset
 * If a ParallelCommunicator is passed, the highest dimension within all
 * procs in the ProcessCommunicator is returned.
 *
 * \param[in]	sh			ISubsetHandler
 * \param[in]	pProcCom	ParallelCommunicator (optinal)
 *
 * \return		dimension					Dimension of Subset
 * 				DIM_SUBSET_EMPTY_GRID		if empty Grid given
 */
inline int DimensionOfSubsets(const ISubsetHandler& sh
#ifdef UG_PARALLEL
                             , const pcl::ProcessCommunicator* pProcCom = NULL
#endif
							);

////////////////////////////////////////////////////////////////////////
/// Returns the dimension of geometric objects, that are contained in the subset
/**
 * This function returns the dimension of the subset. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained in the subset
 * If a ParallelCommunicator is passed, the highest dimension within all
 * procs in the ProcessCommunicator is returned.
 *
 * \param[in]	domain		Domain
 * \param[in]	si			Subset Index
 * \param[in]	pProcCom	ParallelCommunicator (optinal)
 *
 * \return		dimension					Dimension of Subset
 * 				DIM_SUBSET_EMPTY_GRID		if empty Grid given
 */
template <typename TDomain>
inline int DimensionOfDomainSubset(const TDomain& domain, int si
#ifdef UG_PARALLEL
                             , const pcl::ProcessCommunicator* pProcCom = NULL
#endif
							);

////////////////////////////////////////////////////////////////////////
///	Returns the dimension of geometric objects, that are contained in the domain
/**
 * This function returns the dimension of the domain. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained in the domain.
 * If a ParallelCommunicator is passed, the highest dimension within all
 * procs in the ProcessCommunicator is returned.
 *
 * \param[in]	domain		Domain
 * \param[in]	pProcCom	ParallelCommunicator (optinal)
 *
 * \return		dimension					Dimension of Subset
 * 				DIM_SUBSET_EMPTY_GRID		if empty Grid given
 */
template <typename TDomain>
inline int DimensionOfDomain(const TDomain& domain
#ifdef UG_PARALLEL
                             , const pcl::ProcessCommunicator* pProcCom = NULL
#endif
							);

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
#include "domain_util_impl.h"

#endif /* __H__UG__LIB_DISC__DOMAIN_UTIL__ */
