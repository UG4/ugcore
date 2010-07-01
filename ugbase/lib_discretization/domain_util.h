//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d30

#ifndef __H__LIBDISCRETIZATION__DOMAIN_UTIL__
#define __H__LIBDISCRETIZATION__DOMAIN_UTIL__

#include "lib_grid/lg_base.h"
#include "domain.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
///	creates a domain from a grid-file and distributes it over multiple processes.
/** After the grid is loaded it will be distributed to the specified processes.
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
					int autoAssignInnerObjectsToSubset = -2,
					int autoAssignBoundaryObjectsToSubset = -2);

} // end namespace ug


////////////////////////////////
//	include implementation
#ifdef UG_PARALLEL
	#include "domain_util_parallel_impl.hpp"
#else
	#include "domain_util_impl.hpp"
#endif

#endif /* __H__LIBDISCRETIZATION__DOMAIN_UTIL__ */
