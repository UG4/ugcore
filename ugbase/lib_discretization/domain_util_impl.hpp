//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d30

#ifndef __H__LIBDISCRETIZATION__DOMAIN_UTIL_IMPL__
#define __H__LIBDISCRETIZATION__DOMAIN_UTIL_IMPL__

namespace ug{

////////////////////////////////////////////////////////////////////////
template <class TDomain>
bool PrepareDomain(TDomain& domainOut, SubsetHandler& shTopViewOut,
					const char* filename,
					int numProcs,
					bool keepSrcGrid,
					size_t numPreRefinements,
					size_t numPostRefinements,
					bool writeProcessGrids,
					int autoAssignInnerObjectsToSubset,
					int autoAssignBoundaryObjectsToSubset)
{
	UG_LOG("PrepareDomain is currently not implemented for the serial case.\n");
	UG_LOG("  This has to be done in the near future!\n");
	return false;
}

} // end namespace ug


#endif /* __H__LIBDISCRETIZATION__DOMAIN_UTIL_IMPL__ */
