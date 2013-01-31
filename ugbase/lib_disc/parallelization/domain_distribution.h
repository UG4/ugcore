// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 25.05.2011 (m,d,y)

#ifndef __H__UG__domain_distribution__
#define __H__UG__domain_distribution__

#include "lib_disc/domain.h"
#include "lib_grid/tools/partition_map.h"
#include "lib_grid/parallelization/util/partition_weighting_callbacks.h"

namespace ug
{

class PartitionMap;

///	partitions a domain by repeatedly cutting it along the different axis
template <typename TDomain>
static bool PartitionDomain_Bisection(TDomain& domain, PartitionMap& partitionMap,
									  int firstAxisToCut);

///	partitions a domain by sorting all elements into a regular grid
template <typename TDomain>
static bool PartitionDomain_RegularGrid(TDomain& domain, PartitionMap& partitionMap,
										int numCellsX, int numCellsY,
										bool surfaceOnly);

///	partitions a domain by using graph-based partitioning by METIS
template <typename TDomain>
static bool
PartitionDomain_MetisKWay(TDomain& domain, PartitionMap& partitionMap,
						  int numPartitions, size_t baseLevel = 0,
						  int hWeight = 1, int vWeight = 1);

///	partitions a domain by using graph-based partitioning by METIS
/** Weights can be chosen for borders between subsets in order to
 * prevent them from being part of the border of the geometry distribution
 * onto the available processors.
 */
template <typename TDomain>
static bool
PartitionDomain_MetisKWay(TDomain& domain, PartitionMap& partitionMap,
						  int numPartitions, size_t baseLevel,
						  SmartPtr<PartitionWeighting> weightFct);

/// Partitions a domain based on the elements of one level
/**	The elements are thereby weighted by the number of children they have.
 * Partitioning is done using the metis graph library.*/
template <typename TDomain>
static bool
PartitionDomain_LevelBased(TDomain& domain, PartitionMap& partitionMap,
						  	   int numPartitions, size_t level);


///	distributes a already distributed domain onto the specified processes
template <typename TDomain>
static bool DistributeDomain(TDomain& domainOut,
							 PartitionMap& partitionMap,
							 bool createVerticalInterfaces);

}//	end of namespace

////////////////////////////////
//	include implementation
#include "domain_distribution_impl.hpp"

#endif
