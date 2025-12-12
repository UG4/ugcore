/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__domain_distribution__
#define __H__UG__domain_distribution__

//#include "lib_disc/domain.h"
#include "lib_grid/tools/partition_map.h"
#include "lib_grid/parallelization/util/partition_weighting_callbacks.h"

namespace ug {

class PartitionMap;

/*///	partitions a domain by repeatedly cutting it along the different axis
template <typename TDomain>
static bool PartitionDomain_Bisection(TDomain& domain, PartitionMap& partitionMap,
									  int firstAxisToCut);*/

///	partitions a domain by sorting all elements into a regular grid
template <typename TDomain>
static bool PartitionDomain_RegularGrid(TDomain& domain, PartitionMap& partitionMap,
										int numCellsX, int numCellsY, int numCellsZ,
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
