/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Sebastian Reiter
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

#ifndef __H__UG__LIB_DISC__PARALLELIZATION__PARALLELIZATION_UTIL__
#define __H__UG__LIB_DISC__PARALLELIZATION__PARALLELIZATION_UTIL__

//#include "pcl/pcl_util.h"
#include "lib_grid/parallelization/util/compol_interface_status.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#include "lib_algebra/parallelization/communication_policies.h"
//#include "lib_algebra/parallelization/parallel_vector.h"
//#include "lib_algebra/parallelization/parallel_matrix.h"
#include "lib_disc/dof_manager/dof_distribution.h"

namespace ug {

/// creates the index layout for a level given a GridLayoutMap
/**
 * This function creates the Index layout based on a GridLayout. All elements
 * of the GridLayoutMap are loop on grid level and the indices attached to the
 * grid elements are added to the interface. Since the ordering of the grid
 * elements in the interfaces is assumed to be correct, also the ordering in
 * the index layouts are correct.
 *
 * \param[out]		layoutOut		the created index layout
 * \param[in]		dofDistr		the DoF Distribution
 * \param[in]		layoutMap		the grid layout map
 * \param[in]		keyType			key type (e.g. slave or master)
 * \param[in]		level			level, where layouts should be build
 *
 */
bool CreateLevelIndexLayout(	IndexLayout& layoutOut,
                            	DoFDistribution& dofDistr,
                            	GridLayoutMap& layoutMap,
                            	int keyType, int level);

/// creates the index layout for a level given a GridLayoutMap
/**
 * This function creates the Index layout based on a GridLayout. All elements
 * of the GridLayoutMap are loop level by level and the indices attached to the
 * grid elements are added to the interface, if an element does not have a children.
 * Since the ordering of the grid elements in the interfaces is assumed to be
 * correct, also the ordering in the index layouts are correct.
 *
 * \param[out]		layoutOut		the created index layout
 * \param[in]		dofDistr		the DoF Distribution
 * \param[in]		layoutMap		the grid layout map
 * \param[in]		keyType			key type (e.g. slave or master)
 * \param[in]		mg				underlying MultiGrid
 * \param[in]		dGrMgr			distributed Grid Manager
 */
bool CreateSurfaceIndexLayout(	IndexLayout& layoutOut,
                            	DoFDistribution& dofDistr,
                            	GridLayoutMap& layoutMap,
                            	int keyType,
                            	MultiGrid& mg, DistributedGridManager& dGrMgr);


bool CreateIndexLayouts_DomainDecomposition(
						IndexLayout& processLayoutOut,
						IndexLayout& subdomainLayoutOut,
						DoFDistribution& dofDistr,
						GridLayoutMap& layoutMap,
						int keyType, int level,
						pcl::IDomainDecompositionInfo* ddInfoIn);


// returns in a vector all appearencies of an index in a layout
void FindPositionInInterfaces(std::vector<std::pair<int, size_t> >& vIndexInterface,
                                     const IndexLayout& layout, size_t index);

bool AddExtraProcessEntriesToSubdomainLayout(
								size_t numIDs,
								IndexLayout& processMasterLayoutIn,
								IndexLayout& processSlaveLayoutIn,
								IndexLayout& subdomainMasterLayoutInOut,
								IndexLayout& subdomainSlaveLayoutInOut);

/// permutes an IndexLayout for the permutation of indices
/**
 * This Function changes the indices in the layout according to a given
 * permutation of the indices. (The order of the DoFs in the interfaces remains
 * the same, but the DoFs are "renamed")
 * The vector vIndNew must return the new index for each old index,
 * i.e. newIndex = vIndNew[oldIndex].
 *
 * \param[in]	layout		index layout
 * \param[in]	vIndNew		mapping for each index
  */
void PermuteIndicesInIndexLayout(	IndexLayout& layout,
									const std::vector<size_t>& vIndNew);

}//	end of namespace

#endif
