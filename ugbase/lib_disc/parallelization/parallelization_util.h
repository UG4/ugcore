/*
 * parallelization.h
 *
 *  Created on: 21.5.2010
 *      Author: A. Vogel, S.Reiter
 */

#ifndef __H__UG__LIB_DISC__PARALLELIZATION__PARALLELIZATION_UTIL__
#define __H__UG__LIB_DISC__PARALLELIZATION__PARALLELIZATION_UTIL__

#include "pcl/pcl_util.h"
#include "lib_grid/parallelization/parallelization.h"
#include "lib_grid/parallelization/util/compol_interface_status.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#include "lib_algebra/parallelization/communication_policies.h"
#include "lib_algebra/parallelization/parallel_vector.h"
#include "lib_algebra/parallelization/parallel_matrix.h"
#include "lib_disc/dof_manager/dof_distribution.h"

namespace ug
{

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


/// copies all needed parallel informations into a parallel matrix
/**
 * This function copies the layouts and the communicators of a dof distribution
 * into a matrix.
 *
 * \param[in]	dd		DoFDistribution to copy the infos from
 * \param[out]	mat		Matrix filled with infos
 *
 * \tparam	TMatrix 	Sequential Matrix type
 */
template <typename TMatrix>
void CopyLayoutsAndCommunicatorIntoMatrix(ParallelMatrix<TMatrix>& mat,
                                          DoFDistribution& dd)
{
	mat.set_layouts(dd.layouts().master(), dd.layouts().slave());

	mat.set_communicator(dd.layouts().comm());
	mat.set_process_communicator(dd.layouts().proc_comm());
}

/// copies all needed parallel informations into a parallel vector
/**
 * This function copies the layouts and the communicators of a dof distribution
 * into a vector.
 *
 * \param[in]	dd		DoFDistribution to copy the infos from
 * \param[out]	vec		Vector filled with infos
 *
 * \tparam 	TVector		Sequential vector type
 */
template <typename TVector>
void CopyLayoutsAndCommunicatorIntoVector(ParallelVector<TVector>& vec,
                                          DoFDistribution& dd)
{
	//	copy all horizontal layouts (for all domain decomps)
		vec.set_layouts(dd.layouts().master(), dd.layouts().slave());

	//	copy vertical layouts
		vec.set_vertical_layouts(dd.layouts().vertical_master(),
		                         dd.layouts().vertical_slave());

	//	copy communicator
		vec.set_communicator(dd.layouts().comm());
		vec.set_process_communicator(dd.layouts().proc_comm());
}


bool CreateIndexLayouts_DomainDecomposition(
						IndexLayout& processLayoutOut,
						IndexLayout& subdomainLayoutOut,
						DoFDistribution& dofDistr,
						GridLayoutMap& layoutMap,
						int keyType, int level,
						pcl::IDomainDecompositionInfo* ddInfoIn);


// returns in a vector all appearencies of an index in a layout
void FindPositionInInterfaces(std::vector<std::pair<int, size_t> >& vIndexInterface,
                                     IndexLayout& layout, size_t index);

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
 * \param[in]	vIndNew		mapping for each index
 * \returns 	success flag
 */
void PermuteIndicesInIndexLayout(	IndexLayout& layout,
									const std::vector<size_t>& vIndNew);

}//	end of namespace

#endif
