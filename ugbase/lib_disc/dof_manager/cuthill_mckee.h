/*
 * cuthill_mckee.h
 *
 *  Created on: 21.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__CUTHILL_MCKEE__
#define __H__UG__LIB_DISC__DOF_MANAGER__CUTHILL_MCKEE__

#include <vector>
#include "dof_distribution.h"
#include "mg_dof_manager.h"
#include "lib_disc/function_spaces/approximation_space.h"

namespace ug{

/// returns an array describing the needed index mapping for Cuthill-McKee ordering
/**
 * This function computes a index mapping, that transforms a index-graph into
 * Cuthill-McKee ordering. For each index an vector of all adjacent indices
 * must be passed. (If no adjacent index is passed for an index, this index is
 * skipped and not sorted). On exit the index field vNewIndex is filled with
 * the index mapping: newInd = vNewIndex[oldInd]
 *
 * \param[out]	vNewIndex		vector returning new index for old index
 * \param[in]	vvNeighbour		vector of adjacent indices for each index
 * \param[in]	bReverse		flag if "reverse Cuthill-McKee" is used
 * \returns		flag if ordering was successful
 */
bool ComputeCuthillMcKeeOrder(std::vector<size_t>& vNewIndex,
                              std::vector<std::vector<size_t> >& vvNeighbour,
                              bool bReverse = true);

/// orders the dof distribution using Cuthill-McKee
template <typename TDoFImpl>
bool OrderCuthillMcKee(IDoFDistribution<TDoFImpl>& dofDistr,
                       bool bReverse)
{
//	get adjacency graph
	std::vector<std::vector<size_t> > vvConnection;
	if(!dofDistr.get_connections(vvConnection))
		return false;

//	get mapping for cuthill-mckee order
	std::vector<size_t> vNewIndex;
	if(!ComputeCuthillMcKeeOrder(vNewIndex, vvConnection, bReverse))
		return false;

//	reorder indices
	if(!dofDistr.permute_indices(vNewIndex))
		return false;

//	we're done
	return true;
}

/// orders the all DofDistributions of the MGDoFManager using Cuthill-McKee
template <typename TDoFImpl>
bool OrderCuthillMcKee(MGDoFManager<TDoFImpl>& mgDoFManager,
                       bool bReverse)
{
//	order levels
	for(size_t lev = 0; lev < mgDoFManager.num_levels(); ++lev)
		if(!OrderCuthillMcKee(*mgDoFManager.level_dof_distribution(lev), bReverse))
			return false;

//	order surface
	if(!OrderCuthillMcKee(*mgDoFManager.surface_dof_distribution(), bReverse))
		return false;

//	we're done
	return true;
}

/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain, typename TDoFImpl, typename TAlgebra>
bool OrderCuthillMcKee(ApproximationSpace<TDomain, TDoFImpl, TAlgebra>& approxSpace,
                       bool bReverse)
{
//	order levels
	for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev)
		if(!OrderCuthillMcKee(approxSpace.level_dof_distribution(lev), bReverse))
			return false;

//	order surface
	if(!OrderCuthillMcKee(approxSpace.surface_dof_distribution(), bReverse))
		return false;

//	we're done
	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__CUTHILL_MCKEE__ */
