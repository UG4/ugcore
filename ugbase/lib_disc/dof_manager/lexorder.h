/*
 * cuthill_mckee.h
 *
 *  Created on: 21.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__
#define __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__

#include <vector>
#include <utility> // for pair

#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/function_spaces/dof_position_util.h"

namespace ug{

template<int dim>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
                               std::vector<std::pair<MathVector<dim>, size_t> >& vPos);

/// orders the dof distribution using Cuthill-McKee
template <typename TDD, typename TDomain>
void OrderLexForDofDist(SmartPtr<TDD> dd,
                        ConstSmartPtr<TDomain> domain)
{
//	get position attachment
	typedef MathVector<TDomain::dim> vec_type;
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> pos_type;

//	get positions of indices
	std::vector<pos_type> vPositions;
	ExtractPositions<TDomain, TDD>(domain, dd, vPositions);

//	get mapping: old -> new index
	std::vector<size_t> vNewIndex(dd->num_indices());
	ComputeLexicographicOrder<TDomain::dim>(vNewIndex, vPositions);

	//	reorder indices
	dd->permute_indices(vNewIndex);
}


/// orders the all DofDistributions of the ApproximationSpace using lexicographic order
template <typename TDomain>
void OrderLex(ApproximationSpace<TDomain>& approxSpace,
              const char *order)
{
	// TODO: decode order input

	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev)
			OrderLexForDofDist<LevelDoFDistribution, TDomain>(approxSpace.level_dof_distribution(lev), approxSpace.domain());

	//	order surface
	if(approxSpace.top_surface_enabled())
		OrderLexForDofDist<SurfaceDoFDistribution, TDomain>(approxSpace.surface_dof_distribution(), approxSpace.domain());
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__ */
