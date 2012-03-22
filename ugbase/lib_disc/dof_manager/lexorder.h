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

namespace ug{

template<int dim>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
					std::vector<std::pair<MathVector<dim>, size_t> >& vPos);

///	writes positions of vertex dofs into a std::vector
template <typename TDD, typename TDomain>
void ExtractPos(const TDD& dd,
				typename TDomain::position_accessor_type& aaPos,
				std::vector<std::pair<MathVector<TDomain::dim>, size_t> >& vPositions)
{
//	resize positions
	vPositions.resize(dd.num_indices());

	typedef typename TDD::template traits<VertexBase>::const_iterator const_iterator;

//	loop all vertices
	const_iterator iter = dd.template begin<VertexBase>();
	const_iterator iterEnd = dd.template end<VertexBase>();

//	algebra indices vector
	std::vector<size_t> ind;

	for(;iter != iterEnd; ++iter)
	{
	//	get vertex
		VertexBase* v = *iter;

	//	load indices associated with vertex
		dd.inner_algebra_indices(v, ind);

	//	write position
		for(size_t i = 0; i < ind.size(); ++i)
		{
			const size_t index = ind[i];
			vPositions[index].first = aaPos[v];
			vPositions[index].second = index;
		}
	}
}

/// orders the dof distribution using Cuthill-McKee
template <typename TDD, typename TDomain>
void OrderLexForDofDist(TDD& dd,
						typename TDomain::position_accessor_type& aaPos)
{
//	get position attachment
	typedef MathVector<TDomain::dim> vec_type;
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> pos_type;

//	get positions of indices
	std::vector<pos_type> vPositions;
	ExtractPos<TDD, TDomain>(dd, aaPos, vPositions);

//	get mapping: old -> new index
	std::vector<size_t> vNewIndex(dd.num_indices());
	ComputeLexicographicOrder<TDomain::dim>(vNewIndex, vPositions);

	//	reorder indices
	dd.permute_indices(vNewIndex);
}


/// orders the all DofDistributions of the ApproximationSpace using lexicographic order
template <typename TDomain>
void OrderLex(ApproximationSpace<TDomain>& approxSpace,
              const char *order)
{
	// TODO: decode order input

	//	get position attachment
	typedef TDomain domain_type;
	typename domain_type::position_accessor_type& aaPos
			= approxSpace.domain()->position_accessor();

	//	order levels
	if(approxSpace.levels_enabled())
		for(size_t lev = 0; lev < approxSpace.num_levels(); ++lev)
			OrderLexForDofDist<LevelDoFDistribution, TDomain>(*approxSpace.level_dof_distribution(lev), aaPos);

	//	order surface
	if(approxSpace.top_surface_enabled())
		OrderLexForDofDist<SurfaceDoFDistribution, TDomain>(*approxSpace.surface_dof_distribution(), aaPos);
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__ */
