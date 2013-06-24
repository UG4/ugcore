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
#include "lib_disc/reference_element/reference_element_util.h"

namespace ug{

template<int dim>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
                               std::vector<std::pair<MathVector<dim>, size_t> >& vPos);

/// orders the dof distribution using Cuthill-McKee
template <typename TDomain>
void OrderLexForDofDist(SmartPtr<DoFDistribution> dd,
                        ConstSmartPtr<TDomain> domain)
{
//	Lex Ordering is only possible in this cases:
//	b) Same number of DoFs on each geometric object (or no DoFs on object)
//		--> in this case we can order all dofs
//	a) different trial spaces, but DoFs for each trial spaces only on separate
//	   geometric objects. (e.g. one space only vertices, one space only on edges)
//		--> in this case we can order all geometric objects separately

//	a) check for same number of DoFs on every geometric object
	bool bEqualNumDoFOnEachGeomObj = true;
	int numDoFOnGeomObj = -1;
	for(int si = 0; si < dd->num_subsets(); ++si){
		for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
			const int numDoF = dd->num_dofs((ReferenceObjectID)roid, si);

			if(numDoF == 0) continue;

			if(numDoFOnGeomObj == -1){
				numDoFOnGeomObj = numDoF;
			}
			else{
				if(numDoFOnGeomObj != numDoF)
					bEqualNumDoFOnEachGeomObj = false;
			}
		}
	}

//	b) check for non-mixed spaces
	std::vector<bool> bSingleSpaceUsage(NUM_REFERENCE_OBJECTS, true);
	std::vector<bool> vHasDoFs(NUM_REFERENCE_OBJECTS, false);
	for(size_t fct = 0; fct < dd->num_fct(); ++fct){

		LFEID lfeID = dd->local_finite_element_id(fct);
		const CommonLocalDoFSet& locDoF = LocalDoFSetProvider::get(lfeID);

		for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
			const int numDoF = locDoF.num_dof((ReferenceObjectID)roid);

			if(numDoF <= 0) continue;

			if(vHasDoFs[roid] == false){
				vHasDoFs[roid] = true;
			}
			else{
				bSingleSpaceUsage[roid] = false;
			}
		}
	}
	std::vector<bool> bSortableComp(dd->num_fct(), true);
	for(size_t fct = 0; fct < dd->num_fct(); ++fct){

		LFEID lfeID = dd->local_finite_element_id(fct);
		const CommonLocalDoFSet& locDoF = LocalDoFSetProvider::get(lfeID);

		for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
			if(locDoF.num_dof((ReferenceObjectID)roid) != 0){
				if(bSingleSpaceUsage[roid] == false)
					bSortableComp[fct] = false;
			}
		}
	}

//	get position attachment
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> pos_type;

//	get positions of indices
	std::vector<pos_type> vPositions;

//	a) we can order globally
	if(bEqualNumDoFOnEachGeomObj)
	{
		ExtractPositions(domain, dd, vPositions);

	//	get mapping: old -> new index
		std::vector<size_t> vNewIndex(dd->num_indices());
		ComputeLexicographicOrder<TDomain::dim>(vNewIndex, vPositions);

	//	reorder indices
		dd->permute_indices(vNewIndex);
	}
//	b) we can only order some spaces
	else
	{
		UG_LOG("OrderLex: Cannot order globally, trying to order some components:\n");
		for(size_t fct = 0; fct < dd->num_fct(); ++fct){
			if(bSortableComp[fct] == false){
				UG_LOG("OrderLex: '"<<dd->name(fct)<<" NOT SORTED.\n");
				continue;
			}

			ExtractPositions(domain, dd, fct, vPositions);

		//	get mapping: old -> new index
			std::vector<size_t> vNewIndex(dd->num_indices());
			ComputeLexicographicOrder<TDomain::dim>(vNewIndex, vPositions);

		//	reorder indices
			dd->permute_indices(vNewIndex);

			UG_LOG("OrderLex: '"<<dd->name(fct)<<" SORTED.\n");
		}
	}
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
			OrderLexForDofDist<TDomain>(approxSpace.level_dof_distribution(lev), approxSpace.domain());

	//	order surface
	if(approxSpace.top_surface_enabled())
		OrderLexForDofDist<TDomain>(approxSpace.surface_dof_distribution(), approxSpace.domain());
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__LEXORDER__ */
