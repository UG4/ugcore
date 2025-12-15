/*
 * Copyright (c) 2011-2022:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel, Lukas Larisch
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

#include "lexorder.h"

#include <algorithm>
#include <vector>
#include <queue>
#include <utility>

#include "common/common.h"
#include "lib_disc/function_spaces/dof_position_util.h"
//ø #include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/domain.h"

#include "lib_disc/ordering_strategies/algorithms/lexorder_comparators.h"


namespace ug {

template<int dim>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
                               std::vector<std::pair<MathVector<dim>, size_t> >& vPos,
							   size_t orderDim, bool increasing)
{
//	a) order all indices
	if(vNewIndex.size() == vPos.size()){
	//  sort indices based on their position
		if(increasing){
			if (orderDim == 0)
				std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim, 0>);
			else if (orderDim == 1)
				std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim, 1>);
			else if (orderDim == 2)
				std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim, 2>);
			else UG_THROW("Invalid sorting direction.");
		}
		else{
			if (orderDim == 0)
				std::sort(vPos.begin(), vPos.end(), ComparePosDimDec<dim, 0>);
			else if (orderDim == 1)
				std::sort(vPos.begin(), vPos.end(), ComparePosDimDec<dim, 1>);
			else if (orderDim == 2)
				std::sort(vPos.begin(), vPos.end(), ComparePosDimDec<dim, 2>);
			else UG_THROW("Invalid sorting direction.");
		}

	//	write mapping
		for (size_t i=0; i < vPos.size(); ++i)
			vNewIndex[vPos[i].second] = i;
	}
//	b) only some indices to order
	else{
	//	create copy of pos
		std::vector<std::pair<MathVector<dim>, size_t> > vPosOrig(vPos);

	//  sort indices based on their position
		if(increasing){
			if (orderDim == 0){
				std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim, 0>);
			}
			else if (orderDim == 1)
				std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim, 1>);
			else if (orderDim == 2)
				std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim, 2>);
			else UG_THROW("Invalid sorting direction.");
		}
		else{
			if (orderDim == 0){
				std::sort(vPos.begin(), vPos.end(), ComparePosDimDec<dim, 0>);
			}
			else if (orderDim == 1)
				std::sort(vPos.begin(), vPos.end(), ComparePosDimDec<dim, 1>);
			else if (orderDim == 2)
				std::sort(vPos.begin(), vPos.end(), ComparePosDimDec<dim, 2>);
			else UG_THROW("Invalid sorting direction.");
		}

	//	write mapping
		for (size_t i=0; i < vNewIndex.size(); ++i)
			vNewIndex[i] = i;
		for (size_t i=0; i < vPos.size(); ++i)
			vNewIndex[vPos[i].second] = vPosOrig[i].second;
	}
}

/// orders the dof distribution using Cuthill-McKee
template <typename TDomain>
void OrderLexForDofDist(SmartPtr<DoFDistribution> dd, ConstSmartPtr<TDomain> domain, size_t orderDim, bool increasing)
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
			const int numDoF = dd->num_dofs(static_cast<ReferenceObjectID>(roid), si);

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
		const CommonLocalDoFSet& locDoF = LocalFiniteElementProvider::get_dofs(lfeID);

		for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
			const int numDoF = locDoF.num_dof(static_cast<ReferenceObjectID>(roid));

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
		const CommonLocalDoFSet& locDoF = LocalFiniteElementProvider::get_dofs(lfeID);

		for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid){
			if(locDoF.num_dof(static_cast<ReferenceObjectID>(roid)) != 0){
				if(bSingleSpaceUsage[roid] == false)
					bSortableComp[fct] = false;
			}
		}
	}

//	get position attachment
	using pos_type = std::pair<MathVector<TDomain::dim>, size_t>;

//	get positions of indices
	std::vector<pos_type> vPositions;

//	a) we can order globally
	if(bEqualNumDoFOnEachGeomObj)
	{
		ExtractPositions(domain, dd, vPositions);

	//	get mapping: old -> new index
		std::vector<size_t> vNewIndex(dd->num_indices());
		ComputeLexicographicOrder<TDomain::dim>(vNewIndex, vPositions, orderDim, increasing);

/*
		std::vector<bool> vCheck(dd->num_indices(), false);
		for (size_t i = 0; i < vNewIndex.size(); ++i)
		{
			UG_COND_THROW(vCheck.at(vNewIndex[i]), "Double mapping to index " << vNewIndex[i] << ".");
			vCheck.at(vNewIndex[i]) = true;
		}
		for (size_t i = 0; i < vCheck.size(); ++i)
			UG_COND_THROW(!vCheck[i], "Nothing maps to index " << i << ".");
*/

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
			ComputeLexicographicOrder<TDomain::dim>(vNewIndex, vPositions, orderDim, increasing);

		//	reorder indices
			dd->permute_indices(vNewIndex);

			UG_LOG("OrderLex: '"<<dd->name(fct)<<" SORTED.\n");
		}
	}
}


/// orders the all DofDistributions of the ApproximationSpace using lexicographic order
template <typename TDomain>
void OrderLex(ApproximationSpace<TDomain>& approxSpace, const char *order)
{
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();

	size_t len = strlen(order);

	if(len == 0){
		UG_THROW("OrderLex: Empty direction!");
	}

	size_t pos = 0;
	bool increasing = true;
	char sign;
	while(pos < len){
		if(increasing){
			sign = '+';
		}
		else{
			sign = '-';
		}

		switch(order[pos]){
			case '+':
				++pos;
				increasing = true;
				break;
			case '-':
				++pos;
				increasing = false;
				break;
			case 'x':
				UG_LOG("OrderLex: LexOrdering in " << sign << "x direction.\n")
				for (size_t i = 0; i < vDD.size(); ++i)
					OrderLexForDofDist<TDomain>(vDD[i], approxSpace.domain(), 0, increasing);
				++pos;
				increasing = true;
				break;
			case 'y':
				UG_LOG("OrderLex: LexOrdering in " << sign << "y direction.\n")
				for (size_t i = 0; i < vDD.size(); ++i)
					OrderLexForDofDist<TDomain>(vDD[i], approxSpace.domain(), 1, increasing);
				++pos;
				increasing = true;
				break;
			case 'z':
				UG_LOG("OrderLex: LexOrdering in " << sign << "z direction.\n")
				for (size_t i = 0; i < vDD.size(); ++i)
					OrderLexForDofDist<TDomain>(vDD[i], approxSpace.domain(), 2, increasing);
				++pos;
				increasing = true;
				break;
			default:
				UG_THROW("OrderLex: Invalid token in direction string, valid tokens: +, -, x, y, z");
		}
	}
}

#ifdef UG_DIM_1
template void ComputeLexicographicOrder<1>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<1>, size_t> >& vPos, size_t, bool);
template void OrderLexForDofDist<Domain1d>(SmartPtr<DoFDistribution> dd, ConstSmartPtr<Domain1d> domain, size_t, bool);
template void OrderLex<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, const char *order);
#endif
#ifdef UG_DIM_2
template void ComputeLexicographicOrder<2>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<2>, size_t> >& vPos, size_t, bool);
template void OrderLexForDofDist<Domain2d>(SmartPtr<DoFDistribution> dd, ConstSmartPtr<Domain2d> domain, size_t, bool);
template void OrderLex<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, const char *order);
#endif
#ifdef UG_DIM_3
template void ComputeLexicographicOrder<3>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<3>, size_t> >& vPos, size_t, bool);
template void OrderLexForDofDist<Domain3d>(SmartPtr<DoFDistribution> dd, ConstSmartPtr<Domain3d> domain, size_t, bool);
template void OrderLex<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, const char *order);
#endif

}
