/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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
#include "common/common.h"
#include "lib_disc/function_spaces/dof_position_util.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/domain.h"
#include <algorithm>
#include <vector>
#include <queue>
#include <utility>

namespace ug{


// Order for 1D
template<int dim>
bool ComparePosDim(const std::pair<MathVector<dim>, size_t> &p1,
                   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

template<>
bool ComparePosDim<1>(const std::pair<MathVector<1>, size_t> &p1,
                      const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDim<2>(const std::pair<MathVector<2>, size_t> &p1,
                      const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]) {
		return (p1.first[0] < p2.first[0]);
	}
	else {
		return (p1.first[1] < p2.first[1]);
	}
};

// Order for 3D
template<>
bool ComparePosDim<3>(const std::pair<MathVector<3>, size_t> &p1,
                      const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[2]==p2.first[2]){
		if (p1.first[1]==p2.first[1]) {
			return p1.first[0] < p2.first[0];
		}
		else {
			return (p1.first[1] < p2.first[1]);
		}
	}
	else{
		return (p1.first[2] < p2.first[2]);
	}
};


template<int dim>
void ComputeLexicographicOrder(std::vector<size_t>& vNewIndex,
                               std::vector<std::pair<MathVector<dim>, size_t> >& vPos)
{
//	a) order all indices
	if(vNewIndex.size() == vPos.size()){
	//  sort indices based on their position
		std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim>);

	//	write mapping
		for (size_t i=0; i < vPos.size(); ++i)
			vNewIndex[vPos[i].second] = i;
	}
//	b) only some indices to order
	else{
	//	create copy of pos
		std::vector<std::pair<MathVector<dim>, size_t> > vPosOrig(vPos);

	//  sort indices based on their position
		std::sort(vPos.begin(), vPos.end(), ComparePosDim<dim>);

	//	write mapping
		for (size_t i=0; i < vNewIndex.size(); ++i)
			vNewIndex[i] = i;
		for (size_t i=0; i < vPos.size(); ++i)
			vNewIndex[vPos[i].second] = vPosOrig[i].second;
	}
}

/// orders the dof distribution using Cuthill-McKee
template <typename TDomain>
void OrderLexForDofDist(SmartPtr<DoFDistribution> dd, ConstSmartPtr<TDomain> domain)
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
		const CommonLocalDoFSet& locDoF = LocalFiniteElementProvider::get_dofs(lfeID);

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
		const CommonLocalDoFSet& locDoF = LocalFiniteElementProvider::get_dofs(lfeID);

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
void OrderLex(ApproximationSpace<TDomain>& approxSpace, const char *order)
{
	// TODO: decode order input
	if(strcmp(order, "x")){
		UG_THROW("OrderLex: Currently only lexicographic in x-Direction implemented.");
	}

	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();

	for(size_t i = 0; i < vDD.size(); ++i)
		OrderLexForDofDist<TDomain>(vDD[i], approxSpace.domain());
}

#ifdef UG_DIM_1
template void ComputeLexicographicOrder<1>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<1>, size_t> >& vPos);
template void OrderLexForDofDist<Domain1d>(SmartPtr<DoFDistribution> dd, ConstSmartPtr<Domain1d> domain);
template void OrderLex<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, const char *order);
#endif
#ifdef UG_DIM_2
template void ComputeLexicographicOrder<2>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<2>, size_t> >& vPos);
template void OrderLexForDofDist<Domain2d>(SmartPtr<DoFDistribution> dd, ConstSmartPtr<Domain2d> domain);
template void OrderLex<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, const char *order);
#endif
#ifdef UG_DIM_3
template void ComputeLexicographicOrder<3>(std::vector<size_t>& vNewIndex, std::vector<std::pair<MathVector<3>, size_t> >& vPos);
template void OrderLexForDofDist<Domain3d>(SmartPtr<DoFDistribution> dd, ConstSmartPtr<Domain3d> domain);
template void OrderLex<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, const char *order);
#endif

}
