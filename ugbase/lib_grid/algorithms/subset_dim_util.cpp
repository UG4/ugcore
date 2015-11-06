/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#include "subset_dim_util.h"
#include "lib_grid/lg_base.h"

namespace ug{

/// returns if a subset is a regular grid
bool SubsetIsRegularGrid(const SubsetHandler& sh, int si)
{
//	check for constraining/constrained elements
	if(sh.num<ConstrainedVertex>(si) > 0) return false;
	if(sh.num<ConstrainedEdge>(si) > 0) return false;
	if(sh.num<ConstrainingEdge>(si) > 0) return false;
	if(sh.num<ConstrainedTriangle>(si) > 0) return false;
	if(sh.num<ConstrainingTriangle>(si) > 0) return false;
	if(sh.num<ConstrainedQuadrilateral>(si) > 0) return false;
	if(sh.num<ConstrainingQuadrilateral>(si) > 0) return false;

//	if not found, subset describes a regular grid
	return true;
}

/// returns if a subset is a regular grid
bool SubsetIsRegularGrid(const MGSubsetHandler& sh, int si)
{
//	check for constraining/constrained elements
	if(sh.num<ConstrainedVertex>(si) > 0) return false;
	if(sh.num<ConstrainedEdge>(si) > 0) return false;
	if(sh.num<ConstrainingEdge>(si) > 0) return false;
	if(sh.num<ConstrainedTriangle>(si) > 0) return false;
	if(sh.num<ConstrainingTriangle>(si) > 0) return false;
	if(sh.num<ConstrainedQuadrilateral>(si) > 0) return false;
	if(sh.num<ConstrainingQuadrilateral>(si) > 0) return false;

//	if not found, subset describes a regular grid
	return true;
}

/// returns if a subset is a regular grid
bool SubsetIsRegularGrid(const ISubsetHandler& ish, int si)
{
//	test SubsetHandler
	const SubsetHandler* sh = dynamic_cast<const SubsetHandler*>(&ish);
	if(sh != NULL)
		return SubsetIsRegularGrid(*sh, si);

//	test MGSubsetHandler
	const MGSubsetHandler* mgsh = dynamic_cast<const MGSubsetHandler*>(&ish);
	if(mgsh != NULL)
		return SubsetIsRegularGrid(*mgsh, si);

//	unknown type of subset handler
	UG_THROW("Unknown SubsetHandler type.");
	return false;
}

///	returns the current dimension of the subset
int DimensionOfSubset(const ISubsetHandler& sh, int si)
{
	try{
		return sh.subset_info(si).get_property("dim").to_int();
	}
	UG_CATCH_THROW("Make sure to properly set the dim-property of each subset "
					"before calling DimensionOfSubset! Use e.g. "
					"Domain::update_local_subset_dim_property, "
					"Domain::update_global_subset_dim_property, "
					"UpdateMaxDimensionOfSubset or UpdateGlobalMaxDimensionOfSubset.");
}

int DimensionOfSubsets(const ISubsetHandler& sh)
{
//	dimension to be computed
	int dim = DIM_SUBSET_EMPTY_GRID;

//	loop subsets
	for(int si = 0; si < sh.num_subsets(); ++si)
	{
	//	get dimension of subset
		int siDim = DimensionOfSubset(sh, si);

	//	if empty grid given, skip
		if(siDim == DIM_SUBSET_EMPTY_GRID) continue;

	//	check if dimension is higher than already checked subsets
		if(dim < siDim)
			dim = siDim;
	}

//	return computed domain
	return dim;
}

} // end namespace ug
