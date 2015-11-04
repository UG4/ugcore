/*
 * subset_util.cpp
 *
 *  Created on: 05.03.2012
 *      Author: andreasvogel
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
