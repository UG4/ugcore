/*
 * finite_volume_output.h
 *
 *  Created on: 06.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_OUTPUT__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_OUTPUT__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"

// finite volume geometry
#include "./finite_volume_geometry.h"

namespace ug{

//////////////////////////
// Output of Dual Grid
//////////////////////////

template <typename TElem, typename TDomain>
bool CreateDualElements(const TElem& elem, FV1Geometry<TElem, TDomain::dim>& geo, TDomain& domainOut, const TDomain& domain)
{
	// extract dimensions
	static const int worldDim = TDomain::dim;
	static const int refDim = FV1Geometry<TElem, worldDim>::dim;

	// extract grid
	typename TDomain::grid_type& grid = domain.get_grid();

	// extract subset handler
	const typename TDomain::subset_handler_type& sh = domain.get_subset_handler();

	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.get_position_accessor();

	// tmp vector for vertices
	std::vector<VertexBase*> vVert;

	// loop all scv of the element
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		const typename FV1Geometry<TElem, worldDim>::SCV& scv = geo.scv(i);

		// clear vertices
		vVert.clear();

		// loop corners of scv
		for(size_t co = 0; co < scv.num_corners(); ++co)
		{
			//	create a new vertex
				Vertex* vrt = grid.template create<Vertex>();
				vVert.push_back(vrt);

			//	set the coordinates
				aaPos[vrt] = scv.global_corner(co);
		}

		// edge
		if(refDim == 1)
		{
			grid.template create<Edge>(EdgeDescriptor(vVert[0], vVert[1]));
		}
		// face
		if(refDim == 2)
		{
			grid.template create<Quadrilateral>(QuadrilateralDescriptor(vVert[0], vVert[1],
																		vVert[2], vVert[3]));
		}
	}
	return true;
}


template <typename TElem, typename TDomain>
bool ConstructDualDomain(TDomain& domainOut, const TDomain& domain, int si)
{
	// extract world dimension
	static const int dim = TDomain::dim;

	// extract subset handler
	const typename TDomain::subset_handler_type& sh = domain.get_subset_handler();

	// extract grid
	typename TDomain::grid_type& grid = domain.get_grid();

	// Create Geometry
	FV1Geometry<TElem, dim> geo;

	// iterators for primary grid
	typename geometry_traits<TElem>::iterator iter, iterBegin, iterEnd;
	iterBegin = sh.template begin<TElem>(si);
	iterEnd = sh.template end<TElem>(si);

	// corners of element
	std::vector<typename TDomain::position_type> vCornerCoords;

	// iterate over primary grid
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get element
		TElem* elem = *iter;

		// get corner coordinates
		CollectCornerCoordinates(vCornerCoords, *elem, domain);

		// update finite volume geometry
		geo.update(elem, grid, &vCornerCoords[0]);

		// Create dual grid
		CreateDualElement(elem, geo, domainOut, domain);
	}

	return true;
}

template <typename TDomain>
bool ConstructDualDomain(TDomain& domainOut, const TDomain& domain, int si)
{
	const int siDim = DimensionOfSubset(domain, si);
	switch(siDim)
	{
		case 2: if(ConstructDualDomain<Triangle, TDomain>(domainOut, domain, si))
					{UG_LOG("CreateDualGrid: Error while processing Triangles.\n"); return false;}
				if(ConstructDualDomain<Quadrilateral, TDomain>(domainOut, domain, si))
					{UG_LOG("CreateDualGrid: Error while processing Quadrilaterals.\n"); return false;}
				break;
		default: UG_LOG("CreateDualGrid: Dimension " << siDim << " not supported.\n");
					return false;
	}
	return true;
}

template <typename TDomain>
bool WriteDualGridToFile(const char* filename, const TDomain& domain)
{
	// extract subset handler
	const typename TDomain::subset_handler_type& sh = domain.get_subset_handler();

	// create dual domain
	typename TDomain::grid_type dualGrid;
	typename TDomain::subset_handler_type dualSH;
	TDomain dualDomain(dualGrid, dualSH);

	// Construct dual domain
	for(int si = 0; si < sh.num_subsets(); ++si)
	{
		// write only elements with refDim == worldDim
		if(DimensionOfSubset(domain, si) != TDomain::dim) continue;

		if(!ConstructDualDomain(dualDomain, domain, si))
			{UG_LOG("WriteDualGridToFile: Error while writing subset "<<si<<".\n"); return false;}
	}

	if(!WriteDomainToUGX(filename, dualDomain))
		{UG_LOG("WriteDualGridToFile: Cannot write dual grid to file. Aborting.\n"); return false;}

	return true;
}


} // end namespace ug


#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_OUTPUT__ */
