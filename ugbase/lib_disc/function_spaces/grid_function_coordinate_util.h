/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_grid_function_coordinate_util
#define __H__UG_grid_function_coordinate_util

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/domain_util.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/reference_element/reference_mapping.h"

namespace ug{

template <typename TGridFunction>
void AddFunctionValuesToGridCoordinatesP1(
		SmartPtr<TGridFunction> spGridFct,
        size_t fct,
        size_t coordInd,
        const SubsetGroup& ssGrp,
        number timestep)
{
//	check if fast P1 interpolation may be used
	UG_COND_THROW(
		spGridFct->local_finite_element_id(fct).type() != LFEID::LAGRANGE
		|| spGridFct->local_finite_element_id(fct).order() != 1,
		"Fast P1 interpolation may only be used for Lagrange P1 functions.");

//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;

// get position accessor
	typename domain_type::position_accessor_type& aaPos
										= spGridFct->domain()->position_accessor();

	std::vector<DoFIndex> ind;
	typename TGridFunction::template dim_traits<0>::const_iterator iterEnd, iter;
				
	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];
	//	skip if function is not defined in subset
		if(!spGridFct->is_def_in_subset(fct, si)) continue;
	// 	iterate over all elements
		iterEnd = spGridFct->template end<Vertex>(si);
		iter = spGridFct->template begin<Vertex>(si);
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			Vertex* vrt = *iter;

		//	get multiindices of element
			spGridFct->dof_indices(vrt, fct, ind);

		// 	loop all dofs
			for(size_t i = 0; i < ind.size(); ++i)
			{
			//	set value#
				aaPos[vrt][coordInd] += DoFRef(*spGridFct, ind[i]) * timestep;
			}
		}
	}
}

template <typename TGridFunction>
void AddFunctionValuesToGridCoordinatesP1(
		SmartPtr<TGridFunction> spGridFct,
        const char* cmp,
        size_t coordInd)
{
		AddFunctionValuesToGridCoordinatesP1(spGridFct, cmp, coordInd, 1.0);
}


template <typename TGridFunction>
void AddFunctionValuesToGridCoordinatesP1(
		SmartPtr<TGridFunction> spGridFct,
        const char* cmp,
        size_t coordInd,
        number timestep)
{
	//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function found
	if(fct > spGridFct->num_fct())
		UG_THROW("Interpolate: Name of component '"<<cmp<<"' not found.");

	const bool bAllowManyfoldInterpolation =
			(spGridFct->local_finite_element_id(fct).type() == LFEID::LAGRANGE);

//	create subset group
	SubsetGroup ssGrp(spGridFct->domain()->subset_handler());
	//	add all subsets and remove lower dim subsets afterwards
	ssGrp.add_all();
	if(!bAllowManyfoldInterpolation)
		RemoveLowerDimSubsets(ssGrp);

	AddFunctionValuesToGridCoordinatesP1(spGridFct, fct, coordInd, ssGrp, timestep);
}

}//	end of namespace

#endif	//__H__UG_grid_function_coordinate_util
