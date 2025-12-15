/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include "function_pattern.h"

#include "lib_grid/algorithms/subset_dim_util.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif

namespace ug {

void FunctionPattern::set_subset_handler(ConstSmartPtr<ISubsetHandler> spSH){
	if(m_bLocked)
		UG_THROW("FunctionPattern: already locked, but trying to set"
				" new subset handler.");

	m_spSH = spSH;
	clear();
}

void FunctionPattern::add(const std::vector<std::string>& vName, LFEID lfeID)
{
//	add all names
	for(size_t i = 0; i < vName.size(); ++i)
	{
	//	get string
		const char* name = vName[i].c_str();

	// 	if already locked, return false
		if(m_bLocked)
			UG_THROW("FunctionPattern: Already fixed. Cannot change.");

	//	check that space type has been passed
		if(lfeID.type() == LFEID::NONE || lfeID.order() < LFEID::ADAPTIV)
			UG_THROW("FunctionPattern: Specified Local Finite Element Space "
							<<lfeID<< " is not a valid space. "
							"[use e.g. (Lagrange, p), (DG, p), ...].");

	//	check dimension
		if(DimensionOfSubsets(*m_spSH) != lfeID.dim())
			UG_THROW("FunctionPattern: Adding "<<lfeID<<" to whole grid of"
			         " dimension "<<DimensionOfSubsets(*m_spSH)<<".")

	//	create temporary subset group
		SubsetGroup tmpSSGrp;
		tmpSSGrp.set_subset_handler(m_spSH);
		tmpSSGrp.add_all();

	// 	add to function list, everywhere = true, copy SubsetGroup
		m_vFunction.emplace_back(name, lfeID, true, tmpSSGrp);
	}
}

void FunctionPattern::add(const std::vector<std::string>& vName,
                          LFEID lfeID,
                          const SubsetGroup& ssGrp)
{
//	add all names
	for(size_t i = 0; i < vName.size(); ++i)
	{
	//	get string
		const char* name = vName[i].c_str();

	// 	if already locked, return false
		if(m_bLocked)
			UG_THROW("FunctionPattern: Already fixed. Cannot change.");

	//	check that space type has been passed
		if(lfeID.type() == LFEID::NONE || lfeID.order() < LFEID::ADAPTIV)
			UG_THROW("FunctionPattern: "
					" Specified Local Finite Element Space "<<lfeID<< " is not "
					" a valid space. [use e.g. (Lagrange, dim, p), (DG, dim, p), ...].");

	//	check that subset handler are equal
		if(m_spSH.get() != ssGrp.subset_handler().get())
			UG_THROW("FunctionPattern: "
					"SubsetHandler of SubsetGroup does "
					"not match SubsetHandler of FunctionPattern.");

	//	check dimension
		if(ssGrp.get_highest_subset_dimension() != lfeID.dim())
			UG_THROW("FunctionPattern: Adding "<<lfeID<<" to subsets with"
					 " highest dimension "<<ssGrp.get_highest_subset_dimension()<<".")

	// 	add to function list, everywhere = false, copy SubsetGroup as given
		m_vFunction.push_back(Function(name, lfeID, false, ssGrp));
	}
}

void FunctionPattern::add(const std::vector<std::string>& vName,
						  LFEID lfeID,
                          const std::vector<std::string>& vSubset)
{
	add(vName, lfeID, SubsetGroup(m_spSH, vSubset));
}


size_t FunctionPattern::fct_id_by_name(const char* name) const
{
	for(size_t i = 0; i < m_vFunction.size(); ++i)
	{
		if(m_vFunction[i].name == name)
			return i;
	}

	UG_THROW("Function name "<<name<<" not found in pattern.");
}

int FunctionPattern::subset_id_by_name(const char* name) const
{
	for(int i = 0; i < m_spSH->num_subsets(); ++i)
	{
		if(m_spSH->subset_info(i).name == name)
			return i;
	}

	UG_THROW("Subset name "<<name<<" not found in Subset Handler.");
}


} // end namespace ug
