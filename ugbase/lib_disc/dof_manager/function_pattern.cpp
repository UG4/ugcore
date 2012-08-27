/*
 * function_pattern.cpp
 *
 *  Created on: 03.12.2010
 *      Author: andreasvogel
 */

#include "function_pattern.h"
#include "lib_disc/common/subset_util.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif

namespace ug{

void FunctionPattern::
add_fct(const char* names, LFEID lfeID, int dim)
{
//	tokenize names
	std::vector<std::string> vName;
	TokenizeString(names, vName, ',');

//	add all names
	for(size_t i = 0; i < vName.size(); ++i)
	{
	//	trim string
		RemoveWhitespaceFromString(vName[i]);
		const char* name = vName[i].c_str();

	// 	if already locked, return false
		if(m_bLocked)
			UG_THROW("FunctionPattern: Already fixed. Cannot change.\n");

	//	check that space type has been passed
		if(lfeID.type() == LFEID::NONE || lfeID.order() < LFEID::ADAPTIV)
			UG_THROW("FunctionPattern: Specified Local Finite Element Space "
							<<lfeID<< " is not a valid space. "
							"[use e.g. (Lagrange, p), (DG, p), ...].\n");

	//	if no dimension passed, try to get dimension
		if(dim == -1) dim = DimensionOfSubsets(*m_spSH);

	#ifdef UG_PARALLEL
	//	some processes may not have an element of the subset at all. They can not
	//	know the dimension of the subset. Therefore we have to broadcast it.
		pcl::ProcessCommunicator ProcCom; dim = ProcCom.allreduce(dim, PCL_RO_MAX);
	#endif

	//	if still no dimension available, return false
		if(dim == -1)
			UG_THROW("FunctionPattern: Cannot find dimension for new function.\n");

	//	create temporary subset group
		SubsetGroup tmpSSGrp;
		tmpSSGrp.set_subset_handler(m_spSH);
		tmpSSGrp.add_all();

	// 	add to function list, everywhere = true, copy SubsetGroup
		m_vFunction.push_back(Function(name, dim, lfeID, true, tmpSSGrp));

	//	write info
//		UG_LOG("Info: Added discrete function with symbolic name '"<<name<<"' and "
//				"local finite element type " <<lfeID<<" on whole domain (dim="<<dim<<").\n");
	}
}

void FunctionPattern::add_fct(const char* names, LFEID lfeID,
                              const SubsetGroup& ssGrp, int dim)
{
//	tokenize names
	std::vector<std::string> vName;
	TokenizeString(names, vName, ',');

//	add all names
	for(size_t i = 0; i < vName.size(); ++i)
	{
	//	trim string
		RemoveWhitespaceFromString(vName[i]);
		const char* name = vName[i].c_str();

	// 	if already locked, return false
		if(m_bLocked)
			UG_THROW("FunctionPattern: Already fixed. Cannot change.\n");

	//	check that space type has been passed
		if(lfeID.type() == LFEID::NONE || lfeID.order() < LFEID::ADAPTIV)
			UG_THROW("FunctionPattern: "
					" Specified Local Finite Element Space "<<lfeID<< " is not "
					" a valid space. [use e.g. (Lagrange, p), (DG, p), ...].\n");

	//	check that subset handler are equal
		if(m_spSH.get() != ssGrp.subset_handler().get())
			UG_THROW("FunctionPattern: "
					"SubsetHandler of SubsetGroup does "
					"not match SubsetHandler of FunctionPattern.\n");

	//	if no dimension passed, try to get dimension
		if(dim == -1) dim = ssGrp.get_local_highest_subset_dimension();

	#ifdef UG_PARALLEL
	//	some processes may not have an element of the subset at all. They can not
	//	know the dimension of the subset. Therefore we have to broadcast it.
		pcl::ProcessCommunicator ProcCom; dim = ProcCom.allreduce(dim, PCL_RO_MAX);
	#endif

	//	if still no dimension available, return false
		if(dim == -1)
			UG_THROW("FunctionPattern: Cannot find dimension for new function.\n");

	// 	add to function list, everywhere = false, copy SubsetGroup as given
		m_vFunction.push_back(Function(name, dim, lfeID, false, ssGrp));

	//	write info
//		UG_LOG("Info: Added discrete function with symbolic name '"<<name<<"' and "
//				"local finite element type " <<lfeID<<" defined on subsets [");

		for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
		{
			if(i>0) UG_LOG(", ");
			UG_LOG(ssGrp.name(i));
		}
		UG_LOG("](dim="<<dim<<").\n");
	}
}

void FunctionPattern::add_fct(const char* name, LFEID lfeID,
                              const char* subsets, int dim)
{
//	create Function Group
	SubsetGroup ssGrp;

//	convert string to subset group
	try{
		ConvertStringToSubsetGroup(ssGrp, m_spSH, subsets);
	}UG_CATCH_THROW("FunctionPattern: ERROR while parsing Subsets.");

//	forward request
	add_fct(name, lfeID, ssGrp, dim);
}

void FunctionPattern::add_fct(const char* name, const char* fetype, int order)
{
//	convert type to LFEID and forward
	add_fct(name, ConvertStringToLFEID(fetype, order));
}

void FunctionPattern::add_fct(const char* name, const char* fetype)
{
//	convert type to LFEID and forward
	add_fct(name, ConvertStringToLFEID(fetype));
}

void FunctionPattern::add_fct(const char* name, const char* fetype,
                                        int order, const char* subsets)
{
//	convert type to LFEID and forward
	add_fct(name, ConvertStringToLFEID(fetype, order), subsets);
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

///	returns function group by name
FunctionGroup FunctionPattern::fct_grp_by_name(const char* names) const
{
	FunctionGroup fctGrp(*this);
	ConvertStringToFunctionGroup(fctGrp, *this, names);
	return fctGrp;
}


} // end namespace ug
