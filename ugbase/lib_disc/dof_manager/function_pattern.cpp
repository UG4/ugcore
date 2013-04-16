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

void FunctionPattern::set_subset_handler(ConstSmartPtr<ISubsetHandler> spSH){
	if(m_bLocked)
		UG_THROW("FunctionPattern: already locked, but trying to set"
				" new subset handler.");

	m_spSH = spSH;
	clear();
}

void FunctionPattern::add(const std::vector<std::string>& vName, LFEID lfeID, int dim)
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

	//	if no dimension passed, try to get dimension
		if(dim == -1) dim = DimensionOfSubsets(*m_spSH);

	#ifdef UG_PARALLEL
	//	some processes may not have an element of the subset at all. They can not
	//	know the dimension of the subset. Therefore we have to broadcast it.
		pcl::ProcessCommunicator ProcCom; dim = ProcCom.allreduce(dim, PCL_RO_MAX);
	#endif

	//	if still no dimension available, return false
		if(dim == -1)
			UG_THROW("FunctionPattern: Cannot find dimension for new function.");

	//	create temporary subset group
		SubsetGroup tmpSSGrp;
		tmpSSGrp.set_subset_handler(m_spSH);
		tmpSSGrp.add_all();

	// 	add to function list, everywhere = true, copy SubsetGroup
		m_vFunction.push_back(Function(name, dim, lfeID, true, tmpSSGrp));
	}
}

void FunctionPattern::add(const std::vector<std::string>& vName, LFEID lfeID,
                          const SubsetGroup& ssGrp, int dim)
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
					" a valid space. [use e.g. (Lagrange, p), (DG, p), ...].");

	//	check that subset handler are equal
		if(m_spSH.get() != ssGrp.subset_handler().get())
			UG_THROW("FunctionPattern: "
					"SubsetHandler of SubsetGroup does "
					"not match SubsetHandler of FunctionPattern.");

	//	if no dimension passed, try to get dimension
		if(dim == -1) dim = ssGrp.get_local_highest_subset_dimension();

	#ifdef UG_PARALLEL
	//	some processes may not have an element of the subset at all. They can not
	//	know the dimension of the subset. Therefore we have to broadcast it.
		pcl::ProcessCommunicator ProcCom; dim = ProcCom.allreduce(dim, PCL_RO_MAX);
	#endif

	//	if still no dimension available, return false
		if(dim == -1)
			UG_THROW("FunctionPattern: Cannot find dimension for new function.");

	// 	add to function list, everywhere = false, copy SubsetGroup as given
		m_vFunction.push_back(Function(name, dim, lfeID, false, ssGrp));

		for(size_t i = 0; i < ssGrp.size(); ++i)
		{
			if(i>0) UG_LOG(", ");
			UG_LOG(ssGrp.name(i));
		}
		UG_LOG("](dim="<<dim<<").\n");
	}
}

void FunctionPattern::add(const std::vector<std::string>& vName, LFEID lfeID,
                          const std::vector<std::string>& vSubset, int dim)
{
	add(vName, lfeID, SubsetGroup(m_spSH, vSubset), dim);
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

} // end namespace ug
