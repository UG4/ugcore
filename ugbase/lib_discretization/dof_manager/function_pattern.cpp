/*
 * function_pattern.cpp
 *
 *  Created on: 03.12.2010
 *      Author: andreasvogel
 */

#include "function_pattern.h"
#include "lib_discretization/domain_util.h"
#include "lib_discretization/common/groups_util.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif

namespace ug{

bool
FunctionPattern::
add_fct(const char* name, LFEID lfeID, int dim)
{
// 	if already locked, return false
	if(m_bLocked)
	{
		UG_LOG("Already fixed. Cannot change Distributor.\n");
		return false;
	}

//	check that space type has been passed
	if(lfeID.type() == LFEID::NONE || lfeID.order() < LFEID::ADAPTIV)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				" Specified Local Finite Element Space "<<lfeID<< " is not "
				" a valid space. [use e.g. (Lagrange, p), (DG, p), ...].\n");
		return false;
	}
	if(!supports_trial_space(lfeID))
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				" Specified Local Finite Element Space "<<lfeID<< " is not "
				" a supported.\n");
		return false;
	}

//	check that subset handler exists
	if(m_pSH == NULL)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"SubsetHandler not set.\n");
		return false;
	}

//	if no dimension passed, try to get dimension
	if(dim == -1) dim = DimensionOfSubsets(*m_pSH);

#ifdef UG_PARALLEL
//	some processes may not have an element of the subset at all. They can not
//	know the dimension of the subset. Therefore we have to broadcast it.
	pcl::ProcessCommunicator ProcCom; dim = ProcCom.allreduce(dim, PCL_RO_MAX);
#endif

//	if still no dimension available, return false
	if(dim == -1)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"Cannot find dimension for new function.\n");
		return false;
	}

//	create temporary subset group
	SubsetGroup tmpSSGrp;
	tmpSSGrp.set_subset_handler(*m_pSH);
	tmpSSGrp.add_all();

// 	add to function list, everywhere = true, copy SubsetGroup
	m_vFunction.push_back(Function(name, dim, lfeID, true, tmpSSGrp));

//	write info
	UG_LOG("Info: Added discrete function with symbolic name '"<<name<<"' and "
			"local finite element type " <<lfeID<<" on whole domain (dim="<<dim<<").\n");

//	we're done
	return true;
}

bool FunctionPattern::add_fct(const char* name, LFEID lfeID,
                              const SubsetGroup& ssGrp, int dim)
{
// 	if already locked, return false
	if(m_bLocked)
	{
		UG_LOG("Already fixed. Cannot change Distributor.\n");
		return false;
	}

//	check that space type has been passed
	if(lfeID.type() == LFEID::NONE || lfeID.order() < LFEID::ADAPTIV)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				" Specified Local Finite Element Space "<<lfeID<< " is not "
				" a valid space. [use e.g. (Lagrange, p), (DG, p), ...].\n");
		return false;
	}
	if(!supports_trial_space(lfeID))
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				" Specified Local Finite Element Space "<<lfeID<< " is not "
				" a supported.\n");
		return false;
	}

//	check that subset handler exists
	if(m_pSH == NULL)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"SubsetHandler not set.\n");
		return false;
	}

//	check that subset handler are equal
	if(m_pSH != ssGrp.get_subset_handler())
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"SubsetHandler of SubsetGroup does "
				"not match SubsetHandler of FunctionPattern.\n");
		return false;
	}

//	if no dimension passed, try to get dimension
	if(dim == -1) dim = ssGrp.get_local_highest_subset_dimension();

#ifdef UG_PARALLEL
//	some processes may not have an element of the subset at all. They can not
//	know the dimension of the subset. Therefore we have to broadcast it.
	pcl::ProcessCommunicator ProcCom; dim = ProcCom.allreduce(dim, PCL_RO_MAX);
#endif

//	if still no dimension available, return false
	if(dim == -1)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"Cannot find dimension for new function.\n");
		return false;
	}

// 	add to function list, everywhere = false, copy SubsetGroup as given
	m_vFunction.push_back(Function(name, dim, lfeID, false, ssGrp));

//	write info
	UG_LOG("Info: Added discrete function with symbolic name '"<<name<<"' and "
			"local finite element type " <<lfeID<<" defined on subsets [");
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
		if(i>0) UG_LOG(", ");
		UG_LOG(ssGrp.name(i));
	}
	UG_LOG("](dim="<<dim<<").\n");

//	done
	return true;
}

bool FunctionPattern::add_fct(const char* name, LFEID lfeID,
                              const char* subsets, int dim)
{
//	check that subset handler exists
	if(m_pSH == NULL)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"SubsetHandler not set.\n");
		return false;
	}

//	create Function Group
	SubsetGroup ssGrp;

//	convert string to subset group
	if(!ConvertStringToSubsetGroup(ssGrp, *m_pSH, subsets))
	{
		UG_LOG("ERROR while parsing Subsets.\n");
		return false;
	}

//	forward request
	return add_fct(name, lfeID, ssGrp, dim);
}

bool FunctionPattern::add_fct(const char* name, const char* fetype, int order)
{
//	convert type to LFEID and forward
	return add_fct(name, ConvertStringToLFEID(fetype, order));
}

bool FunctionPattern::add_fct(const char* name, const char* fetype,
                                        int order, const char* subsets)
{
//	convert type to LFEID and forward
	return add_fct(name, ConvertStringToLFEID(fetype, order), subsets);
}


} // end namespace ug
