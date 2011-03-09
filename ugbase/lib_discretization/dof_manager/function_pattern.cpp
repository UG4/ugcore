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
add_fct(const char* name, LocalShapeFunctionSetID id, int dim)
{
// 	if already locked, return false
	if(m_bLocked)
	{
		UG_LOG("Already fixed. Cannot change Distributor.\n");
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
	if(dim == -1)
		dim = DimensionOfSubsets(*m_pSH);

#ifdef UG_PARALLEL
//	some processes may not have an element of the subset at all. They can not
//	out the dimension of the subset. Therefore we have to broadcast it.
	pcl::ProcessCommunicator ProcCom;
	int dimGlob;
	ProcCom.allreduce(&dim, &dimGlob, 1, PCL_DT_INT, PCL_RO_MAX);
	dim = dimGlob;
#endif

//	if still no dimension available, return false
	if(dim == -1)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"Cannot find dimension for new function.\n");
		return false;
	}

//	create temporary subset group
	SubsetGroup tmpSubsetGroup;
	tmpSubsetGroup.set_subset_handler(*m_pSH);
	tmpSubsetGroup.add_all();

// 	add to function list, everywhere = true, copy SubsetGroup
	m_vFunction.push_back(Function(name, dim, id, true, tmpSubsetGroup));

//	we're done
	return true;
}

bool
FunctionPattern::
add_fct(const char* name,
                      LocalShapeFunctionSetID id,
                      const SubsetGroup& SubsetIndices,
                      int dim)
{
// 	if already locked, return false
	if(m_bLocked)
	{
		UG_LOG("Already fixed. Cannot change Distributor.\n");
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
	if(m_pSH != SubsetIndices.get_subset_handler())
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"SubsetHandler of SubsetGroup does "
				"not match SubsetHandler of FunctionPattern.\n");
		return false;
	}

//	if no dimension passed, try to get dimension
	if(dim == -1)
	{
		dim = SubsetIndices.get_highest_subset_dimension();
	}

#ifdef UG_PARALLEL
//	some processes may not have an element of the subset at all. They can not
//	out the dimension of the subset. Therefore we have to broadcast it.
	pcl::ProcessCommunicator ProcCom;
	int dimGlob;
	ProcCom.allreduce(&dim, &dimGlob, 1, PCL_DT_INT, PCL_RO_MAX);
	dim = dimGlob;
#endif

//	if still no dimension available, return false
	if(dim == -1)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"Cannot find dimension for new function.\n");
		return false;
	}


// 	add to function list, everywhere = false, copy SubsetGroup as given
	m_vFunction.push_back(Function(name, dim, id, false, SubsetIndices));

	return true;
}

bool
FunctionPattern::
add_fct(const char* name,
                      LocalShapeFunctionSetID id,
                      const char* subsets,
                      int dim)
{
// 	if already locked, return false
	if(m_bLocked)
	{
		UG_LOG("Already fixed. Cannot change Distributor.\n");
		return false;
	}

//	check that subset handler exists
	if(m_pSH == NULL)
	{
		UG_LOG("ERROR in 'FunctionPattern::add_discrete_function': "
				"SubsetHandler not set.\n");
		return false;
	}

//	create Function Group
	SubsetGroup subsetGroup;

//	convert string
	if(!ConvertStringToSubsetGroup(subsetGroup, *m_pSH, subsets))
	{
		UG_LOG("ERROR while parsing Subsets.\n");
		return false;
	}

//	forward request
	return add_fct(name, id, subsetGroup, dim);
}


} // end namespace ug
