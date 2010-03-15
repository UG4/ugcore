/*
 * dofpattern_impl.h
 *
 *  Created on: 05.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DOF_MANAGER__DOF_MANAGER_IMPL__
#define __H__LIBDISCRETIZATION__DOF_MANAGER__DOF_MANAGER_IMPL__

namespace ug{


template <typename TElem>
std::size_t
P1ConformDoFManager::
num_multi_indices(TElem* elem, uint fct)
{
	assert (fct < m_num_single_discrete_functions);
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element;
	return reference_element_traits<reference_element>::num_corners;
};

template<typename TElem>
inline
std::size_t
P1ConformDoFManager::
get_multi_indices(TElem* elem, uint fct, local_index_type& ind)
{
	assert (fct < m_num_single_discrete_functions);
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element;

	for(uint i = 0; i < reference_element_traits<reference_element>::num_corners; ++i)
	{
		VertexBase* vrt = elem->vertex(i);
		ind[i][0] = m_aaIndex[vrt] + (single_index_type) fct;
	}

	return reference_element_traits<reference_element>::num_corners;
}

template<typename TGeomObj>
inline
std::size_t
P1ConformDoFManager::
get_multi_indices_of_geom_obj(TGeomObj* obj, uint fct, local_index_type& ind)
{
	return 0;
}

template<>
inline
std::size_t
P1ConformDoFManager::
get_multi_indices_of_geom_obj(VertexBase* vrt, uint fct, local_index_type& ind)
{
	assert (fct < m_num_single_discrete_functions);

	ind[0][0] = m_aaIndex[vrt] + fct;

	return 1;
}


// groups all dofs of given dofgroups associated with one Geom Object type
template <typename TGeomObj>
bool
P1ConformDoFManager::
group_dof_groups(std::vector<uint>& selected_dofgroup_ids)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager has only one DoF group in a Vertex. You can not alter this pattern." << std::endl;

	return true;
}

template<typename TGeomObj>
bool
P1ConformDoFManager::
group_discrete_functions(std::vector<uint>& selected_functions)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager automatically groups solutions. You can not alter this pattern." << std::endl;

	return true;
}

template <typename TGeomObj>
bool
P1ConformDoFManager::
reset_grouping()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager automatically groups solutions. You can not alter this pattern." << std::endl;
	return true;
}

}

#endif /* __H__LIBDISCRETIZATION__DOF_MANAGER__DOF_MANAGER_IMPL__ */
