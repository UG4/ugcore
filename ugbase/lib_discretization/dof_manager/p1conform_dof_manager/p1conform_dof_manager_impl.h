/*
 * dofpattern_impl.h
 *
 *  Created on: 05.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__P1_DOF_MANAGER__DOF_MANAGER_IMPL__
#define __H__LIBDISCRETIZATION__P1_DOF_MANAGER__DOF_MANAGER_IMPL__

namespace ug{



template <typename TGeomObjContainer>
P1ConformDoFManager<TGeomObjContainer>::
P1ConformDoFManager(geom_obj_container_type& grid) : m_bLocked(false), m_objContainer(grid)
{}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim)
{
	// for a P1 dof manager only Lagrange P1 function space is permitted
	if(id != LSFS_LAGRANGEP1)
	{
		UG_LOG("P1ConformDoFManager: Only LSFS_LAGRANGEP1 functions are supported. Can not add function.\n");
		return false;
	}

	// TODO: Check if dim has been chosen correctly

	// create a new function
	FunctionInfo info;
	info.name = name;
	info.dim = dim;

	m_vFunctionInfo.push_back(info);
	return true;
}


template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
add_discrete_function(std::string name, LocalShapeFunctionSetID id, std::vector<int>& SubsetIndices, int dim)
{
	if(m_bLocked == true)
	{
		UG_LOG("DoF Pattern is locked. Can not alter pattern.\n");
		return false;
	}

	UG_LOG("P1ConformDoFManager: Subset support not implemented yet. Currently only global function assignment possible.\n");
	return false;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
group_discrete_functions(std::vector<size_t>& selected_solutions)
{
	if(m_bLocked == true)
	{
		UG_LOG("DoF Pattern is locked. Can not alter pattern.\n");
		return false;
	}

	UG_LOG("P1ConformDoFManager groups solutions automatically.\n");
	return true;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
reset_grouping()
{
	if(m_bLocked == true)
	{
		UG_LOG("DoF Pattern is locked. Can not alter pattern.\n");
		return false;
	}

	UG_LOG("P1ConformDoFManager groups solutions automatically.\n");
	return true;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
print_info()
{
	UG_LOG("\n================\n");
	UG_LOG("P1 Conform DoF-Pattern INFO: \n");
	UG_LOG("================\n");

	// print discrete functions
	UG_LOG("Number of discrete functions: " << num_fct() << "\n");
	for(size_t fct = 0; fct < num_fct(); ++fct)
	{
		UG_LOG("# " << fct << ":  '" << get_name(fct) << "\n");
	}

	UG_LOG("================\n");
	return true;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
finalize()
{
	if(m_bLocked == true)
	{
		UG_LOG("DoF Pattern is already locked.\n");
		return true;
	}

	// TODO: Generalize for more than one subset
	int si = 0;

	// Attach indices
	m_vSubsetInfo.resize(num_subsets());
	for(si = 0; si < num_subsets(); ++si)
	{
		m_objContainer.enable_subset_attachments(true);
		m_objContainer.template attach_to<VertexBase>(m_vSubsetInfo[si].aDoF, si);
		m_vSubsetInfo[si].aaDoF.access(m_objContainer, m_vSubsetInfo[si].aDoF, si);
	}

	// iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// loop levels
	m_vLevelInfo.resize(num_levels());
	for(size_t level = 0; level < num_levels(); ++level)
	{
		size_t i = 0;

		for(si = 0; si < num_subsets(); ++si)
		{
			iterBegin =  begin<VertexBase>(si, level);
			iterEnd = end<VertexBase>(si, level);

			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
				VertexBase* vrt = *iter;
				(m_vSubsetInfo[si].aaDoF)[vrt] = i;
				i += num_fct();
			}
		}
		m_vLevelInfo[level].numDoFIndex = i;

		UG_LOG( std::setw(8) << i << " DoF indices distributed" <<
					" on level " << std::setw(2) << level << " [each carrying " << 1 << " component(s)]" << std::endl);
	}

	m_bLocked = true;

	///////////////////
	// update IndexInfo

	// loop all levels
	for(size_t level = 0; level < num_levels(); ++level)
	{
		m_vLevelInfo[level].indexInfo.set_leaf(true);
		m_vLevelInfo[level].indexInfo.set_num_index(num_dofs(level));
		m_vLevelInfo[level].indexInfo.set_num_comp(1);
	}
	return true;
}

template <typename TGeomObjContainer>
LocalShapeFunctionSetID
P1ConformDoFManager<TGeomObjContainer>::
get_local_shape_function_set_id(size_t fct) const
{
	return LSFS_LAGRANGEP1;
}

template <typename TGeomObjContainer>
typename P1ConformDoFManager<TGeomObjContainer>::geom_obj_container_type&
P1ConformDoFManager<TGeomObjContainer>::
get_assigned_geom_object_container()
{
	return m_objContainer;
}

template <typename TGeomObjContainer>
P1ConformDoFManager<TGeomObjContainer>::
~P1ConformDoFManager()
{}

template <typename TGeomObjContainer>
template <typename TElem>
size_t
P1ConformDoFManager<TGeomObjContainer>::
num_multi_indices(TElem* elem, size_t fct) const
{
	UG_ASSERT(fct < num_fct(), "Function not defined in DoF Manager");
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element;
	return reference_element_traits<reference_element>::num_corners;
};

template <typename TGeomObjContainer>
template<typename TElem>
inline
size_t
P1ConformDoFManager<TGeomObjContainer>::
get_multi_indices(TElem* elem, size_t fct, local_index_type& ind, size_t offset) const
{
	UG_ASSERT(fct < num_fct(), "Function not defined in DoF Manager");
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element;

	for(int i = 0; i < reference_element_traits<reference_element>::num_corners; ++i)
	{
		VertexBase* vrt = elem->vertex(i);
		int si = m_objContainer.get_subset_index(vrt);

		ind[offset + i][0] = (m_vSubsetInfo[si].aaDoF)[vrt] + fct;
	}

	return reference_element_traits<reference_element>::num_corners;
}

template <typename TGeomObjContainer>
template<typename TGeomObj>
inline
size_t
P1ConformDoFManager<TGeomObjContainer>::
get_multi_indices_of_geom_obj(TGeomObj* obj, size_t fct, local_index_type& ind, size_t offset) const
{
	return get_multi_indices_of_geom_obj_helper(obj, fct, ind, offset);
}

template <typename TGeomObjContainer>
template<typename TGeomObj>
inline
size_t
P1ConformDoFManager<TGeomObjContainer>::
get_multi_indices_of_geom_obj_helper(TGeomObj* obj, size_t fct, local_index_type& ind, size_t offset) const
{
	return 0;
}

template <typename TGeomObjContainer>
inline
size_t
P1ConformDoFManager<TGeomObjContainer>::
get_multi_indices_of_geom_obj_helper(VertexBase* vrt, size_t fct, local_index_type& ind, size_t offset) const
{
	UG_ASSERT(fct < num_fct(), "Function not defined in DoF Manager");

	int si = m_objContainer.get_subset_index(vrt);

	ind[offset + 0][0] = (m_vSubsetInfo[si].aaDoF)[vrt] + fct;

	return 1;
}


// groups all dofs of given dofgroups associated with one Geom Object type
template <typename TGeomObjContainer>
template <typename TGeomObj>
bool
P1ConformDoFManager<TGeomObjContainer>::
group_dof_groups(std::vector<size_t>& selected_dofgroup_ids)
{
	if(m_bLocked == true)
	{
		UG_LOG("DoF Pattern is locked. Can not alter pattern.\n");
		return false;
	}

	UG_LOG("P1ConformDoFManager has only one DoF group in a Vertex. You can not alter this pattern.\n");
	return true;
}

template <typename TGeomObjContainer>
template<typename TGeomObj>
bool
P1ConformDoFManager<TGeomObjContainer>::
group_discrete_functions(std::vector<size_t>& selected_functions)
{
	if(m_bLocked == true)
	{
		UG_LOG("DoF Pattern is locked. Can not alter pattern.\n");
		return false;
	}

	UG_LOG("P1ConformDoFManager automatically groups solutions. You can not alter this pattern.\n");
	return true;
}

template <typename TGeomObjContainer>
template <typename TGeomObj>
bool
P1ConformDoFManager<TGeomObjContainer>::
reset_grouping()
{
	if(m_bLocked == true)
	{
		UG_LOG("DoF Pattern is locked. Can not alter pattern.\n");
		return false;
	}

	UG_LOG("P1ConformDoFManager automatically groups solutions. You can not alter this pattern.\n");
	return true;
}

}

#endif /* __H__LIBDISCRETIZATION__DOF_MANAGER__DOF_MANAGER_IMPL__ */
