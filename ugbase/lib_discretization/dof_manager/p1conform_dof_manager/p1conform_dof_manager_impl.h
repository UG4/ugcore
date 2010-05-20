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
P1ConformDoFManager(geom_obj_container_type& grid) : m_objContainer(grid)
{
	m_num_levels = grid.num_levels();
	m_bLocked = false;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim)
{
	if(id != LSFS_LAGRANGEP1)
	{
		std::cout << "P1ConformDoFManager: Only LSFS_LAGRANGEP1 functions are supported. Can not add function." << std::endl;
		return false;
	}

	m_SingleDiscreteFunctionNames.push_back(name);

	// TODO: Check if dim has been chosen correctly
	m_dim.push_back(dim);
	return true;
}


template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
add_discrete_function(std::string name, LocalShapeFunctionSetID id, std::vector<int>& SubsetIndices, int dim)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager currently only works on a grid. No Subset support." << std::endl;
	return false;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
group_discrete_functions(std::vector<uint>& selected_solutions)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager groups solutions automatically.\n";

	return true;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
reset_grouping()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	std::cout << "P1ConformDoFManager groups solutions automatically.\n";

	return true;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
print_info()
{
	using namespace std;

	cout << "\n================\n";
	cout << "P1 Conform DoF-Pattern INFO: " << endl;
	cout << "================\n";

	// print discrete functions
	cout << "Number of discrete functions: " << m_SingleDiscreteFunctionNames.size() << endl;
	for(uint nr_fct = 0; nr_fct < m_SingleDiscreteFunctionNames.size(); ++nr_fct)
	{
		cout << "# " << nr_fct << ":  '" << m_SingleDiscreteFunctionNames[nr_fct] << endl;
	}

	cout << "================\n";
	return true;
}

template <typename TGeomObjContainer>
bool
P1ConformDoFManager<TGeomObjContainer>::
finalize()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is already locked.\n";
		return true;
	}

	// TODO: Generalize for more than one subset
	int subsetIndex = 0;

	// Attach indices
	m_aIndex.resize(num_subsets());
	m_aaIndex.resize(num_subsets());
	for(subsetIndex = 0; subsetIndex < (int)m_objContainer.num_subsets(); ++subsetIndex)
	{
		m_objContainer.enable_subset_attachments(true);
		m_objContainer.template attach_to<VertexBase>(m_aIndex[subsetIndex], subsetIndex);
		m_aaIndex[subsetIndex].access(m_objContainer, m_aIndex[subsetIndex], subsetIndex);
	}

	// iterators
	geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

	// get number of single grid functions
	m_num_single_discrete_functions = m_SingleDiscreteFunctionNames.size();

	// loop levels
	m_num_dof_index.resize(m_objContainer.num_levels());
	for(uint level = 0; level < m_objContainer.num_levels(); ++level)
	{
		uint i = 0;

		for(subsetIndex = 0; subsetIndex < (int)m_objContainer.num_subsets(); ++subsetIndex)
		{
			iterBegin =  m_objContainer.begin<VertexBase>(subsetIndex, level);
			iterEnd = m_objContainer.end<VertexBase>(subsetIndex, level);

			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
				VertexBase* vrt = *iter;
				m_aaIndex[subsetIndex][vrt] = i;
				i += m_num_single_discrete_functions;
			}
		}
		m_num_dof_index[level] = i;

		UG_LOG( std::setw(8) << i << " DoF indices distributed" <<
					" on level " << std::setw(2) << level << " [each carrying " << 1 << " component(s)]" << std::endl);
	}

	m_bLocked = true;

	///////////////////
	// update IndexInfo

	// loop all levels
	m_IndexInfo.resize(m_num_levels);
	for(uint level = 0; level < m_num_levels; ++level)
	{
		m_IndexInfo[level].set_leaf(true);
		m_IndexInfo[level].set_num_index(m_num_dof_index[level]);
		m_IndexInfo[level].set_num_comp(1);
	}
	return true;
}

template <typename TGeomObjContainer>
LocalShapeFunctionSetID
P1ConformDoFManager<TGeomObjContainer>::
get_local_shape_function_set_id(uint nr_fct) const
{
	return LSFS_LAGRANGEP1;
}

template <typename TGeomObjContainer>
std::string
P1ConformDoFManager<TGeomObjContainer>::
get_name(uint fct) const
{
	return m_SingleDiscreteFunctionNames[fct];
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
std::size_t
P1ConformDoFManager<TGeomObjContainer>::
num_multi_indices(TElem* elem, uint fct) const
{
	assert (fct < m_num_single_discrete_functions);
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element;
	return reference_element_traits<reference_element>::num_corners;
};

template <typename TGeomObjContainer>
template<typename TElem>
inline
std::size_t
P1ConformDoFManager<TGeomObjContainer>::
get_multi_indices(TElem* elem, uint fct, local_index_type& ind, std::size_t offset) const
{
	assert (fct < m_num_single_discrete_functions);
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element;

	for(int i = 0; i < reference_element_traits<reference_element>::num_corners; ++i)
	{
		VertexBase* vrt = elem->vertex(i);
		int subsetIndex = m_objContainer.get_subset_index(vrt);

		ind[offset + i][0] = m_aaIndex[subsetIndex][vrt] + fct;
	}

	return reference_element_traits<reference_element>::num_corners;
}

template <typename TGeomObjContainer>
template<typename TGeomObj>
inline
std::size_t
P1ConformDoFManager<TGeomObjContainer>::
get_multi_indices_of_geom_obj(TGeomObj* obj, uint fct, local_index_type& ind, std::size_t offset) const
{
	return get_multi_indices_of_geom_obj_helper(obj, fct, ind, offset);
}

template <typename TGeomObjContainer>
template<typename TGeomObj>
inline
std::size_t
P1ConformDoFManager<TGeomObjContainer>::
get_multi_indices_of_geom_obj_helper(TGeomObj* obj, uint fct, local_index_type& ind, std::size_t offset) const
{
	return 0;
}

template <typename TGeomObjContainer>
inline
std::size_t
P1ConformDoFManager<TGeomObjContainer>::
get_multi_indices_of_geom_obj_helper(VertexBase* vrt, uint fct, local_index_type& ind, std::size_t offset) const
{
	assert (fct < m_num_single_discrete_functions);

	int subsetIndex = m_objContainer.get_subset_index(vrt);

	ind[offset + 0][0] = m_aaIndex[subsetIndex][vrt] + fct;

	return 1;
}


// groups all dofs of given dofgroups associated with one Geom Object type
template <typename TGeomObjContainer>
template <typename TGeomObj>
bool
P1ConformDoFManager<TGeomObjContainer>::
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

template <typename TGeomObjContainer>
template<typename TGeomObj>
bool
P1ConformDoFManager<TGeomObjContainer>::
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

template <typename TGeomObjContainer>
template <typename TGeomObj>
bool
P1ConformDoFManager<TGeomObjContainer>::
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
