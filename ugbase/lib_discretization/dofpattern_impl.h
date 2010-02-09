/*
 * dofpattern_impl.h
 *
 *  Created on: 05.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__DOFPATTERN_IMPL__
#define __H__LIBDISCRETIZATION__DOFPATTERN_IMPL__

namespace ug{


template <typename TElem>
DoFPattern::size_type
DoFPattern::
num_multi_indices(TElem* elem, uint nr_fct)
{
	return m_SingleDiscreteFunctionVec[nr_fct].get_element_dof_pattern().total_num_dofs();
}

template<typename TElem>
DoFPattern::size_type
DoFPattern::
get_multi_indices(TElem* elem, uint nr_fct, multi_index_type ind[])
{
	assert(nr_fct < m_SingleDiscreteFunctionVec.size());

	DoFPattern::size_type num_dofs = 0;

	const ElementDoFPattern& elemPattern = m_SingleDiscreteFunctionVec[nr_fct].get_element_dof_pattern();

	if(elemPattern.num_dofs<VertexBase>() > 0)
	{
		for(uint i = 0; i < reference_element_traits<TElem>::num_corners; ++i)
		{
			VertexBase* vert = elem->vertex(i);
			num_dofs += get_indices_of_geom_obj<VertexBase>(elem, nr_fct, ind + num_dofs);
		}
	}
	// other Geom Objects
	// [ ... ]

	return num_dofs;
}

template<typename TGeomObj>
DoFPattern::size_type
DoFPattern::
get_indices_of_geom_obj(TGeomObj* obj, uint nr_fct, multi_index_type ind[])
{
	assert(nr_fct < m_SingleDiscreteFunctionVec.size());

	return m_SingleDiscreteFunctionVec[nr_fct].get_multi_indices_of_geom_obj<TGeomObj>(obj, ind);
}

// groups all dofs of given dofgroups associated with one Geom Object type
template <typename TGeomObj>
bool
DoFPattern::
group_dof_groups(std::vector<uint>& selected_dofgroup_ids)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	typename dof_group_traits<TGeomObj>::container_type& container = get_dof_group_container<TGeomObj>();

	// if dofgroup - selection is void, do nothing
	if(selected_dofgroup_ids.size() == 0) return true;

	// sort dofgroups (just for convinience)
	sort(selected_dofgroup_ids.begin(), selected_dofgroup_ids.end());

	// delete ids appearing more than once
	for(uint i = 1; i < selected_dofgroup_ids.size(); ++i)
	{
		if(selected_dofgroup_ids[i] == selected_dofgroup_ids[i-1])
		{
			selected_dofgroup_ids.erase(selected_dofgroup_ids.begin() + i);
			i--;
		}
	}

	// remember subset_handler and subset_index of first DoFGroup
	DoFGroup<TGeomObj>& first_group = container[selected_dofgroup_ids[0]];
	subset_handler_type* sh = first_group.get_subset_handler();
	int subsetIndex = first_group.get_subset_index();

	// check that all dof groups in 'dofgroup' have the same subsethandler & subsetindex
	for(uint i = 0; i < selected_dofgroup_ids.size(); ++i)
	{
		// get id of group
		uint id = selected_dofgroup_ids[i];
		if(container[id].get_subset_handler() != sh)
		{
			std::cout << "ERROR in group_dof_groups: DoFGroups do not have same subset_handler.\n";
			return false;
		}
		if(container[id].get_subset_index() != subsetIndex)
		{
			std::cout << "ERROR in group_dof_groups: DoFGroups do not have same subset_index.\n";
			std::cout << "DoFGroups are: ";
			for(uint j = 0; j < selected_dofgroup_ids.size() - 1; ++j)
			{
				std::cout << selected_dofgroup_ids[j] << ", ";
			}
			std::cout << selected_dofgroup_ids[selected_dofgroup_ids.size()-1] << std::endl;
			return false;
		}
	}

	// create new dof group for subset
	DoFGroup<TGeomObj> group_new(*sh, subsetIndex);

	// add all selected dof groups: copy dofs into new one and erase old group
	for(uint i = 0; i < selected_dofgroup_ids.size(); ++i)
	{
		const uint id_old = selected_dofgroup_ids[i];
		const DoFGroup<TGeomObj>& group_old = container[id_old];

		// add old dofs to new group
		uint nr_group_new = container.size() - 1;
		for(uint comp_old = 0; comp_old < group_old.num_dofs(); ++comp_old)
		{
			const uint nr_fct = group_old.get_nr_fct(comp_old);
			const uint comp_fct = group_old.get_comp_fct(comp_old);
			const uint comp_new = group_new.add_dof(nr_fct, comp_fct);
			SingleDiscreteFunction& fct = m_SingleDiscreteFunctionVec[nr_fct];

			// adjust dofgroup_map for solution 'nr_fct'
			fct.set_nr_dofgroup<TGeomObj>(subsetIndex, comp_fct, nr_group_new);
			fct.set_comp_in_dofgroup<TGeomObj>(subsetIndex, comp_fct, comp_new);
		}

		// remove old group
		container.erase(container.begin() + id_old);

		// adjust dofgroup_maps, since ids (for id >= id_old) have been decreased by 1
		for(uint nr_fct = 0; nr_fct < m_SingleDiscreteFunctionVec.size(); ++nr_fct)
		{
			SingleDiscreteFunction& fct = m_SingleDiscreteFunctionVec[nr_fct];
			for(uint subsetIndex2 = 0; subsetIndex2 < m_num_subsets; ++subsetIndex2)
			{
				for(uint comp = 0; comp < fct.num_dofs<TGeomObj>(subsetIndex2); comp++)
				{
					if(fct.nr_dofgroup<TGeomObj>(subsetIndex2, comp) >= id_old)
					{
						uint comp_old = fct.nr_dofgroup<TGeomObj>(subsetIndex2, comp);
						fct.set_nr_dofgroup<TGeomObj>(subsetIndex2, comp, comp_old-1);
					}
				}
			}
		}

		// adjust selected_dofgroup_ids
		for(uint j = 0; j < selected_dofgroup_ids.size(); ++j)
		{
			selected_dofgroup_ids[j] -= 1;
		}
	}

	// add new grouping to container
	container.push_back(group_new);

	return true;
}

template<typename TGeomObj>
bool
DoFPattern::
group_discrete_functions(std::vector<uint>& selected_functions)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	//check if all solution ids are correct
	for(uint i = 0; i < selected_functions.size(); ++i)
	{
		if(selected_functions[i] >= num_fct())
		{
			std::cout << "Invalid discrete function id.\n";
			return false;
		}
	}

	//loop over all subsets
	const int numSubsets = m_sh.num_subsets();
	for(int subsetIndex = 0; subsetIndex < numSubsets; ++subsetIndex)
	{
		std::vector<uint> selected_dof_groups;
		// loop all selected Discrete Functions
		for(uint i=0; i < selected_functions.size(); ++i)
		{
			const uint nr_fct = selected_functions[i];
			const SingleDiscreteFunction& fct = m_SingleDiscreteFunctionVec[nr_fct];
			for(uint dof = 0; dof < fct.num_dofs<TGeomObj>(subsetIndex); ++dof)
			{
				uint group_id = fct.nr_dofgroup<TGeomObj>(subsetIndex, dof);
				selected_dof_groups.push_back(group_id);
			}
		}
		group_dof_groups<TGeomObj>(selected_dof_groups);
		selected_dof_groups.clear();
	}

	return true;
}

template <typename TGeomObj>
bool
DoFPattern::
reset_grouping()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	// clear Grouping vector in Geom Object
	typename dof_group_traits<TGeomObj>::container_type& container = get_dof_group_container<TGeomObj>();
	container.clear();

	// loop all subsets
	const int numSubsets = m_sh.num_subsets();
	for(int subsetIndex = 0; subsetIndex < numSubsets; ++subsetIndex)
	{
		// loop all discrete functions and add DoF groups again
		for(uint nr_fct = 0; nr_fct < m_SingleDiscreteFunctionVec.size(); ++nr_fct)
		{
			const SingleDiscreteFunction& fct = m_SingleDiscreteFunctionVec[nr_fct];
			uint num_dofs = fct.num_dofs<TGeomObj>(subsetIndex);

			if(add_dof_group<TGeomObj>(nr_fct, subsetIndex, num_dofs) != true) return false;
		}
	}
	return true;
}

// private help function
// adds num_dofs dof groups with one component on TGeomObj for discrete fct nr_fct, on subsetIndex
template <typename TGeomObj>
bool
DoFPattern::
add_dof_group(uint nr_fct, int subsetIndex, uint num_dofs)
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is locked. Can not alter pattern.\n";
		return false;
	}

	typename dof_group_traits<TGeomObj>::container_type& container = get_dof_group_container<TGeomObj>();

	m_SingleDiscreteFunctionVec[nr_fct].set_num_dofs<TGeomObj>(subsetIndex, num_dofs);
	for(uint dof = 0; dof < num_dofs; ++dof)
	{
		container.push_back(DoFGroup<VertexBase>(m_sh, subsetIndex));
		dof_comp_type comp = container.back().add_dof(nr_fct, dof);

		m_SingleDiscreteFunctionVec[nr_fct].set_nr_dofgroup<TGeomObj>(subsetIndex, dof, container.size() - 1);
		m_SingleDiscreteFunctionVec[nr_fct].set_comp_in_dofgroup<TGeomObj>(subsetIndex, dof, comp);
	}

	return true;
}

template <typename TGeomObj>
bool
DoFPattern::
finalize()
{
	if(m_bLocked == true)
	{
		std::cout << "DoF Pattern is already locked.\n";
		return true;
	}

	// clear Grouping vector in Geom Object
	typename dof_group_traits<TGeomObj>::container_type& container = get_dof_group_container<TGeomObj>();

	// loop all dof groups
	for(uint i = 0; i < container.size(); ++i)
	{
		for(uint level = 0; level < m_num_levels; ++level)
			{
				if(container[i].activate(level) != true) return false;
			}
	}

	return true;
}


///////////////////////////////
// DoFGroup Map
///////////////////////////////

inline
uint
DoFPattern::DoFGroupMap::
num_dofs() const
{
	return m_num_dofs;
}

inline
uint
DoFPattern::DoFGroupMap::
nr_dofgroup(uint i) const
{
	return m_nr_dofgroup[i];
}

inline
uint
DoFPattern::DoFGroupMap::
comp_in_dofgroup(uint i) const
{
	return m_comp_in_dofgroup[i];
}



///////////////////////////////
// DoFGroup access
///////////////////////////////


template <> inline std::vector<DoFPattern::DoFGroup<VertexBase> >& DoFPattern::get_dof_group_container() {return m_VertexBaseDoFGroups;}
template <> class DoFPattern::dof_group_traits<VertexBase>{ public: typedef std::vector<DoFPattern::DoFGroup<VertexBase> > container_type;};

///////////////////////////////
// DoFGroup
///////////////////////////////


template <typename TGeomObj>
DoFPattern::DoFGroup<TGeomObj>::
DoFGroup() : m_sh(NULL), m_subsetIndex(-1), m_pActivated(false), m_num_dof(0)
{};


template <typename TGeomObj>
DoFPattern::DoFGroup<TGeomObj>::
DoFGroup(subset_handler_type& sh, int s) : m_sh(&sh), m_subsetIndex(s), m_pActivated(false), m_num_dof(0)
{
	m_nr_fct.clear();
	m_nr_comp.clear();
};

template <typename TGeomObj>
DoFPattern::DoFGroup<TGeomObj>::
~DoFGroup()
{
	if(m_pActivated == true) deactivate();
};

template <typename TGeomObj>
bool
DoFPattern::DoFGroup<TGeomObj>::
activate(uint level)
{
	assert(m_sh != NULL);
	if(!m_aaIndex.valid())
	{
		m_sh->enable_subset_attachments(true);
		m_sh->attach_to<TGeomObj>(m_aIndex, m_subsetIndex);
		m_aaIndex.access(*m_sh, m_aIndex, m_subsetIndex);
	}

	typename geometry_traits<TGeomObj>::iterator iter, iterBegin, iterEnd;
	iterBegin =  m_sh->begin<TGeomObj>(m_subsetIndex, level);
	iterEnd = m_sh->end<TGeomObj>(m_subsetIndex, level);

	dof_index_type i = 0;
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
		TGeomObj* elem = *iter;
		m_aaIndex[elem] = i++;
	}
	m_num_dof_index[level] = i;

	std::cout << i << " DoF indices distributed on subset " << m_subsetIndex <<
				" on level " << level << " [each carrying " << m_num_dof << " component(s)]" << std::endl;

	m_pActivated = true;
	return true;
};

template <typename TGeomObj>
bool
DoFPattern::DoFGroup<TGeomObj>::
deactivate()
{
	m_sh->detach_from<TGeomObj>(m_aIndex, m_subsetIndex);
	m_num_dof_index.clear();
	m_pActivated = false;
}

template <typename TGeomObj>
inline
DoFPattern::dof_index_type
DoFPattern::DoFGroup<TGeomObj>::
index(TGeomObj* elem) const
{
	assert(m_aaIndex.valid());
	return m_aaIndex[elem];
}

template <typename TGeomObj>
inline
DoFPattern::dof_comp_type
DoFPattern::DoFGroup<TGeomObj>::
num_dofs() const
{
	return m_num_dof;
}

template <typename TGeomObj>
DoFPattern::dof_comp_type
DoFPattern::DoFGroup<TGeomObj>::
add_dof(uint nr_fct, uint fct_comp)
{
	m_nr_fct.push_back(nr_fct);
	m_nr_comp.push_back(fct_comp);
	dof_comp_type comp = m_num_dof++;
	return comp;
}

template <typename TGeomObj>
DoFPattern::subset_handler_type*
DoFPattern::DoFGroup<TGeomObj>::
get_subset_handler() const
{
	return m_sh;
}

template <typename TGeomObj>
int
DoFPattern::DoFGroup<TGeomObj>::
get_subset_index() const
{
	return m_subsetIndex;
}

template <typename TGeomObj>
uint
DoFPattern::DoFGroup<TGeomObj>::
get_nr_fct(uint comp) const
{
	assert(comp < m_nr_fct.size());
	return m_nr_fct[comp];
}

template <typename TGeomObj>
uint
DoFPattern::DoFGroup<TGeomObj>::
get_comp_fct(uint comp) const
{
	assert(comp < m_nr_comp.size());
	return m_nr_comp[comp];
}


template <typename TGeomObj>
inline
DoFPattern::size_type
DoFPattern::SingleDiscreteFunction::
get_multi_indices_of_geom_obj(TGeomObj* obj, multi_index_type ind[])
{
	int subsetIndex = m_sh->get_subset_index(obj);
	uint i = geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID;
	size_type dof;
	for(dof = 0; dof < m_dofGroupMap[i][subsetIndex].num_dofs(); ++dof)
	{
		typename dof_group_traits<TGeomObj>::container_type& container = get_dof_group_container<TGeomObj>();
		const uint nr_dofgroup = m_dofGroupMap[i][subsetIndex].nr_dofgroup(dof);
		const DoFGroup<TGeomObj>& group = container[nr_dofgroup];
		ind[dof][3] = subsetIndex;
		ind[dof][2] = nr_dofgroup;
		ind[dof][1] = group.index(obj);
		ind[dof][0] = m_dofGroupMap[i][subsetIndex].comp_in_dofgroup(dof);
	}
	return dof;
}

template <typename TElem>
bool
DoFPattern::SingleDiscreteFunction::
set_num_dofs(int subsetIndex, uint num)
{
	uint i = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	return m_dofGroupMap[i][subsetIndex].set_num_dofs(num);
}

template <typename TElem>
inline
uint
DoFPattern::SingleDiscreteFunction::
num_dofs(int subsetIndex) const
{
	const uint i = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	return m_dofGroupMap[i][subsetIndex].num_dofs();
}

template <typename TElem>
inline
bool
DoFPattern::SingleDiscreteFunction::
set_nr_dofgroup(int subsetIndex, uint dof, uint group)
{
	const uint i = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	return m_dofGroupMap[i][subsetIndex].set_nr_dofgroup(dof, group);
}

template <typename TElem>
inline
uint
DoFPattern::SingleDiscreteFunction::
nr_dofgroup(int subsetIndex, uint dof) const
{
	const uint i = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	return m_dofGroupMap[i][subsetIndex].nr_dofgroup(dof);
}

template <typename TElem>
inline
bool
DoFPattern::SingleDiscreteFunction::
set_comp_in_dofgroup(int subsetIndex, uint dof, uint comp)
{
	const uint i = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	return m_dofGroupMap[i][subsetIndex].set_comp_in_dofgroup(dof, comp);
}

template <typename TElem>
inline
uint
DoFPattern::SingleDiscreteFunction::
comp_in_dofgroup(int subsetIndex, uint dof) const
{
	const uint i = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	return m_dofGroupMap[i][subsetIndex].comp_in_dofgroup(dof);
}


}

#endif /* __H__LIBDISCRETIZATION__DOFPATTERN_IMPL__ */
