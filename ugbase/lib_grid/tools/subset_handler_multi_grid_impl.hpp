/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIBGRID__SUBSET_HANDLER_MULTI_GRID_IMPL__
#define __H__LIBGRID__SUBSET_HANDLER_MULTI_GRID_IMPL__

#include "subset_handler_multi_grid.h"

#include <cassert>

namespace ug {

template <typename TElem>
typename geometry_traits<TElem>::iterator
MultiGridSubsetHandler::
begin(int subsetIndex, int level)
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets_in_list()) &&
			"ERROR in SubsetHandler::begin(): bad subset index.");

	level_required(level);

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				section_container<TElem>(subsetIndex, level).begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				section_container<TElem>(subsetIndex, level).section_begin(sectionInd));
}

template <typename TElem>
typename geometry_traits<TElem>::iterator
MultiGridSubsetHandler::
end(int subsetIndex, int level)
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets_in_list()) &&
			"ERROR in SubsetHandler::end(): bad subset index.");

	level_required(level);

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				section_container<TElem>(subsetIndex, level).end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				section_container<TElem>(subsetIndex, level).section_end(sectionInd));
}

template <typename TElem>
typename geometry_traits<TElem>::const_iterator
MultiGridSubsetHandler::
begin(int subsetIndex, int level) const
{
	UG_ASSERT(subsetIndex >= 0, "-1 is not a valid subset index when accessing iterators!");
	if(subsetIndex >= static_cast<int>(num_subsets_in_list()) || level >= static_cast<int>(num_levels()))
		return end<TElem>(subsetIndex, level);

	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				section_container<TElem>(subsetIndex, level).begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				section_container<TElem>(subsetIndex, level).section_begin(sectionInd));
}

template <typename TElem>
typename geometry_traits<TElem>::const_iterator
MultiGridSubsetHandler::
end(int subsetIndex, int level) const
{
	UG_ASSERT(subsetIndex >= 0, "-1 is not a valid subset index when accessing iterators!");
	if(subsetIndex >= static_cast<int>(num_subsets_in_list()) || level >= static_cast<int>(num_levels())){
		static typename geometry_traits<TElem>::const_iterator dummyEndIter;
		return dummyEndIter;
	}

	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				section_container<TElem>(subsetIndex, level).end());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				section_container<TElem>(subsetIndex, level).section_end(sectionInd));
}
template <typename TElem>
void
MultiGridSubsetHandler::
clear_subset_elements(int subsetIndex)
{
	for(int i = 0; i < static_cast<int>(num_levels()); ++i)
		clear_subset_elements<TElem>(subsetIndex, i);
}

template <typename TElem>
void MultiGridSubsetHandler::
clear_subset_elements(int subsetIndex, int level)
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets_in_list()) &&
			"ERROR in SubsetHandler::clear_subsets_elements(): bad subset index.");

//	iterate through the elements of type TElem and erase them from the subsets list.
	if(m_pGrid != nullptr)
	{
		typename Grid::traits<TElem>::SectionContainer& secCon =
									section_container<TElem>(subsetIndex, level);

		typename geometry_traits<TElem>::iterator iter = begin<TElem>(subsetIndex, level);
		while(iter != end<TElem>(subsetIndex, level))
		{
			typename geometry_traits<TElem>::iterator iterErase = iter++;
			alter_subset_index(*iterErase) = -1;
			secCon.erase(iterErase, sectionInd);
		}
	}
}

template <typename TElem>
uint
MultiGridSubsetHandler::
num(int subsetIndex, int level) const
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets_in_list()) &&
			"ERROR in SubsetHandler::num_elements(): bad subset index.");

	// TODO: Maybe the passed level should be of uint-type
	// if level does not exist, there are no elements in it
	if(static_cast<uint>(level) >= num_levels()) return 0;

	if(sectionInd < 0)
		return section_container<TElem>(subsetIndex, level).num_elements();
	else
		return section_container<TElem>(subsetIndex, level).num_elements(sectionInd);
}

template <typename TElem>
uint
MultiGridSubsetHandler::
num(int subsetIndex) const
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets_in_list()) &&
			"ERROR in SubsetHandler::num_elements(): bad subset index.");

	uint numElems = 0;
	if(sectionInd < 0)
	{
		for(size_t i = 0; i < m_levels.size(); ++i)
			numElems += section_container<TElem>(subsetIndex, static_cast<int>(i)).num_elements();
	}
	else
	{
		for(size_t i = 0; i < m_levels.size(); ++i)
			numElems += section_container<TElem>(subsetIndex, static_cast<int>(i)).num_elements(sectionInd);
	}

	return numElems;
}

template <typename TElem>
uint
MultiGridSubsetHandler::
num() const
{
	uint n = 0;
	for(size_t i = 0; i < num_subsets_in_list(); ++i)
		n += num<TElem>(i);

	return n;
}

template <typename TElem>
void MultiGridSubsetHandler::
change_elem_subset_indices(int indOld, int indNew)
{
	using iterator = typename geometry_traits<TElem>::iterator;

	for(size_t i = 0; i < m_levels.size(); ++i)
	{
		for(iterator iter = begin<TElem>(indOld, i);
			iter != end<TElem>(indOld, i); ++iter)
			alter_subset_index(*iter, indNew);
	}
}

inline void MultiGridSubsetHandler::
level_required(int level)
{
	while(static_cast<int>(m_levels.size()) <= level) {
		add_level();
	}
}

inline void MultiGridSubsetHandler::
level_required(int level) const
{
	if(level >= static_cast<int>(num_levels())){
		UG_THROW("Can't create additional levels in const MGSubsetHandler. "
						<< "num current levels: " << num_levels()
						<< " required level: " << level);
	}
}

template <typename TElem>
typename Grid::traits<TElem>::SectionContainer&
MultiGridSubsetHandler::
section_container(int si, int lvl)
{
	Subset* sub = subset(si, lvl);
	return SectionContainerSelector<typename geometry_traits<TElem>::grid_base_object>::
			section_container(sub->m_vertices, sub->m_edges, sub->m_faces, sub->m_volumes);
}


template <typename TElem>
const typename Grid::traits<TElem>::SectionContainer&
MultiGridSubsetHandler::
section_container(int si, int lvl) const
{
	const Subset* sub = subset(si, lvl);
	return SectionContainerSelector<typename geometry_traits<TElem>::grid_base_object>::
			section_container(sub->m_vertices, sub->m_edges, sub->m_faces, sub->m_volumes);
}

}//	end of namespace

#endif

