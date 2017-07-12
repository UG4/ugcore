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

#ifndef __H__LIBGRID__SUBSET_HANDLER_GRID_IMPL__
#define __H__LIBGRID__SUBSET_HANDLER_GRID_IMPL__

#include <cassert>

namespace ug
{

template <class TElem>
typename geometry_traits<TElem>::iterator
GridSubsetHandler::
begin(int subsetIndex)
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
					typename Grid::traits<TElem>::AttachedElementList::iterator());

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				section_container<TElem>(subsetIndex).begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				section_container<TElem>(subsetIndex).section_begin(sectionInd));
}

template <class TElem>
typename geometry_traits<TElem>::iterator
GridSubsetHandler::
end(int subsetIndex)
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
					typename Grid::traits<TElem>::AttachedElementList::iterator());

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				section_container<TElem>(subsetIndex).end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				section_container<TElem>(subsetIndex).section_end(sectionInd));
}

template <class TElem>
typename geometry_traits<TElem>::const_iterator
GridSubsetHandler::
begin(int subsetIndex) const
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
					typename Grid::traits<TElem>::AttachedElementList::iterator());

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				section_container<TElem>(subsetIndex).begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				section_container<TElem>(subsetIndex).section_begin(sectionInd));
}

template <class TElem>
typename geometry_traits<TElem>::const_iterator
GridSubsetHandler::
end(int subsetIndex) const
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
					typename Grid::traits<TElem>::AttachedElementList::iterator());

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				section_container<TElem>(subsetIndex).end());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				section_container<TElem>(subsetIndex).section_end(sectionInd));
}


template <class TElem>
void
GridSubsetHandler::
clear_subset_elements(int subsetIndex)
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return;

//	iterate through the elements of type TElem and erase them from the subsets list.
	if(m_pGrid != NULL)
	{
		typename Grid::traits<TElem>::SectionContainer& secCon =
											section_container<TElem>(subsetIndex);

		typename geometry_traits<TElem>::iterator iter = begin<TElem>(subsetIndex);
		while(iter != end<TElem>(subsetIndex))
		{
			typename geometry_traits<TElem>::iterator iterErase = iter++;
			alter_subset_index(*iterErase) = -1;
			secCon.erase(iterErase, sectionInd);
		}
	}
}

template <class TElem>
uint
GridSubsetHandler::
num_elements(int subsetIndex) const
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if((subsetIndex < 0) || (subsetIndex >= (int)num_subsets_in_list()))
		return 0;
		
	if(sectionInd < 0)
		return section_container<TElem>(subsetIndex).num_elements();
	else
		return section_container<TElem>(subsetIndex).num_elements(sectionInd);
}

template <class TElem>
uint
GridSubsetHandler::
num(int subsetIndex) const
{
	const int sectionInd = geometry_traits<TElem>::CONTAINER_SECTION;

	if((subsetIndex < 0) || (subsetIndex >= (int)num_subsets_in_list()))
		return 0;

	if(sectionInd < 0)
		return section_container<TElem>(subsetIndex).num_elements();
	else
		return section_container<TElem>(subsetIndex).num_elements(sectionInd);
}

template <class TElem>
uint
GridSubsetHandler::
num() const
{
	uint n = 0;
	for(uint i = 0; i < num_subsets_in_list(); ++i)
		n += num<TElem>(i);
		
	return n;
}

template <class TElem> inline
bool GridSubsetHandler::
empty() const
{
	return num<TElem>() == 0;
}

inline bool GridSubsetHandler::
empty() const
{
	return empty<Vertex>() && empty<Edge>()
		   && empty<Face>() && empty<Volume>();
}

template <class TElem> inline
bool GridSubsetHandler::
empty(int subsetIndex) const
{
	return num<TElem>(subsetIndex) == 0;
}

inline bool GridSubsetHandler::
empty(int subsetIndex) const
{
	return empty<Vertex>(subsetIndex) && empty<Edge>(subsetIndex)
		   && empty<Face>(subsetIndex) && empty<Volume>(subsetIndex);
}

template<class TElem>
void GridSubsetHandler::
change_elem_subset_indices(int indOld, int indNew)
{
	typedef typename geometry_traits<TElem>::iterator iterator;
	for(iterator iter = begin<TElem>(indOld);
		iter != end<TElem>(indOld); iter++){
		ISubsetHandler::alter_subset_index(*iter, indNew);
	}
}

template <class TElem>
bool GridSubsetHandler::perform_self_tests()
{
	typedef typename geometry_traits<TElem>::iterator iterator;
	
	bool bSuccess = true;
	
	LOG("performing self tests on GridSubsetHandler\n");
	LOG("  num subets: " << num_subsets_in_list() << std::endl);
	
//	iterate through the subsets and check whether the assigned
//	elements have the correct subset index
	LOG("  checking subset indices\n");
	for(size_t i = 0; i < num_subsets_in_list(); ++i){
		LOG("  checking subset " << i);
		for(iterator iter = begin<TElem>(i); iter != end<TElem>(i); ++iter)
		{
			if(get_subset_index(*iter) != i){
				LOG(" bad element subset index: "
					<< get_subset_index(*iter));
				bSuccess = false;
				break;
			}
		}
		LOG("\n");
	}
	
	return bSuccess;
}

template <class TElem>
typename Grid::traits<TElem>::SectionContainer&
GridSubsetHandler::
section_container(int si)
{
	Subset* sub = m_subsets[si];
	return SectionContainerSelector<typename geometry_traits<TElem>::grid_base_object>::
			section_container(sub->m_vertices, sub->m_edges, sub->m_faces, sub->m_volumes);
}


template <class TElem>
const typename Grid::traits<TElem>::SectionContainer&
GridSubsetHandler::
section_container(int si) const
{
	const Subset* sub = m_subsets[si];
	return SectionContainerSelector<typename geometry_traits<TElem>::grid_base_object>::
			section_container(sub->m_vertices, sub->m_edges, sub->m_faces, sub->m_volumes);
}

}//	end of namespace

#endif

