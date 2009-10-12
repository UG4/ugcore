//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24

#ifndef __H__LIBGRID__SUBSET_HANDLER__IMPL__
#define __H__LIBGRID__SUBSET_HANDLER__IMPL__

#include <cassert>

namespace ug
{

template <class TElem>
void
SubsetHandler::
clear_subsets_elements(int subsetIndex)
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::clear_subsets_elements(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::clear_subsets_elements(): bad element type.");

//	iterate through the elements of type TElem and erase them from the subsets list.
	if(m_pGrid != NULL)
	{
		Grid::AttachmentAccessor<TElem, ASubsetIndex> aaSI(*m_pGrid, m_aSubsetIndex);
		typename geometry_traits<TElem>::iterator iter = begin<TElem>(subsetIndex);
		while(iter != end<TElem>(subsetIndex))
		{
			typename geometry_traits<TElem>::iterator iterErase = iter++;
			aaSI[*iterErase] = -1;
			m_subsets[subsetIndex]->m_elements[baseObjID].erase(iterErase, sectionInd);
		}
	}
}

template <class TElem>
uint
SubsetHandler::
num_elements(int subsetIndex)
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::num_elements(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::num_elements(): bad element type.");

	if(sectionInd < 0)
		return m_subsets[subsetIndex]->m_elements[baseObjID].num_elements();
	else
		return m_subsets[subsetIndex]->m_elements[baseObjID].num_elements(sectionInd);
}

template <class TElem>
uint
SubsetHandler::
num(int subsetIndex)
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::num_elements(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::num_elements(): bad element type.");

	if(sectionInd < 0)
		return m_subsets[subsetIndex]->m_elements[baseObjID].num_elements();
	else
		return m_subsets[subsetIndex]->m_elements[baseObjID].num_elements(sectionInd);
}

template <class TIterator>
void
SubsetHandler::
assign_subset(TIterator iterBegin, TIterator iterEnd, int subsetIndex)
{
	if(subsetIndex >= (int)num_subsets())
		resize_subset_vec(subsetIndex + 1);

	for(TIterator iter = iterBegin; iter != iterEnd; iter++)
		assign_subset(*iter, subsetIndex);
}

template <class TElem>
typename geometry_traits<TElem>::iterator
SubsetHandler::
begin(int subsetIndex)
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::begin(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::begin(): bad element type.");

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				m_subsets[subsetIndex]->m_elements[baseObjID].begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				m_subsets[subsetIndex]->m_elements[baseObjID].section_begin(sectionInd));
}

template <class TElem>
typename geometry_traits<TElem>::iterator
SubsetHandler::
end(int subsetIndex)
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::end(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::end(): bad element type.");

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				m_subsets[subsetIndex]->m_elements[baseObjID].end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				m_subsets[subsetIndex]->m_elements[baseObjID].section_end(sectionInd));
}

template <class TElem>
void
SubsetHandler::
reset_subset_indices(typename geometry_traits<TElem>::iterator iterBegin,
									typename geometry_traits<TElem>::iterator iterEnd)
{
	if(m_pGrid != NULL)
	{
		Grid::AttachmentAccessor<TElem, ASubsetIndex> aaSI(*m_pGrid, m_aSubsetIndex);
	//	set all subset indices to -1
		for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
			aaSI[*iter] = -1;
	}
}

template <class TElem>
void
SubsetHandler::
set_subset_indices(typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd, int subsetIndex)
{
//	does not alter any iterators
	if(m_pGrid != NULL)
	{
		Grid::AttachmentAccessor<TElem, ASubsetIndex> aaSI(*m_pGrid, m_aSubsetIndex);
	//	set all subset indices to -1
		for(typename geometry_traits<TElem>::iterator iter = iterBegin; iter != iterEnd; iter++)
			aaSI[*iter] = subsetIndex;
	}
}

}//	end of namespace

#endif

