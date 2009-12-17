//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24   (reworked y09 m12 d15)

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
GridSubsetHandler::
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
GridSubsetHandler::
clear_subset_elements(int subsetIndex)
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
		typename geometry_traits<TElem>::iterator iter = begin<TElem>(subsetIndex);
		while(iter != end<TElem>(subsetIndex))
		{
			typename geometry_traits<TElem>::iterator iterErase = iter++;
			alter_subset_index(*iterErase) = -1;
			m_subsets[subsetIndex]->m_elements[baseObjID].erase(iterErase, sectionInd);
		}
	}
}

template <class TElem>
uint
GridSubsetHandler::
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
GridSubsetHandler::
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

template<class TElem>
void GridSubsetHandler::
change_elem_subset_indices(int indOld, int indNew)
{
	typedef typename geometry_traits<TElem>::iterator iterator;
	for(iterator iter = begin<TElem>(indOld);
		iter != end<TElem>(indOld); iter++)
		ISubsetHandler::alter_subset_index(*iter, indNew);
}

}//	end of namespace

#endif

