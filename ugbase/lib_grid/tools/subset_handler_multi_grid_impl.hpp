//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24   (reworked y09 m12 d15)

#ifndef __H__LIBGRID__SUBSET_HANDLER_MULTI_GRID_IMPL__
#define __H__LIBGRID__SUBSET_HANDLER_MULTI_GRID_IMPL__

#include <cassert>

namespace ug
{

template <class TElem>
typename geometry_traits<TElem>::iterator
MultiGridSubsetHandler::
begin(int subsetIndex, int level)
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::begin(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::begin(): bad element type.");

	level_required(level);

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				subset(subsetIndex, level)->m_elements[baseObjID].begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				subset(subsetIndex, level)->m_elements[baseObjID].section_begin(sectionInd));
}

template <class TElem>
typename geometry_traits<TElem>::iterator
MultiGridSubsetHandler::
end(int subsetIndex, int level)
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::end(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::end(): bad element type.");

	level_required(level);

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				subset(subsetIndex, level)->m_elements[baseObjID].end());
	else
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
				subset(subsetIndex, level)->m_elements[baseObjID].section_end(sectionInd));
}

template <class TElem>
void
MultiGridSubsetHandler::
clear_subset_elements(int subsetIndex)
{
	for(int i = 0; i < (int)num_levels(); ++i)
		clear_subset_elements<TElem>(subsetIndex, i);
}

template <class TElem>
void MultiGridSubsetHandler::
clear_subset_elements(int subsetIndex, int level)
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
		Subset* pSubset = subset(subsetIndex, level);
		typename geometry_traits<TElem>::iterator iter = begin<TElem>(subsetIndex, level);
		while(iter != end<TElem>(subsetIndex, level))
		{
			typename geometry_traits<TElem>::iterator iterErase = iter++;
			alter_subset_index(*iterErase) = -1;
			pSubset->m_elements[baseObjID].erase(iterErase, sectionInd);
		}
	}
}

template <class TElem>
uint
MultiGridSubsetHandler::
num(int subsetIndex, int level) const
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::num_elements(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::num_elements(): bad element type.");

	// TODO: Maybe the passed level should be of uint-type
	// if level does not exist, there are no elements in it
	if((uint)level >= num_levels()) return 0;

	if(sectionInd < 0)
		return subset(subsetIndex, level)->m_elements[baseObjID].num_elements();
	else
		return subset(subsetIndex, level)->m_elements[baseObjID].num_elements(sectionInd);
}

template <class TElem>
uint
MultiGridSubsetHandler::
num(int subsetIndex) const
{
	int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	assert((subsetIndex >= 0) && (subsetIndex < (int)num_subsets()) &&
			"ERROR in SubsetHandler::num_elements(): bad subset index.");
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::num_elements(): bad element type.");

	uint numElems = 0;
	if(sectionInd < 0)
	{
		for(size_t i = 0; i < m_levels.size(); ++i)
			numElems += subset(subsetIndex, i)->m_elements[baseObjID].num_elements();
	}
	else
	{
		for(size_t i = 0; i < m_levels.size(); ++i)
			numElems += subset(subsetIndex, i)->m_elements[baseObjID].num_elements(sectionInd);
	}

	return numElems;
}

template <class TElem>
uint
MultiGridSubsetHandler::
num() const
{
	uint n = 0;
	for(size_t i = 0; i < num_subsets(); ++i)
		n += num<TElem>(i);

	return n;
}

template<class TElem>
void MultiGridSubsetHandler::
change_elem_subset_indices(int indOld, int indNew)
{
	typedef typename geometry_traits<TElem>::iterator iterator;

	for(size_t i = 0; i < m_levels.size(); ++i)
	{
		for(iterator iter = begin<TElem>(indOld, i);
			iter != end<TElem>(indOld, i); iter++)
			alter_subset_index(*iter, indNew);
	}
}

}//	end of namespace

#endif

