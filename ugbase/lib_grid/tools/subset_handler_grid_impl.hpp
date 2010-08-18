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
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	const int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
												m_invalidContainer.end());
			
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
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	const int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return iterator_cast<typename geometry_traits<TElem>::iterator>(
												m_invalidContainer.end());

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
typename geometry_traits<TElem>::const_iterator
GridSubsetHandler::
begin(int subsetIndex) const
{
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	const int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
												m_invalidContainer.end());
			
	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::begin(): bad element type.");

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				m_subsets[subsetIndex]->m_elements[baseObjID].begin());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				m_subsets[subsetIndex]->m_elements[baseObjID].section_begin(sectionInd));
}

template <class TElem>
typename geometry_traits<TElem>::const_iterator
GridSubsetHandler::
end(int subsetIndex) const
{
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	const int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
												m_invalidContainer.end());

	assert((baseObjID >= 0) && (baseObjID < NUM_GEOMETRIC_BASE_OBJECTS) &&
			"ERROR in SubsetHandler::end(): bad element type.");

	if(sectionInd < 0)
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				m_subsets[subsetIndex]->m_elements[baseObjID].end());
	else
		return iterator_cast<typename geometry_traits<TElem>::const_iterator>(
				m_subsets[subsetIndex]->m_elements[baseObjID].section_end(sectionInd));
}


template <class TElem>
void
GridSubsetHandler::
clear_subset_elements(int subsetIndex)
{
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	const int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	if(subsetIndex < 0 || subsetIndex >= (int)num_subsets_in_list())
		return;

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
num_elements(int subsetIndex) const
{
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	const int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	if((subsetIndex < 0) || (subsetIndex >= (int)num_subsets_in_list()))
		return 0;
		
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
num(int subsetIndex) const
{
	const int baseObjID = geometry_traits<TElem>::BASE_OBJECT_TYPE_ID;
	const int sectionInd = geometry_traits<TElem>::SHARED_PIPE_SECTION;

	if((subsetIndex < 0) || (subsetIndex >= (int)num_subsets_in_list()))
		return 0;
		
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
num() const
{
	uint n = 0;
	for(size_t i = 0; i < num_subsets_in_list(); ++i)
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
	return empty<VertexBase>() && empty<EdgeBase>()
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
	return empty<VertexBase>(subsetIndex) && empty<EdgeBase>(subsetIndex)
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
	for(int i = 0; i < num_subsets_in_list(); ++i){
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

}//	end of namespace

#endif

