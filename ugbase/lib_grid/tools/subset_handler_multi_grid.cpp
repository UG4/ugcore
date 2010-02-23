//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24

#include <cassert>
#include "subset_handler_multi_grid.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	MultiGridSubsetHandler implementation
MultiGridSubsetHandler::
MultiGridSubsetHandler(uint supportedElements) : ISubsetHandler(supportedElements)
{
	m_numSubsets = 0;
	m_pMG = NULL;
}

MultiGridSubsetHandler::
MultiGridSubsetHandler(MultiGrid& mg, uint supportedElements) : ISubsetHandler(supportedElements)
{
	m_numSubsets = 0;
	m_pMG = &mg;
	assign_grid(mg);
}

MultiGridSubsetHandler::
MultiGridSubsetHandler(const MultiGridSubsetHandler& sh) :
	ISubsetHandler(sh.m_supportedElements)
{
	MultiGrid* pGrid = sh.m_pMG;

	if(pGrid){
		assign_grid(*pGrid);		
		assign_subset_handler(sh);
	}
}

MultiGridSubsetHandler::~MultiGridSubsetHandler()
{
	if(m_pGrid != NULL)
		m_pGrid->unregister_observer(this);
}


void MultiGridSubsetHandler::erase_subset_lists()
{
	for(size_t level = 0; level < m_levels.size(); ++level)
	{
		for(size_t i = 0; i < m_levels[level].size(); ++i)
			delete m_levels[level][i];

		m_levels[level].clear();
	}

	m_numSubsets = 0;
}

void MultiGridSubsetHandler::clear_subset_lists(int index)
{
	if(m_pGrid)
	{
		for(int i = 0; i < NUM_GEOMETRIC_BASE_OBJECTS; ++i)
		{
			for(size_t level = 0; level < m_levels.size(); ++level)
				m_levels[level][index]->m_elements[i].clear();
		}
	}
}

template<class TElemPtr>
void
MultiGridSubsetHandler::
assign_subset(TElemPtr elem, int subsetIndex, int elemType)
{
	assert((m_pGrid != NULL) && "ERROR in SubsetHandler::assign_subset(): No grid assigned to SubsetHandler.");

	int level = m_pMG->get_level(elem);
	subset_required(subsetIndex);
	level_required(level);

//	check if we have to remove elem from a subset.
	int oldIndex = get_subset_index(elem);
	
	if(oldIndex != -1)
		m_levels[level][oldIndex]->m_elements[elemType].erase(get_list_iterator(elem), elem->shared_pipe_section());

//	add the element to the subset.
	if(subsetIndex != -1)
	{
		ISubsetHandler::iterator iter = m_levels[level][subsetIndex]->m_elements[elemType].insert(elem, elem->shared_pipe_section());
		subset_assigned(elem, iter, subsetIndex);
	}
}

void MultiGridSubsetHandler::assign_subset(VertexBase* elem, int subsetIndex)
{
	assign_subset(elem, subsetIndex, VERTEX);
}

void MultiGridSubsetHandler::assign_subset(EdgeBase* elem, int subsetIndex)
{
	assign_subset(elem, subsetIndex, EDGE);
}

void MultiGridSubsetHandler::assign_subset(Face* elem, int subsetIndex)
{
	assign_subset(elem, subsetIndex, FACE);
}

void MultiGridSubsetHandler::assign_subset(Volume* elem, int subsetIndex)
{
	assign_subset(elem, subsetIndex, VOLUME);
}

void MultiGridSubsetHandler::
change_subset_indices(int indOld, int indNew)
{
	if(m_pGrid)
	{
		if(elements_are_supported(VERTEX))
			change_elem_subset_indices<VertexBase>(indOld, indNew);
		if(elements_are_supported(EDGE))
			change_elem_subset_indices<EdgeBase>(indOld, indNew);
		if(elements_are_supported(FACE))
			change_elem_subset_indices<Face>(indOld, indNew);
		if(elements_are_supported(VOLUME))
			change_elem_subset_indices<Volume>(indOld, indNew);
	}
}

void MultiGridSubsetHandler::add_required_subset_lists(int maxIndex)
{
	while((int)num_subsets() <= maxIndex)
		add_subset_to_all_levels();
}

void MultiGridSubsetHandler::erase_subset_lists(int index)
{
	for(size_t level = 0; level < m_levels.size(); ++level)
	{
		SubsetVec& subsets = m_levels[level];
		delete subsets[index];
		for(uint i = index + 1; i < num_subsets(); ++i)
			subsets[i-1] = subsets[i];
		subsets.resize(subsets.size() - 1);
	}
	m_numSubsets--;
}

void MultiGridSubsetHandler::swap_subset_lists(int ind1, int ind2)
{
	for(size_t level = 0; level < m_levels.size(); ++level)
	{
		SubsetVec& subsets = m_levels[level];
		Subset* pTmp = subsets[ind1];
		subsets[ind1] = subsets[ind2];
		subsets[ind2] = pTmp;
	}
}

void MultiGridSubsetHandler::move_subset_lists(int indexFrom, int indexTo)
{
	int moveDir = 0;
//	we have to distinguish two cases
	if(indexFrom < indexTo)
		moveDir = 1;
	else if(indexTo < indexFrom)
		moveDir = -1;

	if(moveDir != 0)
	{
		for(size_t level = 0; level < m_levels.size(); ++level)
		{
			SubsetVec& subsets = m_levels[level];

		//	store pointer to the from-subset
			Subset* pFrom = subsets[indexFrom];

			for(int i = indexFrom; i != indexTo; i+= moveDir)
			{
				int iNext = i+moveDir;

			//	move pointer
				subsets[i] = subsets[iNext];
			}

		//	assign stored pointer
			subsets[indexTo] = pFrom;
		}
	}
}

void MultiGridSubsetHandler::
register_subset_elements_at_pipe()
{
	for(int l = 0; l < num_levels(); ++l)
	{
		for(int i = 0; i < num_subsets(); ++i)
		{
		//	register vertices
			for(VertexBaseIterator iter = begin<VertexBase>(i, l);
				iter != end<VertexBase>(i, l); ++iter)
				register_at_pipe(*iter);

		//	register edges
			for(EdgeBaseIterator iter = begin<EdgeBase>(i, l);
				iter != end<EdgeBase>(i, l); ++iter)
				register_at_pipe(*iter);

		//	register faces
			for(FaceIterator iter = begin<Face>(i, l);
				iter != end<Face>(i, l); ++iter)
				register_at_pipe(*iter);

		//	register volumes
			for(VolumeIterator iter = begin<Volume>(i, l);
				iter != end<Volume>(i, l); ++iter)
				register_at_pipe(*iter);
		}
	}
}


GeometricObjectCollection
MultiGridSubsetHandler::
get_goc(int subsetIndex, int level)
{
	subset_required(subsetIndex);
	level_required(level);

	return GeometricObjectCollection(&m_levels[level][subsetIndex]->m_elements[VERTEX],
									 &m_levels[level][subsetIndex]->m_elements[EDGE],
									 &m_levels[level][subsetIndex]->m_elements[FACE],
									 &m_levels[level][subsetIndex]->m_elements[VOLUME]);
}

MultiLevelGeometricObjectCollection
MultiGridSubsetHandler::
get_mlgoc_by_subset(int subsetIndex)
{
	subset_required(subsetIndex);
	MultiLevelGeometricObjectCollection mgoc(m_levels.size());
	for(size_t i = 0; i < m_levels.size(); ++i)
	{
		mgoc.add_level(	&m_levels[i][subsetIndex]->m_elements[VERTEX],
						&m_levels[i][subsetIndex]->m_elements[EDGE],
						&m_levels[i][subsetIndex]->m_elements[FACE],
						&m_levels[i][subsetIndex]->m_elements[VOLUME]);
	}
	
	return mgoc;
}

MultiLevelGeometricObjectCollection
MultiGridSubsetHandler::
get_mlgoc_by_level(int level)
{
	level_required(level);
	uint numSubsets = num_subsets();
	MultiLevelGeometricObjectCollection mgoc(numSubsets);
	for(uint i = 0; i < numSubsets; ++i)
	{
		mgoc.add_level(	&m_levels[level][i]->m_elements[VERTEX],
						&m_levels[level][i]->m_elements[EDGE],
						&m_levels[level][i]->m_elements[FACE],
						&m_levels[level][i]->m_elements[VOLUME]);
	}
	
	return mgoc;
}

void MultiGridSubsetHandler::add_level()
{
	m_levels.push_back(SubsetVec());
	int topLevel = m_levels.size() - 1;
	for(uint i = 0; i < num_subsets(); ++i)
		m_levels[topLevel].push_back(new Subset);
}

void MultiGridSubsetHandler::add_subset_to_all_levels()
{
	for(size_t i = 0; i < m_levels.size(); ++i)
		m_levels[i].push_back(new Subset);

	m_numSubsets++;
}

void MultiGridSubsetHandler::
unregistered_from_grid(Grid* grid)
{
	assert(m_pGrid && "this method should only be called while the handler is registered at a grid.");
	
	if(m_pGrid)
	{
		erase_subset_lists();
	}

	ISubsetHandler::unregistered_from_grid(grid);
}

}//	end of namespace
