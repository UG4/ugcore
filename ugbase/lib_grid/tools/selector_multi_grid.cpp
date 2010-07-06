// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

#include <cassert>
#include "selector_multi_grid.h"
#include "lib_grid/multi_grid.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
MGSelector::MGSelector(uint supportedElements) :
	ISelector(supportedElements)
{
	m_pMultiGrid = NULL;
}

MGSelector::MGSelector(MultiGrid& grid, uint supportedElements) :
	ISelector(grid, supportedElements)
{
	m_pMultiGrid = &grid;
}

MGSelector::~MGSelector()
{
}

void MGSelector::assign_grid(MultiGrid& grid)
{
	m_pMultiGrid = &grid;
	BaseClass::set_grid(&grid);
}

void MGSelector::assign_grid(MultiGrid* grid)
{
	m_pMultiGrid = grid;
	BaseClass::set_grid(grid);
}

void MGSelector::clear_lists()
{
	for(size_t i = 0; i < m_levels.size(); ++i)
		delete m_levels[i];
	m_levels.clear();
}

void MGSelector::clear()
{
	clear<VertexBase>();
	clear<EdgeBase>();
	clear<Face>();
	clear<Volume>();
}

void MGSelector::clear(int level)
{
	clear<VertexBase>(level);
	clear<EdgeBase>(level);
	clear<Face>(level);
	clear<Volume>(level);
}

MGSelector::iterator
MGSelector::add_to_list(VertexBase* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	MGSelector::iterator iter = get_section_container<VertexBase>(level).insert(elem,
								elem->shared_pipe_section());

	return iter;
}

MGSelector::iterator 
MGSelector::add_to_list(EdgeBase* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	return get_section_container<EdgeBase>(level).insert(elem,
								elem->shared_pipe_section());
}

MGSelector::iterator 
MGSelector::add_to_list(Face* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	return get_section_container<Face>(level).insert(elem,
								elem->shared_pipe_section());
}

MGSelector::iterator 
MGSelector::add_to_list(Volume* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	return get_section_container<Volume>(level).insert(elem,
								elem->shared_pipe_section());
}	

void MGSelector::erase_from_list(VertexBase* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	get_section_container<VertexBase>(level).erase(get_iterator(elem),
						elem->shared_pipe_section());
}
void MGSelector::erase_from_list(EdgeBase* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	get_section_container<EdgeBase>(level).erase(get_iterator(elem),
						elem->shared_pipe_section());
}
void MGSelector::erase_from_list(Face* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	get_section_container<Face>(level).erase(get_iterator(elem),
						elem->shared_pipe_section());
}

void MGSelector::erase_from_list(Volume* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	get_section_container<Volume>(level).erase(get_iterator(elem),
						elem->shared_pipe_section());
}

//	geometric-object-collection
GeometricObjectCollection 
MGSelector::get_geometric_object_collection()
{
	uint numLevels = num_levels();
	GeometricObjectCollection goc(numLevels);
	
	for(uint i = 0; i < numLevels; ++i)
	{
		goc.add_level(	&m_levels[i]->m_elements[VERTEX],
						&m_levels[i]->m_elements[EDGE],
						&m_levels[i]->m_elements[FACE],
						&m_levels[i]->m_elements[VOLUME]);
	}
	
	return goc;
}

/*
void MGSelector::registered_at_grid(Grid* grid)
{
	m_pMultiGrid = dynamic_cast<MultiGrid*>(grid);
	assert(m_pMultiGrid && "bad grid type!");
	ISelector::registered_at_grid(grid);
}

void MGSelector::unregistered_from_grid(Grid* grid)
{
	ISelector::unregistered_from_grid(grid);
	clear_lists();
}
*/
void MGSelector::grid_to_be_destroyed(Grid* grid)
{
	ISelector::grid_to_be_destroyed(grid);
	m_pMultiGrid = NULL;
	clear_lists();
}

}//	end of namespace
