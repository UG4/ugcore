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
	ISelector(supportedElements),
	m_aSharedEntryVRT("MGSelector_SharedListEntryVRT"),
	m_aSharedEntryEDGE("MGSelector_SharedListEntryEDGE"),
	m_aSharedEntryFACE("MGSelector_SharedListEntryFACE"),
	m_aSharedEntryVOL("MGSelector_SharedListEntryVOL")
{
	m_pMultiGrid = NULL;
}

MGSelector::MGSelector(MultiGrid& grid, uint supportedElements) :
	ISelector(supportedElements), m_pMultiGrid(NULL),
	m_aSharedEntryVRT("MGSelector_SharedListEntryVRT"),
	m_aSharedEntryEDGE("MGSelector_SharedListEntryEDGE"),
	m_aSharedEntryFACE("MGSelector_SharedListEntryFACE"),
	m_aSharedEntryVOL("MGSelector_SharedListEntryVOL")
{
	assign_grid(&grid);
}

MGSelector::~MGSelector()
{
	cleanup();
}

void MGSelector::cleanup()
{
//	delete all levels
	for(size_t i = 0; i < m_levels.size(); ++i)
		delete m_levels[i];

	m_levels.clear();

//	detach shared entry-attachments of section containers
	if(m_pMultiGrid){
		disable_element_support(m_supportedElements);
	//	unregister the previously registered callback
		m_callbackId = MessageHub::SPCallbackId(NULL);

		m_pMultiGrid = NULL;
	}

	BaseClass::set_grid(NULL);
}

void MGSelector::assign_grid(MultiGrid& grid)
{
	assign_grid(&grid);
}

void MGSelector::assign_grid(MultiGrid* grid)
{
	if(grid != m_pMultiGrid){
		uint elementSupport = m_supportedElements;
	//	if a grid already exists, we'll perform cleanup
		if(m_pMultiGrid)
			cleanup();

		m_supportedElements = SE_NONE;
		m_pMultiGrid = grid;
		BaseClass::set_grid(grid);

	//	attach shared entry-attachments to section containers
		if(m_pMultiGrid){
			enable_element_support(elementSupport);
		//	register the callback
			m_callbackId = m_pMultiGrid->message_hub()->register_class_callback(
								this, &ug::MGSelector::multigrid_changed);
			level_required(m_pMultiGrid->num_levels());
		}
		m_supportedElements = elementSupport;
	}
}

void MGSelector::set_supported_elements(uint shElements)
{
//	do this in two steps:
//	1: disable the element-support that is no longer required.
//	2: enable the element-support that was not already enabled.
//	disable the elements that shall be disabled.

//	(the ones which shall not be set, but are currently active.)
	disable_element_support((~shElements) & m_supportedElements);

//	enable the elements that are not already enabled
	enable_element_support(shElements & (~m_supportedElements));
}

void MGSelector::enable_element_support(uint shElements)
{
	if((shElements & SE_VERTEX) && (!elements_are_supported(SE_VERTEX)))
		m_pMultiGrid->attach_to_vertices(m_aSharedEntryVRT);

	if((shElements & SE_EDGE) && (!elements_are_supported(SE_EDGE)))
		m_pMultiGrid->attach_to_edges(m_aSharedEntryEDGE);

	if((shElements & SE_FACE) && (!elements_are_supported(SE_FACE)))
		m_pMultiGrid->attach_to_faces(m_aSharedEntryFACE);

	if((shElements & SE_VOLUME) && (!elements_are_supported(SE_VOLUME)))
		m_pMultiGrid->attach_to_volumes(m_aSharedEntryVOL);

	for(size_t i = 0; i < m_levels.size(); ++i){
		Level& lvl = *m_levels[i];
		if((shElements & SE_VERTEX) && (!elements_are_supported(SE_VERTEX)))
			lvl.m_vertices.get_container().set_pipe(
				&m_pGrid->get_attachment_pipe<Vertex>(), m_aSharedEntryVRT);

		if((shElements & SE_EDGE) && (!elements_are_supported(SE_EDGE)))
			lvl.m_edges.get_container().set_pipe(
				&m_pGrid->get_attachment_pipe<EdgeBase>(), m_aSharedEntryEDGE);

		if((shElements & SE_FACE) && (!elements_are_supported(SE_FACE)))
			lvl.m_faces.get_container().set_pipe(
				&m_pGrid->get_attachment_pipe<Face>(), m_aSharedEntryFACE);

		if((shElements & SE_VOLUME) && (!elements_are_supported(SE_VOLUME)))
			lvl.m_volumes.get_container().set_pipe(
				&m_pGrid->get_attachment_pipe<Volume>(), m_aSharedEntryVOL);
	}
	ISelector::enable_element_support(shElements);
}

void MGSelector::disable_element_support(uint shElements)
{
	//	release the attachments in the current grid
	for(size_t i = 0; i < m_levels.size(); ++i){
		Level& lvl = *m_levels[i];
		if((shElements & SE_VERTEX) && elements_are_supported(SE_VERTEX))
			lvl.m_vertices.get_container().set_pipe(NULL);

		if((shElements & SE_EDGE) && elements_are_supported(SE_EDGE))
			lvl.m_edges.get_container().set_pipe(NULL);

		if((shElements & SE_FACE) && elements_are_supported(SE_FACE))
			lvl.m_faces.get_container().set_pipe(NULL);

		if((shElements & SE_VOLUME) && elements_are_supported(SE_VOLUME))
			lvl.m_volumes.get_container().set_pipe(NULL);
	}

	if((shElements & SE_VERTEX) && elements_are_supported(SE_VERTEX))
		m_pMultiGrid->detach_from_vertices(m_aSharedEntryVRT);

	if((shElements & SE_EDGE) && elements_are_supported(SE_EDGE))
		m_pMultiGrid->detach_from_edges(m_aSharedEntryEDGE);

	if((shElements & SE_FACE) && elements_are_supported(SE_FACE))
		m_pMultiGrid->detach_from_faces(m_aSharedEntryFACE);

	if((shElements & SE_VOLUME) && elements_are_supported(SE_VOLUME))
		m_pMultiGrid->detach_from_volumes(m_aSharedEntryVOL);

	ISelector::disable_element_support(shElements);
}

void MGSelector::add_level()
{
//	adds a level and and initializes associated section containers.
	Level* pLvl = new Level;
	if(elements_are_supported(SE_VERTEX))
		pLvl->m_vertices.get_container().set_pipe(
				&m_pGrid->get_attachment_pipe<Vertex>(), m_aSharedEntryVRT);
	if(elements_are_supported(SE_EDGE))
		pLvl->m_edges.get_container().set_pipe(
				&m_pGrid->get_attachment_pipe<EdgeBase>(), m_aSharedEntryEDGE);
	if(elements_are_supported(SE_FACE))
		pLvl->m_faces.get_container().set_pipe(
				&m_pGrid->get_attachment_pipe<Face>(), m_aSharedEntryFACE);
	if(elements_are_supported(SE_VOLUME))
		pLvl->m_volumes.get_container().set_pipe(
				&m_pGrid->get_attachment_pipe<Volume>(), m_aSharedEntryVOL);

	m_levels.push_back(pLvl);
}

void MGSelector::clear_lists()
{
	for(size_t i = 0; i < m_levels.size(); ++i)
		delete m_levels[i];
	m_levels.clear();
}

void MGSelector::clear()
{
	clear<Vertex>();
	clear<EdgeBase>();
	clear<Face>();
	clear<Volume>();
}

void MGSelector::clear(int level)
{
	clear<Vertex>(level);
	clear<EdgeBase>(level);
	clear<Face>(level);
	clear<Volume>(level);
}

void MGSelector::add_to_list(Vertex* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	section_container<Vertex>(level).insert(elem,
								elem->container_section());
}

void MGSelector::add_to_list(EdgeBase* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	section_container<EdgeBase>(level).insert(elem,
								elem->container_section());
}

void MGSelector::add_to_list(Face* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	section_container<Face>(level).insert(elem,
								elem->container_section());
}

void MGSelector::add_to_list(Volume* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	section_container<Volume>(level).insert(elem,
								elem->container_section());
}	

void MGSelector::erase_from_list(Vertex* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	section_container<Vertex>(level).erase(get_level_iterator(elem),
						elem->container_section());
}
void MGSelector::erase_from_list(EdgeBase* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	section_container<EdgeBase>(level).erase(get_level_iterator(elem),
						elem->container_section());
}
void MGSelector::erase_from_list(Face* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	section_container<Face>(level).erase(get_level_iterator(elem),
						elem->container_section());
}

void MGSelector::erase_from_list(Volume* elem)
{
	const int level = m_pMultiGrid->get_level(elem);
	section_container<Volume>(level).erase(get_level_iterator(elem),
						elem->container_section());
}

//	geometric-object-collection
GridObjectCollection MGSelector::get_grid_objects() const
{
	uint numLevels = num_levels();
	GridObjectCollection goc(numLevels);
	
	for(uint i = 0; i < numLevels; ++i)
	{
		goc.add_level(	&m_levels[i]->m_vertices,
						&m_levels[i]->m_edges,
						&m_levels[i]->m_faces,
						&m_levels[i]->m_volumes);
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

void MGSelector::
multigrid_changed(const GridMessage_MultiGridChanged& gm)
{
	if(gm.message_type() == GMMGCT_LEVEL_ADDED)
		level_required(gm.num_levels_in_grid() - 1);
}

}//	end of namespace
