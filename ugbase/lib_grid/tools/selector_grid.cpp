// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d15

#include <cassert>
#include "selector_grid.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
Selector::Selector(uint supportedElements) :
	ISelector(supportedElements)
{
}

Selector::Selector(Grid& grid, uint supportedElements) :
	ISelector(supportedElements)
{
	assign_grid(&grid);
}

Selector::~Selector()
{
	if(m_pGrid){
	//	release the attachments in the current grid
		if(elements_are_supported(SE_VERTEX))
			section_container<VertexBase>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_EDGE))
			section_container<EdgeBase>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_FACE))
			section_container<Face>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_VOLUME))
			section_container<Volume>().get_container().set_pipe(NULL);
	}
}

void Selector::assign_grid(Grid& grid)
{
	assign_grid(&grid);
}

void Selector::assign_grid(Grid* grid)
{
	if(grid != m_pGrid){
		uint elementSupport = m_supportedElements;

		if(m_pGrid){
		//	release the attachments in the current grid
			disable_element_support(elementSupport);
		}

		m_supportedElements = SE_NONE;
		BaseClass::set_grid(grid);

		if(m_pGrid){
		//	initialize attachment lists
			enable_element_support(elementSupport);
		}
	}
}

void Selector::set_supported_elements(uint shElements)
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

void Selector::enable_element_support(uint shElements)
{
	if((shElements & SE_VERTEX) && (!elements_are_supported(SE_VERTEX)))
		section_container<VertexBase>().get_container().
				set_pipe(&m_pGrid->get_attachment_pipe<VertexBase>());

	if((shElements & SE_EDGE) && (!elements_are_supported(SE_EDGE)))
		section_container<EdgeBase>().get_container().
				set_pipe(&m_pGrid->get_attachment_pipe<EdgeBase>());

	if((shElements & SE_FACE) && (!elements_are_supported(SE_FACE)))
		section_container<Face>().get_container().
				set_pipe(&m_pGrid->get_attachment_pipe<Face>());

	if((shElements & SE_VOLUME) && (!elements_are_supported(SE_VOLUME)))
		section_container<Volume>().get_container().
				set_pipe(&m_pGrid->get_attachment_pipe<Volume>());

	ISelector::enable_element_support(shElements);
}

void Selector::disable_element_support(uint shElements)
{
	//	release the attachments in the current grid
	if((shElements & SE_VERTEX) && elements_are_supported(SE_VERTEX))
		section_container<VertexBase>().get_container().set_pipe(NULL);

	if((shElements & SE_EDGE) && elements_are_supported(SE_EDGE))
		section_container<EdgeBase>().get_container().set_pipe(NULL);

	if((shElements & SE_FACE) && elements_are_supported(SE_FACE))
		section_container<Face>().get_container().set_pipe(NULL);

	if((shElements & SE_VOLUME) && elements_are_supported(SE_VOLUME))
		section_container<Volume>().get_container().set_pipe(NULL);

	ISelector::disable_element_support(shElements);
}

void Selector::clear_lists()
{
	section_container<VertexBase>().clear();
	section_container<EdgeBase>().clear();
	section_container<Face>().clear();
	section_container<Volume>().clear();
}

void Selector::clear()
{
	clear<VertexBase>();
	clear<EdgeBase>();
	clear<Face>();
	clear<Volume>();
}

void Selector::add_to_list(VertexBase* elem)
{
	section_container<VertexBase>().insert(elem,
								elem->container_section());
}

void Selector::add_to_list(EdgeBase* elem)
{
	section_container<EdgeBase>().insert(elem,
								elem->container_section());
}

void Selector::add_to_list(Face* elem)
{
	section_container<Face>().insert(elem,
								elem->container_section());
}

void Selector::add_to_list(Volume* elem)
{
	section_container<Volume>().insert(elem,
								elem->container_section());
}	

void Selector::erase_from_list(VertexBase* elem)
{
	section_container<VertexBase>().erase(get_iterator(elem),
						elem->container_section());
}
void Selector::erase_from_list(EdgeBase* elem)
{
	section_container<EdgeBase>().erase(get_iterator(elem),
						elem->container_section());
}
void Selector::erase_from_list(Face* elem)
{
	section_container<Face>().erase(get_iterator(elem),
						elem->container_section());
}

void Selector::erase_from_list(Volume* elem)
{
	section_container<Volume>().erase(get_iterator(elem),
						elem->container_section());
}

//	geometric-object-collection
GridObjectCollection Selector::get_grid_objects() const
{
//TODO: ugly casts! GenericElementSelector should store its selected elements
//		in a GridObjectSectionContainer!
	return GridObjectCollection(const_cast<VertexSectionContainer*>(&m_vertices),
									 const_cast<EdgeSectionContainer*>(&m_edges),
									 const_cast<FaceSectionContainer*>(&m_faces),
									 const_cast<VolumeSectionContainer*>(&m_volumes));
}

/*
void Selector::unregistered_from_grid(Grid* grid)
{
	clear_lists();
}
*/
void Selector::grid_to_be_destroyed(Grid* grid)
{
	ISelector::grid_to_be_destroyed(grid);
	clear_lists();
}

}//	end of namespace
