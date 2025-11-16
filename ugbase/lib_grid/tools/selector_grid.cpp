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
			section_container<Vertex>().get_container().set_pipe(nullptr);

		if(elements_are_supported(SE_EDGE))
			section_container<Edge>().get_container().set_pipe(nullptr);

		if(elements_are_supported(SE_FACE))
			section_container<Face>().get_container().set_pipe(nullptr);

		if(elements_are_supported(SE_VOLUME))
			section_container<Volume>().get_container().set_pipe(nullptr);
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
		section_container<Vertex>().get_container().
				set_pipe(&m_pGrid->get_attachment_pipe<Vertex>());

	if((shElements & SE_EDGE) && (!elements_are_supported(SE_EDGE)))
		section_container<Edge>().get_container().
				set_pipe(&m_pGrid->get_attachment_pipe<Edge>());

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
		section_container<Vertex>().get_container().set_pipe(nullptr);

	if((shElements & SE_EDGE) && elements_are_supported(SE_EDGE))
		section_container<Edge>().get_container().set_pipe(nullptr);

	if((shElements & SE_FACE) && elements_are_supported(SE_FACE))
		section_container<Face>().get_container().set_pipe(nullptr);

	if((shElements & SE_VOLUME) && elements_are_supported(SE_VOLUME))
		section_container<Volume>().get_container().set_pipe(nullptr);

	ISelector::disable_element_support(shElements);
}

void Selector::clear_lists()
{
	section_container<Vertex>().clear();
	section_container<Edge>().clear();
	section_container<Face>().clear();
	section_container<Volume>().clear();
}

void Selector::clear()
{
	clear<Vertex>();
	clear<Edge>();
	clear<Face>();
	clear<Volume>();
}

void Selector::add_to_list(Vertex* elem)
{
	section_container<Vertex>().insert(elem,
								elem->container_section());
}

void Selector::add_to_list(Edge* elem)
{
	section_container<Edge>().insert(elem,
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

void Selector::erase_from_list(Vertex* elem)
{
	section_container<Vertex>().erase(get_iterator(elem),
						elem->container_section());
}
void Selector::erase_from_list(Edge* elem)
{
	section_container<Edge>().erase(get_iterator(elem),
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
