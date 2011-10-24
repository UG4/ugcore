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
			get_section_container<VertexBase>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_EDGE))
			get_section_container<EdgeBase>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_FACE))
			get_section_container<Face>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_VOLUME))
			get_section_container<Volume>().get_container().set_pipe(NULL);
	}
}

void Selector::assign_grid(Grid& grid)
{
	assign_grid(&grid);
}

void Selector::assign_grid(Grid* grid)
{
	if(m_pGrid){
	//	release the attachments in the current grid
		if(elements_are_supported(SE_VERTEX))
			get_section_container<VertexBase>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_EDGE))
			get_section_container<EdgeBase>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_FACE))
			get_section_container<Face>().get_container().set_pipe(NULL);

		if(elements_are_supported(SE_VOLUME))
			get_section_container<Volume>().get_container().set_pipe(NULL);
	}

	BaseClass::set_grid(grid);

	if(m_pGrid){
	//	initialize attachment lists
		if(elements_are_supported(SE_VERTEX))
			get_section_container<VertexBase>().get_container().
					set_pipe(&m_pGrid->get_attachment_pipe<VertexBase>());

		if(elements_are_supported(SE_EDGE))
			get_section_container<EdgeBase>().get_container().
					set_pipe(&m_pGrid->get_attachment_pipe<EdgeBase>());

		if(elements_are_supported(SE_FACE))
			get_section_container<Face>().get_container().
					set_pipe(&m_pGrid->get_attachment_pipe<Face>());

		if(elements_are_supported(SE_VOLUME))
			get_section_container<Volume>().get_container().
					set_pipe(&m_pGrid->get_attachment_pipe<Volume>());
	}
}

void Selector::clear_lists()
{
	get_section_container<VertexBase>().clear();
	get_section_container<EdgeBase>().clear();
	get_section_container<Face>().clear();
	get_section_container<Volume>().clear();
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
	get_section_container<VertexBase>().insert(elem,
								elem->shared_pipe_section());
}

void Selector::add_to_list(EdgeBase* elem)
{
	get_section_container<EdgeBase>().insert(elem,
								elem->shared_pipe_section());
}

void Selector::add_to_list(Face* elem)
{
	get_section_container<Face>().insert(elem,
								elem->shared_pipe_section());
}

void Selector::add_to_list(Volume* elem)
{
	get_section_container<Volume>().insert(elem,
								elem->shared_pipe_section());
}	

void Selector::erase_from_list(VertexBase* elem)
{
	get_section_container<VertexBase>().erase(get_iterator(elem),
						elem->shared_pipe_section());
}
void Selector::erase_from_list(EdgeBase* elem)
{
	get_section_container<EdgeBase>().erase(get_iterator(elem),
						elem->shared_pipe_section());
}
void Selector::erase_from_list(Face* elem)
{
	get_section_container<Face>().erase(get_iterator(elem),
						elem->shared_pipe_section());
}

void Selector::erase_from_list(Volume* elem)
{
	get_section_container<Volume>().erase(get_iterator(elem),
						elem->shared_pipe_section());
}

//	geometric-object-collection
GeometricObjectCollection Selector::get_geometric_objects()
{
//TODO: ugly casts! GenericElementSelector should store its selected elements
//		in a GeometricObjectSectionContainer!
	return GeometricObjectCollection(&m_elements[VERTEX],
									&m_elements[EDGE],
									&m_elements[FACE],
									&m_elements[VOLUME]);
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
