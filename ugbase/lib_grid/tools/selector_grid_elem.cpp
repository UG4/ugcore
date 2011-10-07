// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m02 d16

#include <cassert>
#include "selector_grid_elem.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class BaseElem>
TElemSelector<BaseElem>::TElemSelector() :
	ISelector(SE_NONE)
{
//	enable element support
//TODO:	unify constants in ISubsetHandler and ISelector and add them to the traits.
	int objID = geometry_traits<BaseElem>::BASE_OBJECT_TYPE_ID;
	switch(objID){
		case VERTEX:
			set_supported_elements(SE_VERTEX);
			break;
		case EDGE:
			set_supported_elements(SE_EDGE);
			break;
		case FACE:
			set_supported_elements(SE_FACE);
			break;
		case VOLUME:
			set_supported_elements(SE_VOLUME);
			break;
	}
}

template <class BaseElem>
TElemSelector<BaseElem>::TElemSelector(Grid& grid) :
	ISelector(SE_NONE)
{
//	enable element support
//TODO:	unify constants in ISubsetHandler and ISelector and add them to the traits.
	int objID = geometry_traits<BaseElem>::BASE_OBJECT_TYPE_ID;
	switch(objID){
		case VERTEX:
			set_supported_elements(SE_VERTEX);
			break;
		case EDGE:
			set_supported_elements(SE_EDGE);
			break;
		case FACE:
			set_supported_elements(SE_FACE);
			break;
		case VOLUME:
			set_supported_elements(SE_VOLUME);
			break;
	}

	assign_grid(grid);
}

template <class BaseElem>
void TElemSelector<BaseElem>::assign_grid(Grid& grid)
{
	BaseClass::set_grid(&grid);

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

template <class BaseElem>
void TElemSelector<BaseElem>::clear_lists()
{
	if(elements_are_supported(SE_VERTEX))
		get_section_container<VertexBase>().clear();
	if(elements_are_supported(SE_EDGE))
		get_section_container<EdgeBase>().clear();
	if(elements_are_supported(SE_FACE))
		get_section_container<Face>().clear();
	if(elements_are_supported(SE_VOLUME))
		get_section_container<Volume>().clear();
}

template <class BaseElem>
void TElemSelector<BaseElem>::clear()
{
	if(elements_are_supported(SE_VERTEX))
		clear<VertexBase>();
	if(elements_are_supported(SE_EDGE))
		clear<EdgeBase>();
	if(elements_are_supported(SE_FACE))
		clear<Face>();
	if(elements_are_supported(SE_VOLUME))
		clear<Volume>();
}

template <class BaseElem>
void
TElemSelector<BaseElem>::add_to_list(VertexBase* elem)
{
	get_section_container<VertexBase>().insert(elem,
								elem->shared_pipe_section());
}

template <class BaseElem>
void
TElemSelector<BaseElem>::add_to_list(EdgeBase* elem)
{
	get_section_container<EdgeBase>().insert(elem,
								elem->shared_pipe_section());
}

template <class BaseElem>
void
TElemSelector<BaseElem>::add_to_list(Face* elem)
{
	get_section_container<Face>().insert(elem,
								elem->shared_pipe_section());
}

template <class BaseElem>
void
TElemSelector<BaseElem>::add_to_list(Volume* elem)
{
	get_section_container<Volume>().insert(elem,
								elem->shared_pipe_section());
}	

template <class BaseElem>
void TElemSelector<BaseElem>::erase_from_list(VertexBase* elem)
{
	get_section_container<VertexBase>().erase(get_iterator(elem),
						elem->shared_pipe_section());
}

template <class BaseElem>
void TElemSelector<BaseElem>::erase_from_list(EdgeBase* elem)
{
	get_section_container<EdgeBase>().erase(get_iterator(elem),
						elem->shared_pipe_section());
}

template <class BaseElem>
void TElemSelector<BaseElem>::erase_from_list(Face* elem)
{
	get_section_container<Face>().erase(get_iterator(elem),
						elem->shared_pipe_section());
}

template <class BaseElem>
void TElemSelector<BaseElem>::erase_from_list(Volume* elem)
{
	get_section_container<Volume>().erase(get_iterator(elem),
						elem->shared_pipe_section());
}

//	geometric-object-collection
template <class BaseElem>
GeometricObjectCollection 
TElemSelector<BaseElem>::get_geometric_objects()
{
//TODO: ugly casts! GenericElementSelector should store its selected elements
//		in a GeometricObjectSectionContainer!
	return GeometricObjectCollection(NULL,
									&m_elements,
									NULL,
									NULL);
}

template <class BaseElem>
void TElemSelector<BaseElem>::unregistered_from_grid(Grid* grid)
{
	clear_lists();
}


//	explicit instatiation for the four main geometric objects
template class TElemSelector<VertexBase>;
template class TElemSelector<EdgeBase>;
template class TElemSelector<Face>;
template class TElemSelector<Volume>;

}//	end of namespace
