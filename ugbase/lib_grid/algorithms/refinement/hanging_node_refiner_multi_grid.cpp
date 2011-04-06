// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 11.01.2011 (m,d,y)
 
#include "common/assert.h"
#include "hanging_node_refiner_multi_grid.h"

using namespace std;

namespace ug{

HangingNodeRefiner_MultiGrid::
HangingNodeRefiner_MultiGrid(IRefinementCallback* refCallback) :
	HangingNodeRefinerBase(refCallback),
	m_pMG(NULL)
{

}

HangingNodeRefiner_MultiGrid::
HangingNodeRefiner_MultiGrid(MultiGrid& mg,
						IRefinementCallback* refCallback) :
	HangingNodeRefinerBase(refCallback),
	m_pMG(NULL)
{
	set_grid(&mg);
}

HangingNodeRefiner_MultiGrid::
~HangingNodeRefiner_MultiGrid()
{
	set_grid(NULL);
}

void HangingNodeRefiner_MultiGrid::
assign_grid(MultiGrid& mg)
{
	set_grid(&mg);
}

void HangingNodeRefiner_MultiGrid::
set_grid(MultiGrid* mg)
{
//	call base class implementation
	HangingNodeRefinerBase::set_grid(mg);

	m_pMG = mg;
}

void HangingNodeRefiner_MultiGrid::
pre_refine()
{
	if(!m_pMG)
		throw(UGError("HangingNodeRefiner_MultiGrid::refine: No grid assigned."));

	MultiGrid& mg = *m_pMG;

//	create a child vertex for each marked vertex
	for(VertexBaseIterator iter = m_selMarkedElements.begin<VertexBase>();
		iter != m_selMarkedElements.end<VertexBase>(); ++iter)
	{
		UG_ASSERT(!mg.has_children(*iter), "Only vertices without children should be selected.");
		if(refinement_is_allowed(*iter)){
			VertexBase* vrt = *mg.create<Vertex>(*iter);
			if(m_refCallback)
				m_refCallback->new_vertex(vrt, *iter);
		}
	}
}

void HangingNodeRefiner_MultiGrid::
collect_objects_for_refine()
{
	if(!m_pMG)
		throw(UGError("HangingNodeRefiner_MultiGrid::collect_objects_for_refine: No grid assigned."));

	MultiGrid& mg = *m_pMG;

//	make sure that no element which has children is selected
//	(except constraining edges and faces)
	for(EdgeBaseIterator iter = m_selMarkedElements.begin<EdgeBase>();
		iter != m_selMarkedElements.end<EdgeBase>();)
	{
		EdgeBase* e = *iter;
		++iter;

		if(mg.has_children(e)){
			if(!ConstrainingEdge::type_match(e))
				m_selMarkedElements.deselect(e);
		}
	}

	for(FaceIterator iter = m_selMarkedElements.begin<Face>();
		iter != m_selMarkedElements.end<Face>();)
	{
		Face* f = *iter;
		++iter;

		if(mg.has_children(f)){
			if(!ConstrainingFace::type_match(f))
				m_selMarkedElements.deselect(f);
		}
	}

	for(VolumeIterator iter = m_selMarkedElements.begin<Volume>();
		iter != m_selMarkedElements.end<Volume>();)
	{
		Volume* v = *iter;
		++iter;

		if(mg.has_children(v)){
			m_selMarkedElements.deselect(v);
		}
	}

//	first call base implementation
	HangingNodeRefinerBase::collect_objects_for_refine();

//	make sure that no element which has children is selected
//	(except constraining edges and faces)
	for(VertexBaseIterator iter = m_selMarkedElements.begin<VertexBase>();
		iter != m_selMarkedElements.end<VertexBase>();)
	{
		VertexBase* v = *iter;
		++iter;

		if(mg.has_children(v)){
			m_selMarkedElements.deselect(v);
		}
	}

	for(EdgeBaseIterator iter = m_selMarkedElements.begin<EdgeBase>();
		iter != m_selMarkedElements.end<EdgeBase>();)
	{
		EdgeBase* e = *iter;
		++iter;

		if(mg.has_children(e)){
			if(!ConstrainingEdge::type_match(e))
				m_selMarkedElements.deselect(e);
		}
	}

	for(FaceIterator iter = m_selMarkedElements.begin<Face>();
		iter != m_selMarkedElements.end<Face>();)
	{
		Face* f = *iter;
		++iter;

		if(mg.has_children(f)){
			if(!ConstrainingFace::type_match(f))
				m_selMarkedElements.deselect(f);
		}
	}

	for(VolumeIterator iter = m_selMarkedElements.begin<Volume>();
		iter != m_selMarkedElements.end<Volume>();)
	{
		Volume* v = *iter;
		++iter;

		if(mg.has_children(v)){
			m_selMarkedElements.deselect(v);
		}
	}

//	finally we have to select all associated vertices of marked objects,
//	since we have to create new vertices in the next levels of the hierarchies.
//	only vertices which do not already have child vertices are selected.
	for(EdgeBaseIterator iter = m_selMarkedElements.begin<EdgeBase>();
		iter != m_selMarkedElements.end<EdgeBase>(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!mg.has_children((*iter)->vertex(i)))
				mark((*iter)->vertex(i));
		}
	}

	for(FaceIterator iter = m_selMarkedElements.begin<Face>();
		iter != m_selMarkedElements.end<Face>(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!mg.has_children((*iter)->vertex(i)))
				mark((*iter)->vertex(i));
		}
	}

	for(VolumeIterator iter = m_selMarkedElements.begin<Volume>();
		iter != m_selMarkedElements.end<Volume>(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!mg.has_children((*iter)->vertex(i)))
				mark((*iter)->vertex(i));
		}
	}
}

void HangingNodeRefiner_MultiGrid::
refine_edge_with_normal_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, e);

//	call original implementation
	HangingNodeRefinerBase::refine_edge_with_normal_vertex(e, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_edge_with_hanging_vertex(EdgeBase* e, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, e);

//	call original implementation
	HangingNodeRefinerBase::refine_edge_with_hanging_vertex(e, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_face_with_normal_vertex(Face* f, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, f);

//	call original implementation
	HangingNodeRefinerBase::refine_face_with_normal_vertex(f, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_face_with_hanging_vertex(Face* f, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, f);

//	call original implementation
	HangingNodeRefinerBase::refine_face_with_hanging_vertex(f, &childCorners.front());
}

void HangingNodeRefiner_MultiGrid::
refine_volume_with_normal_vertex(Volume* v, VertexBase** newCornerVrts)
{
//	collect child corners
	std::vector<VertexBase*> childCorners;
	collect_child_corners(childCorners, v);

//	call original implementation
	HangingNodeRefinerBase::refine_volume_with_normal_vertex(v, &childCorners.front());
}

VertexBase* HangingNodeRefiner_MultiGrid::
get_center_vertex(EdgeBase* e)
{
	return m_pMG->get_child_vertex(e);
}


void HangingNodeRefiner_MultiGrid::
set_center_vertex(EdgeBase* e, VertexBase* v)
{
//	the child vertex should already have been associated since the
//	multigrid is an observer itself.
	UG_ASSERT(m_pMG->get_child_vertex(e) == v, "child/vertex mismatch.");
}


VertexBase* HangingNodeRefiner_MultiGrid::
get_center_vertex(Face* f)
{
	return m_pMG->get_child_vertex(f);
}


void HangingNodeRefiner_MultiGrid::
set_center_vertex(Face* f, VertexBase* v)
{
//	the child vertex should already have been associated since the
//	multigrid is an observer itself.
	UG_ASSERT(m_pMG->get_child_vertex(f) == v, "child/vertex mismatch.");
}

void HangingNodeRefiner_MultiGrid::
refine_constraining_edge(ConstrainingEdge* cge)
{
//	call base implementation to perform refinement
	HangingNodeRefinerBase::refine_constraining_edge(cge);

//	the constrained edge is now a normal edge
	m_pMG->create_and_replace<Edge>(cge);
}

}// end of namespace
