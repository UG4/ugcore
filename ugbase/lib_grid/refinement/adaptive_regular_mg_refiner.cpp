/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "lib_grid/algorithms/selection_util.h"
#include "adaptive_regular_mg_refiner.h"

using namespace std;

namespace ug{

AdaptiveRegularRefiner_MultiGrid::
AdaptiveRegularRefiner_MultiGrid(SPRefinementProjector projector) :
	HangingNodeRefiner_MultiGrid(projector)
{
}

AdaptiveRegularRefiner_MultiGrid::
AdaptiveRegularRefiner_MultiGrid(MultiGrid& mg,
								 SPRefinementProjector projector) :
	HangingNodeRefiner_MultiGrid(projector)
{
	set_grid(&mg);
}

AdaptiveRegularRefiner_MultiGrid::
~AdaptiveRegularRefiner_MultiGrid()
{
}

void AdaptiveRegularRefiner_MultiGrid::
assign_grid(MultiGrid& mg)
{
	set_grid(&mg);
}

void AdaptiveRegularRefiner_MultiGrid::
set_grid(MultiGrid* mg)
{
	HangingNodeRefiner_MultiGrid::set_grid(mg);

	m_closureElems.enable_autoselection(false);
	m_closureElems.enable_selection_inheritance(false);
	m_closureElems.enable_strict_inheritance(false);

	m_closureElems.assign_grid(mg);
}


void AdaptiveRegularRefiner_MultiGrid::
remove_closure_elements()
{
//	remove all closure elements from the current grid

	if(!multi_grid()){
		UG_THROW("AdaptiveRegularRefiner_MultiGrid has to be associated with a grid!");
	}

	EraseSelectedObjects(m_closureElems);
}


void AdaptiveRegularRefiner_MultiGrid::
create_closure_elements()
{
	if(!multi_grid()){
		UG_THROW("AdaptiveRegularRefiner_MultiGrid has to be associated with a grid!");
	}

	MultiGrid& mg = *multi_grid();

	if(mg.num<Volume>() > 0)
		create_closure_elements_3d();
	else if(mg.num<Face>() > 0)
		create_closure_elements_2d();
}


void AdaptiveRegularRefiner_MultiGrid::
create_closure_elements_2d()
{
//todo:	This method currently only works for one type of elements. No manifolds
//		are currently supported.

//	for each surface element (faces in 2d, volumes in 3d) adjacent to a constraining
//	element, we generate closure elements in the level above.

	if(!multi_grid()){
		UG_THROW("AdaptiveRegularRefiner_MultiGrid has to be associated with a grid!");
	}

//	remove all closure elements. This is currently required, since the selector
//	m_closureElems is currently also used to store non-closer elements, which are
//	to be refined.
	remove_closure_elements();

	MultiGrid& mg = *multi_grid();

//	iterate over all constraining edges and collect associated surface faces
//	Those then have to be refined to generate a closure.

	Grid::face_traits::secure_container	assElems;
	Grid::edge_traits::secure_container assEdges;
	Face::ConstVertexArray vrts;
	std::vector<Vertex*> newVrtVrts;
	std::vector<Vertex*> newEdgeVrts;
	std::vector<Face*>	newFaces;
	EdgeDescriptor ed;

//	we'll select all new elements on the fly
	m_closureElems.enable_autoselection(true);

	for(Grid::traits<ConstrainingEdge>::iterator i_edge = mg.begin<ConstrainingEdge>();
		i_edge != mg.end<ConstrainingEdge>(); ++i_edge)
	{
	//	check all associated elements of i_edge, whether one is a surface element.
	//	If so, create the closure.
		mg.associated_elements(assElems, *i_edge);

		for(size_t i_elem = 0; i_elem < assElems.size(); ++i_elem){
			Face* elem = assElems[i_elem];
			if(mg.has_children(elem))
				continue;

			newVrtVrts.clear();
			newEdgeVrts.clear();

		//	copy associated vertices and edges to the next level
			vrts = elem->vertices();
			size_t numVrts = elem->num_vertices();
			for(size_t i_vrt = 0; i_vrt < numVrts; ++i_vrt){
				Vertex* vrt = vrts[i_vrt];
				if(!mg.has_children(vrt)){
					newVrtVrts.push_back(*mg.create<RegularVertex>(vrt));
					if(m_projector.valid())
						m_projector->new_vertex(newVrtVrts.back(), vrt);
				}
				else
					newVrtVrts.push_back(mg.get_child_vertex(vrt));
			}

			mg.associated_elements_sorted(assEdges, elem);
			for(size_t i = 0; i < assEdges.size(); ++i){
				Edge* e = assEdges[i];
				if(!mg.has_children(e)){
					ed.set_vertices(mg.get_child_vertex(e->vertex(0)),
									mg.get_child_vertex(e->vertex(1)));
					mg.create<RegularEdge>(ed, e);
					newEdgeVrts.push_back(nullptr);
				}
				else
					newEdgeVrts.push_back(mg.get_child_vertex(e));
			}

		//	refine the element
			Vertex* newFaceVrt;
			if(elem->refine(newFaces, &newFaceVrt, &newEdgeVrts.front(),
							nullptr, &newVrtVrts.front()))
			{
				if(newFaceVrt){
					mg.register_element(newFaceVrt, elem);
					if(m_projector.valid())
						m_projector->new_vertex(newFaceVrt, elem);
				}
				for(size_t i = 0; i < newFaces.size(); ++i)
					mg.register_element(newFaces[i], elem);
			}
		}
	}

//	stop selecting new elements
	m_closureElems.enable_autoselection(false);
}

void AdaptiveRegularRefiner_MultiGrid::
create_closure_elements_3d()
{
//todo:	This method currently only works for one type of elements. No manifolds
//		are currently supported.

//	for each surface element (faces in 2d, volumes in 3d) adjacent to a constraining
//	element, we generate closure elements in the level above.

	if(!multi_grid()){
		UG_THROW("AdaptiveRegularRefiner_MultiGrid has to be associated with a grid!");
	}

//	remove all closure elements. This is currently required, since the selector
//	m_closureElems is currently also used to store non-closer elements, which are
//	to be refined.
	remove_closure_elements();

	MultiGrid& mg = *multi_grid();

//	iterate over all constraining edges and collect associated surface faces / volumes.
//	Those then have to be refined to generate a closure.

	Grid::volume_traits::secure_container	assElems;
	Grid::edge_traits::secure_container assEdges;
	Grid::face_traits::secure_container assFaces;
	Volume::ConstVertexArray vrts;
	std::vector<Vertex*> newVolVrtVrts;
	std::vector<Vertex*> newVolEdgeVrts;
	std::vector<Vertex*> newVolFaceVrts;
	std::vector<Volume*> newVols;
	EdgeDescriptor ed;
	FaceDescriptor fd;
	
//	when refining the associated faces, we need some structs, too
	Grid::edge_traits::secure_container assFaceEdges;
	std::vector<Vertex*> newFaceVrtVrts;
	std::vector<Vertex*> newFaceEdgeVrts;
	std::vector<Face*> newFaces;

//	we'll select all new elements on the fly
	m_closureElems.enable_autoselection(true);

	for(Grid::traits<ConstrainingEdge>::iterator i_edge = mg.begin<ConstrainingEdge>();
		i_edge != mg.end<ConstrainingEdge>(); ++i_edge)
	{
	//	check all associated elements of i_edge, whether one is a surface element.
	//	If so, select it into m_closureElems. This is only temporary, since it
	//	isn't a real closure element.
		mg.associated_elements(assElems, *i_edge);

		for(size_t i_elem = 0; i_elem < assElems.size(); ++i_elem){
			Volume* elem = assElems[i_elem];
			if(mg.has_children(elem))
				continue;

			newVolVrtVrts.clear();
			newVolEdgeVrts.clear();
			newVolFaceVrts.clear();

		//	copy associated vertices and edges to the next level
			vrts = elem->vertices();
			size_t numVrts = elem->num_vertices();
			for(size_t i_vrt = 0; i_vrt < numVrts; ++i_vrt){
				Vertex* vrt = vrts[i_vrt];
				if(!mg.has_children(vrt)){
					newVolVrtVrts.push_back(*mg.create<RegularVertex>(vrt));
					if(m_projector.valid())
						m_projector->new_vertex(newVolVrtVrts.back(), vrt);
				}
				else
					newVolVrtVrts.push_back(mg.get_child_vertex(vrt));
			}

			mg.associated_elements_sorted(assEdges, elem);
			for(size_t i = 0; i < assEdges.size(); ++i){
				Edge* e = assEdges[i];
				if(!mg.has_children(e)){
					ed.set_vertices(mg.get_child_vertex(e->vertex(0)),
									mg.get_child_vertex(e->vertex(1)));
					mg.create<RegularEdge>(ed, e);
					newVolEdgeVrts.push_back(nullptr);
				}
				else
					newVolEdgeVrts.push_back(mg.get_child_vertex(e));
			}

		//	we have to either refine or copy associated faces
			mg.associated_elements_sorted(assFaces, elem);
			for(size_t i_face = 0; i_face < assFaces.size(); ++i_face){
				
				Face* f = assFaces[i_face];
				if(mg.has_children(f)){
					newVolFaceVrts.push_back(mg.get_child_vertex(f));
					continue;
				}
				
			//	collect child vertices of associated edges
				mg.associated_elements_sorted(assFaceEdges, f);
				newFaceEdgeVrts.clear();
				
				bool faceRefinement = false;
				for(size_t i = 0; i < assFaceEdges.size(); ++i){
					Vertex* child = mg.get_child_vertex(assFaceEdges[i]);
					newFaceEdgeVrts.push_back(child);
					faceRefinement |= (child != nullptr);
				}				
				
				if(faceRefinement){
					newFaceVrtVrts.clear();
					for(size_t i = 0; i < f->num_vertices(); ++i)
						newFaceVrtVrts.push_back(mg.get_child_vertex(f->vertex(i)));
					
					Vertex* newFaceVrt = nullptr;
					if(f->refine(newFaces, &newFaceVrt, &newFaceEdgeVrts.front(),
								 nullptr, &newFaceVrtVrts.front()))
					{
						if(newFaceVrt){
							mg.register_element(newFaceVrt, f);
							if(m_projector.valid())
								m_projector->new_vertex(newFaceVrt, f);
						}
						for(size_t i = 0; i < newFaces.size(); ++i)
							mg.register_element(newFaces[i], f);
					}
					newVolFaceVrts.push_back(newFaceVrt);
				}	
				else{
					Face::ConstVertexArray fvrts = f->vertices();
					if(f->num_vertices() == 3)
						mg.create<Triangle>(TriangleDescriptor(
												mg.get_child_vertex(fvrts[0]),
												mg.get_child_vertex(fvrts[1]),
												mg.get_child_vertex(fvrts[2])),
											f);
					else if(f->num_vertices() == 4)
						mg.create<Quadrilateral>(QuadrilateralDescriptor(
													mg.get_child_vertex(fvrts[0]),
													mg.get_child_vertex(fvrts[1]),
													mg.get_child_vertex(fvrts[2]),
													mg.get_child_vertex(fvrts[3])),
												f);
					newVolFaceVrts.push_back(nullptr);
				}
			}

		//	refine the element
		//	if we're performing tetrahedral refinement, we have to collect
		//	the corner coordinates, so that the refinement algorithm may choose
		//	the best interior diagonal.
			vector3 corners[4];
			vector3* pCorners = nullptr;
			if((elem->num_vertices() == 4) && m_projector.valid()){
				for(size_t i = 0; i < 4; ++i){
					corners[i] = m_projector->geometry()->pos(vrts[i]);
				}
				pCorners = corners;
			}

			Vertex* newVolVrt;
			if(elem->refine(newVols, &newVolVrt, &newVolEdgeVrts.front(),
							&newVolFaceVrts.front(), nullptr, RegularVertex(),
							&newVolVrtVrts.front(), pCorners))
			{
				if(newVolVrt){
					mg.register_element(newVolVrt, elem);
					if(m_projector.valid())
						m_projector->new_vertex(newVolVrt, elem);
				}

				for(size_t i = 0; i < newVols.size(); ++i)
					mg.register_element(newVols[i], elem);
			}
		}
	}

//	stop selecting new elements
	m_closureElems.enable_autoselection(false);
}

template <class TElem>
void AdaptiveRegularRefiner_MultiGrid::
get_parents_of_marked_closure_elements(std::vector<GridObject*>& parents,
									   Selector::status_t mark)
{
	UG_ASSERT(multi_grid(), "A multi grid has to be assigned to the refiner.");
	MultiGrid& mg = *multi_grid();

	using TIter = typename selector_t::traits<TElem>::iterator;
	for(TIter iter = m_selMarkedElements.begin<TElem>();
		iter != m_selMarkedElements.end<TElem>(); ++iter)
	{
		TElem* e = *iter;
		if(!m_closureElems.is_selected(e))
			continue;

		if(get_mark(e) & mark){
			GridObject* parent = mg.get_parent(e);
			if(parent && !m_closureElems.is_selected(parent))
				parents.push_back(parent);
		}
	}
}

void AdaptiveRegularRefiner_MultiGrid::
perform_refinement()
{
//	todo: copy refinement marks from closure elements to their parents
	vector<GridObject*> parents;
	Selector::status_t refMark = RM_REFINE | RM_ANISOTROPIC;
	get_parents_of_marked_closure_elements<Vertex>(parents, refMark);
	get_parents_of_marked_closure_elements<Edge>(parents, refMark);
	get_parents_of_marked_closure_elements<Face>(parents, refMark);
	get_parents_of_marked_closure_elements<Volume>(parents, refMark);

	remove_closure_elements();

//	mark parents of formerly marked closure elements for refinement
	mark(parents.begin(), parents.end(), RM_REFINE);

	HangingNodeRefiner_MultiGrid::perform_refinement();
	create_closure_elements();
}

bool AdaptiveRegularRefiner_MultiGrid::
perform_coarsening()
{
//	todo: copy coarsen marks from closure elements to their parents
	remove_closure_elements();
	bool bSuccess = HangingNodeRefiner_MultiGrid::perform_coarsening();
	create_closure_elements();
	return bSuccess;
}

}// end of namespace
