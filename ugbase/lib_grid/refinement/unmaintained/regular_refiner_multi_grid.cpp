øunused
/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#include "regular_refiner_multi_grid.h"
#include "common/util/vec_for_each.h"
#include "lib_grid/iterators/lg_for_each.h"


namespace ug{

RegularRefiner_MultiGrid::
RegularRefiner_MultiGrid () :
	m_pMG(nullptr)
{}

RegularRefiner_MultiGrid::
RegularRefiner_MultiGrid (MultiGrid* pmg) :
	m_marks(*pmg),
	m_pMG(pmg)
{}


void RegularRefiner_MultiGrid::
set_grid (MultiGrid* pmg) {
	m_pMG = pmg;
	m_marks.assign_grid(pmg);
}


MultiGrid* RegularRefiner_MultiGrid::
multi_grid () const {
	return m_pMG;
}


Grid* RegularRefiner_MultiGrid::
get_associated_grid () {
	return m_pMG;
}


Grid* RegularRefiner_MultiGrid::
grid () {
	return m_pMG;
}


bool RegularRefiner_MultiGrid::
adaptivity_supported () const {
	return true;
}


bool RegularRefiner_MultiGrid::
coarsening_supported () const {
	return false;
}



bool RegularRefiner_MultiGrid::
mark (Vertex* v, RefinementMark refMark) {
	m_marks.select(v, refMark);
	return true;
}


bool RegularRefiner_MultiGrid::
mark (Edge* e, RefinementMark refMark) {
	m_marks.select(e, refMark);
	return true;
}


bool RegularRefiner_MultiGrid::
mark (Face* f, RefinementMark refMark) {
	m_marks.select(f, refMark);
	return true;
}


bool RegularRefiner_MultiGrid::
mark (Volume* v, RefinementMark refMark) {
	m_marks.select(v, refMark);
	return true;
}



void RegularRefiner_MultiGrid::
mark_neighborhood (size_t numIterations, RefinementMark refMark, bool sideNbrsOnly) {
	UG_THROW("NOT YET IMPLEMENTED!");
}


RefinementMark RegularRefiner_MultiGrid::
get_mark (Vertex* v) {
	return static_cast<RefinementMark>( m_marks.get_selection_status(v) );
}


RefinementMark RegularRefiner_MultiGrid::
get_mark (Edge* e) {
	return static_cast<RefinementMark>( m_marks.get_selection_status(e) );
}


RefinementMark RegularRefiner_MultiGrid::
get_mark (Face* f) {
	return static_cast<RefinementMark>( m_marks.get_selection_status(f) );
}


RefinementMark RegularRefiner_MultiGrid::
get_mark (Volume* v) {
	return static_cast<RefinementMark>( m_marks.get_selection_status(v) );
}


bool RegularRefiner_MultiGrid::
save_marks_to_file (const char* filename) {
	UG_THROW("NOT YET IMPLEMENTED!");
}


void RegularRefiner_MultiGrid::
perform_refinement () {

	collect_objects_for_refine();

//	emit refinement-begins message

//	perform refinement

//	emit refinement-ends message

//	clean-up
}

template <typename TElem>
void RegularRefiner_MultiGrid::
adjust_side_states (
	size_t lvl,
	uint considerElemMarks,
	uint ignoreSideMarks,
	RefinementMark newSideMark,
	bool closure)
{
	using side_t = typename TElem::side;

	MultiGrid& mg = *m_pMG;
	MGSelector& sel = m_marks;
	typename Grid::traits<side_t>::secure_container sides;

	lg_for_each_in_lvl_template(TElem, elem, sel, lvl){
		uint m = get_mark(elem);
		if(m & considerElemMarks){
			mg.associated_elements(sides, elem);
			for_each_in_vec(side_t* s, sides){
			//	unless the side is already marked for regular refinement,
			//	we'll copy the state to the side
				if((get_mark(s) & ignoreSideMarks) == 0){
				//	don't use mark here, since we also want to mark shadow-rims
					if(!refinement_is_allowed(s))
						sel.select(s, RM_DUMMY);
					else
						sel.select(s, newSideMark);
					if(closure)
						m_aaClosure[s] = true;
				}
			}end_for;
		}
	}lg_end_for;
}

template <typename TElem>
void RegularRefiner_MultiGrid::
copy_state_to_sides (
	size_t lvl,
	uint considerElemMarks,
	bool closure)
{
	using side_t = typename TElem::side;

	MultiGrid& mg = *m_pMG;
	MGSelector& sel = m_marks;
	typename Grid::traits<side_t>::secure_container sides;

	lg_for_each_in_lvl_template(TElem, elem, sel, lvl){
		RefinementMark m = get_mark(elem);
		if(m & considerElemMarks){
			mg.associated_elements(sides, elem);
			for_each_in_vec(side_t* s, sides){
				if(get_mark(s) < m){
					if(!refinement_is_allowed(s))
						sel.select(s, RM_DUMMY);
					else
						sel.select(s, m);
					if(closure)
						m_aaClosure[s] = true;
				}
			}end_for;
		}
	}lg_end_for;
}

template <typename TSide>
void RegularRefiner_MultiGrid::
adjust_side_of_states (
	size_t lvl,
	uint considerSideMarks,
	uint ignoreElemMarks,
	RefinementMark newElemMark,
	bool closure)
{
	using elem_t = typename TSide::sideof ;

	MultiGrid& mg = *m_pMG;
	MGSelector& sel = m_marks;
	typename Grid::traits<elem_t>::secure_container elems;

	lg_for_each_in_lvl_template(TSide, s, sel, lvl){
		uint m = get_mark(s);
		if(m & considerSideMarks){
			mg.associated_elements(elems, s);
			for_each_in_vec(elem_t* elem, elems){
				if((get_mark(elem) & ignoreElemMarks) != 0){
					if(!refinement_is_allowed(s))
						sel.select(s, RM_DUMMY);
					else
						sel.select(s, newElemMark);
					if(closure)
						m_aaClosure[s] = true;
				}
			}end_for;
		}
	}lg_end_for;
}

template <typename TElem>
void RegularRefiner_MultiGrid::
clear_dummies () {
	for(size_t lvl = 0; lvl < m_marks.num_levels(); ++lvl){
		lg_for_each_in_lvl_template(TElem, e, m_marks, lvl){
			if(get_mark(e) == RM_DUMMY)
				mark(e, RM_NONE);
		}lg_end_for;
	}
}

template <typename TElem>
void RegularRefiner_MultiGrid::
mark_by_level_discrepancy (int lvl, Grid::VertexAttachmentAccessor<AInt> aaLvl)
{
	MultiGrid& mg = *m_pMG;
	MGSelector& sel = m_marks;
	typename Grid::traits<TElem>::secure_container elems;

	lg_for_each_in_lvl(Vertex, vrt, sel, lvl){
		if(aaLvl[vrt] > lvl + 1){
			mg.associated_elements(elems, vrt);
			for_each_in_vec(TElem* e, elems){
				if(((get_mark(e) & LIFT) == 0) && refinement_is_allowed(e)){
					mark(e, RM_ANISOTROPIC);
				}
			}end_for;
		}
	}lg_end_for;
}

void RegularRefiner_MultiGrid::
collect_objects_for_refine ()
{
	MultiGrid& mg = *m_pMG;
	MGSelector& sel = m_marks;

	Grid::edge_traits::secure_container assEdges;
	Grid::face_traits::secure_container assFaces;
	Grid::face_traits::secure_container assVols;

//	first make sure that no element is marked with RM_DUMMY
	clear_dummies<Vertex>();
	clear_dummies<Edge>();
	clear_dummies<Face>();
	clear_dummies<Volume>();

	AInt aLvl;
	mg.attach_to_vertices(aLvl);
	Grid::VertexAttachmentAccessor<AInt> aaLvl(mg, aLvl);
	SetAttachmentValues(aaLvl, mg.begin<Vertex>(), mg.end<Vertex>(), -1);

//	iterate over all levels from top to bottom and adjust marks
	for(int lvl = (int)mg.num_levels(); lvl >= 0; --lvl){
	//	mark connected elements which have a level discrepancy > 1 with any
	//	associated vertex. Note that only marked vertices can cause such a
	//	discrepancy.
		mark_by_level_discrepancy<Edge>(lvl, aaLvl);
		mark_by_level_discrepancy<Face>(lvl, aaLvl);
		mark_by_level_discrepancy<Volume>(lvl, aaLvl);

	//	sides of marked elements have to be marked too
		copy_state_to_sides<Volume>(lvl, LIFT, false);
		copy_state_to_sides<Face>(lvl, LIFT | RM_DUMMY, false);
		copy_state_to_sides<Edge>(lvl, LIFT | RM_DUMMY, false);

	//todo: communicate marks horizontally for vertices, edges, faces

	//	unmarked elements with marked edges will be used to create a closure
		adjust_side_of_states<Edge>(lvl, LIFT, LIFT, RM_ANISOTROPIC, true);
		adjust_side_of_states<Face>(lvl, LIFT, LIFT, RM_ANISOTROPIC, true);

	//	sides of closure elements which are not yet marked have to be marked for closure, too
		adjust_side_states<Volume>(lvl, LIFT, LIFT, RM_COPY, true);
		adjust_side_states<Face>(lvl, LIFT | RM_DUMMY, LIFT, RM_COPY, true);
		adjust_side_states<Edge>(lvl, LIFT | RM_DUMMY, LIFT, RM_COPY, true);

	//todo: communicate marks horizontally for vertices, edges, faces
	
	//	visit all marked vertices in the current level and assign their target-levels
		lg_for_each_in_lvl(Vertex, vrt, sel, lvl){
			if(aaLvl[vrt] == -1)
				aaLvl[vrt] = lvl + 1;
		}lg_end_for;

	//todo: notify parent vertices. Don't use mark here, since we'll also select
	//		shadow vertices.
	//note: Only shadow-vertices which are connected to elements that are allowed
	//		for refinement are of interest here. Those are guaranteed not to be ghosts.
	//		We thus don't have to communicate over vertical interfaces here.
		if(lvl > 0){
			lg_for_each_in_lvl(Vertex, vrt, sel, lvl){
				Vertex* parent = dynamic_cast<Vertex*> (mg.get_parent(vrt));
				if(parent){
					sel.select(parent, RM_DUMMY);
					aaLvl[parent] = aaLvl[vrt];
				}
			}lg_end_for;
		}
	}
}


bool RegularRefiner_MultiGrid::
perform_coarsening () {
	return false;
}


void RegularRefiner_MultiGrid::
num_marked_edges_local (std::vector<int>& numMarkedEdgesOut) {
	numMarkedEdgesOut.resize( m_marks.num_levels() );
	for(size_t lvl = 0; lvl < m_marks.num_levels(); ++lvl){
		numMarkedEdgesOut[lvl] = m_marks.num<Edge>(lvl);
	}
}


void RegularRefiner_MultiGrid::
num_marked_faces_local (std::vector<int>& numMarkedFacesOut) {
	numMarkedFacesOut.resize( m_marks.num_levels() );
	for(size_t lvl = 0; lvl < m_marks.num_levels(); ++lvl){
		numMarkedFacesOut[lvl] = m_marks.num<Face>(lvl);
	}
}


void RegularRefiner_MultiGrid::
num_marked_volumes_local (std::vector<int>& numMarkedVolsOut) {
	numMarkedVolsOut.resize( m_marks.num_levels() );
	for(size_t lvl = 0; lvl < m_marks.num_levels(); ++lvl){
		numMarkedVolsOut[lvl] = m_marks.num<Volume>(lvl);
	}
}


bool RegularRefiner_MultiGrid::
refinement_is_allowed (Vertex* v) {
//todo: consider parallelization and ghost elements
	return !m_pMG->has_children(v);
}


bool RegularRefiner_MultiGrid::
refinement_is_allowed (Edge* e) {
//todo: consider parallelization and ghost elements
	return !m_pMG->has_children(e);
}


bool RegularRefiner_MultiGrid::
refinement_is_allowed (Face* f) {
//todo: consider parallelization and ghost elements
	return !m_pMG->has_children(f);
}


bool RegularRefiner_MultiGrid::
refinement_is_allowed (Volume* v) {
//todo: consider parallelization and ghost elements
	return !m_pMG->has_children(v);
}



}//	end of namespace
