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

#include "bool_marker.h"

namespace ug{

BoolMarker::BoolMarker() :
	m_pGrid(nullptr),
	m_defaultMark(false),
	m_markInheritanceEnabled(true),
	m_strictInheritanceEnabled(false)
{
}

BoolMarker::BoolMarker(Grid& g) :
	m_pGrid(nullptr),
	m_defaultMark(false),
	m_markInheritanceEnabled(true),
	m_strictInheritanceEnabled(false)
{
	assign_grid(&g);
}

BoolMarker::~BoolMarker()
{
	assign_grid(nullptr);
}

void BoolMarker::assign_grid(Grid* g)
{
	if(m_pGrid == g)
		return;

	if(m_pGrid){
		m_pGrid->detach_from_all(m_aBool);
		m_pGrid->unregister_observer(this);
		m_aaMarkVRT.invalidate();
		m_aaMarkEDGE.invalidate();
		m_aaMarkFACE.invalidate();
		m_aaMarkVOL.invalidate();
	}

	m_pGrid = g;
	if(g){
		g->register_observer(this, OT_GRID_OBSERVER | OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER |
									OT_FACE_OBSERVER | OT_VOLUME_OBSERVER);
		g->attach_to_all(m_aBool);
		m_aaMarkVRT.access(*g, m_aBool);
		m_aaMarkEDGE.access(*g, m_aBool);
		m_aaMarkFACE.access(*g, m_aBool);
		m_aaMarkVOL.access(*g, m_aBool);
	}
}


bool BoolMarker::is_marked(GridObject* e) const
{
	switch(e->base_object_id()){
		case VERTEX: return is_marked(static_cast<Vertex*>(e));
		case EDGE: return is_marked(static_cast<Edge*>(e));
		case FACE: return is_marked(static_cast<Face*>(e));
		case VOLUME: return is_marked(static_cast<Volume*>(e));
		default: return false;
	}
}




void BoolMarker::grid_to_be_destroyed(Grid* grid)
{
	assign_grid(nullptr);
}

void BoolMarker::clear()
{
	assert(m_pGrid);
	unmark(m_pGrid->begin<Vertex>(), m_pGrid->end<Vertex>());
	unmark(m_pGrid->begin<Edge>(), m_pGrid->end<Edge>());
	unmark(m_pGrid->begin<Face>(), m_pGrid->end<Face>());
	unmark(m_pGrid->begin<Volume>(), m_pGrid->end<Volume>());
}


void BoolMarker::
vertex_created(Grid* grid, Vertex* vrt, GridObject* pParent,
				bool replacesParent)
{
	if(!pParent){
		mark(vrt, default_mark());
		return;
	}

	if(strict_inheritance_enabled() && (pParent->base_object_id() != VERTEX)){
		mark(vrt, default_mark());
		return;
	}

	if(mark_inheritance_enabeld())
		mark(vrt, is_marked(pParent));

	mark(vrt, default_mark());
}

void BoolMarker::
edge_created(Grid* grid, Edge* e, GridObject* pParent,
			 bool replacesParent)
{
	if(!pParent){
		mark(e, default_mark());
		return;
	}

	if(strict_inheritance_enabled() && (pParent->base_object_id() != EDGE)){
		mark(e, default_mark());
		return;
	}

	if(mark_inheritance_enabeld())
		mark(e, is_marked(pParent));

	mark(e, default_mark());
}

void BoolMarker::
face_created(Grid* grid, Face* f, GridObject* pParent,
			 bool replacesParent)
{
	if(!pParent){
		mark(f, default_mark());
		return;
	}

	if(strict_inheritance_enabled() && (pParent->base_object_id() != FACE)){
		mark(f, default_mark());
		return;
	}

	if(mark_inheritance_enabeld())
		mark(f, is_marked(pParent));

	mark(f, default_mark());
}

void BoolMarker::
volume_created(Grid* grid, Volume* vol, GridObject* pParent,
			   bool replacesParent)
{
	if(!pParent){
		mark(vol, default_mark());
		return;
	}

	if(strict_inheritance_enabled() && (pParent->base_object_id() != VOLUME)){
		mark(vol, default_mark());
		return;
	}

	if(mark_inheritance_enabeld())
		mark(vol, is_marked(pParent));

	mark(vol, default_mark());
}

void BoolMarker::
vertices_to_be_merged(Grid* grid, Vertex* target,
					  Vertex* elem1, Vertex* elem2)
{
	mark(target, is_marked(elem1) || is_marked(elem2));
}

void BoolMarker::
edges_to_be_merged(Grid* grid, Edge* target,
				   Edge* elem1, Edge* elem2)
{
	mark(target, is_marked(elem1) || is_marked(elem2));
}

void BoolMarker::
faces_to_be_merged(Grid* grid, Face* target,
				   Face* elem1, Face* elem2)
{
	mark(target, is_marked(elem1) || is_marked(elem2));
}

void BoolMarker::
volumes_to_be_merged(Grid* grid, Volume* target,
					 Volume* elem1, Volume* elem2)
{
	mark(target, is_marked(elem1) || is_marked(elem2));
}

}// end of namespace
