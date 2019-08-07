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

#include "mg_hnode_adjuster.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug{

void MGHNodeAdjuster::
ref_marks_changed(IRefiner& ref,
			   	  const std::vector<Vertex*>& vrts,
			   	  const std::vector<Edge*>& edges,
			   	  const std::vector<Face*>& faces,
			   	  const std::vector<Volume*>& vols)
{
	UG_ASSERT(ref.grid(), "A refiner has to operate on a grid, before marks can be adjusted!");
	if(!ref.grid()){
		return;
	}
	
	MultiGrid* pmg = dynamic_cast<MultiGrid*>(ref.grid());
	UG_ASSERT(pmg, "MGHNodeAdjuster can only operate on multi-grids, not on standard grids.");
	if(!pmg)
		return;
	MultiGrid& mg = *pmg;

//	select all associated vertices of marked objects,
//	since we have to create new vertices in the next levels of the hierarchies.
//	only vertices which do not already have child vertices are selected.
	for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
		Edge* e = edges[i_edge];
		for(size_t i = 0; i < e->num_vertices(); ++i){
			if(!mg.has_children(e->vertex(i)))
				ref.mark(e->vertex(i));
		}
	}

	for(size_t i_face = 0; i_face < faces.size(); ++i_face){
		Face* f = faces[i_face];
		for(size_t i = 0; i < f->num_vertices(); ++i){
			if(!mg.has_children(f->vertex(i)))
				ref.mark(f->vertex(i));
		}
	}

	for(size_t i_vol = 0; i_vol < vols.size(); ++i_vol){
		Volume* v = vols[i_vol];
		for(size_t i = 0; i < v->num_vertices(); ++i){
			if(!mg.has_children(v->vertex(i)))
				ref.mark(v->vertex(i));
		}
	}


//	since we have to make sure that surface elements which meet at a vertex have
//	a level-distance of at most 1, we now mark parent vertices of marked vertices
//	which are connected to a constrained edge.
//	Those marked parents are then used to mark associated edges, faces and volumes
	Grid::edge_traits::secure_container assEdges;
	Grid::face_traits::secure_container assFaces;
	Grid::volume_traits::secure_container assVols;
	for(size_t i_vrt = 0; i_vrt < vrts.size(); ++i_vrt){
		Vertex* vrt = vrts[i_vrt];
		if(vrt->is_constrained())
			continue;

		if(mg.num_child_vertices(vrt) > 0){
		//	we have to mark all associated edges, faces and volumes
			mg.associated_elements(assEdges, vrt);
			for(size_t i = 0; i < assEdges.size(); ++i)
				if (ref.get_mark(assEdges[i]) == RM_NONE) // might already be marked RM_CLOSURE, do not overwrite
					ref.mark(assEdges[i]);
			mg.associated_elements(assFaces, vrt);
			for(size_t i = 0; i < assFaces.size(); ++i)
				if (ref.get_mark(assFaces[i]) == RM_NONE)
					ref.mark(assFaces[i]);
			mg.associated_elements(assVols, vrt);
			for(size_t i = 0; i < assVols.size(); ++i)
				if (ref.get_mark(assVols[i]) == RM_NONE)
					ref.mark(assVols[i]);
		}
		else if(ref.get_mark(vrt) != RM_DUMMY){
		//	we don't have to select parents of dummy vertices, since we assume
		//	that the maximum level-distance is 1
			Vertex* parent = dynamic_cast<Vertex*>(mg.get_parent(vrt));
			if(parent)
				ref.mark(parent, RM_DUMMY);
		}
	}
}
}// end of namespace
