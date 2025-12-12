/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__duplicate__
#define __H__UG__duplicate__

#include <vector>

#include "lib_grid/lg_base.h"
#include "selection_util.h"

namespace ug {

///	Duplicates the selected part of a grid and translates it by an offset
/**	NOTE: This method is not suited for MultiGrids.
 * Duplicates all selected elements and associated lower dimensional elements.
 * New vertices are shifted by the given offset.
 * You may specify whether the old selection shall be cleared and whether the
 * new elements shall be selected. Note that this overrules the settings
 * regarding autoselection and selection-inheritance in the given selector.
 */
template <typename TAPos>
bool Duplicate(Grid& grid, Selector& sel, const typename TAPos::ValueType& offset,
			   TAPos& aPos, bool deselectOld = true, bool selectNew = true)
{
	using namespace std;

	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG("ERROR in Duplicate. Given position attachment not attached to vertices.\n");
	}

	Grid::VertexAttachmentAccessor<TAPos> aaPos(grid, aPos);


//	first of all we have to select all associated elements
	SelectAssociatedGridObjects(sel);

//	store the currently selected elements in arrays
	vector<Vertex*>	oldVrts;
	vector<Edge*> oldEdges;
	vector<Face*> oldFaces;
	vector<Volume*> oldVols;
	oldVrts.reserve(sel.num<Vertex>());
	oldVrts.assign(sel.vertices_begin(), sel.vertices_end());
	oldEdges.reserve(sel.num<Edge>());
	oldEdges.assign(sel.edges_begin(), sel.edges_end());
	oldFaces.reserve(sel.num<Face>());
	oldFaces.assign(sel.faces_begin(), sel.faces_end());
	oldVols.reserve(sel.num<Volume>());
	oldVols.assign(sel.volumes_begin(), sel.volumes_end());

//	depending on the options we'll now clear the selector and enable autoselection
	bool autoselectionEnabled = sel.autoselection_enabled();
	bool selectionInheritanceEnabled = sel.selection_inheritance_enabled();
	if(deselectOld)
		sel.clear();
	if(selectNew)
		sel.enable_autoselection(true);
	else{
		sel.enable_autoselection(false);
		sel.enable_selection_inheritance(false);
	}

//	new vertices have to be stored in an array, so that we can access them later on
	vector<Vertex*> vrts;
	vrts.reserve(sel.num<Vertex>());

//	attach an integer to the vertices. This is required, since
//	we have to access them by index later on.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);

//	create new vertices and push them to the vrts array.
//	Store the new index in the aInt attachment of old ones.
	for(vector<Vertex*>::iterator iter = oldVrts.begin();
		iter != oldVrts.end(); ++iter)
	{
		Vertex* oldVrt = *iter;
		Vertex* vrt = *grid.create_by_cloning(oldVrt, oldVrt);
		VecAdd(aaPos[vrt], aaPos[oldVrt], offset);
		aaInt[oldVrt] = (int)vrts.size();
		vrts.push_back(vrt);
	}

//	create new elements for all vector-entries
//	edges
	for(size_t i = 0; i < oldEdges.size(); ++i){
		Edge* oldElem = oldEdges[i];
		grid.create_by_cloning(oldElem, EdgeDescriptor(vrts[aaInt[oldElem->vertex(0)]],
											   	   	   vrts[aaInt[oldElem->vertex(1)]]),
							   oldElem);
	}

//	faces
	{
		FaceDescriptor desc;
		for(size_t i = 0; i < oldFaces.size(); ++i){
			Face* oldElem = oldFaces[i];

			desc.set_num_vertices(oldElem->num_vertices());
			for(size_t j = 0; j < desc.num_vertices(); ++j)
				desc.set_vertex(j, vrts[aaInt[oldElem->vertex(j)]]);

			grid.create_by_cloning(oldElem, desc, oldElem);
		}
	}

//	volumes
	{
		VolumeDescriptor desc;
		for(size_t i = 0; i < oldVols.size(); ++i){
			Volume* oldElem = oldVols[i];

			desc.set_num_vertices(oldElem->num_vertices());
			for(size_t j = 0; j < desc.num_vertices(); ++j)
				desc.set_vertex(j, vrts[aaInt[oldElem->vertex(j)]]);

			grid.create_by_cloning(oldElem, desc, oldElem);
		}
	}

//	restore autoselection property of the selector
	sel.enable_autoselection(autoselectionEnabled);
	sel.enable_selection_inheritance(selectionInheritanceEnabled);
	return true;
}

}//	end of namespace

#endif
