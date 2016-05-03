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

#include <vector>
#include "boost/container/vector.hpp"	// for bool-vectors
#include "regular_refinement_new.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/selection_util.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/callbacks/selection_callbacks.h"

using namespace std;

namespace ug
{

bool RefineNew(Grid& grid, Selector& sel,
			RefinementProjector* projector,
			bool useSnapPoints)
{
	AInt aInt;
	if(grid.num<Face>() > 0)
		grid.attach_to_edges(aInt);
	if(grid.num<Volume>() > 0)
		grid.attach_to_faces(aInt);
		
	bool bSuccess = RefineNew(grid, sel, aInt, projector, useSnapPoints);

	if(grid.num<Face>() > 0)
		grid.detach_from_edges(aInt);
	if(grid.num<Volume>() > 0)
		grid.detach_from_faces(aInt);

	return bSuccess;
}

static void AdjustSelection(Grid& grid, Selector& sel)
{
//	select all edges of selected faces
	vector<Edge*> vEdges;
	for(FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter){
		CollectEdges(vEdges, grid, *iter);
		for(size_t i = 0; i < vEdges.size(); ++i)
			sel.select(vEdges[i]);
	}
	
//	select all edges of selected volumes
	for(VolumeIterator iter = sel.begin<Volume>();
		iter != sel.end<Volume>(); ++iter)
	{
		CollectEdges(vEdges, grid, *iter);
		for(size_t i = 0; i < vEdges.size(); ++i)
			sel.select(vEdges[i]);
	}
	
//	select all faces and volumes which are adjacent to selected edges
	vector<Face*> vFaces;
	vector<Volume*> vVols;
	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
	{
		CollectFaces(vFaces, grid, *iter);
		for(size_t i = 0; i < vFaces.size(); ++i)
			sel.select(vFaces[i]);
			
		CollectVolumes(vVols, grid, *iter);
		for(size_t i = 0; i < vVols.size(); ++i)
			sel.select(vVols[i]);
	}
}


bool RefineNew(Grid& grid, Selector& sel, AInt& aInt,
			RefinementProjector* projector,
			bool useSnapPoints)
{
//	position data is required
	if(!grid.has_vertex_attachment(aPosition)){
		LOG("  WARNING in Refine: aPosition is not attached to the vertices of the grid. Aborting.\n");
		return false;
	}

//	if the refinement-callback is empty, use a linear one.
	RefinementProjector defaultProjector;
	if(!projector){
		defaultProjector.set_geometry(make_sp(new Geometry<3, 3>(grid, aPosition)));
		projector = &defaultProjector;
	}
		
//	make sure that GRIDOPT_VERTEXCENTRIC_INTERCONNECTION is enabled
	if(grid.num_edges() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_EDGES))){
		LOG("  INFO in Refine: autoenabling VRTOPT_STORE_ASSOCIATED_EDGES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES);
	}
	if(grid.num_faces() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_FACES))){
		LOG("  INFO in Refine: autoenabling VRTOPT_STORE_ASSOCIATED_FACES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_FACES);
	}
	if(grid.num_volumes() && (!grid.option_is_enabled(VRTOPT_STORE_ASSOCIATED_VOLUMES))){
		LOG("  INFO in Refine: autoenabling VRTOPT_STORE_ASSOCIATED_VOLUMES\n");
		grid.enable_options(VRTOPT_STORE_ASSOCIATED_VOLUMES);
	}

//	make sure that FACEOPT_AUTOGENERATE_EDGES is enabled
	if(grid.num<Face>() > 0 && (!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))){
		LOG("  INFO in Refine: autoenabling FACEOPT_AUTOGENERATE_EDGES\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}
	
//	if there are volumes, make sure that VOLOPT_AUTOGENERATE_FACES is enabled.
	if(grid.num<Volume>() > 0 && (!grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES)))
	{
		LOG("  INFO in Refine: autoenabling VOLOPT_AUTOGENERATE_FACES\n");
		grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}
	
//	snap-points are recorded here, since we alter the vertex selection and have
//	to restore them later on.
	vector<Vertex*> snapPoints;
	if(useSnapPoints)
		snapPoints.insert(snapPoints.end(), sel.begin<Vertex>(), sel.end<Vertex>());

//	adjust selection
	AdjustSelection(grid, sel);

//	we will select associated vertices, too, since we have to
//	notify the refinement-callback, that they are involved in refinement.
	sel.clear<Vertex>();
	SelectAssociatedVertices(sel, sel.begin<Edge>(), sel.end<Edge>(), ISelector::SELECTED);
	SelectAssociatedVertices(sel, sel.begin<Face>(), sel.end<Face>(), ISelector::SELECTED);
	SelectAssociatedVertices(sel, sel.begin<Volume>(), sel.end<Volume>(), ISelector::SELECTED);
	
//	select snap-vertices with a special mark
	const ISelector::status_t snapSelVal = ISelector::SELECTED + 1;
	if(useSnapPoints){
		for(vector<Vertex*>::iterator iter = snapPoints.begin();
			iter != snapPoints.end(); ++iter)
		{
			sel.select(*iter, snapSelVal);
		}
	}

//	aInt has to be attached to the edges of the grid
	if(sel.num<Face>() > 0 && (!grid.has_edge_attachment(aInt))){
		LOG("  WARNING in Refine: aInt is not attached to the edges of the grid. Aborting.\n");
		return false;
	}

//	if there are selected volumes,
//	aInt has to be attached to the faces of the grid
	if(sel.num<Volume>() && (!grid.has_face_attachment(aInt))){
		LOG("  WARNING in Refine: aInt is not attached to the faces of the grid. Aborting.\n");
		return false;
	}

//	number of edges, faces and volumes that will be refined
	const size_t numRefEdges = sel.num<Edge>();
	const size_t numRefFaces = sel.num<Face>();
	const size_t numRefVols = sel.num<Volume>();
	
//	we need several arrays.
//	this one stores pointers to all edges that shall be refined
	vector<Edge*> edges(numRefEdges);
//	one that stores the new vertex for each edge that shall be refined.
	vector<RegularVertex*>	edgeVrts(numRefEdges);
//	one that stores the selected faces
	vector<Face*> faces(numRefFaces);
//	one that stores vertices which are created on faces
	vector<Vertex*> faceVrts;
	if(numRefVols > 0)
		faceVrts.resize(numRefFaces);
//	one that stores selected volumes
	vector<Volume*> vols(sel.num<Volume>());
	
//	acces the int-attachment
	Grid::EdgeAttachmentAccessor<AInt> aaIntEDGE;
	Grid::FaceAttachmentAccessor<AInt> aaIntFACE;
	
	if(numRefFaces > 0)
		aaIntEDGE.access(grid, aInt);
	if(numRefVols > 0)
		aaIntFACE.access(grid, aInt);


//	notify refinement projector about the start of refinement
	projector->refinement_begins(
				SubGrid<IsSelected>(
					sel.get_grid_objects(), IsSelected(sel)));

//	store the geometry
	IGeometry3d& geom = *projector->geometry();

////////////////////////////////
//	fill the edges- and edgeVrts-array and assign indices to selected edges
	{
		EdgeIterator edgesEnd = sel.end<Edge>();
		int i = 0;
		for(EdgeIterator iter = sel.begin<Edge>();
			iter != edgesEnd; ++iter, ++i)
		{
		//	store the edge
			edges[i] = *iter;
			if(numRefFaces > 0)		
				aaIntEDGE[*iter] = i;
				
		//	create the new vertex
			edgeVrts[i] = *grid.create<RegularVertex>(*iter);
			sel.select(edgeVrts[i]);
		//	calculate new position
			projector->new_vertex(edgeVrts[i], *iter);
		}
	}

////////////////////////////////
//	refine the selected edges
	vector<Edge*> newEdges;
	newEdges.reserve(2);
	for(size_t i = 0; i < edges.size(); ++i){
		Edge* e = edges[i];
		if(e->refine(newEdges, edgeVrts[i])){
			for(size_t j = 0; j < newEdges.size(); ++j)
				grid.register_element(newEdges[j], e);
		}
		else{
			LOG("  WARNING in Refine: could not refine edge.\n");
		}
	}
	
////////////////////////////////
//	set up face arrays
	{		
		FaceIterator facesEnd = sel.end<Face>();
		int i = 0;
		for(FaceIterator iter = sel.begin<Face>();
			iter != facesEnd; ++iter, ++i)
		{
			Face* f = *iter;
			faces[i] = f;
			if(numRefVols > 0)
				aaIntFACE[f] = i;
		}
	}

////////////////////////////////
//	refine the selected faces
	vector<Face*> newFaces;
	newFaces.reserve(4);
//	we need a container that stores the vertex for each edge of a face
//	entries will be set to NULL if the associated edge will not be refined
	vector<Vertex*> faceEdgeVrts;
	faceEdgeVrts.reserve(4);
	
	for(size_t i = 0; i < faces.size(); ++i){
		Face* f = faces[i];
		Vertex* newVrt;
		
	//	collect vertices of associated edges
		faceEdgeVrts.clear();
		for(uint j = 0; j < f->num_edges(); ++j){
			Edge* e = grid.get_edge(f, j);
			if(sel.is_selected(e))
				faceEdgeVrts.push_back(edgeVrts[aaIntEDGE[e]]);
			else
				faceEdgeVrts.push_back(NULL);
		}

		int snapPointIndex = -1;
		if(useSnapPoints){
			Face::ConstVertexArray vrts = f->vertices();
			const size_t numVrts = f->num_vertices();
			for(size_t s = 0; s < numVrts; ++s){
				if(sel.get_selection_status(vrts[s]) == snapSelVal){
					if(snapPointIndex != -1){
						UG_LOG("WARNING: Only one snap-point per face is allowed, "
							"but more are present. Ignoring snap-points for this face.\n");
						snapPointIndex = -1;
						break;
					}
					snapPointIndex = static_cast<int>(s);
				}
			}
		}

		if(f->refine(newFaces, &newVrt, &faceEdgeVrts.front(), NULL, NULL, snapPointIndex)){
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				grid.register_element(newVrt, f);
				projector->new_vertex(newVrt, f);
				sel.select(newVrt);
			}
			
		//	if volumes are refined too, we have to store the vertex
			if(numRefVols > 0)
				faceVrts[i] = newVrt;
				
		//	register the new faces
			for(size_t j = 0; j < newFaces.size(); ++j)
				grid.register_element(newFaces[j], f);
		}
		else{
			LOG("  WARNING in Refine: could not refine face.\n");
		}
	}

////////////////////////////////
//	set up volume arrays
	{							
		VolumeIterator volsEnd = sel.end<Volume>();
		int i = 0;
		for(VolumeIterator iter = sel.begin<Volume>();
			iter != volsEnd; ++iter, ++i)
		{
			Volume* v = *iter;
			vols[i] = v;
		}
	}
	
////////////////////////////////
//	refine the selected volumes
	vector<Volume*> newVols;
	newVols.reserve(8);
//	we need a container that stores the vertex for each edge of a volume
//	entries will be set to NULL if the associated edge will not be refined
	vector<Vertex*> volEdgeVrts;
	volEdgeVrts.reserve(12);
//	we need a container that stores the vertex for each face of a volume
//	entries will be set to NULL if the associated face will not be refined
	vector<Vertex*> volFaceVrts;
	volFaceVrts.reserve(6);
	
//	only used for tetrahedron refinement
	vector<vector3> corners(4, vector3(0, 0, 0));

// //	DEBUG
// 	UG_LOG("> VOL-REF-BEGIN\n");
// 	UG_LOG("> DEBUG-ACCESSOR...\n");
// 	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	boost::container::vector<bool> isSnapPoint;
	if(useSnapPoints)
		isSnapPoint.reserve(8);

	for(size_t i = 0; i < vols.size(); ++i){
		Volume* v = vols[i];
		Vertex* newVrt;
		
	//	collect vertices of associated edges
		volEdgeVrts.clear();
		for(uint j = 0; j < v->num_edges(); ++j){
			Edge* e = grid.get_edge(v, j);
			if(sel.is_selected(e))
				volEdgeVrts.push_back(edgeVrts[aaIntEDGE[e]]);
			else
				volEdgeVrts.push_back(NULL);
		}

	//	collect vertices of associated faces
		volFaceVrts.clear();
		for(uint j = 0; j < v->num_faces(); ++j){
			Face* f = grid.get_face(v, j);
			if(sel.is_selected(f))
				volFaceVrts.push_back(faceVrts[aaIntFACE[f]]);
			else
				volFaceVrts.push_back(NULL);
		}
		
	//	if we're performing tetrahedral refinement, we have to collect
	//	the corner coordinates, so that the refinement algorithm may choose
	//	the best interior diagonal.
		vector3* pCorners = NULL;
		if((v->num_vertices() == 4)){
			for(size_t i = 0; i < 4; ++i){
				corners[i] = geom.pos(v->vertex(i));
			}
			pCorners = &corners.front();
		}

		bool* pIsSnapPoint = NULL;
		if(useSnapPoints){
			Volume::ConstVertexArray vrts = v->vertices();
			const size_t numVrts = v->num_vertices();
			isSnapPoint.clear();

			bool gotOne = false;
			for(size_t s = 0; s < numVrts; ++s){
				const bool val = sel.get_selection_status(vrts[s]) == snapSelVal;
				isSnapPoint.push_back(val);
				gotOne = gotOne | val;
			}

			if(gotOne)
				pIsSnapPoint = &isSnapPoint.front();
		}

		if(v->refine(newVols, &newVrt, &volEdgeVrts.front(),
					&volFaceVrts.front(), NULL, RegularVertex(), NULL,
					pCorners, pIsSnapPoint))
		{
		//	if a new vertex was generated, we have to register it
			if(newVrt){
				grid.register_element(newVrt, v);
				projector->new_vertex(newVrt, v);
				sel.select(newVrt);
			}
				
		//	register the new volumes
			for(size_t j = 0; j < newVols.size(); ++j)
				grid.register_element(newVols[j], v);
		}
		else{
			LOG("  WARNING in Refine: could not refine volume.\n");
		}
	}
// UG_LOG("> VOL-REF-END\n");

//	erase old volumes
	grid.erase(vols.begin(), vols.end());
//	erase old faces
	grid.erase(faces.begin(), faces.end());
//	erase old edges
	grid.erase(edges.begin(), edges.end());

//	notify refinement projector about the end of refinement
	projector->refinement_ends(
				SubGrid<IsSelected>(
					sel.get_grid_objects(), IsSelected(sel)));
	return true;
}

}//	end of namespace
