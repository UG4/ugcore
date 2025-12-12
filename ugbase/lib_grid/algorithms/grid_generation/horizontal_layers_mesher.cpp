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

#include "horizontal_layers_mesher.h"

#include <algorithm>
#include <vector>

#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/algorithms/extrusion/extrude.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "lib_grid/file_io/file_io_asc.h"
#include "lib_grid/iterators/associated_elements_iterator.h"
//ø #include "lib_grid/iterators/lg_for_each.h"

using namespace std;

namespace ug {

////////////////////////////////////////////////////////////////////////////////
void MeshLayerBoundaries(Grid& grid, const RasterLayers& layers,
						 Grid::VertexAttachmentAccessor<AVector3> aaPos,
						 ISubsetHandler* pSH)
{
	int defSubInd = -1;
	if(pSH)
		defSubInd = pSH->get_default_subset_index();

	for(size_t lvl = 0; lvl < layers.size(); ++lvl){
		if(pSH)
			pSH->set_default_subset_index((int)lvl);
		CreateGridFromFieldBoundary(grid, layers.heightfield(lvl), aaPos);
	}

	if(pSH)
		pSH->set_default_subset_index(defSubInd);
}

////////////////////////////////////////////////////////////////////////////////
void MeshLayers(Grid& grid, const RasterLayers& layers,
				Grid::VertexAttachmentAccessor<AVector3> aaPos,
				ISubsetHandler* pSH)
{
	int defSubInd = -1;
	if(pSH)
		defSubInd = pSH->get_default_subset_index();

	for(size_t lvl = 0; lvl < layers.size(); ++lvl){
		if(pSH)
			pSH->set_default_subset_index((int)lvl);
		CreateGridFromField(grid, layers.heightfield(lvl), aaPos);
	}

	if(pSH)
		pSH->set_default_subset_index(defSubInd);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	EXTRUDE LAYERS
struct ConnectedToOneMarkedVrt{
	ConnectedToOneMarkedVrt(Grid& grid) : m_grid(grid) {}
	bool operator () (Edge* e) const{
		return	(m_grid.is_marked(e->vertex(0)) || m_grid.is_marked(e->vertex(1)))
			&& !(m_grid.is_marked(e->vertex(0)) && m_grid.is_marked(e->vertex(1)));
	}
	Grid& m_grid;
};

void ExtrudeLayers (
		Grid& grid, 
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler& sh,
		bool allowForTetsAndPyras,
		const ANumber* aRelZOut)
{
	if(allowForTetsAndPyras){
		ExtrudeLayersMixed(grid, layers, aaPos, sh, aRelZOut);
		return;
	}

	UG_COND_THROW(layers.size() < 2, "At least 2 layers are required to perform extrusion!");

	grid.begin_marking();

	vector<Vertex*> curVrts;	// list of vertices that are considered for extrusion
	vector<Vertex*> tmpVrts;	// used to determine the set of vertices that can be extruded
	vector<Vertex*> smoothVrts;
	vector<Face*> curFaces;		// list of faces that are considered for extrusion
	vector<Face*> tmpFaces;		// used to determine the set of faces that can be extruded
	vector<number> vrtHeightVals;	// here we'll record height-values at which new vertices will be placed
	vector<number> volHeightVals;	// here we'll record height-values at the center of new volumes
	vector<number> curUpperLayerHeight;
	vector<int> volSubsetInds;
	vector<Volume*> newVols;
	queue<Volume*> volCandidates;

	Grid::edge_traits::callback cbIsMarked = IsMarked(grid);
	Grid::edge_traits::callback cbConnectedToOneMarkedVrt = ConnectedToOneMarkedVrt(grid);
	AssocElemIter<Vertex, Edge> assocVrtEdgeIterMarkedEdge(cbIsMarked);
	AssocElemIter<Vertex, Edge> assocVrtEdgeIterOneMarked(cbConnectedToOneMarkedVrt);
	AssocElemIter<Face, Edge> assocFaceEdgeIter;
	AssocElemIter<Volume, Edge> assocVolEdgeIter(cbConnectedToOneMarkedVrt);
	
	grid.reserve<Vertex>(grid.num_vertices() * layers.size());
	grid.reserve<Volume>(grid.num_faces() * layers.size() - 1);

//	todo: this accessor is primarily used during smoothing. Maybe it can be removed?
	ANumber aHeight;
	grid.attach_to_vertices(aHeight);
	Grid::VertexAttachmentAccessor<ANumber> aaHeight(grid, aHeight);

	Grid::VertexAttachmentAccessor<ANumber> aaRelZOut;
	if(aRelZOut){
		ANumber aRelZ = *aRelZOut;
		if(!grid.has_vertex_attachment(aRelZ))
			grid.attach_to_vertices(aRelZ);
		aaRelZOut.access(grid, aRelZ);
	}

//	we have to determine the vertices that can be projected onto the top of the
//	given layers-stack. Only those will be used during extrusion
//	all considered vertices will be marked.
	const RasterLayers::layer_t& top = layers.top();
	const int topLayerInd = (int)layers.num_layers() - 1;
	for(VertexIterator i = grid.begin<Vertex>(); i != grid.end<Vertex>(); ++i){
		Vertex* v = *i;
		vector2 c(aaPos[v].x(), aaPos[v].y());
		// number val = top.heightfield.interpolate(c);
		number val = layers.relative_to_global_height(c, topLayerInd);
		if(val != top.heightfield.no_data_value()){
			aaPos[v].z() = val;
			curVrts.push_back(v);
			curUpperLayerHeight.push_back(val);
			grid.mark(v);
			sh.assign_subset(v, topLayerInd);
		}
	}

	if(curVrts.size() < 3){
		UG_LOG("Not enough vertices lie in the region of the surface layer\n");
		return;
	}

	if(aaRelZOut.valid()){
		for(size_t ivrt = 0; ivrt < curVrts.size(); ++ivrt)
			aaRelZOut[curVrts[ivrt]] = topLayerInd;
	}

//	all faces of the initial triangulation that are only connected to marked
//	vertices are considered for extrusion
	for(FaceIterator fi = grid.begin<Face>(); fi != grid.end<Face>(); ++fi){
		Face* f = *fi;
		bool allMarked = true;
		Face::ConstVertexArray vrts = f->vertices();
		const size_t numVrts = f->num_vertices();
		for(size_t i = 0; i < numVrts; ++i){
			if(!grid.is_marked(vrts[i])){
				allMarked = false;
				break;
			}
		}
		if(allMarked)
			curFaces.push_back(f);
	}

	if(curFaces.empty()){
		UG_LOG("Not enough faces lie in the region of the surface layer\n");
		return;
	}


	tmpVrts.reserve(curVrts.size());
	smoothVrts.reserve(curVrts.size());
	curUpperLayerHeight.reserve(curVrts.size());
	tmpFaces.reserve(curFaces.size());
	newVols.reserve(curFaces.size());
	volHeightVals.reserve(curFaces.size());
	volSubsetInds.reserve(curFaces.size());
	const int invalidSub = max<int>(sh.num_subsets(), layers.size() + 1);
	vector<Volume*> invalidVols;

	for(int ilayer = (int)layers.size() - 2; ilayer >= 0; --ilayer){

		tmpVrts.clear();
		vrtHeightVals.clear();
		tmpFaces.clear();
		newVols.clear();
		volHeightVals.clear();
		volSubsetInds.clear();
		grid.clear_marks();
	//	trace rays from the current vertices down through the layers until the
	//	next valid entry is found. If none is found, the vertex will either be ignored
	//	from then on (allowForTetsAndPyras == false) or a dummy vertex will be inserted.
		for(size_t icur = 0; icur < curVrts.size(); ++icur){
			Vertex* v = curVrts[icur];
			vector2 c(aaPos[v].x(), aaPos[v].y());
			number upperVal = curUpperLayerHeight[icur];
			number height = layers.relative_to_global_height(c, (number) ilayer);

			// bool searching = true;
			int curLayer = ilayer;

			// while(searching && curLayer >= 0){
				// searching = false;
				pair<int, number> val = layers.trace_line_down(c, curLayer);
				number lowerVal = 0;
				if(val.first >= 0)
					lowerVal = layers.relative_to_global_height(c, (number) val.first);

				// if(curLayer == 0 && val.first >= 0
				// 	&& upperVal - lowerVal < layers.min_height(curLayer))
				// {
				// 	val.first = -1;
				// }
				
				if(val.first >= 0){
					// if(upperVal - lowerVal < layers.min_height(curLayer)){
					// 	searching = true;
					// 	--curLayer;
					// 	continue;
					// }
				//	if val.first == curLayer height will equal lowerVal. If not,
				//	a linear interpolation is performed, considering the height-val
				//	of the current vertex, the layer distance and the target value.
					tmpVrts.push_back(v);
					vrtHeightVals.push_back(height);
					aaHeight[v] = upperVal - lowerVal;//total height of curLayer
					sh.assign_subset(v, val.first);
					grid.mark(v);
					if(val.first == ilayer){
						// UG_LOG("DBG: setting upper value: " << lowerVal << "(old: " << upperVal << ")\n");
						curUpperLayerHeight[icur] = lowerVal;
					}
				}
				else if(allowForTetsAndPyras){
				//	we insert a dummy-vertex which will later on allow for easier
				//	edge-collapses of inner vertical rim edges
					tmpVrts.push_back(v);
					vrtHeightVals.push_back(height);
					aaHeight[v] = upperVal - val.second;//total height of ilayer
					sh.assign_subset(v, invalidSub);
					grid.mark(v);
				}
			// }
		}
	//	now find the faces which connect those vertices
		for(size_t iface = 0; iface < curFaces.size(); ++iface){
			Face* f = curFaces[iface];

			vector3 center = CalculateCenter(f, aaPos);
			vector2 c(center.x(), center.y());
			pair<int, number> val = layers.trace_line_down(c, ilayer);
			if((val.first < 0) && !allowForTetsAndPyras)
				continue;

		//	now check whether all vertices are marked
			bool allMarked = true;
			Face::ConstVertexArray vrts = f->vertices();
			size_t numVrts = f->num_vertices();
			for(size_t i = 0; i < numVrts; ++i){
				if(!grid.is_marked(vrts[i])){
					allMarked = false;
					break;
				}
			}

			if(allMarked){
				int subsetIndex = invalidSub;

				// if(val.first >= 0) {
				// //	subset from center trace down
				// 	subsetIndex = val.first;
				// }
				// else {
				//	the highest valid side of the volume determines its subset.
					number maxHeight = 0;
					for(size_t i = 0; i < numVrts; ++i){
						int si = sh.get_subset_index(vrts[i]);
						if(si == invalidSub)
							continue;

						number vh = aaHeight[vrts[i]];
						if(vh > maxHeight){
							maxHeight = vh;
							subsetIndex = si;
						}
					}
				// }



				pair<int, number> upVal = layers.trace_line_up(c, ilayer+1);
				if(val.first < 0 || upVal.first < 0){
					if(allowForTetsAndPyras){
						tmpFaces.push_back(f);
						// volSubsetInds.push_back(invalidSub);
						volSubsetInds.push_back(subsetIndex);
						volHeightVals.push_back(layers.min_height(ilayer));
					}
				}
				else{
					tmpFaces.push_back(f);
					// volSubsetInds.push_back(val.first);
					volSubsetInds.push_back(subsetIndex);
					volHeightVals.push_back(upVal.second - val.second);
				}
			}
		}

		if(tmpFaces.empty())
			break;

	//	swap containers and perform extrusion
		curVrts.swap(tmpVrts);
		curFaces.swap(tmpFaces);
		Extrude(grid, &curVrts, nullptr, &curFaces, vector3(0, 0, 0), aaPos,
		        ExtrusionOptions::EO_DEFAULT, &newVols);

	// assign pre-determined subsets
		for(size_t ivol = 0; ivol < newVols.size(); ++ivol){
			sh.assign_subset(newVols[ivol], volSubsetInds[ivol]);
			if(volSubsetInds[ivol] == invalidSub)
				invalidVols.push_back(newVols[ivol]);
		}

	//	set the precalculated height of new vertices
		for(size_t ivrt = 0; ivrt < curVrts.size(); ++ivrt){
			aaPos[curVrts[ivrt]].z() = vrtHeightVals[ivrt];
		}

		if(aaRelZOut.valid()){
			for(size_t ivrt = 0; ivrt < curVrts.size(); ++ivrt)
				aaRelZOut[curVrts[ivrt]] = ilayer;
		}
	}
	
//	remove unnecessary prisms through edge-collapses and thus introduce pyramids and tetrahedra
	if(allowForTetsAndPyras){
		ABool aInterface;
		grid.attach_to_vertices_dv(aInterface, false, true);
		Grid::VertexAttachmentAccessor<ABool> aaInterface(grid, aInterface);

		Grid::volume_traits::secure_container	assVols;

	//	all triangle-interface-elements shall store 'true' in aaIsInterface
		for(Grid::traits<Face>::iterator _feI = grid.begin<Face>(); _feI != grid.end<Face>(); ++_feI){ Face* f = *_feI;{
			if(f->num_vertices() != 3)
				continue;

			grid.associated_elements(assVols, f);
			if(assVols.size() == 1){
				for(size_t _feI = 0; _feI < f->num_vertices(); ++_feI){ Vertex* vrt = f->vertex(_feI);{
					aaInterface[vrt] = true;
				}};
			}
			else{
				int si = -1;
				for(size_t _vfeI = 0; _vfeI < assVols.size(); ++_vfeI){ Volume* v = assVols[_vfeI];{
					if(si == -1)
						si = sh.get_subset_index(v);
					else if(sh.get_subset_index(v) != si){
						for(size_t _feI = 0; _feI < f->num_vertices(); ++_feI){ Vertex* vrt = f->vertex(_feI);{
							aaInterface[vrt] = true;
						}};
						break;
					}
				}};
			}
		}};

	//	all lower vertices of elements in invalidSub are collapse candidates and
	//	are thus not considered to be interface elements
		Grid::edge_traits::secure_container assEdges;
		for(size_t _vfeI = 0; _vfeI < invalidVols.size(); ++_vfeI){ Volume* vol = invalidVols[_vfeI];{
			grid.associated_elements(assEdges, vol);
			for(size_t _vfeI = 0; _vfeI < assEdges.size(); ++_vfeI){ Edge* e = assEdges[_vfeI];{
				Vertex* v0 = e->vertex(0);
				Vertex* v1 = e->vertex(1);

				vector3 dir;
				VecSubtract(dir, aaPos[v1], aaPos[v0]);
				if((dir.z() > SMALL) && (fabs(dir.x()) < SMALL) && (fabs(dir.y()) < SMALL)){
					aaInterface[v0] = false;
				}
				else if((-dir.z() > SMALL) && (fabs(dir.x()) < SMALL) && (fabs(dir.y()) < SMALL)){
					aaInterface[v1] = false;
				}
			}};
		}};

	//	all unmarked vertices are collapse candidates
		vector<Vertex*> candidates;
		for(Grid::traits<Vertex>::iterator _feI = grid.begin<Vertex>(); _feI != grid.end<Vertex>(); ++_feI){ Vertex* vrt = *_feI;{
			if(!aaInterface[vrt])
				candidates.push_back(vrt);
		}};

	//	merge each candidate with the next upper vertex.
		for(size_t _vfeI = 0; _vfeI < candidates.size(); ++_vfeI){ Vertex* vrt = candidates[_vfeI];{
			grid.associated_elements(assEdges, vrt);
			for(size_t _vfeI = 0; _vfeI < assEdges.size(); ++_vfeI){ Edge* e = assEdges[_vfeI];{
				Vertex* cv = GetConnectedVertex(e, vrt);
				if(cv == vrt)
					continue;

				vector3 dir;
				VecSubtract(dir, aaPos[cv], aaPos[vrt]);
				if((dir.z() > SMALL) && (fabs(dir.x()) < SMALL) && (fabs(dir.y()) < SMALL)){
					CollapseEdge(grid, e, cv);
					break;
				}
			}};
		}};

	//	clean up
		grid.detach_from_vertices(aInterface);

	//	delete invalidSub from the subset handler
		if(invalidSub < sh.num_subsets())
			sh.erase_subset(invalidSub);
	}
	grid.end_marking();
	grid.detach_from_vertices(aHeight);
}


void ExtrudeLayersMixed (
		Grid& grid, 
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler& sh,
		const ANumber* aRelZOut)
{
	UG_COND_THROW(layers.size() < 2, "At least 2 layers are required to perform extrusion!");

	vector<Vertex*> curVrts;	// list of vertices that are considered for extrusion
	vector<Vertex*> nextVrts;	// used to determine the set of vertices that can be extruded
	vector<Face*> curFaces;		// list of faces that are considered for extrusion
	vector<Face*> nextFaces;		// used to determine the set of faces that can be extruded

	grid.reserve<Vertex>(grid.num_vertices() * layers.size());
	grid.reserve<Volume>(grid.num_faces() * layers.size() - 1);
	
	AInt aVrtInd;
	grid.attach_to_vertices(aVrtInd);
	Grid::VertexAttachmentAccessor<AInt> aaVrtInd(grid, aVrtInd);

	Grid::VertexAttachmentAccessor<ANumber> aaRelZOut;
	if(aRelZOut){
		ANumber aRelZ = *aRelZOut;
		if(!grid.has_vertex_attachment(aRelZ))
			grid.attach_to_vertices(aRelZ);
		aaRelZOut.access(grid, aRelZ);
	}

//	we have to determine the vertices that can be projected onto the top of the
//	given layers-stack. Only those will be used during extrusion
//	all considered vertices will be marked.
	const RasterLayers::layer_t& top = layers.top();
	const int topLayerInd = (int)layers.num_layers() - 1;
	for(VertexIterator i = grid.begin<Vertex>(); i != grid.end<Vertex>(); ++i){
		Vertex* v = *i;
		vector2 c(aaPos[v].x(), aaPos[v].y());
		// number val = layers.relative_to_global_height(c, topLayerInd);
		number val = layers.heightfield(topLayerInd).interpolate(c, 1);
		if(val != top.heightfield.no_data_value()){
			aaPos[v].z() = val;
			aaVrtInd[v] = (int)curVrts.size();
			curVrts.push_back(v);
			sh.assign_subset(v, topLayerInd);
		}
		else
			aaVrtInd[v] = -1;
	}

	if(curVrts.size() < 3){
		UG_LOG("Not enough vertices lie in the region of the surface layer\n");
		return;
	}

	if(aaRelZOut.valid()){
		for(size_t ivrt = 0; ivrt < curVrts.size(); ++ivrt)
			aaRelZOut[curVrts[ivrt]] = topLayerInd;
	}

//	make sure that only triangles are used
	if(grid.num<Quadrilateral>() > 0){
		UG_LOG("WARNING: Converting all quadrilaterals to triangles\n");
		Triangulate(grid, grid.begin<Quadrilateral>(), grid.end<Quadrilateral>(), &aaPos);
	}

//	all faces of the initial triangulation that are only connected to marked
//	vertices are considered for extrusion
	for(FaceIterator fi = grid.begin<Face>(); fi != grid.end<Face>(); ++fi){
		Face* f = *fi;
		bool allConsidered = true;
		Face::ConstVertexArray vrts = f->vertices();
		const size_t numVrts = f->num_vertices();
		for(size_t i = 0; i < numVrts; ++i){
			if(aaVrtInd[vrts[i]] == -1){
				allConsidered = false;
				break;
			}
		}
		if(allConsidered)
			curFaces.push_back(f);
	}

	if(curFaces.empty()){
		UG_LOG("Not enough faces lie in the region of the surface layer\n");
		return;
	}


	nextVrts.reserve(curVrts.size());
	nextFaces.reserve(curFaces.size());

	
	for(int ilayer = (int)layers.size() - 1; ilayer > 0; --ilayer){

		nextVrts.clear();
		nextFaces.clear();

		const int nextLayer = ilayer - 1;

	//	trace rays from the current vertices down through the layers until the
	//	next valid entry is found. If none is found, the vertex will either be ignored
	//	from then on (allowForTetsAndPyras == false) or a dummy vertex will be inserted.
		for(size_t icur = 0; icur < curVrts.size(); ++icur){
			Vertex* v = curVrts[icur];
			vector2 c(aaPos[v].x(), aaPos[v].y());
			number curH = aaPos[v].z();
			// number nextH = layers.relative_to_global_height(c, (number) nextLayer);
			number nextH = layers.heightfield(nextLayer).interpolate(c, 1);
			if(curH - nextH < layers.min_height(nextLayer)){
				nextVrts.push_back(v);
			}
			else{
				Vertex* nextVrt = *grid.create<RegularVertex>(v);
				aaVrtInd[nextVrt] = aaVrtInd[v];
				aaPos[nextVrt] = aaPos[v];
				aaPos[nextVrt].z() = nextH;
				sh.assign_subset(nextVrt, nextLayer);
				nextVrts.push_back(nextVrt);

				if(aaRelZOut.valid()){
					aaRelZOut[nextVrt] = nextLayer;
				}
			}
		}

	//	now find the faces which connect those vertices
		for(size_t iface = 0; iface < curFaces.size(); ++iface){
			Face* f = curFaces[iface];
			Face::ConstVertexArray vrts = f->vertices();
			int numNew = 0;
			int firstSame = -1;
			int firstNew = -1;
			TriangleDescriptor nextTriDesc;
			for(int i = 0; i < 3; ++i){
				int vrtInd = aaVrtInd[vrts[i]];
				nextTriDesc.set_vertex(i, nextVrts[vrtInd]);
				if(nextVrts[vrtInd] == curVrts[vrtInd]){
					if(firstSame == -1)
						firstSame = i;
				}
				else{
					++numNew;
					if(firstNew == -1)
						firstNew = i;
				}
			}

			Face* nextFace = f;
			if(numNew > 0){
				nextFace = *grid.create<Triangle>(nextTriDesc, f);
				sh.assign_subset (nextFace, nextLayer);
			}
			nextFaces.push_back(nextFace);

			Volume* newVol = nullptr;

			switch(numNew){
				case 0:// no new volume to be built.
					break;
				case 1: {// build a tetrahedron
					TetrahedronDescriptor 	tetDesc (nextTriDesc.vertex(0),
					                                 nextTriDesc.vertex(1),
													 nextTriDesc.vertex(2),
													 vrts[firstNew]);
					newVol = *grid.create<Tetrahedron>(tetDesc);
				} break;

				case 2: {// build a pyramid
					int i0 = (firstSame + 1) % 3;
					int i1 = (firstSame + 2) % 3;
					PyramidDescriptor	pyrDesc (nextTriDesc.vertex(i0),
					                             nextTriDesc.vertex(i1),
					                             vrts[i1],
					                             vrts[i0],
					                             vrts[firstSame]);
					newVol = *grid.create<Pyramid>(pyrDesc);
				} break;

				case 3: {// build a prism
					PrismDescriptor		prismDesc (nextTriDesc.vertex(0),
												   nextTriDesc.vertex(1),
												   nextTriDesc.vertex(2),
												   vrts[0],
												   vrts[1],
												   vrts[2]);
					newVol = *grid.create<Prism>(prismDesc);
				} break;
			}

			if(newVol)
				sh.assign_subset(newVol, nextLayer);
		}

	//	swap containers and perform extrusion
		curVrts.swap(nextVrts);
		curFaces.swap(nextFaces);
	}
	
//	clean up
	grid.detach_from_vertices(aVrtInd);
}

///	projects the given (surface-) grid to the specified raster
void ProjectToLayer(
		Grid& grid,
		const RasterLayers& layers,
		int layerIndex,
		Grid::VertexAttachmentAccessor<AVector3> aaPos)
{
	UG_COND_THROW(layerIndex < 0 || layerIndex >= (int)layers.size(),
				  "Bad layerIndex in ProjectToLayer");

	const RasterLayers::layer_t& layer = layers.layer(layerIndex);
	for(VertexIterator i = grid.begin<Vertex>(); i != grid.end<Vertex>(); ++i){
		Vertex* v = *i;
		vector2 c(aaPos[v].x(), aaPos[v].y());
		number val = layers.relative_to_global_height(c, layerIndex);

		if(val != layer.heightfield.no_data_value()){
			aaPos[v].z() = val;
		}
	}
}


void SnapToHorizontalRaster(
		Grid& grid,
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos)
{
	UG_COND_THROW(layers.empty(), "Can't snap to empty raster!");

	const Heightfield& field = layers.top().heightfield;

	for(VertexIterator ivrt = grid.vertices_begin();
		ivrt != grid.vertices_end(); ++ivrt)
	{
		Vertex* vrt = *ivrt;
		vector3& v = aaPos[vrt];
		pair<int, int> ci = field.coordinate_to_index(v.x(), v.y());
		vector2 nv = field.index_to_coordinate(ci.first, ci.second);
		v.x() = nv.x();
		v.y() = nv.y();
	}
}

}//	end of namespace
