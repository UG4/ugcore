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

//#include <boost/function.hpp>
#include <stack>
#include <vector>
#include "expand_layers.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/grid/grid_util.h"
//#include "lib_grid/util/simple_algebra/least_squares_solver.h"

using namespace std;

namespace ug{


///	This class can be used in Element callbacks.
/**	Returns true, if the attachmed value in the given element matches a predefined value.*/
template <class TElem, class TAttachmentAccessor>
class AttachmentUnequal{
	public:
		AttachmentUnequal(const TAttachmentAccessor& aa,
						const typename TAttachmentAccessor::ValueType& val) :
			m_aa(aa), m_val(val) {}

		bool operator() (TElem* e) {return m_aa[e] != m_val;}

	private:
		TAttachmentAccessor						m_aa;
		typename TAttachmentAccessor::ValueType	m_val;

};

///	returns true if the vertex lies on a surface.
/**	This method uses Grid::mark.
 *
 *	This method tries to find a closed surface around the vertex.
 *	Note that it returns true, even if the vertex lies on a
 *	boundary edge at the same time (this can happen if surfaces
 *	intersect each other).
 *
 *	Requires the option FACEOPT_AUTOGENERATE_EDGES.
 */
/*
static bool VertexLiesOnSurface(Grid& grid, Vertex* vrt,
						 CB_ConsiderFace funcIsSurfFace)
{
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		UG_LOG("WARNING in VertexLiesOnSurface: autoenabling FACEOPT_AUTOGENERATE_EDGES.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

	if(grid.associated_edges_begin(vrt) == grid.associated_edges_end(vrt))
	{
		return  false;
	}

	grid.begin_marking();

	stack<Edge*> stk;
	vector<Edge*> edges;
	vector<Face*> faces;

//	collect associated faces of vrt, which lie on the surface
	for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt);
		iter != grid.associated_faces_end(vrt); ++iter)
	{
		if(funcIsSurfFace(*iter))
			faces.push_back(*iter);
	}

//	start with an associated edge of vrt
	stk.push(*grid.associated_edges_begin(vrt));
	grid.mark(stk.top());

	while(!stk.empty()){
		Edge* e = stk.top();

	//	find an unmarked associated face
		Face* f = NULL;
		for(size_t i = 0; i < faces.size(); ++i){
			if(!grid.is_marked(faces[i])){
				if(FaceContains(faces[i], e)){
					f = faces[i];
					break;
				}
			}
		}

		if(!f){
			stk.pop();
			continue;
		}

		grid.mark(f);

	//	collect associated edges of f
		CollectEdges(edges, grid, f);

	//	find the edge that is not e an which is connected to vrt
		for(size_t i = 0; i < edges.size(); ++i){
			Edge* ne = edges[i];
			if(ne != e){
				if(EdgeContains(ne, vrt)){
				//	if the edge is marked, then we found a surface
					if(grid.is_marked(ne)){
						grid.end_marking();
						return true;
					}
					else{
						grid.mark(ne);
						stk.push(ne);
						break;
					}
				}
			}
		}
	}

	grid.end_marking();
	return false;
}
*/
///	calculates the normal of the crease vertex vrt on the side of f
/**
 *	This algorithm uses grid::mark.
 *	f has to contain vrt. vrt has to have at least two associated
 *	crease edges.
 *
 *	Note that the resulting normal is not normalized. This is important,
 *	since rounding errors could else lead to problems with normals which
 *	have a length of nearly 0.
 *
 *	This algorithm requires the option FACEOPT_AUTOGENERATE_EDGES.
 *	The option is automatically enabled if required.
 */
template <class TAAPosVRT>
typename TAAPosVRT::ValueType
CalculateCreaseNormal(Grid& grid, Face* f, Vertex* vrt,
						Grid::edge_traits::callback funcIsCreaseEdge,
						TAAPosVRT& aaPos)
{
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		UG_LOG("WARNING in CalculateCreaseNormal: grid option FACEOPT_AUTOGENERATE_EDGES autoenabled.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

	typedef typename TAAPosVRT::ValueType vector_t;
	vector_t n;
	VecSet(n, 0);

	grid.begin_marking();

//	we'll use a stack to find all faces on this side of the crease.
//	each face in the stack has to be marked.
	stack<Face*> stk;
	stk.push(f);
	grid.mark(f);

//	objects for temporary results
	vector<Edge*> edges;
	vector<Face*> faces;

//	we'll loop while there are faces in the stack
	while(!stk.empty()){
		Face* face = stk.top();
		stk.pop();

	//	the center might be required later on
		vector_t center = CalculateCenter(face, aaPos);

	//	iterate over the edges of face
		CollectEdges(edges, grid, face);
		for(size_t i_e = 0; i_e < edges.size(); ++i_e){
			Edge* e = edges[i_e];
			if(EdgeContains(e, vrt)){
			//	check whether e is a crease
				if(funcIsCreaseEdge(e)){
				//	we have to add the edges normal to n.
				//	to make sure that the algorithm works for manifolds too,
				//	we'll perform a more complicated calculation.

				//	project the center onto the edge
					vector_t p;
					DropAPerpendicular(p, center, aaPos[e->vertex(0)], aaPos[e->vertex(1)]);

				//	vector from projection to center is the unnormalized normal
					vector_t tmpN;
					VecSubtract(tmpN, center, p);
					VecNormalize(tmpN, tmpN);
					VecAdd(n, n, tmpN);
				}
				else{
				//	since the edge is not a crease edge, we have to add associated unmarked
				//	faces to the stack
					CollectFaces(faces, grid, e);
					for(size_t i = 0; i < faces.size(); ++i){
						if(!grid.is_marked(faces[i])){
							grid.mark(faces[i]);
							stk.push(faces[i]);
						}
					}
				}
			}
		}
	}

	grid.end_marking();

	//VecNormalize(n, n);
	return n;
}

///	calculates the normal of the crease vertex vrt on the side of vol
/**
 *	This algorithm uses grid::mark.
 *	vol has to contain vrt. vrt has to have at least two associated
 *	crease faces.
 *
 *	Note that the resulting normal is not normalized. This is important,
 *	since rounding errors could else lead to problems with normals which
 *	have a length of nearly 0.
 *
 *	This algorithm requires the option VOLOPT_AUTOGENERATE_FACES.
 *	The option is automatically enabled if required.
 */
template <class TAAPosVRT>
typename TAAPosVRT::ValueType
CalculateCreaseNormal(Grid& grid, Volume* vol, Vertex* vrt,
						Grid::face_traits::callback funcIsCreaseFace,
						TAAPosVRT& aaPos)
{
	if(!grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES)){
		UG_LOG("WARNING in CalculateCreaseNormal: grid option VOLOPT_AUTOGENERATE_FACES autoenabled.\n");
		grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}

	typedef typename TAAPosVRT::ValueType vector_t;
	vector_t n;
	VecSet(n, 0);

	grid.begin_marking();

//	we'll use a stack to find all volumes on this side of the crease.
//	each volume in the stack has to be marked.
	stack<Volume*> stk;
	stk.push(vol);
	grid.mark(vol);

//	objects for temporary results
	FaceDescriptor fd;
	vector<Volume*> vols;

//	we'll loop while there are faces in the stack
	while(!stk.empty()){
		Volume* curVol = stk.top();
		stk.pop();

	//	iterate over all sides of curVol
		for(size_t i_side = 0; i_side < curVol->num_sides(); ++i_side){
			Face* f = grid.get_side(curVol, i_side);
			if(f){
			//	only proceed if the face contains vrt
				if(FaceContains(f, vrt)){
				//	check whether f is a crease
					if(funcIsCreaseFace(f)){
					//	calculate the normal of the side
						vector_t tmpN;
						curVol->face_desc(i_side, fd);
						CalculateNormal(tmpN, &fd, aaPos);

					//	the normal points away from the volume, so
					//	we have to subtract it from n
						VecSubtract(n, n, tmpN);
					}
					else{
					//	since the face is not a crease face, we have to add associated unmarked
					//	volumes to the stack
						CollectVolumes(vols, grid, f);
						for(size_t i = 0; i < vols.size(); ++i){
							if(!grid.is_marked(vols[i])){
								grid.mark(vols[i]);
								stk.push(vols[i]);
							}
						}
					}
				}
			}
		}
	}

	grid.end_marking();

	//VecNormalize(n, n);
	return n;
}

/**
 * This algorithm indirectly uses Grid::mark.
 *
 * 1 dimensional fractures specified in fracInfos are expanded to 2 dimensional subsets.
 * the resulting fractures will then consist of 2 layers of quadrilaterals. On the
 * boundaries triangles are inserted.
 *
 * Through expandFracBoundaries you can tell the algorithm whether inner fracture
 * boundaries shall be expanded. Note that this means that an additional node is
 * introduced at each inner fracture boundary vertex and that the associated
 * fracture elements are connected at two sides.
 * Note that fractures are always expanded at boundaries which lie on the geometries
 * boundary.
 *
 *	This algorithm requires the option FACEOPT_AUTOGENERATE_EDGES.
 *	The option is automatically enabled if required.
 */
bool ExpandFractures2d(Grid& grid, SubsetHandler& sh, const vector<FractureInfo>& fracInfos,
						bool expandInnerFracBnds, bool expandOuterFracBnds)
{
//	access position attachment
	if(!grid.has_vertex_attachment(aPosition)){
		UG_LOG("Error in ExpandFractures: Missing position attachment");
		return false;
	}
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		UG_LOG("WARNING in CalculateCreaseNormal: grid option FACEOPT_AUTOGENERATE_EDGES autoenabled.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

//	objects for temporary results
	FaceDescriptor fd;
	vector<Edge*> edges; // used for temporary results.
	vector<Face*> faces; // used for temporary results.

//	vectors that allow to access fracture properties by subset index
	vector<FractureInfo> fracInfosBySubset(sh.num_subsets(), FractureInfo(-1, -1, 0));
	for(size_t i = 0; i < fracInfos.size(); ++i){
		if(fracInfos[i].subsetIndex >= sh.num_subsets()){
			throw(UGError("Bad subsetIndex in given fracInfos."));
		}

		fracInfosBySubset[fracInfos[i].subsetIndex] = fracInfos[i];
	}

////////////////////////////////
//	Collect surrounding faces of all fractures in a selector
//	and select edges and vertices too.
	Selector sel(grid);
	sel.enable_autoselection(false);
	sel.enable_selection_inheritance(false);

	AInt aAdjMarker;	// used to mark how many adjacent fractures a vertex has.
						// 0: no frac, 1: frac-boundary, >1: inner frac vertex
	grid.attach_to_vertices_dv(aAdjMarker, 0);
	Grid::VertexAttachmentAccessor<AInt> aaMarkVRT(grid, aAdjMarker);
	grid.attach_to_edges_dv(aAdjMarker, 0);
	Grid::EdgeAttachmentAccessor<AInt> aaMarkEDGE(grid, aAdjMarker);

//	iterate over the given fracture infos and select all fracture edges
//	and fracture vertices.
	for(size_t i_fi = 0; i_fi < fracInfos.size(); ++i_fi){
		int fracInd = fracInfos[i_fi].subsetIndex;
		for(EdgeIterator iter = sh.begin<Edge>(fracInd);
			iter != sh.end<Edge>(fracInd); ++iter)
		{
		//	mark edge and vertices
			sel.select(*iter);
			aaMarkEDGE[*iter] = 1;

		//	select associated vertices
			for(size_t i = 0; i < 2; ++i){
				Vertex* v = (*iter)->vertex(i);
				sel.select(v);

			//	if fracture boundaries are expanded, we'll regard all fracture vertices
			//	as inner vertices
				if(expandInnerFracBnds){
					if(!expandOuterFracBnds){
						if(IsBoundaryVertex2D(grid, v))
							aaMarkVRT[v]++;
						else
							aaMarkVRT[v] = 2;
					}
					else
						aaMarkVRT[v] = 2;
				}
				else
					aaMarkVRT[v]++;
			}
		}
	}

//	Make sure that selected vertices that lie on the boundary of the geometry
//	are treated as inner fracture vertices.
//	This is only required if frac-boundaries are not expanded anyways.
	if(expandOuterFracBnds && !expandInnerFracBnds){
		for(VertexIterator iter = sel.vertices_begin();
			iter != sel.vertices_end(); ++iter)
		{
			Vertex* v = *iter;
			if(aaMarkVRT[v] == 1){
				if(IsBoundaryVertex2D(grid, v))
					aaMarkVRT[v] = 2;
			}
		}
	}

//	Select all edges and faces which are connected to inner fracture vertices
	for(VertexIterator iter = sel.begin<Vertex>();
		iter != sel.end<Vertex>(); ++iter)
	{
		if(aaMarkVRT[*iter] > 1){
			sel.select(grid.associated_edges_begin(*iter),
						grid.associated_edges_end(*iter));
			sel.select(grid.associated_faces_begin(*iter),
						grid.associated_faces_end(*iter));
		}
	}

////////////////////////////////
//	create new vertices

//	we have to associate a vector of vertices with each node in the fracture.
//	since an empty vector is quite small, we can associate one with each vertex in
//	the whole grid. This could be optimized if required, by using subset attachments.
	typedef Attachment<vector<Vertex*> > AVrtVec;
	AVrtVec aVrtVec;
	grid.attach_to_vertices(aVrtVec);
	Grid::VertexAttachmentAccessor<AVrtVec> aaVrtVecVRT(grid, aVrtVec);

//	we also have to associate a vector of vertices for each face adjacent to the frac.
//	it will store the a second set of vertices. An entry contains the new vertex, if the
//	corresponding vertex is an inner fracture vertex, and NULL if not.
	grid.attach_to_faces(aVrtVec);
	Grid::FaceAttachmentAccessor<AVrtVec> aaVrtVecFACE(grid, aVrtVec);

//	a callback that returns true if the edge is a fracture edge
	AttachmentUnequal<Edge, Grid::EdgeAttachmentAccessor<AInt> > isFracEdge(aaMarkEDGE, 0);

//	iterate over all surrounding faces and create new vertices.
	for(FaceIterator iter_sf = sel.faces_begin(); iter_sf != sel.faces_end(); ++iter_sf)
	{
		Face* sf = *iter_sf;

		vector<Vertex*>& newVrts = aaVrtVecFACE[sf];
		newVrts.resize(sf->num_vertices());

	//	check for each vertex whether it lies in the fracture
	//	(aaMarkVRT > 1 in this case)
	//	if so, we have to copy or create a vertex from/in aaVrtVec[vrt] which is
	//	associated with the crease normal on the side of sf.
		for(size_t i_vrt = 0; i_vrt < sf->num_vertices(); ++i_vrt)
		{
			newVrts[i_vrt] = NULL;
			Vertex* vrt = sf->vertex(i_vrt);
			if(aaMarkVRT[vrt] > 1){
			//	calculate the normal on this side of the frac
				vector3 n = CalculateCreaseNormal(grid, sf, vrt, isFracEdge, aaPos);
				//UG_LOG("calculated crease normal: " << n << endl);

			//	check whether aaVrtVecs already contains a vertex associated with n.
			//	the normal of new vrts is stored in their position attachment
				vector<Vertex*>& vrtVec = aaVrtVecVRT[vrt];
				for(size_t i = 0; i < vrtVec.size(); ++i){
					//UG_LOG("comparing to: " << aaPos[vrtVec[i]] << endl);
					if(VecDistanceSq(aaPos[vrtVec[i]], n) < SMALL){
					//	got one
						newVrts[i_vrt] = vrtVec[i];
						break;
					}
				}

			//	if we didn't find one then create and associate one.
			//	store the normal in the position attachment of the new vertex
				if(!newVrts[i_vrt]){
					newVrts[i_vrt] = *grid.create<RegularVertex>();
					aaPos[newVrts[i_vrt]] = n;
					aaVrtVecVRT[vrt].push_back(newVrts[i_vrt]);
				}
			}
		}
	}

//	assign the new positions
	for(VertexIterator iter = sel.vertices_begin();
		iter != sel.vertices_end(); ++iter)
	{
		Vertex* vrt = *iter;

	//	calculate the width as the maximum of associated fracture widths
		CollectEdges(edges, grid, vrt);

		number width = 0;
		for(size_t i = 0; i < edges.size(); ++i){
			if(aaMarkEDGE[edges[i]])
				width = max<number>(width, fracInfosBySubset.at(sh.get_subset_index(edges[i])).width);
		}

	//	iterate over associated vertices
		vector<Vertex*>& vrtVec = aaVrtVecVRT[vrt];

	//	note that the position attachment of new vertices holds their normal.
		for(size_t i = 0; i < vrtVec.size(); ++i){
			Vertex* nVrt = vrtVec[i];
			if(width > 0){
				vector3 n = aaPos[nVrt];
				if(VecLengthSq(n) > SMALL)
					VecNormalize(n, n);

				VecScale(n, n, width / 2.);

				UG_LOG("n: " << n << endl);

			//	n now holds the offset for nVrt relative to vrt.
			//	if width is higher than 0, we'll have to adjust the offset at
			//	boundary vertices.
				if(IsBoundaryVertex2D(grid, vrt)){
				//	First determine the normal pointing outwards
					vector3 nOut;
					CalculateBoundaryVertexNormal2D(nOut, grid, vrt, aaPos);

				//	flip it by 90 degrees
					number tmp = nOut.x();
					nOut.x() = -nOut.y();
					nOut.y() = tmp;

					UG_LOG("nOut: " << nOut << endl);

				//	now project the offset onto this vector
					VecScale(nOut, nOut, VecDot(nOut, n));

				//	and now scale the new offset so that we receive the final offset.
					number dot = VecDot(n, nOut);
					if(dot > SMALL)
						VecScale(n, nOut, VecLengthSq(n) / dot);
				}

				UG_LOG("nFinal: " << n << endl);
				VecAdd(aaPos[nVrt], n, aaPos[vrt]);
				UG_LOG("\n");
			}
			else
				aaPos[nVrt] = aaPos[vrt];
		}

	//	the current position is only a guess. Especially vertices where
	//	fractures cross, this is not yet optimal.
	//todo: create an iterative spring system to find the new position.
	}

////////////////////////////////
//	create new elements

//	first we create new edges from selected ones which are connected to
//	inner vertices. This allows to preserve old subsets.
//	Since we have to make sure that we use the right vertices,
//	we have to iterate over the selected faces and perform all actions on the edges
//	of those faces.
	for(FaceIterator iter_sf = sel.faces_begin(); iter_sf != sel.faces_end(); ++iter_sf)
	{
		Face* sf = *iter_sf;
	//	check for each edge whether it has to be copied.
		for(size_t i_edge = 0; i_edge < sf->num_edges(); ++i_edge){
			Edge* e = grid.get_edge(sf, i_edge);
			if(sel.is_selected(e)){
			//	check the associated vertices through the volumes aaVrtVecVol attachment.
			//	If at least one has an associated new vertex and if no edge between the
			//	new vertices already exists, we'll create the new edge.
				size_t ind0 = i_edge;
				size_t ind1 = (i_edge + 1) % sf->num_edges();

				Vertex* nv0 = (aaVrtVecFACE[sf])[ind0];
				Vertex* nv1 = (aaVrtVecFACE[sf])[ind1];

				if(nv0 || nv1){
				//	if one vertex has no associated new one, then we use the vertex itself
					if(!nv0)
						nv0 = sf->vertex(ind0);
					if(!nv1)
						nv1 = sf->vertex(ind1);

				//	create the new edge if it not already exists.
					if(!grid.get_edge(nv0, nv1))
						grid.create_by_cloning(e, EdgeDescriptor(nv0, nv1), e);
				}
			}
		}
	}

//	iterate over all surrounding faces and create new vertices.
//	Since faces are replaced on the fly, we have to take care with the iterator.
	for(FaceIterator iter_sf = sel.faces_begin(); iter_sf != sel.faces_end();)
	{
		Face* sf = *iter_sf;
		++iter_sf;

		vector<Vertex*> newVrts = aaVrtVecFACE[sf];

	//	all new vertices have been assigned to newVrts.
	//	Note that if newVrts[i] == NULL, then we have to take the
	//	old vertex sf->vertex(i).
	//	now expand the fracture edges of sf to faces.
		for(size_t i_vrt = 0; i_vrt < sf->num_vertices(); ++i_vrt){
			size_t iv1 = i_vrt;
			size_t iv2 = (i_vrt + 1) % sf->num_vertices();
			Edge* tEdge = grid.get_edge(sf->vertex(iv1), sf->vertex(iv2));
			if(tEdge){
				if(aaMarkEDGE[tEdge]){
					Face* expFace = NULL;
					if(newVrts[iv1] && newVrts[iv2]){
					//	create a new quadrilateral
						expFace = *grid.create<Quadrilateral>(
									QuadrilateralDescriptor(sf->vertex(iv1), sf->vertex(iv2),
															newVrts[iv2], newVrts[iv1]));
					}
					else if(newVrts[iv1]){
					//	create a new triangle
						expFace = *grid.create<Triangle>(
									TriangleDescriptor(sf->vertex(iv1), sf->vertex(iv2), newVrts[iv1]));
					}
					else if(newVrts[iv2]){
					//	create a new triangle
						expFace = *grid.create<Triangle>(
									TriangleDescriptor(sf->vertex(iv1), sf->vertex(iv2), newVrts[iv2]));
					}
					else{
					//	this code-block should never be entered. If it is entered then
					//	we selected the wrong faces. This is would be a BUG!!!
					//	remove the temporary attachments and throw an error
						grid.detach_from_vertices(aVrtVec);
						grid.detach_from_faces(aVrtVec);
						grid.detach_from_vertices(aAdjMarker);
						grid.detach_from_edges(aAdjMarker);
						throw(UGError("Implementation error in ExpandFractures2d."));
					}

					sh.assign_subset(expFace, fracInfosBySubset.at(sh.get_subset_index(tEdge)).newSubsetIndex);
				}
			}
		}


	//	now set up a new face descriptor and replace the face.
		if(fd.num_vertices() != sf->num_vertices())
			fd.set_num_vertices(sf->num_vertices());

		for(size_t i_vrt = 0; i_vrt < sf->num_vertices(); ++i_vrt){
			if(newVrts[i_vrt])
				fd.set_vertex(i_vrt, newVrts[i_vrt]);
			else
				fd.set_vertex(i_vrt, sf->vertex(i_vrt));
		}
		grid.create_by_cloning(sf, fd, sf);
		grid.erase(sf);
	}

//	we have to clean up unused edges.
//	All selected edges with mark 0 have to be deleted.
	for(EdgeIterator iter = sel.edges_begin(); iter != sel.edges_end();)
	{
	//	be careful with the iterator
		Edge* e = *iter;
		++iter;

		if(!aaMarkEDGE[e])
			grid.erase(e);
	}

//	remove the temporary attachments
	grid.detach_from_vertices(aVrtVec);
	grid.detach_from_faces(aVrtVec);
	grid.detach_from_vertices(aAdjMarker);
	grid.detach_from_edges(aAdjMarker);

	return true;
}

/**	Selects all involved geometic objects and assigns marks to them.
 * If required, som edges may be split, so that we always operate
 * on a fully expandable fracture.
 *
 * Make sure that selection_inheritance is enabled and that
 * strict_inheritance is disabled in sel.*/
static void DistributeExpansionMarks3D(Grid& grid, SubsetHandler& sh, Selector& sel,
									const vector<FractureInfo>& fracInfos,
									bool expandInnerFracBnds,
									bool expandOuterFracBnds,
									Grid::VertexAttachmentAccessor<AInt>& aaMarkVRT,
									Grid::EdgeAttachmentAccessor<AInt>& aaMarkEDGE,
									Grid::FaceAttachmentAccessor<AInt>& aaMarkFACE)
{
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	objects for temporary results
	//VolumeDescriptor vd;
	vector<Edge*> edges; // used for temporary results.
	vector<Face*> faces; // used for temporary results.
	vector<Volume*> vols; // used for temporary results.

//	iterate over the given fracture infos and select all fracture faces,
//	fracture edges and fracture vertices
	for(size_t i_fi = 0; i_fi < fracInfos.size(); ++i_fi){
		int fracInd = fracInfos[i_fi].subsetIndex;
		for(FaceIterator iter = sh.begin<Face>(fracInd);
			iter != sh.end<Face>(fracInd); ++iter)
		{
			Face* f = *iter;

		//	mark face and vertices
			sel.select(f);

		//	collect associated volumes of the vertices and add them to surroundingFaces
			for(size_t i = 0; i < f->num_vertices(); ++i){
				Vertex* v = f->vertex(i);
				sel.select(v);
//				CollectVolumes(vols, grid, v);
//				sel.select(vols.begin(), vols.end());
			}

		//	collect associated edges of the faces and
		//	increase their adjacency counter too, since this helps
		//	to identify whether they lie on the selection boundary.
			CollectEdges(edges, grid, f);

			for(size_t i = 0; i < edges.size(); ++i){
				aaMarkEDGE[edges[i]]++;
				sel.select(edges[i]);
			}
		}
	}

//	all edges that lie on the geometries boundary have to be regarded as inner edges
	if(expandOuterFracBnds){
		for(EdgeIterator iter = sel.edges_begin(); iter != sel.edges_end(); ++iter)
		{
			Edge* e = *iter;
			if(aaMarkEDGE[e] == 1){
				if(IsBoundaryEdge3D(grid, e))
					aaMarkEDGE[e] = 2;
			}
		}
	}

//	iterate over selected vertices and check the adjacency status of associated edges.
//	the vertex can either be a boundary vertex (1) or a surface vertex(2) or both (3).
//	Note that vertices lying on the geometries boundary are regarded as surface vertices.
//	At this point it is important to mark all vertices that lie on any boundary as
//	a boundary vertex, since we have to split inner edges connecting boundary-vertices
//	in the next step.
//	Note that this is important even for the degenerate case, since otherwise identical
//	elements may be created.
//todo:	Currently only surface or boundary vertices are regarded. The mix is
//		treated as a boundary vertex for the non-degenerated case.
//		This is because regarding the mix of both would result in problematic element shapes.
//		For the degenerated case we regard the mix as an inner vertex.
	for(VertexIterator iter = sel.vertices_begin();
		iter != sel.vertices_end(); ++iter)
	{
		CollectEdges(edges, grid, *iter);
		aaMarkVRT[*iter] = 2;
		for(size_t i = 0; i < edges.size(); ++i){
			if(aaMarkEDGE[edges[i]] == 1){
//				if(VertexLiesOnSurface(grid, *iter, IsSelected(sel))){
//					UG_LOG("found mixed\n");
//					aaMarkVRT[*iter] = 3;
//				}
				aaMarkVRT[*iter] = 1;
				break;
			}
		}
	}

// todo:	Quadrilaterals which have more than 2 boundary vertices or have two boundary
//			vertices, which are not adjacent to each other, have to be transformed to
//			triangles.


//	now make sure that no inner edge is associated with two
//	boundary vertices (referring to the selection)
	edges.clear();
	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
	{
		Edge* e = *iter;
		if(aaMarkVRT[e->vertex(0)] != 2 &&
			aaMarkVRT[e->vertex(1)] != 2 &&
			aaMarkEDGE[e] > 1)
		{
			edges.push_back(e);
		}
	}

	for(size_t i = 0; i < edges.size(); ++i){
		vector3 center = CalculateCenter(edges[i], aaPos);
		RegularVertex* v =	SplitEdge<RegularVertex>(grid, edges[i], false);
		aaPos[v] = center;
		aaMarkVRT[v] = 2;
		sel.select(v);
	//	assign adjacency values for associated selected edges (2 to each)
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(v);
			iter != grid.associated_edges_end(v); ++iter)
		{
			if(sel.is_selected(*iter))
				aaMarkEDGE[*iter] = 2;
		}
	}

//	if fracture boundaries shall be extended, then we have to regard all
//	boundary vertices as inner vertices
	if(expandInnerFracBnds && expandOuterFracBnds){
		for(VertexIterator iter = sel.vertices_begin(); iter != sel.vertices_end(); ++iter)
		{
			Vertex* v = *iter;
			if(aaMarkVRT[v] == 1)
				aaMarkVRT[v] = 2;
		}
	}
	else if(expandInnerFracBnds){
		for(VertexIterator iter = sel.vertices_begin(); iter != sel.vertices_end(); ++iter)
		{
			Vertex* v = *iter;
			if(!IsBoundaryVertex3D(grid, v)){
				if(aaMarkVRT[v] == 1)
					aaMarkVRT[v] = 2;
			}
		}
	}
	else if(expandOuterFracBnds){
		for(VertexIterator iter = sel.vertices_begin(); iter != sel.vertices_end(); ++iter)
		{
			Vertex* v = *iter;
			if(IsBoundaryVertex3D(grid, v)){
			//	get state from marked associated boundary edges
				aaMarkVRT[v] = 0;
				CollectAssociated(edges, grid, v);

				for(size_t i = 0; i < edges.size(); ++i){
					Edge* e = edges[i];
					if(aaMarkEDGE[e] > 1){
						if(IsBoundaryEdge3D(grid, e))
							aaMarkVRT[v]++;
					}
				}
			}
		}
	}


//	All fracture faces shall be marked with 1. We do this here, since new
//	faces may have been created during the edge-splits.
	for(FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter)
		aaMarkFACE[*iter] = 1;

//	select all non-fracture edges, faces and volumes, which are connected to an
//	inner fracture vertex.
	for(VertexIterator iter = sel.begin<Vertex>();
		iter != sel.end<Vertex>(); ++iter)
	{
		Vertex* vrt = *iter;
		if(aaMarkVRT[vrt] > 1){
		//	select all associated edges, faces and volumes
			sel.select(grid.associated_edges_begin(vrt),
						grid.associated_edges_end(vrt));
			sel.select(grid.associated_faces_begin(vrt),
						grid.associated_faces_end(vrt));
			sel.select(grid.associated_volumes_begin(vrt),
						grid.associated_volumes_end(vrt));
		}
	}
}


/**
 * This algorithm indirectly uses Grid::mark.
 *
 * 2 dimensional fractures specified in fracInfos are expanded to 3 dimensional subsets.
 * the resulting fractures will then consist of 2 layers of hexahedrons. On the
 * boundaries tetrahedrons, prisms and pyramids are inserted.
 *
 * Through expandFracBoundaries you can tell the algorithm whether inner fracture
 * boundaries shall be expanded. Note that this means that an additional node is
 * introduced at each inner fracture boundary vertex and that the associated
 * fracture elements are connected at two sides.
 * Note that fractures are always expanded at boundaries which lie on the geometries
 * boundary.
 *
 *	This algorithm requires the option FACEOPT_AUTOGENERATE_EDGES.
 *	The option is automatically enabled if required.
 *
 *	This algorithm requires the option VOLOPT_AUTOGENERATE_FACES.
 *	The option is automatically enabled if required.
 */
bool ExpandFractures3d(Grid& grid, SubsetHandler& sh, const vector<FractureInfo>& fracInfos,
						bool expandInnerFracBnds, bool expandOuterFracBnds)
{
//	access position attachment
	if(!grid.has_vertex_attachment(aPosition)){
		UG_LOG("Error in ExpandFractures: Missing position attachment");
		return false;
	}
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

//	make sure that the required options are enabled.
	if(!grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES)){
		UG_LOG("WARNING in CalculateCreaseNormal: grid option VOLOPT_AUTOGENERATE_FACES autoenabled.\n");
		grid.enable_options(VOLOPT_AUTOGENERATE_FACES);
	}

	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		UG_LOG("WARNING in CalculateCreaseNormal: grid option FACEOPT_AUTOGENERATE_EDGES autoenabled.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

//	objects for temporary results
	FaceDescriptor fd;
	VolumeDescriptor vd;
	vector<Edge*> edges; // used for temporary results.
	vector<Face*> faces; // used for temporary results.
	vector<Volume*> vols; // used for temporary results.

////////////////////////////////
//	Collect surrounding volumes, faces and edges of all fractures in a selector
//	and select fracture faces, edges and vertices too.
	Selector sel(grid);
	//Selector& sel = tmpSel;
	sel.enable_autoselection(false);
	sel.enable_selection_inheritance(true);	//required for DistributeExpansionMarks3D. disabled later on.
	sel.enable_strict_inheritance(false);

	AInt aAdjMarker;	// used to mark where in a fracture an edge or a vertex lies.
						// 0: no frac, 1: frac-boundary, >1: inner frac
	grid.attach_to_vertices_dv(aAdjMarker, 0);
	Grid::VertexAttachmentAccessor<AInt> aaMarkVRT(grid, aAdjMarker);
	grid.attach_to_edges_dv(aAdjMarker, 0);
	Grid::EdgeAttachmentAccessor<AInt> aaMarkEDGE(grid, aAdjMarker);
	grid.attach_to_faces_dv(aAdjMarker, 0);
	Grid::FaceAttachmentAccessor<AInt> aaMarkFACE(grid, aAdjMarker);

//	Distribute marks and select elements.

	DistributeExpansionMarks3D(grid, sh, sel, fracInfos, expandInnerFracBnds,
						expandOuterFracBnds, aaMarkVRT, aaMarkEDGE, aaMarkFACE);

//	We'll now store all fracture faces in a vector.
//	They will help to adjust positions of new vertices later on.
	std::vector<Face*> originalFractureFaces;
	for(FaceIterator iter = sel.begin<Face>(); iter != sel.end<Face>(); ++iter){
		if(aaMarkFACE[*iter] == 1)
			originalFractureFaces.push_back(*iter);
	}

//	vectors that allow to access fracture properties by subset index
	vector<FractureInfo> fracInfosBySubset(sh.num_subsets(), FractureInfo(-1, -1, 0));
	for(size_t i = 0; i < fracInfos.size(); ++i){
		if(fracInfos[i].subsetIndex >= sh.num_subsets()){
			throw(UGError("Bad subsetIndex in given fracInfos."));
		}

		fracInfosBySubset[fracInfos[i].subsetIndex] = fracInfos[i];
	}

//	disable selection inheritance to avoid infinite recursion.
	sel.enable_selection_inheritance(false);
//	clear buffers for later use
	edges.clear();

////////////////////////////////
//	create new vertices
//	we have to associate a vector of vertices with each node in the fracture.
//	since an empty vector is quite small, we can associate one with each vertex in
//	the whole grid. This could be optimized if required, by using subset attachments.
	typedef Attachment<vector<Vertex*> > AVrtVec;
	AVrtVec aVrtVec;
	grid.attach_to_vertices(aVrtVec);
	Grid::VertexAttachmentAccessor<AVrtVec> aaVrtVecVRT(grid, aVrtVec);

//	we also have to associate a vector of vertices for each volume adjacent to the frac.
//	it will store a second set of vertices. An entry contains the new vertex, if the
//	corresponding vertex is an inner fracture vertex, and NULL if not.
	grid.attach_to_volumes(aVrtVec);
	Grid::VolumeAttachmentAccessor<AVrtVec> aaVrtVecVOL(grid, aVrtVec);

//	a callback which tells whether a face is inside the fracture or not
	AttachmentUnequal<Face, Grid::FaceAttachmentAccessor<AInt> > faceIsInFrac(aaMarkFACE, 0);

//	iterate over all surrounding volumes and create new vertices.
	for(VolumeIterator iter_sv = sel.volumes_begin(); iter_sv != sel.volumes_end(); ++iter_sv)
	{
		Volume* sv = *iter_sv;

		vector<Vertex*>& newVrts = aaVrtVecVOL[sv];
		newVrts.resize(sv->num_vertices());

	//	check for each vertex whether it lies in the fracture
	//	(aaMarkVRT > 1 in this case)
	//	if so, we have to copy or create a vertex from/in aaVrtVec[vrt] which is
	//	associated with the crease normal on the side of sf.
		for(size_t i_vrt = 0; i_vrt < sv->num_vertices(); ++i_vrt)
		{
			newVrts[i_vrt] = NULL;
			Vertex* vrt = sv->vertex(i_vrt);
			if(aaMarkVRT[vrt] > 1){
			//	calculate the normal on this side of the frac
				vector3 n = CalculateCreaseNormal(grid, sv, vrt, faceIsInFrac, aaPos);
				//UG_LOG("calculated crease normal: " << n << endl);

			//	check whether aaVrtVecs already contains a vertex associated with n.
			//	the normal of new vrts is stored in their position attachment
				vector<Vertex*>& vrtVec = aaVrtVecVRT[vrt];
				for(size_t i = 0; i < vrtVec.size(); ++i){
					//UG_LOG("comparing to: " << aaPos[vrtVec[i]] << endl);
					if(VecDistanceSq(aaPos[vrtVec[i]], n) < SMALL){
					//	got one
						newVrts[i_vrt] = vrtVec[i];
						break;
					}
				}

			//	if we didn't find one then create and associate one.
			//	store the normal in the position attachment of the new vertex
				if(!newVrts[i_vrt]){
					newVrts[i_vrt] = *grid.create<RegularVertex>();
					aaPos[newVrts[i_vrt]] = n;
					aaVrtVecVRT[vrt].push_back(newVrts[i_vrt]);
				}
			}
		}
	}

////////////////////////////////
//	assign positions
	for(VertexIterator iter = sel.vertices_begin();
		iter != sel.vertices_end(); ++iter)
	{
		Vertex* vrt = *iter;

	//	calculate the width as the maximum of associated fracture widths
		CollectFaces(faces, grid, vrt);

		number width = 0;
		for(size_t i = 0; i < faces.size(); ++i){
			if(aaMarkFACE[faces[i]])
				width = max<number>(width, fracInfosBySubset.at(sh.get_subset_index(faces[i])).width);
		}

	//	iterate over associated vertices
		vector<Vertex*>& vrtVec = aaVrtVecVRT[vrt];

	//	note that the position attachment of new vertices holds their normal.
		for(size_t i = 0; i < vrtVec.size(); ++i){
			Vertex* nVrt = vrtVec[i];
			if(width > 0){
				vector3 n = aaPos[nVrt];
				if(VecLengthSq(n) > SMALL)
					VecNormalize(n, n);

				VecScale(n, n, width / 2.);

				if(IsBoundaryVertex3D(grid, vrt)){
				//	First determine the normal pointing outwards
					vector3 nOut;
					CalculateBoundaryVertexNormal3D(nOut, grid, vrt, aaPos);

				//	project the normal into the plane with the normal nOut
					vector3 nNew;
					ProjectPointToPlane(nNew, n, vector3(0, 0, 0), nOut);

				//	and now scale the new offset so that we receive the final offset.
					number dot = VecDot(n, nNew);
					if(dot > SMALL)
						VecScale(n, nNew, VecLengthSq(n) / dot);

				}

				VecAdd(aaPos[nVrt], n, aaPos[vrt]);
			}
			else
				aaPos[nVrt] = aaPos[vrt];
		}

	//	the current position is only a guess. Especially at vertices where
	//	fractures cross, this is not yet optimal.
	//todo: create an iterative spring system to find the new position.
	}

////////////////////////////////
//	create new elements

//	holds local side vertex indices
	vector<size_t>	locVrtInds;
//	all new vertices have been assigned to newVrts.
//	Note that if newVrts[i] == NULL, then we have to take the
//	old vertex sf->vertex(i).

//	first we create new edges from selected ones which are connected to
//	inner vertices. This allows to preserve old subsets.
//	Since we have to make sure that we use the right vertices,
//	we have to iterate over the selected volumes and perform all actions on the edges
//	of those volumes.
	for(VolumeIterator iter_sv = sel.volumes_begin(); iter_sv != sel.volumes_end(); ++iter_sv)
	{
		Volume* sv = *iter_sv;
	//	check for each edge whether it has to be copied.
		for(size_t i_edge = 0; i_edge < sv->num_edges(); ++i_edge){
			Edge* e = grid.get_edge(sv, i_edge);
			if(sel.is_selected(e)){
			//	check the associated vertices through the volumes aaVrtVecVol attachment.
			//	If at least one has an associated new vertex and if no edge between the
			//	new vertices already exists, we'll create the new edge.
				size_t ind0, ind1;
				sv->get_vertex_indices_of_edge(ind0, ind1, i_edge);
				Vertex* nv0 = (aaVrtVecVOL[sv])[ind0];
				Vertex* nv1 = (aaVrtVecVOL[sv])[ind1];

				if(nv0 || nv1){
				//	if one vertex has no associated new one, then we use the vertex itself
					if(!nv0)
						nv0 = sv->vertex(ind0);
					if(!nv1)
						nv1 = sv->vertex(ind1);

				//	create the new edge if it not already exists.
					if(!grid.get_edge(nv0, nv1))
						grid.create_by_cloning(e, EdgeDescriptor(nv0, nv1), e);
				}
			}
		}
	}

//	now we create new faces from selected ones which are connected to
//	inner vertices. This allows to preserve old subsets.
//	Since we have to make sure that we use the right vertices,
//	we have to iterate over the selected volumes and perform all actions on the side-faces
//	of those volumes.
	for(VolumeIterator iter_sv = sel.volumes_begin(); iter_sv != sel.volumes_end(); ++iter_sv)
	{
		Volume* sv = *iter_sv;
	//	check for each face whether it has to be copied.
		for(size_t i_face = 0; i_face < sv->num_faces(); ++i_face){
			Face* sf = grid.get_face(sv, i_face);
			if(sel.is_selected(sf)){
			//	check the associated vertices through the volumes aaVrtVecVol attachment.
			//	If no face between the new vertices already exists, we'll create the new face.
				sv->get_vertex_indices_of_face(locVrtInds, i_face);
				fd.set_num_vertices(sf->num_vertices());

				for(size_t i = 0; i < sf->num_vertices(); ++i){
					Vertex* nVrt = (aaVrtVecVOL[sv])[locVrtInds[i]];
					if(nVrt)
						fd.set_vertex(i, nVrt);
					else
						fd.set_vertex(i, sv->vertex(locVrtInds[i]));
				}

			//	if the new face does not already exist, we'll create it
				if(!grid.get_face(fd))
					grid.create_by_cloning(sf, fd, sf);
			}
		}
	}

//	Expand all faces.
//	Since volumes are replaced on the fly, we have to take care with the iterator.
//	record all new volumes in a vector. This will help to adjust positions later on.
	vector<Volume*> newFractureVolumes;
	for(VolumeIterator iter_sv = sel.volumes_begin(); iter_sv != sel.volumes_end();)
	{
		Volume* sv = *iter_sv;
		++iter_sv;

	//	now expand the fracture faces of sv to volumes.
		for(size_t i_side = 0; i_side < sv->num_sides(); ++i_side){
		//	get the local vertex indices of the side of the volume
			sv->get_vertex_indices_of_face(locVrtInds, i_side);

			Face* tFace = grid.get_side(sv, i_side);

			if(tFace){
				if(aaMarkFACE[tFace]){
					Volume* expVol = NULL;
					if(locVrtInds.size() == 3){
						size_t iv0 = locVrtInds[0];
						size_t iv1 = locVrtInds[1];
						size_t iv2 = locVrtInds[2];

						if((aaVrtVecVOL[sv])[iv0] && (aaVrtVecVOL[sv])[iv1] && (aaVrtVecVOL[sv])[iv2]){
						//	create a new prism
							expVol = *grid.create<Prism>(
										PrismDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
														(aaVrtVecVOL[sv])[iv2], (aaVrtVecVOL[sv])[iv1], (aaVrtVecVOL[sv])[iv0]));
						}
						else if((aaVrtVecVOL[sv])[iv0] && (aaVrtVecVOL[sv])[iv1]){
						//	create a new Pyramid
							expVol = *grid.create<Pyramid>(
										PyramidDescriptor(sv->vertex(iv0), sv->vertex(iv1),
											(aaVrtVecVOL[sv])[iv1], (aaVrtVecVOL[sv])[iv0], sv->vertex(iv2)));
						}
						else if((aaVrtVecVOL[sv])[iv1] && (aaVrtVecVOL[sv])[iv2]){
						//	create a new Pyramid
							expVol = *grid.create<Pyramid>(
										PyramidDescriptor(sv->vertex(iv1), sv->vertex(iv2),
											(aaVrtVecVOL[sv])[iv2], (aaVrtVecVOL[sv])[iv1], sv->vertex(iv0)));
						}
						else if((aaVrtVecVOL[sv])[iv0] && (aaVrtVecVOL[sv])[iv2]){
						//	create a new Pyramid
							expVol = *grid.create<Pyramid>(
										PyramidDescriptor(sv->vertex(iv2), sv->vertex(iv0),
											(aaVrtVecVOL[sv])[iv0], (aaVrtVecVOL[sv])[iv2], sv->vertex(iv1)));
						}
						else if((aaVrtVecVOL[sv])[iv0]){
						//	create a new Tetrahedron
							expVol = *grid.create<Tetrahedron>(
										TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
															 (aaVrtVecVOL[sv])[iv0]));
						}
						else if((aaVrtVecVOL[sv])[iv1]){
						//	create a new Tetrahedron
							expVol = *grid.create<Tetrahedron>(
										TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
															 (aaVrtVecVOL[sv])[iv1]));
						}
						else if((aaVrtVecVOL[sv])[iv2]){
						//	create a new Tetrahedron
							expVol = *grid.create<Tetrahedron>(
										TetrahedronDescriptor(sv->vertex(iv2), sv->vertex(iv1), sv->vertex(iv0),
															 (aaVrtVecVOL[sv])[iv2]));
						}
						else{
						//	this code-block should never be entered. If it is entered then
						//	we either selected the wrong faces (this shouldn't happen), or there
						//	are selected faces, which have fracture-boundary-vertices only.
						//	This is the same is if inner fracture edges exists, which are
						//	connected to two boundary vertices.
						//	Since we tried to remove those edges above, something went wrong.
						//	remove the temporary attachments and throw an error
							grid.detach_from_vertices(aVrtVec);
							grid.detach_from_volumes(aVrtVec);
							grid.detach_from_vertices(aAdjMarker);
							grid.detach_from_edges(aAdjMarker);
							throw(UGError("Error in ExpandFractures3d. Implementation Error."));
						}
					}
					else{
					//	currently only tetrahedrons are supported. This section thus raises an error
						grid.detach_from_vertices(aVrtVec);
						grid.detach_from_volumes(aVrtVec);
						grid.detach_from_vertices(aAdjMarker);
						grid.detach_from_edges(aAdjMarker);
						throw(UGError("Incomplete implementation error in ExpandFractures3d: Only tetrahedrons are supported in the current implementation."));
					}
					if(expVol){
						sh.assign_subset(expVol, fracInfosBySubset.at(sh.get_subset_index(tFace)).newSubsetIndex);
						newFractureVolumes.push_back(expVol);
					}
				}
			}
		}

	//	now set up a new volume descriptor and replace the volume.
		if(vd.num_vertices() != sv->num_vertices())
			vd.set_num_vertices(sv->num_vertices());

		for(size_t i_vrt = 0; i_vrt < sv->num_vertices(); ++i_vrt){
			if((aaVrtVecVOL[sv])[i_vrt])
				vd.set_vertex(i_vrt, (aaVrtVecVOL[sv])[i_vrt]);
			else
				vd.set_vertex(i_vrt, sv->vertex(i_vrt));
		}

		grid.create_by_cloning(sv, vd, sv);
		grid.erase(sv);
	}

//	we have to clean up unused faces and edges.
//	note that all selected edges with mark 0 may safley be deleted.
	for(EdgeIterator iter = sel.begin<Edge>();
		iter!= sel.end<Edge>();)
	{
	//	take care of the iterator
		Edge* e = *iter;
		++iter;

		if(aaMarkEDGE[e] == 0)
			grid.erase(e);
	}

//	make sure that no unused faces linger around (This should never happen!)
	bool foundUnusedFaces = false;
	for(FaceIterator iter = sel.begin<Face>();
		iter != sel.end<Face>();)
	{
		Face* f = *iter;
		++iter;

		if(aaMarkFACE[f] == 0){
			foundUnusedFaces = true;
			grid.erase(f);
		}
	}

	if(foundUnusedFaces){
		UG_LOG("WARNING in ExpandFractures3D: Unused faces encountered during cleanup. Removing them...\n");
	}
/*
//	finally assign the new positions
//	find the maximal fracture width. If it is 0, we only have to copy positions.
	number maxFracWidth = 0;
	for(size_t i = 0; i < fracInfos.size(); ++i)
		maxFracWidth = max(fracInfos[i].width, maxFracWidth);

//	in this case equality with 0 is desired.
	if(maxFracWidth == 0){
	//	set all positions of new vertices to the positions of their parents
		for(VertexIterator iter = sel.vertices_begin();
			iter != sel.vertices_end(); ++iter)
		{
			Vertex* vrt = *iter;
		//	iterate over associated vertices
			vector<Vertex*>& vrtVec = aaVrtVecVRT[vrt];
			for(size_t i = 0; i < vrtVec.size(); ++i){
				Vertex* nVrt = vrtVec[i];
				aaPos[nVrt] = aaPos[vrt];
			}
		}
	}
	else{
	//	we will find the new positions by solving a minimal least squares problem.
	//todo: Special treatment has to be to boundary vertices.

	//	currently all original fracture vertices are selected.
	}*/

//	remove the temporary attachments
	grid.detach_from_vertices(aVrtVec);
	grid.detach_from_faces(aVrtVec);
	grid.detach_from_vertices(aAdjMarker);
	grid.detach_from_edges(aAdjMarker);
	return true;
}

}// end of namespace


//	This method is unused. And I can't really remeber what its use was.
//	It looks quite complicated for a rather simple task. Probably some
//	deeper thought was involved...
///	this method returns true if the face has to be treated as
///	an inner fracture boundary face.
/**	This method uses Grid::mark.
 *
 *	This method starts at all open boundary edges and tries
 *	to find f by checking faces which are connected to one
 *	of the boundary-edges endpoints and are reachable from
 *	the initial selection by traversing faces over regular surface
 *	edges.
 *	If it can be found, the face is a fracture boundary face,
 *	if not, it is a inner fracture face.
 *
 *	This method only makes sense, if funcIsSurfFace(f) evaluates to true.
 */
/*
bool FractureBoundaryFace(Grid& grid, Face* f,
						  CB_ConsiderFace funcIsSurfFace)
{
	if(!funcIsSurfFace(f))
		return false;

	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		UG_LOG("WARNING in FractureBoundaryFace: autoenabling FACEOPT_AUTOGENERATE_EDGES.\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

	grid.begin_marking();

	stack<Edge*> stk;
	vector<Edge*> edges;
	vector<Face*> allFaces;
	vector<Face*> faces;

//	collect all associated fracture-boundary-edges of vertices of e
	for(size_t i = 0; i < f->num_vertices(); ++i){
		for(Grid::AssociatedEdgeIterator iter = grid.associated_edges_begin(f->vertex(i));
			iter != grid.associated_edges_end(f->vertex(i)); ++iter)
		{
			Edge* e = *iter;
			if(!grid.is_marked(e)){
				if(IsBoundaryEdge(grid, e, funcIsSurfFace)){
					grid.mark(e);
					grid.mark(e->vertex(0));
					grid.mark(e->vertex(1));
					stk.push(e);
				}
			}
		}
	}

	while(!stk.empty()){
		Edge* e = stk.top();
	//	we can safely pop e from the stack
		stk.pop();

	//	collect unmarked surface-faces associated with e
		CollectFaces(allFaces, grid, e);
		faces.clear();
		for(size_t i = 0; i < allFaces.size(); ++i){
			if(!grid.is_marked(allFaces[i])){
				if(funcIsSurfFace(allFaces[i]))
					faces.push_back(allFaces[i]);
			}
		}

	//	make sure that we found exactly one face. Otherwise
	//	we would not cross a regular surface edge
		if(faces.size() != 1)
			continue;

		Face* nf = faces[0];

	//	if we found f, it is clear that f is a boundary face.
		if(nf == f){
			grid.end_marking();
			return true;
		}

		grid.mark(nf);

	//	add unmarked edges to the stack and mark them, if they
	//	are connected to a marked vertex
		CollectEdges(edges, grid, nf);

		for(size_t i = 0; i < edges.size(); ++i){
			Edge* e = edges[i];
			if(!grid.is_marked(e)){
				grid.mark(e);
				if(grid.is_marked(e->vertex(0)) || grid.is_marked(e->vertex(1))){
					stk.push(e);
				}
			}
		}
	}

	grid.end_marking();
	return false;
}
*/
