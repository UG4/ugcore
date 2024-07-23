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

#include <boost/function.hpp>
#include <stack>
#include <vector>
#include "expand_layers.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/grid/grid_util.h"
//#include "lib_grid/util/simple_algebra/least_squares_solver.h"

#include <utility>
#include <vector>
#include <type_traits>
#include <limits>
#include <atomic>
#include <cstddef>
#include <bitset>
#include <string>

#include "support.h"


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

using VertFracTrip = VertexFractureTriple<Edge*, Face*, vector3>;

//using VecVertFracTrip = std::vector<VertFracTrip>;

//using VvftIterator = VecVertFracTrip::iterator;

using AttVrtVec = Attachment<vector<Vertex*> >;

//using VertexOfFaceInfo = VertexFractureTriple< std::pair<Edge*, Edge*>, Face*, std::pair<vector3,vector3> >;
//
//using SegmentsFractExtrus = std::vector<VecVertexOfFaceInfo>;

using IndexType = unsigned short;

//using ShiftInfoBasis = std::pair<Edge*, vector3>;
//
//using ShiftInfoSegment = std::pair<ShiftInfoBasis,ShiftInfoBasis>;

//using CrossVertInf = CrossingVertexInfo<Vertex*, IndexType, Edge*, Face* >;
//using CrossVertInf = CrossingVertexInfo<Vertex*, IndexType, Edge*, ShiftInfoSegment  >;
using CrossVertInf = CrossingVertexInfo<Vertex*, IndexType >; //, Edge* >;


// for cases with one fracture and no crossing points
template <typename ASOF >
bool expandSingleFractureAtGivenSide( vector3 const & nOne, vector3 const & nTwo,
 									  Edge * edgeOne, Edge * edgeTwo,
									  Face * facOne, Face * facTwo,
									  vector<FractureInfo> const & fracInfosBySubset,
									  vector3 const & posOldVrt,
									  Grid::VertexAttachmentAccessor<APosition> & aaPos,
									  Grid & grid, SubsetHandler & sh,
									  ASOF const & assoFaces,
									  std::vector<Vertex *> const & nextFracVrt,
									  Grid::FaceAttachmentAccessor<AttVrtVec> & aaVrtVecFace,
									  int & dbg_flachen_passiert,
									  Vertex * iterV									  )
{

	CrossVertInf cvi( nullptr, 0 );

	return expandSingleFractureAtGivenSide( nOne, nTwo,
			  edgeOne, edgeTwo,
			  facOne, facTwo,
			  fracInfosBySubset,
			  posOldVrt,
			  aaPos,
			  grid, sh,
			  assoFaces,
			  nextFracVrt,
			  aaVrtVecFace,
			  dbg_flachen_passiert,
			  iterV,
			  cvi,
			  false
			  );
}

// for the case of crossing points
template <typename ASOF >
bool expandSingleFractureAtGivenSide( vector3 const & nOne, vector3 const & nTwo,
 									  Edge * edgeOne, Edge * edgeTwo,
									  Face * facOne, Face * facTwo,
									  vector<FractureInfo> const & fracInfosBySubset,
									  vector3 const & posOldVrt,
									  Grid::VertexAttachmentAccessor<APosition> & aaPos,
									  Grid & grid, SubsetHandler & sh,
									  ASOF const & assoFaces,
									  std::vector<Vertex *> const & nextFracVrt,
									  Grid::FaceAttachmentAccessor<AttVrtVec> & aaVrtVecFace,
									  int & dbg_flachen_passiert,
									  Vertex * iterV,
									  CrossVertInf & crossVrtInf,
									  bool insertCrossVrtInf = true
									  )
{

#if 1
	// gleiche Seite vermutet oder gegeben

	// average the normals

	vector3 normSum;

	VecAdd( normSum, nOne, nTwo );

	vector3 normSumNormed;

	VecNormalize(normSumNormed, normSum);

	UG_LOG("averaged normal " << normSumNormed << std::endl);

	std::vector<Edge * > attEdg;
	std::vector<Face * > attFac;

	attEdg.push_back( edgeOne );
	attEdg.push_back( edgeTwo );

	attFac.push_back( facOne );
	attFac.push_back( facTwo );

	// jetzt neuen Vertex erzeugen in Richtung der Normalen
	// sonst ist das attachment Schwachsinn!

	vector3 posNewVrt;

	vector3 moveVrt;

	auto subsIndEdgOne = sh.get_subset_index(edgeOne);

	auto subsIndEdgTwo = sh.get_subset_index(edgeTwo);

	if( subsIndEdgOne != subsIndEdgTwo )
	{
		UG_THROW("subsets passen nicht Vereinheitlichung" << std::endl );
	}

	number width = fracInfosBySubset.at(subsIndEdgOne).width;

	// der Faktor ist Käse und muss noch aus den Eingaben übernommen werden
	VecScale(moveVrt, normSumNormed, width/2. );

	VecAdd(posNewVrt, posOldVrt, moveVrt );

	UG_LOG("neuer Vertex " << posNewVrt << std::endl );

	// TODO FIXME hier ist das PROBLEM, SEGFAULT durch create regular vertex

	Vertex * newShiftVrtx = *grid.create<RegularVertex>();
	aaPos[newShiftVrtx] = posNewVrt;

	sh.assign_subset(newShiftVrtx, subsIndEdgOne );

	if( insertCrossVrtInf )
	{
		crossVrtInf.addShiftVrtx(newShiftVrtx, true);
	}
	// only needed in case of crossing vertices



	// alle anhängenden faces müssen noch zu wissen bekommen
	// dass es diesen neuen Vertex gibt, nicht nur die
	// an den edges anhängenden
	// vielleicht gibt es einen Loop über attached faces des
	// Vertex, für die schon bekannten direkt angehängten klar
	// wenn auch dort vermerkt werden muss im Attachment von Seb
	// bei den anderen, die keine Edge haben von der Kluft
	// da muss man die Normale ins Zentrum bestimmen
	// um heraus zu finden, ob sie auf dieser seite sind
	// am besten dann das Attachment der faces für vertizes
	// von Seb recyclen

	// loop über assosciated faces des vertex am besten
	// vermutlich auch noch assosciated edges, um
	// zu markieren, welche weg fallen sollen, wenn
	// nicht von Kluft selber, sondern quasi verschoben
	// und neu erzeugt

	int dbg_FaceIterator = 0;

	for( auto const & ifac : assoFaces )
	{
		bool isFromFrac = false;

		int dbg_innterFacFracIt = 0;

		for( auto const & facFrac : attFac )
		{

			static_assert( std::is_same<  decltype( (facFrac) ), decltype ( ifac ) >::value );

			if( ifac == facFrac )
			{
				isFromFrac = true;

				static_assert( std::is_same< decltype( (facFrac) ), Face * const & >::value  );
				static_assert( std::is_same< decltype( (facFrac) ), decltype( ifac ) >::value  );

			}
		}

		bool atRightSide = false;

		if( isFromFrac )
			atRightSide = true;

		if( !isFromFrac )
		{
			// check if on same side of edge where the normal points to: compute cosinus between vector of face center
			//  perpendicular to the edge

			vector3 facCenter = CalculateCenter( ifac, aaPos );

			vector3 perpendicu;

			if( nextFracVrt.size() != 2 )
			{
				UG_THROW("komische Groesse" << std::endl);
			}

			DropAPerpendicular(perpendicu, facCenter, aaPos[nextFracVrt[0]], aaPos[nextFracVrt[1]]);

			vector3 tmpN;

			VecSubtract(tmpN, facCenter, perpendicu );

			VecNormalize(tmpN, tmpN);

			UG_LOG("Normale zum Face ist " << tmpN << std::endl);

			number cosBetwFracEdgAndDirection2Face = VecDot(tmpN, normSumNormed );

			UG_LOG("Cosinus zur Normalen ist " << cosBetwFracEdgAndDirection2Face << std::endl);

			if( cosBetwFracEdgAndDirection2Face > 0 )
			{
				UG_LOG("assuming face to be on richt side" << std::endl);

				atRightSide = true;

#if ANSCHAULICH_ERZEUGE_SUDOS_ANHANG

				Vertex * otherFacCent = *grid.create<RegularVertex>();
				aaPos[otherFacCent] = facCenter;
				sh.assign_subset(otherFacCent, 5 );

				Vertex * pp = *grid.create<RegularVertex>();
				aaPos[pp] = perpendicu;
				sh.assign_subset(pp, 6 );

				sh.assign_subset(*iterFac,7);
#endif

			}
			else
			{
				UG_LOG("assuming face to be on wrong side" << std::endl);
			}


			dbg_flachen_passiert++;
		}


		if( atRightSide ) // atRightSide ) NOCH FALSCH TODO FIXME muss nur auf richtiger Seite sein
		{

			// ACHTUNG neue Variable Face klein geschrieben im Gegensatz zu Prof. Reiter! nicht später falsche verwenden!
			vector<Vertex*>& newVrts4Fac = aaVrtVecFace[ ifac ];

			IndexType vrtxFnd = 0;

			for(size_t indVrt = 0; indVrt < (ifac)->num_vertices();  indVrt++ )
			{
				Vertex* facVrt = (ifac)->vertex(indVrt);

				if(  facVrt == iterV )
				{
					newVrts4Fac[ indVrt ] = newShiftVrtx;
					vrtxFnd++;
				}
			}

//			if( insertCrossVrtInf )
//			{
////				crossVrtInf.setShiftVrtx(newVrts4Fac);
//				crossVrtInf.addShiftVrtx(newVrts4Fac, true);
//			}
//			// only needed in case of crossing vertices

			if( vrtxFnd <= 0 )
			{
				UG_THROW("vertex not found!" << std::endl);
												}
			else if( vrtxFnd > 1 )
			{
				UG_THROW("vertex zu oft gefunden " << vrtxFnd << std::endl );
			}
			else if ( vrtxFnd == 1 )
			{
			}
			else
			{
				UG_THROW("vertex finden komisch " << std::endl);
			}


		}

		dbg_innterFacFracIt++;



		dbg_FaceIterator++;

	}
#endif

	return true;

}

#ifndef NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX
#define NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX 1
#endif

#ifndef ANSCHAULICH_ERZEUGE_SUDOS_ANHANG
#define ANSCHAULICH_ERZEUGE_SUDOS_ANHANG 0
#endif

#ifndef OLD_PROFREITER_STUFF
#define OLD_PROFREITER_STUFF 0
#endif

bool ExpandFractures2dArte(Grid& grid, SubsetHandler& sh, vector<FractureInfo> const & fracInfos,
						bool expandInnerFracBnds, bool expandOuterFracBnds)
{

//	for(EdgeIterator iter = sh.begin<Edge>(1); iter != sh.end<Edge>(1); ++iter)
//	{
//		size_t sm = sh.num_subsets();
//
//		sh.assign_subset(*iter,sm);
//
//		SplitEdge<RegularVertex>(grid,*iter,true);
//
//		UG_LOG("vertex gesplittet" << std::endl);
//
//		return true;
//
//	}


//	static constexpr bool dehneInnereKluftGrenzpunkteAus = false;

//	expandInnerFracBnds = false;

//	expandOuterFracBnds = true;

//	access position attachment
	if(!grid.has_vertex_attachment(aPosition)){
		UG_LOG("Error in ExpandFractures 2D Arte: Missing position attachment");
		return false;
	}
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);

	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES)){
		UG_LOG("WARNING in Arte 2D CalculateCreaseNormal: grid option FACEOPT_AUTOGENERATE_EDGES autoenabled.\n");
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

#if OLD_PROFREITER_STUFF
	AInt aAdjMarker;	// used to mark how many adjacent fractures a vertex has.
						// 0: no frac, 1: frac-boundary, >1: inner frac vertex
	// TODO FIXME das sieht komisch aus, ist das immer so, wenn z.B. an einer Grenze sich zwei fracs treffen?
	grid.attach_to_vertices_dv(aAdjMarker, 0);
	Grid::VertexAttachmentAccessor<AInt> aaMarkVRT(grid, aAdjMarker);
	grid.attach_to_edges_dv(aAdjMarker, 0);
	Grid::EdgeAttachmentAccessor<AInt> aaMarkEDGE(grid, aAdjMarker);
#endif

//	using IndexType = unsigned short;
	using AttVerFracProp = Attachment<VertexFractureProperties<IndexType> >;
	// attachment pair boundary is fracture, number fractures crossing

	AttVerFracProp aAdjMarkerVFP;

	VertexFractureProperties<IndexType> vfp0( false, 0 );
	// default value: no boundary fracture, no fractures crossing

	grid.attach_to_vertices_dv(aAdjMarkerVFP, vfp0 );
	Grid::VertexAttachmentAccessor<AttVerFracProp> aaMarkVrtVFP(grid, aAdjMarkerVFP);

	ABool aAdjMarkerB; // used to know if an edge is frac edge
	grid.attach_to_edges_dv(aAdjMarkerB, false);
	Grid::EdgeAttachmentAccessor<ABool> aaMarkEdgeB(grid, aAdjMarkerB);


	// die Vertizes, Faces und Edges, die mit einer Kluft zu tun haben

//	using VertFracTrip = VertexFractureTriple<Edge*, Face*, vector3>;
//
	using VecVertFracTrip = std::vector<VertFracTrip>;

	VecVertFracTrip vertexNoInfo;

	using AttVecVertFracTrip = Attachment<VecVertFracTrip>;

	AttVecVertFracTrip aAdjInfoAVVFT;

	grid.attach_to_vertices_dv( aAdjInfoAVVFT, vertexNoInfo );
	Grid::VertexAttachmentAccessor<AttVecVertFracTrip> aaVrtInfoFraTri(grid,  aAdjInfoAVVFT );


	using VecEdge = std::vector<Edge*>;
	using VecFace = std::vector<Face*>;

	using AttVecEdge = Attachment<VecEdge>;
	using AttVecFace = Attachment<VecFace>;

	VecEdge noEdge;
	VecFace noFace;
	AttVecEdge aAdjInfoEdges;
	AttVecFace aAdjInfoFaces;

	grid.attach_to_vertices_dv( aAdjInfoEdges, noEdge );
	Grid::VertexAttachmentAccessor<AttVecEdge> aaVrtInfoAssoEdges( grid, aAdjInfoEdges );

	grid.attach_to_vertices_dv( aAdjInfoFaces, noFace );
	Grid::VertexAttachmentAccessor<AttVecFace> aaVrtInfoAssoFaces( grid, aAdjInfoFaces );


	// das ist Käse, ich brauche für jeden Vertex ein Attachment der Form
	// class VertexTriple, bzw std vektoren von solchen vertex triplen
	// da weiss dann jeder Vertex das triple
	// Kante (damit subdom) - Face - normal (von Kante in Face rein)
	// dann kann man nämlich anhand des Winkels von zwei Normalen
	// von solchen Vertizes bestimmtn, ob sie auf die gleiche Seite der Kante zeigen
	// und dann kann man sie mitteln, sofern die Vertizes keine Kreuzungs Vertizes sind
	// oder keine äusseren Vertizes - ob sie dsa sind, dafür haben wir schon attachments!

	// TODO FIXME diese komischen accessoren sollen jetzt so zugewiesen
	// werden, dass
	/*
	 * jeder Kluftvertex soll wissen, welche Kluftedges ihm anliegen, und welche faces
	 * jede Kluftedge soll wissen, welche Vertizes ihr anliegen, und welche faces
	 * jedes Face, das an eine Kluft anliegt, soll wissen, welche Vertizes und Edges
	 * ihm anliegen
	 * letzteres muss man vielleicht noch ausdehnen, indem man subdomain indizes
	 * dazu bringt
	 * das kann auch notwendig sein fuer alle anderen - wirklich?
	 * die edges und faces kennen ihre subdomain
	 * nur beim Vertex kann in Kreuzungspunkten ein Problem sein
	 * zudem ist die subdomain der faces EGAL, im Zweifel entscheidet die subdomain
	 * der edges, die anliegen bzw der edge, die anliegt!
	 *
	 */

//	iterate over the given fracture infos and select all fracture edges
//	and fracture vertices.
	for(size_t i_fi = 0; i_fi < fracInfos.size(); ++i_fi)
	{
		int fracInd = fracInfos[i_fi].subsetIndex;

		for(EdgeIterator iter = sh.begin<Edge>(fracInd); iter != sh.end<Edge>(fracInd); ++iter)
		{
		//	mark edge and vertices
			sel.select(*iter);

#if OLD_PROFREITER_STUFF
			aaMarkEDGE[*iter] = 1;
#endif

			aaMarkEdgeB[*iter] = true;

		//	select associated vertices
			for(size_t i = 0; i < 2; ++i)
			{
				Vertex* v = (*iter)->vertex(i);
				sel.select(v);

				// wird in jedem Fall inkrimiert, da der Vertex auf jeden Fall mit einer Kante einer frac verbunden sein muss, sonst darf der loop gar nicht darüber gehen
				aaMarkVrtVFP[v]++;

				if( IsBoundaryVertex2D(grid, v) )
					aaMarkVrtVFP[v].setIsBndFracVertex();


#if OLD_PROFREITER_STUFF
				// das ist Sebastians loop, den nicht mehr lassen lassen
				//	if fracture boundaries are expanded, we'll regard all fracture vertices
				//	as inner vertices
				if(expandInnerFracBnds)
				{
					if(!expandOuterFracBnds)
					{
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
#endif

			}
		}
	}


#if OLD_PROFREITER_STUFF
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
	// TODO FIXME was soll das fuer ein Kaese sein?
#endif

	int dbg_flachen_passiert = 0;

	for(VertexIterator iter = sel.begin<Vertex>(); iter != sel.end<Vertex>(); ++iter)
	{

		bool wahl = true;


		// so stimmt es vielleicht, aber ist auch ein komischer Fall, innen expandieren und aussen nicht...... die Frage ist, ob es oonst Sinn macht.....
		if( expandInnerFracBnds && !expandOuterFracBnds && aaMarkVrtVFP[*iter].getIsBndFracVertex() )
			wahl = false;

		static_assert( std::is_same< decltype(*iter), Vertex * >::value );

		bool isBnd = aaMarkVrtVFP[ *iter ].getIsBndFracVertex();
		auto numCrosFrac = aaMarkVrtVFP[ *iter ].getNumberFracEdgesInVertex();

		if( ! isBnd && numCrosFrac == 1 )
		{
			wahl = false;
		}

//		if( ! dehneInnereKluftGrenzpunkteAus )
//		{
//			if( numCrosFrac == 1 ) // inner frac boundary vertex
//			{
//				wahl = false;
//			}
//		}


		if( wahl )
		{
			sel.select(grid.associated_edges_begin(*iter),
						grid.associated_edges_end(*iter));
			sel.select(grid.associated_faces_begin(*iter),
						grid.associated_faces_end(*iter));

			// TODO FIXME hier ein attachment der associated faces und vertizes, am besten als Klasse, die std vertizes davon frisst, an jeden Vertex anhängen
			// so muss man später nicht nochmal über alle Faces und Edges laufen, sondern hat die angehängten schon zur Hand
			// im Fall, dass eine Kreuzung vorliegt, das Ganze irgendwann ordnen, dass nebeneinander liegende Faces und edges in eine verkettete Liste kommen
			// vielleicht das aber nur bei den Schnittvertizes später
			// und vielleicht sollen die Faces Zeiger auf ihre Edges bekommen, und die edges auf die faces, aber wird irgendwann zu wild verzeigert.....
			// vielleicht einfach attachment von faces und edges unabhängig an jeden Fracture Vertex......

			// testen, ob ein Schnittvertex vor liegt, indem die Anzahl der touches getestet wird, anhand einfacher Geometrien testen, was die Anzahl ist

			// mit UG_LOG ausgeben, was die Koordinaten sind, und die Anzahl der hits

			VecFace assFac;
			VecEdge assEdg;

//			UG_LOG("----------" << std::endl);

			for( std::vector<Face *>::iterator iterFac = grid.associated_faces_begin(*iter); iterFac != grid.associated_faces_end(*iter); iterFac++ )
			{
				assFac.push_back(*iterFac);
//				vector3 facCenter = CalculateCenter( *iterFac, aaPos );
//				UG_LOG("fac center " << facCenter << std::endl);

//				sh.assign_subset(*iterFac, 10);
			}

//			UG_LOG("----------" << std::endl);

			for( std::vector<Edge *>::iterator iterEdg = grid.associated_edges_begin(*iter); iterEdg != grid.associated_edges_end(*iter); iterEdg++ )
			{
				assEdg.push_back(*iterEdg);
//				sh.assign_subset(*iterEdg,10);
			}

			aaVrtInfoAssoFaces[*iter] = assFac;
			aaVrtInfoAssoEdges[*iter] = assEdg;

			UG_LOG("marked vertex wahl: " << aaPos[*iter] << " is bnd " << isBnd << " number cross frac " << numCrosFrac << std::endl );

			// fuer Nicht Boundary Vertizes muessen wir durch 2 teilen, damit wir richtige Anzahl der
			// Fracs haben, die durch den spezifischen Vertex durch geht
			// FALSCH war mal so, ist schon abgestellt, es wird angezeigt, wieviele Ecken von Kanten rein kommen

		}
		else
		{
			UG_LOG("marked vertex unwahl: " << aaPos[*iter] << " is bnd " << isBnd << " number cross frac " << numCrosFrac << std::endl );
		}

	}

//	return true;

	using pairIndDbl = std::pair<IndexType,double>;

	std::vector< pairIndDbl > fracSubdom_facePerpendMinVal;

	for( auto const & pf: fracInfos )
	{
		fracSubdom_facePerpendMinVal.push_back( pairIndDbl( pf.subsetIndex, std::numeric_limits<double>::max() ) );
	}

	T_min<double> minDistPerpOverall( std::numeric_limits<double>::max() );

//	for( auto fI : fracInfos )
//	for( size_t fraInfInd = 0; fraInfInd < fracInfos.size(); fraInfInd++ )
	for( auto & fsfpmv : fracSubdom_facePerpendMinVal )
	{
//		int fracInd = fracInfos[fraInfInd].subsetIndex;
//		int fracInd = fI.subsetIndex;

		auto fracInd = fsfpmv.first;

		T_min<double> minDistPerpThisFrac( fsfpmv.second );

		for(EdgeIterator iterEdg = sh.begin<Edge>(fracInd); iterEdg != sh.end<Edge>(fracInd); iterEdg++ )
		{

			// get subdomain of edge

			auto sudoEdg = sh.get_subset_index(*iterEdg);

			static_assert( std::is_same< decltype(sudoEdg), int >::value );

			// get vertices of edge, always 2

			std::vector<Vertex* > verticesEdg;

			static_assert( std::is_same< Vertex*, decltype( (*iterEdg)->vertex(0) ) >::value );

			for( size_t i = 0; i < 2; ++i )
				verticesEdg.push_back( (*iterEdg)->vertex(i) );

			// get attached faces

#if 0
			std::vector<Face* > assFace;

			// TODO FIXME dieser loop kann vielleicht vereinfacht werden, wenn dieser ganze loop umgebaut wird
			// denn die Vertizes kennen dann die assosciated faces schon
			// der Name AssociatedFaceIterator ist sowieso verwirrend, hat mit assosciated nix zu tun, bezieht sich auf alle std Vektoren von Face *
			UG_LOG("XXXXXXXXXXXX" << std::endl);
			for(Grid::AssociatedFaceIterator iterAFI = grid.associated_faces_begin( verticesEdg[0] );
				iterAFI != grid.associated_faces_end( verticesEdg[0] );
				iterAFI++ )
			{

				if(FaceContains( *iterAFI, *iterEdg ))
				{
					assFace.push_back( *iterAFI );


					vector3 facCenter = CalculateCenter( *iterAFI, aaPos );
					UG_LOG("fac center " << facCenter << std::endl);



//					sh.assign_subset( *iterAFI, sh.get_subset_index(*iterEdg));
				}

			}
			UG_LOG("XXXXXXXXX" << std::endl);
#else

			std::vector<Face* > & assFaceVrt0 = aaVrtInfoAssoFaces[verticesEdg[0]];

			std::vector<Face* > assFace;

//			static_assert( std::is_same< decltype( aaVrtInfoAssoFaces[verticesEdg[0]] )[0], std::vector<Face *> >::value );
			//static_assert( std::is_same< decltype( *(aaVrtInfoAssoFaces[verticesEdg[0]]) ), Face * >::value );

//			UG_LOG("XXXXXXXXXXXX" << std::endl);

			for( auto const & ifa : assFaceVrt0 )
			{
				if(FaceContains( ifa, *iterEdg ))
				{
					assFace.push_back( ifa );

//					vector3 facCenter = CalculateCenter( ifa, aaPos );
//					UG_LOG("fac center " << facCenter << std::endl);
				}
			}

			std::vector<Face* > & assFaceVrt1 = aaVrtInfoAssoFaces[verticesEdg[1]];

			for( auto const & ifa : assFaceVrt1 )
			{
				if(FaceContains( ifa, *iterEdg ))
				{
					bool faceContained = false;

					for( auto const & afa : assFace )
					{
						if( afa == ifa )
							faceContained = true;
					}

					if( !faceContained )
						assFace.push_back( ifa );
				}
			}


//			UG_LOG("XXXXXXXXX" << std::endl);

#endif
			// von hier lernen:
			//	VecFace & assoFaces = aaVrtInfoAssoFaces[*iterV ist verticesEdg[0] ];
			//		for( auto const & ifac : assoFaces )
			//		{
			//			static_assert( std::is_same< decltype( ifac ), Face * const & >::value );
			//		}



			// compute normal of edge

			std::vector< vector3 > edgeNormals;

			std::vector<double> perpendDistances;

			for( auto const & fac : assFace )
			{
				vector3 facCenter = CalculateCenter( fac, aaPos );

				vector3 perpendicu;

				DropAPerpendicular(perpendicu, facCenter, aaPos[verticesEdg[0]], aaPos[verticesEdg[1]]);

				double perpendDist = VecLength( perpendicu ); // betrag perpendicu

				perpendDistances.push_back( perpendDist );

				minDistPerpThisFrac( perpendDist );


			//	vector from projection to center is the unnormalized normal
				vector3 tmpN;

				VecSubtract(tmpN, facCenter, perpendicu );

				VecNormalize(tmpN, tmpN);

				edgeNormals.push_back( tmpN );

				static_assert( std::is_same< Edge*, decltype(*iterEdg) >::value );

				static_assert( std::is_same< Face * const &, decltype(fac) >::value );
				static_assert( std::is_same< Face *, decltype( const_cast<Face*>(fac) ) >::value );
				static_assert( std::is_same< vector3, decltype( tmpN ) >::value );

				VertFracTrip infoVertizesThisEdge( *iterEdg, fac, tmpN );

//				UG_LOG("TypE Fac " << typeid(const_cast<Face*>(fac) ).name() << std::endl);
//				UG_LOG("TypE Edg " << typeid( *iterEdg ).name() << std::endl);
//				UG_LOG("TypE Vec " << typeid( tmpN ).name() << std::endl);


				for( auto const & v : verticesEdg )
				{
					static_assert( std::is_same< decltype(v), Vertex * const & >::value );
					static_assert( std::is_same< decltype(const_cast<Vertex*>(v)), Vertex *  >::value );
					aaVrtInfoFraTri[v].push_back( infoVertizesThisEdge );

//					VecVertFracTrip allInfosVrtxThisEdg = aaVrtInfoFraTri[v];

//					static_assert( std::is_same< decltype(  aaVrtInfoFraTri[v] ),  VecVertFracTrip >::value );

//					UG_LOG("type Fac " << typeid( aaVrtInfoFraTri[v][ aaVrtInfoFraTri[v].size() - 1 ].getFace() ).name() << std::endl);
//					UG_LOG("type Edg " << typeid( aaVrtInfoFraTri[v][ aaVrtInfoFraTri[v].size() - 1 ].getEdge() ).name() << std::endl);
//					UG_LOG("type Vec " << typeid( aaVrtInfoFraTri[v][ aaVrtInfoFraTri[v].size() - 1 ].getNormal() ).name() << std::endl);

					static_assert( std::is_same< decltype( aaVrtInfoFraTri[v][ aaVrtInfoFraTri[v].size() - 1 ].getFace() ), Face * >::value );
					static_assert( std::is_same< decltype( aaVrtInfoFraTri[v][ aaVrtInfoFraTri[v].size() - 1 ].getEdge() ), Edge * >::value );
					static_assert( std::is_same< decltype( aaVrtInfoFraTri[v][ aaVrtInfoFraTri[v].size() - 1 ].getNormal() ), vector3 const >::value );
				}

			}

			// damit speichern wir die plus und minus Normale, und das ist alles, und auch
			// gleich wieder weg
			// TODO FIXME besser, wir speichern die gemittelte Normale an den Vertizes
			// vielleihct als attachment pair, das die subdom kennt der frac und die Normale dazu?
			// ziemlich nutzlos, die Normale wie hier gemacht in einen kurzen Vektor zu speichern, der schnell
			// wieder weg ist......
			// wir brauchen alle Normalen zu jeder fracture an jedem fracture Vertex, also die Mittelung vermutlich
			// diese Mittelung kann aber erst stattfinden, wenn wir nachher über die vertizes loopen,
			// hier kennen wir nur die Vertizes, die an derselben Edge anliegen

			UG_LOG("EDGE NORMALS: " << sh.get_subset_index(*iterEdg) << " -> ");

			int j = 0;

			for( auto const & en: edgeNormals )
			{

				for( size_t i = 0; i < 3; i++ )
					UG_LOG( en[i] << ", "  );


				UG_LOG(" --- " << perpendDistances[j] << " ///////// " );

				j++;
			}

			UG_LOG(std::endl);

		}

		fsfpmv.second = minDistPerpThisFrac();

		minDistPerpOverall( fsfpmv.second );

		UG_LOG("first " << fsfpmv.first << " second " << fsfpmv.second << std::endl);
	}

	for( auto const & fsfpmv : fracSubdom_facePerpendMinVal )
	{
		UG_LOG("min dist sd " << fsfpmv.first << " -> " << fsfpmv.second << std::endl  );
	}

	UG_LOG("overall min dist " << minDistPerpOverall() << std::endl);



	// von Sebastian teilweise das Prinzip übernommen, dass die Faces wissen können, was ihre neuen Vertizes sein sollen
	//	create new vertices

	//	we have to associate a vector of vertices with each node in the fracture.
	//	since an empty vector is quite small, we can associate one with each vertex in
	//	the whole grid. This could be optimized if required, by using subset attachments.

	// es reicht vielleicht, hier statt einem Vektor einfach nur einen Vertex * zu storen
//	using AttVrtVec = Attachment<vector<Vertex*> >;
	AttVrtVec attVrtVec;

	//	we  associate a vector of vertices for each face adjacent to the frac.
	//	it will store a set of vertices. An entry contains the new vertex, if the
	//	corresponding vertex is an inner fracture vertex, and NULL if not.
	grid.attach_to_faces(attVrtVec);
	Grid::FaceAttachmentAccessor<AttVrtVec> aaVrtVecFace(grid, attVrtVec);


	//	iterate over all surrounding faces to enable face changes, this loop taken from SR but shortened
	for(FaceIterator iterSurrFac = sel.faces_begin(); iterSurrFac != sel.faces_end(); ++iterSurrFac )
	{
		Face* sf = *iterSurrFac;

		std::vector<Vertex*>& newVrts = aaVrtVecFace[sf];
		newVrts.resize(sf->num_vertices());

		for(size_t i_vrt = 0; i_vrt < sf->num_vertices(); ++i_vrt)
		{
			newVrts[i_vrt] = NULL;
		}
			// erstmal so tun, als ob keine neuen Vertizes erzeugt werden an den alten Vertizes
	}


	// attachment to allow fracture vertizes to know the newly created vertizes
	// due to extrusion which are related to them, in connection with
	// the normals which are an average of the related edges and the faces
	// defining the original normal

	// usage: store edges and associated faces in SAME order in std vectors!
	using ExpandVertexMultiplett = VertexFractureTriple< std::vector<Edge*>, std::vector<Face*>, vector3 >;
	// holds the averaged normal of the related edges and their corresponding faces which give direction

	using VecExpandVertexMultiplett = std::vector<ExpandVertexMultiplett>;


	VecExpandVertexMultiplett vertexMultiplettEmpty;

	using AttVecExpandVertexMultiplett = Attachment<VecExpandVertexMultiplett>;

	AttVecExpandVertexMultiplett aAdjInfoVVEVM;

	grid.attach_to_vertices_dv( aAdjInfoVVEVM, vertexMultiplettEmpty );
	Grid::VertexAttachmentAccessor<AttVecExpandVertexMultiplett> aaVrtExpMP(grid, aAdjInfoVVEVM );


	// am Ende dieser Prozedur sollten alle Vertizes wissen, welche Tripel vom Typ Edge - Face - Normal zum Face hin an ihnen angelagert sind

	// damit weiss, wenn es stimmt, jeder Vertex, der an einer Fracture ist, wieviele Schnittpunkte von Fractures er hat,
	// ob er ein boundary vertex ist, und was für einen Vektor von Tripeln an ihm angehängt sind
	// die subdomain der Fracture muss anhand der subdomain der edge bestimmt werden immer

	UG_LOG("loop over all marked vertizes " << std::endl);

	int dbg_vertizesPassiert = 0;

//	std::vector<Vertex *> crossVrts;
//	std::vector<Vertex *> teeVrts;

//	using CrossVertInf = CrossingVertexInfo<Vertex*, IndexType, Edge*, Face* >;

//	std::vector<CrossingVertexInfo<Vertex*, IndexType> > vecCrossVrtInf;
	std::vector<CrossVertInf > vecCrossVrtInf;

	// needed for crossing points
	using VertexOfFaceInfo = VertexFractureTriple< std::pair<Edge*, Edge*>, Face*, std::pair<vector3,vector3> >;
	// all edges of the attached face - must always be two, the face itself, and the normal vectors of the face in direction of the two edges
	// the size of the normal vector vector also must be two
	// however, if an edge of the face is not a fracture edge, we do not compute the normal, but assign zero as norm
	// for those edges and faces which are Kluft edges, we assign the normal known from the info computed before, vertex fracture triple

	using VecVertexOfFaceInfo = std::vector<VertexOfFaceInfo>;

	using SegmentsFractExtrus = std::vector<VecVertexOfFaceInfo>;


	// jetzt können wir alle Vertizes ablaufen und an ihnen neue Vertizes erzeugen, die anhand der gemittelten Normalen von den Vertizes weg gehen
	// ob zwei anhängende Faces auf der gleichen Seite liegen, wenn es kein Schnittvertex von zwei oder mehr Klüften ist
	// das kann man anhand des Winkels zwischen zwei face Normalen unterscheiden vermutlich
	// dabei müssen die edges sowieso disjunkt sein, sonst ist man sowieso sicher auf verschiedenen Seiten
	// wenn wir es mit einem boundary Vertex zu tun haben, müssen wir weiter überlegen, wie wir die Verschiebung auf die äussere Kante projizieren
	// muss auch mit dem Winkel zu tun haben
	for(VertexIterator iterV = sel.begin<Vertex>(); iterV != sel.end<Vertex>(); ++iterV)
	{

		// POsition dieses Vertex
		vector3 posOldVrt = aaPos[*iterV];

		// vielleicht muss man, wenn die neuen Vertizes da sind, diese auch gleich mit den umliegenden Knoten per neuer Kanten verbinden
		// und die neuen faces erzeugen nach Löschen der alten?
		// oder alle neuen Vertizes wie bei Prof Reiter in einen std Vektor, der als attachment den bisherigen Face Vertizes angehängt wird
		// und Edge Vernichtung und Erzeugung neuer edges und faces wie bei Prof Reiter in Folgeschritten?

		VecVertFracTrip & vecVertFracTrip = aaVrtInfoFraTri[*iterV];

		std::vector<Edge*> & allAssoEdges = aaVrtInfoAssoEdges[*iterV];

		static_assert( std::is_same< decltype( vecVertFracTrip[ vecVertFracTrip.size() - 1 ].getFace() ), Face * >::value );
		static_assert( std::is_same< decltype( vecVertFracTrip[ vecVertFracTrip.size() - 1 ].getEdge() ), Edge * >::value );
		static_assert( std::is_same< decltype( vecVertFracTrip[ vecVertFracTrip.size() - 1 ].getNormal() ), vector3 const >::value );

		for( auto const & vft : vecVertFracTrip )
		{
			static_assert( std::is_same< decltype( vft.getFace() ), Face * >::value );
			static_assert( std::is_same< decltype( vft.getEdge() ), Edge * >::value );
			static_assert( std::is_same< decltype( vft.getNormal() ), vector3 const >::value );

			Face * f = vft.getFace();
			Edge * e = vft.getEdge();
			vector3 n = vft.getNormal();

		}

		using VvftIterator = VecVertFracTrip::iterator;

		VecFace & assoFaces = aaVrtInfoAssoFaces[*iterV];
		// TODO FIXME hier braucht man das nicht zu ordnen
		// aber bei Kreuzpunkten von Klueften muss es so geordnet werden, wie es nebeneinander liegt
		// bei den Edges gibt es auch die benachbarten, und die edges haben das attachment, ob sie Kluftedges sind

//		for( auto const & ifac : assoFaces )
//		{
//			static_assert( std::is_same< decltype( ifac ), Face * const & >::value );
//		}



		// Anzahl der Kreuzungspunkte auslesen und danach unterscheiden, erstmal keine Kreuzung! TODO FIXME

		// irgendwie muessen wir diese Infos jetzt verwerten, um als erstes neue Vertizes zu erzeugen, anfangs für eine Kluft nur
		// und danach die alten Edges und faces löschen und an neuer Stelle neu erzeugen, plus die sowieso neuen,
		// oder Edges verschieben, wenn es möglich ist, die Vertizes zu verschieben, und die Edges und in Folge faces passen sich an,
		// dann müssen nur die neuen edges und faces neu erzeugt werden
		// verschieben der Position des Vertex löst Kaskade aus, dass Edge und Face auch verschoben werden, kann also angewendet werden
		// allerdings Problem, dass die Vertizes dafür verdoppelt werden müssen und die Kanten, sonst kann man sie nicht nach aussen verschieben
		// also doch komplette Neuerzeugung vermutlich..... oder doppeltes Klonen, und das alte bleibt in der Mitte.....

		vector3 posThisVrt =  aaPos[*iterV];

		UG_LOG("vertex at " << posThisVrt << std::endl );

		bool vrtxIsBndVrt = aaMarkVrtVFP[*iterV].getIsBndFracVertex();
		// alternativ wäre möglich: IsBoundaryVertex2D(grid, *iterV)

		UG_LOG("is bndry " << vrtxIsBndVrt << std::endl);

		IndexType numFracsCrossAtVrt = aaMarkVrtVFP[*iterV].getNumberFracEdgesInVertex();

		UG_LOG("number crossing fracs " << numFracsCrossAtVrt << std::endl);

		size_t numbAttTripl = vecVertFracTrip.size();

		UG_LOG("sizes of vft " << numbAttTripl << std::endl );

		if( ! vrtxIsBndVrt )
		{

			if( numFracsCrossAtVrt < 1 )
			{
				UG_THROW("no fracs crossing but marked vertex? << std::endl");
			}
			else if( numFracsCrossAtVrt == 1 )
			{

//				if( numbAttTripl != 0 )
//				{
//					UG_THROW("Anzahl der angehaengten Triples kann nicht stimmen, Vertex einer Kluft ohne Schnittpunkte, nicht am Rand, Kluftende " << std::endl);
//				}



				UG_LOG("END THIS VERTEX NORMAL INNER ENDING CLEFT" << std::endl);


//				if( ! dehneInnereKluftGrenzpunkteAus )
//				{
//					break;
//				}
				// inner vertex where fracture ends
				// TODO FIXME

				// in this case, we have two attached edges, and each of these edges has two attached faces
				// the faces have a naormal, and based on the normal, we can decide which faces belong to the same side of the edges

#if 0

				if( numbAttTripl != 2 )
				{
					UG_THROW("Anzahl der angehaengten Triples kann nicht stimmen, Vertex einer Kluft ohne Schnittpunkte, nicht am Rand, Kluftende " << std::endl);
				}

				// Zuordnung der Edges und Faces, die auf der gleichen Seite der fracture sind

				// und gleich auch Erzeugung der neuen Knoten, die dann
				// in einem Doublett zusammen mit ihren Normalen an die alten Vertizes
				// angehängt werden; der Winkel zur Normalen hilft später, die Seite
				// heraus zu finden, Seite von den Edges

				int dbg_iteratorAblaufen = 0;


#if NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX

				int dbg_laenge = 0;

				for( auto const & vft : vecVertFracTrip )
				{
					dbg_laenge++;

					UG_LOG("VERTEXFRACTRIP" << std::endl);

					vector3 ve = vft.getNormal();

					UG_LOG("NORMAL " << ve << std::endl);

					UG_LOG("laenge " << dbg_laenge << std::endl );
				}

				int dbg_laenge_eins = 0;

#endif



				for( VvftIterator vvftV = vecVertFracTrip.begin();
						vvftV != vecVertFracTrip.end();
						vvftV++
				)
				{

#if NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX
					dbg_laenge_eins++;

					if( dbg_laenge_eins > dbg_laenge )
					{
						break;
					}

#endif

					vector3 nV = vvftV->getNormal();

					Edge * edgeV = vvftV->getEdge();



#if NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX

					UG_LOG("NORMAL " << vvftV->getNormal() << std::endl);
					UG_LOG("LAENGE EINZ " << dbg_laenge_eins << std::endl );
#endif

					Vertex * nextFracVrt;

					IndexType foundThisVrt = 0;

					for( size_t i = 0; i < 2; ++i )
					{
						Vertex * vrtEdgEnd = edgeV->vertex(i);

						if( vrtEdgEnd == *iterV )
						{
							foundThisVrt++;
						}
						else
						{
							nextFracVrt = vrtEdgEnd ;
						}

					}

					if( foundThisVrt != 1 )
					{
						UG_THROW("zu viel zu wenig vertizex one " << std::endl);
					}



					// Klasse schreiben, die als attachment an einen Fracture-Vertex
					// die neuen Vertizes samt ihrer gemittelten Normalen speichert
					// also std::vector von dieser neuen Klasse als Vertex attachment

					std::vector<Edge * > attEdg;
					std::vector<Face * > attFac;

					attEdg.push_back( edgeV );

					Face * facV = vvftV->getFace();

					attFac.push_back( facV );

					// jetzt neuen Vertex erzeugen in Richtung der Normalen
					// sonst ist das attachment Schwachsinn!

					vector3 posNewVrt;

					vector3 moveVrt;

					auto subsIndEdgV = sh.get_subset_index(edgeV);

					number width = fracInfosBySubset.at(subsIndEdgV).width;

//					if( expandInnerFracBnds )
//					{
//						// der Faktor ist Käse und muss noch aus den Eingaben übernommen werden
//						VecScale(moveVrt, nV, width/2. );
//					}
//					else
//					{
//						// auf Annes Wunsch hin werden die Normalen innendrin an einer endenen Kluft zu Null gesetzt
//
//						VecScale(moveVrt, nV, 0. );
//
//					}

					VecScale(moveVrt, nV, width/2. );

					VecAdd(posNewVrt, posOldVrt, moveVrt );

					UG_LOG("neuer Vertex " << posNewVrt << std::endl );

					// TODO FIXME hier ist das PROBLEM, SEGFAULT durch create regular vertex



					Vertex * newShiftVrtx = *grid.create<RegularVertex>();
					aaPos[newShiftVrtx] = posNewVrt;

					sh.assign_subset(newShiftVrtx, subsIndEdgV );



					// fuer was braucheh wir das eigentlich? selber schon vergessen.....

					ExpandVertexMultiplett vrtMtpl( attEdg, attFac, nV );

					aaVrtExpMP[ *iterV ].push_back( vrtMtpl );



					// alle anhängenden faces müssen noch zu wissen bekommen
					// dass es diesen neuen Vertex gibt, nicht nur die
					// an den edges anhängenden
					// vielleicht gibt es einen Loop über attached faces des
					// Vertex, für die schon bekannten direkt angehängten klar
					// wenn auch dort vermerkt werden muss im Attachment von Seb
					// bei den anderen, die keine Edge haben von der Kluft
					// da muss man die Normale ins Zentrum bestimmen
					// um heraus zu finden, ob sie auf dieser seite sind
					// am besten dann das Attachment der faces für vertizes
					// von Seb recyclen

					// loop über assosciated faces des vertex am besten
					// vermutlich auch noch assosciated edges, um
					// zu markieren, welche weg fallen sollen, wenn
					// nicht von Kluft selber, sondern quasi verschoben
					// und neu erzeugt

					int dbg_FaceIterator = 0;



					for( auto const & ifac : assoFaces )
					{
						bool isFromFrac = false;

						for( auto const & facFrac : attFac )
						{

//											static_assert( std::is_same<  decltype( const_cast<Face* & >(facFrac) ), decltype ( ifac ) >::value );
							static_assert( std::is_same<  decltype( (facFrac) ), decltype ( ifac ) >::value );

							if( ifac == facFrac )
							{
								isFromFrac = true;

//													static_assert( std::is_same< decltype( const_cast<Face* & >(facFrac) ), Face * & >::value  );
								static_assert( std::is_same< decltype( (facFrac) ), Face * const & >::value  );
//												static_assert( std::is_same< decltype( const_cast<Face* & >(facFrac) ), decltype( ifac ) >::value  );
								static_assert( std::is_same< decltype( (facFrac) ), decltype( ifac ) >::value  );

							}
						}


						bool atRightSide = false;

						if( isFromFrac )
							atRightSide = true;

						if( !isFromFrac )
						{
							// check if on same side of edge where the normal points to: compute cosinus between vector of face center
							//  perpendicular to the edge




							vector3 facCenter = CalculateCenter( ifac, aaPos );

							vector3 perpendicu;


//							UG_LOG("pos 0 " << aaPos[nextFracVrt[0]] << std::endl);
//							UG_LOG("pos 1 " << aaPos[*iterV] << std::endl);
//							UG_LOG("fac ce " << facCenter << std::endl);

							DropAPerpendicular(perpendicu, facCenter, aaPos[nextFracVrt], aaPos[*iterV]);

//							if( dbg_FaceIterator == 1 )
//							{
//								UG_LOG("huhu a0" << std::endl);
//								return true;
//							}


							vector3 tmpN;

							VecSubtract(tmpN, facCenter, perpendicu );

							VecNormalize(tmpN, tmpN);

							UG_LOG("Normale zum Face ist " << tmpN << std::endl);

							number cosBetwFracEdgAndDirection2Face = VecDot(tmpN, nV );

							UG_LOG("Cosinus zur Normalen ist " << cosBetwFracEdgAndDirection2Face << std::endl);

//							if( dbg_FaceIterator == 1 )
//							{
//								UG_LOG("huhu a" << std::endl);
////								return true;
//							}


							if( cosBetwFracEdgAndDirection2Face > 0 )
							{
								UG_LOG("assuming face to be on richt side" << std::endl);

								atRightSide = true;

#if ANSCHAULICH_ERZEUGE_SUDOS_ANHANG

								Vertex * otherFacCent = *grid.create<RegularVertex>();
								aaPos[otherFacCent] = facCenter;
								sh.assign_subset(otherFacCent, 6 );

								Vertex * pp = *grid.create<RegularVertex>();
								aaPos[pp] = perpendicu;
								sh.assign_subset(pp, 7 );

								sh.assign_subset(ifac,8);

#endif

							}
							else
							{
								UG_LOG("assuming face to be on wrong side" << std::endl);
							}


							dbg_flachen_passiert++;
						}


//						if( dbg_FaceIterator == 1 )
//						{
//							UG_LOG("huhu b" << std::endl);
////							return true;
//						}



						if( atRightSide ) // atRightSide ) NOCH FALSCH TODO FIXME muss nur auf richtiger Seite sein
						{

							// ACHTUNG neue Variable Face klein geschrieben im Gegensatz zu Prof. Reiter! nicht später falsche verwenden!
							vector<Vertex*>& newVrts4Fac = aaVrtVecFace[ ifac ];

							IndexType vrtxFnd = 0;

							for(size_t indVrt = 0; indVrt < (ifac)->num_vertices();  indVrt++ )
							{
								Vertex* facVrt = (ifac)->vertex(indVrt);

								if(  facVrt == *iterV )
								{
									newVrts4Fac[ indVrt ] = newShiftVrtx;
									vrtxFnd++;
								}
							}


							if( vrtxFnd <= 0 )
							{
								UG_THROW("vertex not found!" << std::endl);
							}
							else if( vrtxFnd > 1 )
							{
								UG_THROW("vertex zu oft gefunden " << vrtxFnd << std::endl );
							}
							else if ( vrtxFnd == 1 )
							{
							}
							else
							{
								UG_THROW("vertex finden komisch " << std::endl);
							}


						}

						dbg_FaceIterator++;

					}








				}



				dbg_iteratorAblaufen++;



//				// Ziel: die beiden parallelen Normalen mitteln, und in die jeweiligen beiden Richtungen je einen neuen Vertex erzeugen
//				// irgendwie muss der Vertex oder die Edge besser sogar wissen, dass sie einen neuen Verschiebevertex bekommen hat
//				// denn später müssen neue Edges und neue Faces basierend auf den neuen Vertizes erzeugt werden
//				// vielleicht braucht die edge und das face ein Attachment, das ihnen das sagt, ähnlihc wie VertexTrible std Vektoren?
//
//
//
				UG_LOG("END THIS VERTEX NORMAL INNER ENDING CLEFT" << std::endl);

#endif


//				return true;


			}
			else if( numFracsCrossAtVrt == 2 ) // free line of fracture, no crossing point, not at boundary
			{
				// in this case, we have two attached edges, and each of these edges has two attached faces
				// the faces have a naormal, and based on the normal, we can decide which faces belong to the same side of the edges


				if( numbAttTripl != 4 )
				{

					UG_LOG("NUMBER OF TRIPLETTS " << numbAttTripl << std::endl);

//					return true;

					UG_THROW("Anzahl der angehaengten Triples kann nicht stimmen, Vertex einer Kluft ohne Schnittpunkte, nicht am Rand " << std::endl);
				}

				// Zuordnung der Edges und Faces, die auf der gleichen Seite der fracture sind

				// und gleich auch Erzeugung der neuen Knoten, die dann
				// in einem Doublett zusammen mit ihren Normalen an die alten Vertizes
				// angehängt werden; der Winkel zur Normalen hilft später, die Seite
				// heraus zu finden, Seite von den Edges



				int dbg_iteratorAblaufen = 0;

#if NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX

				int dbg_laenge = 0;

				for( auto const & vft : vecVertFracTrip )
				{
					dbg_laenge++;

					UG_LOG("VERTEXFRACTRIP" << std::endl);

					vector3 ve = vft.getNormal();

					UG_LOG("NORMAL " << ve << std::endl);

					UG_LOG("laenge " << dbg_laenge << std::endl );
				}


				UG_LOG("SINGLE" << std::endl);


				for( VvftIterator vvftOne = vecVertFracTrip.begin();
						vvftOne != vecVertFracTrip.end() - 1;
						vvftOne++
				)
				{


					Edge * edgeOne = vvftOne->getEdge();
					vector3 nOne = vvftOne->getNormal();


					for( VvftIterator vvftTwo = vvftOne + 1;
							vvftTwo != vecVertFracTrip.end();
							vvftTwo++
					)
					{
						Edge * edgeTwo = vvftTwo->getEdge();
						vector3 nTwo = vvftTwo->getNormal();

						number cosinus = VecDot( nOne, nTwo );
						if( edgeOne != edgeTwo )
						{
							UG_LOG("COSI  between " << nOne << " and " << nTwo << " -> " << cosinus << std::endl );
						}

					}
				}

				UG_LOG("SINGLE END" << std::endl);

				int dbg_laenge_eins = 0;

#endif


				for( VvftIterator vvftOne = vecVertFracTrip.begin();
						vvftOne != vecVertFracTrip.end() - 1;
						vvftOne++
				)
				{

#if NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX
					dbg_laenge_eins++;

					if( dbg_laenge_eins > dbg_laenge )
					{
						break;
					}

					int dbg_laenge_zwei = dbg_laenge_eins;
#endif
					int dbg_zweiterIteratorAblaufen = 0;

					vector3 nOne = vvftOne->getNormal();

					Edge * edgeOne = vvftOne->getEdge();



					for( VvftIterator vvftTwo = vvftOne + 1;
							vvftTwo != vecVertFracTrip.end();
							vvftTwo++
					)
					{

#if NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX
						dbg_laenge_zwei++;

						if( dbg_laenge_zwei > dbg_laenge )
						{
							break;
						}

						UG_LOG("NORMAL ONE " << vvftOne->getNormal() << std::endl);
						UG_LOG("NORMAL TWO " << vvftTwo->getNormal() << std::endl);
						UG_LOG("LAENGE EINZ ZWO " << dbg_laenge_eins << " " << dbg_laenge_zwei << std::endl );
#endif

						// dieselben brauchen wir nicht vergleichen
						if( vvftOne == vvftTwo )
						{
							// sollte nie vorkommen!
							UG_THROW("Unsinn " << std::endl);
						}
						else
						{

							Edge * edgeTwo = vvftTwo->getEdge();

							// noch testen, ob nicht die Kante dieselbe ist, geht das?
							// bei der gleichen Ecke ist es unnötig, da es gegensätzlich sein muss



							if( edgeOne != edgeTwo )
							{

								std::vector<Vertex *> nextFracVrt;

								IndexType foundThisVrtOne = 0;

								for( size_t i = 0; i < 2; ++i )
								{
									Vertex * vrtEdgEnd = edgeOne->vertex(i);

									if( vrtEdgEnd == *iterV )
									{
										foundThisVrtOne++;
									}
									else
									{
										nextFracVrt.push_back( vrtEdgEnd );
									}

								}

								if( foundThisVrtOne != 1 )
								{
									UG_THROW("zu viel zu wenig vertizex one " << std::endl);
								}


								IndexType foundThisVrtTwo = 0;

								for( size_t i = 0; i < 2; ++i )
								{
									Vertex * vrtEdgEnd = edgeTwo->vertex(i);

									if( vrtEdgEnd == *iterV )
									{
										foundThisVrtTwo++;
									}
									else
									{
										nextFracVrt.push_back( vrtEdgEnd );
									}

								}

								if( foundThisVrtTwo != 1 )
								{
									UG_THROW("zu viel zu wenig vertizex two " << std::endl);
								}



								vector3 nTwo = vvftTwo->getNormal();

								number cosinus = VecDot( nOne, nTwo );

//								bool vz = ! std::signbit(cosinus);

								UG_LOG("cosinus " << dbg_vertizesPassiert << " between " << nOne << " and " << nTwo << " -> " << cosinus << std::endl );
								//UG_LOG("sign between " << nOne << " and " << nTwo << " -> " << vz << std::endl );



								if( cosinus > 0 )
								{
									// gleiche Seite vermutet

									// sind die edges dieselben? pruefen! gleiche unnoetig - wird oben schon abgefragt

									// Klasse schreiben, die als attachment an einen Fracture-Vertex
									// die neuen Vertizes samt ihrer gemittelten Normalen speichert
									// also std::vector von dieser neuen Klasse als Vertex attachment

#if 1

									Face * facOne = vvftOne->getFace();
									Face * facTwo = vvftTwo->getFace();

									expandSingleFractureAtGivenSide( nOne, nTwo,
																	 edgeOne, edgeTwo,
																	 facOne, facTwo,
																	 fracInfosBySubset,
																	 posOldVrt,
																	 aaPos,
																	 grid, sh,
																	 assoFaces
																	 ,
																	 nextFracVrt,
																	 aaVrtVecFace,
																	 dbg_flachen_passiert,
																	 *iterV
																	);

#else

									// average the normals

									vector3 normSum;

									VecAdd( normSum, nOne, nTwo );

									vector3 normSumNormed;

									VecNormalize(normSumNormed, normSum);

									UG_LOG("averaged normal " << normSumNormed << std::endl);

									std::vector<Edge * > attEdg;
									std::vector<Face * > attFac;

									attEdg.push_back( edgeOne );
									attEdg.push_back( edgeTwo );

//									Face * facOne = vvftOne->getFace();
//									Face * facTwo = vvftTwo->getFace();

									attFac.push_back( facOne );
									attFac.push_back( facTwo );

									// jetzt neuen Vertex erzeugen in Richtung der Normalen
									// sonst ist das attachment Schwachsinn!

									vector3 posNewVrt;

									vector3 moveVrt;

									auto subsIndEdgOne = sh.get_subset_index(edgeOne);

									auto subsIndEdgTwo = sh.get_subset_index(edgeTwo);


									if( subsIndEdgOne != subsIndEdgTwo )
									{
										UG_THROW("subsets passen nicht" << std::endl );
									}




									number width = fracInfosBySubset.at(subsIndEdgOne).width;

									// der Faktor ist Käse und muss noch aus den Eingaben übernommen werden
									VecScale(moveVrt, normSumNormed, width/2. );

									VecAdd(posNewVrt, posOldVrt, moveVrt );

									UG_LOG("neuer Vertex " << posNewVrt << std::endl );

									// TODO FIXME hier ist das PROBLEM, SEGFAULT durch create regular vertex

									Vertex * newShiftVrtx = *grid.create<RegularVertex>();
									aaPos[newShiftVrtx] = posNewVrt;

									sh.assign_subset(newShiftVrtx, subsIndEdgOne );



									// fuer was braucheh wir das eigentlich? selber schon vergessen.....

									ExpandVertexMultiplett vrtMtpl( attEdg, attFac, normSumNormed );

									aaVrtExpMP[ *iterV ].push_back( vrtMtpl );



									// alle anhängenden faces müssen noch zu wissen bekommen
									// dass es diesen neuen Vertex gibt, nicht nur die
									// an den edges anhängenden
									// vielleicht gibt es einen Loop über attached faces des
									// Vertex, für die schon bekannten direkt angehängten klar
									// wenn auch dort vermerkt werden muss im Attachment von Seb
									// bei den anderen, die keine Edge haben von der Kluft
									// da muss man die Normale ins Zentrum bestimmen
									// um heraus zu finden, ob sie auf dieser seite sind
									// am besten dann das Attachment der faces für vertizes
									// von Seb recyclen

									// loop über assosciated faces des vertex am besten
									// vermutlich auch noch assosciated edges, um
									// zu markieren, welche weg fallen sollen, wenn
									// nicht von Kluft selber, sondern quasi verschoben
									// und neu erzeugt

									int dbg_FaceIterator = 0;

#if 0
//									for( auto iterFac = grid.associated_faces_begin(*iterV); iterFac != grid.associated_faces_end(*iterV); iterFac++ )
									for( std::vector<Face *>::iterator iterFac = grid.associated_faces_begin(*iterV); iterFac != grid.associated_faces_end(*iterV); iterFac++ )
									{
										bool isFromFrac = false;

//										for( std::vector<Face *>::iterator iterF2 = attFac.begin(); iterF2 != attFac.end(); iterF2++ )
//										{
//											static_assert( std::is_same< decltype( *iterF2 ), decltype ( *iterFac ) >::value );
//
//										}

										int dbg_innterFacFracIt = 0;

										for( auto const & facFrac : attFac )
										{


//											UG_LOG("type iter facFrac " << typeid( facFrac ).name() << std::endl);
//
//											UG_LOG("type iter Fac " << typeid( *iterFac ).name() << std::endl);

											static_assert( std::is_same<  decltype( const_cast<Face* & >(facFrac) ), decltype ( *iterFac ) >::value );




											if( *iterFac == facFrac )
											{
												isFromFrac = true;

												static_assert( std::is_same< decltype( const_cast<Face* & >(facFrac) ), Face * & >::value  );
												static_assert( std::is_same< decltype( const_cast<Face* & >(facFrac) ), decltype( * iterFac ) >::value  );

											}
										}

										bool atRightSide = false;

										if( isFromFrac )
											atRightSide = true;

										if( !isFromFrac )
										{
											// check if on same side of edge where the normal points to: compute cosinus between vector of face center
											//  perpendicular to the edge
											// TODO FIXME
											// KAESE!!!

											vector3 facCenter = CalculateCenter( *iterFac, aaPos );

											vector3 perpendicu;

											if( nextFracVrt.size() != 2 )
											{
												UG_THROW("komische Groesse" << std::endl);
											}

											DropAPerpendicular(perpendicu, facCenter, aaPos[nextFracVrt[0]], aaPos[nextFracVrt[1]]);

											vector3 tmpN;

											VecSubtract(tmpN, facCenter, perpendicu );

											VecNormalize(tmpN, tmpN);

											UG_LOG("Normale zum Face ist " << tmpN << std::endl);

											number cosBetwFracEdgAndDirection2Face = VecDot(tmpN, normSumNormed );

											UG_LOG("Cosinus zur Normalen ist " << cosBetwFracEdgAndDirection2Face << std::endl);

											if( cosBetwFracEdgAndDirection2Face > 0 )
											{
												UG_LOG("assuming face to be on richt side" << std::endl);

												atRightSide = true;

#if ANSCHAULICH_ERZEUGE_SUDOS_ANHANG

												Vertex * otherFacCent = *grid.create<RegularVertex>();
												aaPos[otherFacCent] = facCenter;
												sh.assign_subset(otherFacCent, 5 );

												Vertex * pp = *grid.create<RegularVertex>();
												aaPos[pp] = perpendicu;
												sh.assign_subset(pp, 6 );

												sh.assign_subset(*iterFac,7);
#endif

											}
											else
											{
												UG_LOG("assuming face to be on wrong side" << std::endl);
											}

//											if( dbg_flachen_passiert == 0 )
//											{
//												UG_LOG("passiert " << dbg_flachen_passiert << std::endl);
//
//												Vertex * otherFacCent = *grid.create<RegularVertex>();
//												aaPos[otherFacCent] = facCenter;
//												sh.assign_subset(otherFacCent, 5 );
//
//												Vertex * pp = *grid.create<RegularVertex>();
//												aaPos[pp] = perpendicu;
//												sh.assign_subset(pp, 6 );
//
//												sh.assign_subset(*iterFac,7);
//
//
//												sh.assign_subset(*iterFac,3);
//
//												UG_LOG("is from frac " << isFromFrac << std::endl);
//
//												return true;
//											}


											dbg_flachen_passiert++;
										}


										if( atRightSide ) // atRightSide ) NOCH FALSCH TODO FIXME muss nur auf richtiger Seite sein
										{


											// ACHTUNG neue Variable Face klein geschrieben im Gegensatz zu Prof. Reiter! nicht später falsche verwenden!
											vector<Vertex*>& newVrts4Fac = aaVrtVecFace[ * iterFac ];

											IndexType vrtxFnd = 0;

											for(size_t indVrt = 0; indVrt < (*iterFac)->num_vertices();  indVrt++ )
											{
												Vertex* facVrt = (*iterFac)->vertex(indVrt);

												if(  facVrt == *iterV )
												{
													newVrts4Fac[ indVrt ] = newShiftVrtx;
	//													UG_LOG("vertex found " <<  indVrt << std::endl );
													vrtxFnd++;
												}
											}


											if( vrtxFnd <= 0 )
											{
												UG_THROW("vertex not found!" << std::endl);
											}
											else if( vrtxFnd > 1 )
											{
												UG_THROW("vertex zu oft gefunden " << vrtxFnd << std::endl );
											}
											else if ( vrtxFnd == 1 )
											{
	//												UG_LOG("vertex found abgeschlossen" << std::endl);
											}
											else
											{
												UG_THROW("vertex finden komisch " << std::endl);
											}


										}

										dbg_innterFacFracIt++;



//
//
//										if( ! isFromFrac )
//										{
//											// Vektor zum Zentrum von KNoten aus berechnen und Winkel zur Normalen bestimmen zur Unterscheidung der Seite
//											// wenn auf richtiger Seite, zuweisen
//										}

										dbg_FaceIterator++;

									}
#else
//									std::vector<Face* > & assFaceVrt = aaVrtInfoAssoFaces[*iterV];

									//									VecFace & assoFaces = aaVrtInfoAssoFaces[*iterV];
																		// TODO FIXME hier braucht man das nicht zu ordnen
																		// aber bei Kreuzpunkten von Klueften muss es so geordnet werden, wie es nebeneinander liegt
																		// bei den Edges gibt es auch die benachbarten, und die edges haben das attachment, ob sie Kluftedges sind

									//									for( auto const & ifac : assoFaces )
									//									{
									//										static_assert( std::is_same< decltype( ifac ), Face * const & >::value );
									//
									//										// TODO FIXME folgenden loop durch diesen ersetzen
									//										// Achtung: Zeigerproblematik, Referenzen, etc.....
									//										// *iterFac ersetzen durch ifac vermutlich, aber wer weiss
									//									}


									//									for( auto iterFac = grid.associated_faces_begin(*iterV); iterFac != grid.associated_faces_end(*iterV); iterFac++ )
									//for( std::vector<Face *>::iterator iterFac = grid.associated_faces_begin(*iterV); iterFac != grid.associated_faces_end(*iterV); iterFac++ )
									for( auto const & ifac : assoFaces )
									{
										bool isFromFrac = false;


										int dbg_innterFacFracIt = 0;

										for( auto const & facFrac : attFac )
										{

//											static_assert( std::is_same<  decltype( const_cast<Face* & >(facFrac) ), decltype ( ifac ) >::value );
											static_assert( std::is_same<  decltype( (facFrac) ), decltype ( ifac ) >::value );

											if( ifac == facFrac )
											{
												isFromFrac = true;

//												static_assert( std::is_same< decltype( const_cast<Face* & >(facFrac) ), Face * & >::value  );
												static_assert( std::is_same< decltype( (facFrac) ), Face * const & >::value  );
//												static_assert( std::is_same< decltype( const_cast<Face* & >(facFrac) ), decltype( ifac ) >::value  );
												static_assert( std::is_same< decltype( (facFrac) ), decltype( ifac ) >::value  );

											}
										}

										bool atRightSide = false;

										if( isFromFrac )
											atRightSide = true;

										if( !isFromFrac )
										{
											// check if on same side of edge where the normal points to: compute cosinus between vector of face center
											//  perpendicular to the edge

											vector3 facCenter = CalculateCenter( ifac, aaPos );

											vector3 perpendicu;

											if( nextFracVrt.size() != 2 )
											{
												UG_THROW("komische Groesse" << std::endl);
											}

											DropAPerpendicular(perpendicu, facCenter, aaPos[nextFracVrt[0]], aaPos[nextFracVrt[1]]);

											vector3 tmpN;

											VecSubtract(tmpN, facCenter, perpendicu );

											VecNormalize(tmpN, tmpN);

											UG_LOG("Normale zum Face ist " << tmpN << std::endl);

											number cosBetwFracEdgAndDirection2Face = VecDot(tmpN, normSumNormed );

											UG_LOG("Cosinus zur Normalen ist " << cosBetwFracEdgAndDirection2Face << std::endl);

											if( cosBetwFracEdgAndDirection2Face > 0 )
											{
												UG_LOG("assuming face to be on richt side" << std::endl);

												atRightSide = true;

#if ANSCHAULICH_ERZEUGE_SUDOS_ANHANG

												Vertex * otherFacCent = *grid.create<RegularVertex>();
												aaPos[otherFacCent] = facCenter;
												sh.assign_subset(otherFacCent, 5 );

												Vertex * pp = *grid.create<RegularVertex>();
												aaPos[pp] = perpendicu;
												sh.assign_subset(pp, 6 );

												sh.assign_subset(*iterFac,7);
#endif

											}
											else
											{
												UG_LOG("assuming face to be on wrong side" << std::endl);
											}


											dbg_flachen_passiert++;
										}


										if( atRightSide ) // atRightSide ) NOCH FALSCH TODO FIXME muss nur auf richtiger Seite sein
										{

											// ACHTUNG neue Variable Face klein geschrieben im Gegensatz zu Prof. Reiter! nicht später falsche verwenden!
											vector<Vertex*>& newVrts4Fac = aaVrtVecFace[ ifac ];

											IndexType vrtxFnd = 0;

											for(size_t indVrt = 0; indVrt < (ifac)->num_vertices();  indVrt++ )
											{
												Vertex* facVrt = (ifac)->vertex(indVrt);

												if(  facVrt == *iterV )
												{
													newVrts4Fac[ indVrt ] = newShiftVrtx;
													vrtxFnd++;
												}
											}


											if( vrtxFnd <= 0 )
											{
												UG_THROW("vertex not found!" << std::endl);
																				}
											else if( vrtxFnd > 1 )
											{
												UG_THROW("vertex zu oft gefunden " << vrtxFnd << std::endl );
											}
											else if ( vrtxFnd == 1 )
											{
											}
											else
											{
												UG_THROW("vertex finden komisch " << std::endl);
											}


										}

										dbg_innterFacFracIt++;



										dbg_FaceIterator++;

									}

#endif

#endif

								}
								else
								{
									// andere Seite vermutet, nichts tun!
								}



							}


						}

						dbg_zweiterIteratorAblaufen++;

					}

					dbg_iteratorAblaufen++;

				}


//				// Ziel: die beiden parallelen Normalen mitteln, und in die jeweiligen beiden Richtungen je einen neuen Vertex erzeugen
//				// irgendwie muss der Vertex oder die Edge besser sogar wissen, dass sie einen neuen Verschiebevertex bekommen hat
//				// denn später müssen neue Edges und neue Faces basierend auf den neuen Vertizes erzeugt werden
//				// vielleicht braucht die edge und das face ein Attachment, das ihnen das sagt, ähnlihc wie VertexTrible std Vektoren?
//
//
//


				UG_LOG("END THIS VERTEX NORMAL COSINE" << std::endl);



			}
			else // two fractures completely crossing, numFracsCrossAtVrt >= 3, i.e. T crossing and two fractures completely crossing
			{

//				CrossingVertexInfo<Vertex*, IndexType> crossVrtInf( *iterV, numFracsCrossAtVrt );

				UG_LOG("number fracs " << numFracsCrossAtVrt << std::endl);

				CrossVertInf crossVrtInf( *iterV, numFracsCrossAtVrt );

				UG_LOG("Nummer vorbei " << std::endl);

//				for( auto const & aae : allAssoEdges )
//				{
//					crossVrtInf.addOriginalFracEdge( aae );
//				}

//				crossVrtInf.setOriginalFracEdge(allAssoEdges);

				// TODO FIXME in case of three fractures, we have to use the method for eine durchgehende fracture
				// auf der Seite, wo die zweite fracture NICHT rein geht

				IndexType countedCrossingFracEdgs = 0;

				// TODO FIXME kreuzende Fractures im Innenraum -> Arte in Reinform implementieren

				// verkettete Liste der anhängenden fractures in Reihenfolge
				// der Anhängung mit INfo, ob eine Kluft vorliegt

				for( auto const & attVFT : vecVertFracTrip )
				{
					Edge * edg = attVFT.getEdge();
					Face * fac = attVFT.getFace();
					vector3 nv = attVFT.getNormal();
				}

//				// hier werden  ALLE attached Faces benötigt, auch die, die zwischen den direkt an den fractures liegenden Faces sind
//
				// copies of all faces and of fractured ones
				auto vVFT = vecVertFracTrip; // caution: COPY, not reference!
				auto aF = assoFaces; // caution: COPY, not reference!

				UG_LOG("Gesamtanzahl faces um Knoten " <<  aF.size() << std::endl );

				//  erstmal die ganzen anhaengenden Faces ordnen, dass wir wissen, in welcher Reihenfolge wir durchlaufen muessen
				// jede Edge hat ein bool attachment schon, das weiss, ob sie Fracture edge ist oder nicht
				// Reihenfolge der faces und die edges auch dazu, vielleicht neues Triple oder dergleiche, dabei zwei edges und zwei normals
				// und wie gesagt, die edges wissen, ob sie fractures sind, dazu keine neuen Variablen notwendig

//				using VertexOfFaceInfo = VertexFractureTriple< std::pair<Edge*, Edge*>, Face*, std::pair<vector3,vector3> >;
//				// all edges of the attached face - must always be two, the face itself, and the normal vectors of the face in direction of the two edges
//				// the size of the normal vector vector also must be two
//				// however, if an edge of the face is not a fracture edge, we do not compute the normal, but assign zero as norm
//				// for those edges and faces which are Kluft edges, we assign the normal known from the info computed before, vertex fracture triple
//
//				using VecVertexOfFaceInfo = std::vector<VertexOfFaceInfo>;

				VecVertexOfFaceInfo orderedFaces;

//				using SegmentsFractExtrus = std::vector<VecVertexOfFaceInfo>;

				SegmentsFractExtrus segments;
				// single components always from one fracture edge to the next one

				VecVertexOfFaceInfo segmentPart;

				// note: we do not attach this info to the vertex, as we only need it local; in principle, in case of further need, it would
				// be usful to establish some sort of attachment

				if( vVFT.size() == 0 )
					UG_THROW("vertex frac triple zu klein an Kreuzung " << std::endl);

				// we start with the first fracture face edge stuff, copy it,  and delete this immidiately
				VertFracTrip startVertFracTrip = vVFT[0];

				vVFT.erase(vVFT.begin());

				bool atFirstTriple = true;

				Face* fracFac = startVertFracTrip.getFace();
				Edge* fracEdg = startVertFracTrip.getEdge();
				vector3 fracNorm = startVertFracTrip.getNormal();

				Edge* originalStartEdge = startVertFracTrip.getEdge();

				if( fracEdg != 0 )
				{
					countedCrossingFracEdgs++;
				}

				// do not change this pointer
				Edge* startEdg = fracEdg;
				Face* startFace = fracFac;

				vector3 startNormal = fracNorm;

				Face* nextFace = NULL;

				UG_LOG("Gesamtanzahl faces um Knoten vor while " <<  aF.size() << std::endl );


				while( aF.size() != 0 )
				{

					UG_LOG("Gesamtanzahl faces um Knoten Anfang while " <<  aF.size() << std::endl );


					Face* face2Append = startFace;
					Edge* startEdg2Append = startEdg;


					IndexType fndCommEdg = 0;
					vector3 nuVe(0,0,0);

					Edge* nextEdge = NULL;

					std::pair<Edge*, Edge *> edge2Append( startEdg2Append, nextEdge );
					std::pair<vector3, vector3 > normal2Append( startNormal, nuVe );


					// if start face and start edge from a triple, then has to be erased this triple, exept for the entire start, as already erased
					if( ! atFirstTriple )
					{
						for( VecVertFracTrip::iterator itAttVFT = vVFT.begin(); itAttVFT !=  vVFT.end(); itAttVFT++ )
						{
							auto vft = *itAttVFT;

							Edge * edgIt = vft.getEdge();

							Face * facIt = vft.getFace();

							if( edgIt == startEdg && facIt == startFace )
							{
								// the first edge if from a fracture and the face is connected to it

								vVFT.erase(itAttVFT);

								normal2Append.first = vft.getNormal();

								if( ! FaceContains( facIt, startEdg ))
								{
									UG_THROW("Face does not contain start edge of its edge" << std::endl);
								}

								break;
							}
						}

					}
					else // we can save the investigation if we have a triple, and we do not need to erase, as already erased.....
					{
						atFirstTriple = false;
					}


					for( auto const & iE : allAssoEdges ) // werden nicht gelöscht, deswegen Zugriff auf attachment direkt
					{
						if( FaceContains(face2Append, iE) )
						{
							fndCommEdg++;

							if( iE != startEdg )
							{
								nextEdge = iE;

								edge2Append.second = iE;

							}
						}


					}

					if( fndCommEdg != 2 )
					{
						UG_THROW("komische Anzahl gemeinsamer Ecke " << fndCommEdg << std::endl);
					}

					if( nextEdge == NULL )
					{
						UG_THROW("wieso keine zweite Ecke gefunden???? " << std::endl);
					}

						if( edge2Append.first == NULL || edge2Append.second == NULL )
					{
						UG_THROW("null immer noch?" << std::endl);
					}

					// erase the face from the list

					IndexType faceFound = 0;

					for( std::vector<Face*>::iterator itFac = aF.begin(); itFac != aF.end(); itFac++ )
					{
						Face * iFa = *itFac;

						if( iFa == startFace &&  FaceContains( iFa, nextEdge ) && FaceContains(iFa, startEdg))
						{
							faceFound++;
						}
					}

					int totalSubsNum = sh.num_subsets();

					int newSubsToAdd = totalSubsNum;

					if( faceFound != 1 )
					{


						sh.assign_subset(startFace,newSubsToAdd++);
						sh.assign_subset(startEdg,newSubsToAdd++);
						sh.assign_subset(nextEdge,newSubsToAdd++);

						int faNum = aF.size();

						UG_LOG("Gesamtzahl faces vor Absturz " << faNum << std::endl);

						UG_LOG("Gesicht in falscher Anztahl gefunden " << faceFound << std::endl);

//						return true;



						UG_THROW("Gesicht in falscher Anztahl gefunden " << faceFound << std::endl);
					}
					else
					{
//						sh.assign_subset(startFace,newSubsToAdd++);
//						sh.assign_subset(startEdg,newSubsToAdd++);
//						sh.assign_subset(nextEdge,newSubsToAdd++);

						int faNum = aF.size();

						UG_LOG("Gesamtzahl faces ohne Absturz " << faNum << std::endl);

					}

					for( std::vector<Face*>::iterator itFac = aF.begin(); itFac != aF.end(); itFac++ )
					{
						Face * iFa = *itFac;

						if( iFa == startFace && FaceContains( iFa, nextEdge ) && FaceContains(iFa, startEdg) )
						{
							aF.erase(itFac);
							break;
						}
					}




					bool sndEdgIsFracEdgeAlso = aaMarkEdgeB[nextEdge];

					bool tripFound = false;

					if( sndEdgIsFracEdgeAlso )
					{

						if( nextEdge != originalStartEdge )
							countedCrossingFracEdgs++;

						// we need to have a look for the next triple

						// check if the next normal is a frac normal which contains the face as well

						for( VecVertFracTrip::iterator itAttVFT = vVFT.begin(); itAttVFT !=  vVFT.end(); itAttVFT++ )
						{
							auto vft = *itAttVFT;

							Edge * edgIt = vft.getEdge();

							Face * facIt = vft.getFace();

							if( edgIt == nextEdge && facIt == face2Append )
							{
								// the second edge if from a fracture and the face is connected to it

								tripFound = true;

								vVFT.erase(itAttVFT);

								normal2Append.second = vft.getNormal();

								if( ! FaceContains( facIt, nextEdge ))
								{
									UG_THROW("Face does not contain edge of its edge" << std::endl);
								}

								break;
							}
						}

					}

					if( ! tripFound && sndEdgIsFracEdgeAlso )
					{
						UG_THROW("Triple nicht gefunden trotz markierter Edge" << std::endl);
					}


					// check if aF or vVFT still contain the former or the next face - must not be the case!

					VertexOfFaceInfo vOFI( edge2Append, face2Append, normal2Append );

					orderedFaces.push_back( vOFI );

					segmentPart.push_back( vOFI );

					if( sndEdgIsFracEdgeAlso )
					{
						segments.push_back( segmentPart );

						segmentPart.clear();
					}


					// what is next face, what is next edge?
					// wie kriegen wir es hin, auch das nächste Triple zu erasen, wenn es jetzt kommt als nächstes?


					startNormal = nuVe;
					startEdg = nextEdge;

					if( aF.size() == 0 )
					{
						if( nextEdge != originalStartEdge )
						{
							UG_THROW("Gesichter leer, aber keine Anfangsecke gefunden" << std::endl);
						}
						else
						{
							break; // while loop zu Ende, raus aus dem while loop, den Rest nicht mehr machen, würde schief gehen zwingendermassen
						}

					}


					// bleibt noch das nächste Gesicht heraus zu finden, dafür kommt eigentlich nur noch eines in Frage, da das zweite Gesicht vom edge
					// geloescht sein muss in aF, es muss das einzig übrige face sein, das die jetzt start edge enthält, davon darf es nur eines geben, wir löschen aber noch nicht

					IndexType nextFaceFound = 0;

					for( std::vector<Face*>::iterator itFac = aF.begin(); itFac != aF.end(); itFac++ )
					{
						Face * iFa = *itFac;

						if( FaceContains(iFa, startEdg ) )
						{
							nextFaceFound++;
						}
					}

					if( nextFaceFound != 1 )
					{
						UG_THROW("folgendes Gesicht in falscher Anztahl gefunden " << nextFaceFound << std::endl);
					}

					for( std::vector<Face*>::iterator itFac = aF.begin(); itFac != aF.end(); itFac++ )
					{
						Face * iFa = *itFac;

						if( FaceContains(iFa, startEdg ) )
						{
							startFace = iFa;
							break;
						}
					}


				}

				if( vVFT.size() != 0 )
				{
					UG_THROW("not all triples found! " << std::endl);
				}

				if( aF.size() != 0 )
					UG_THROW("not all faces found " << std::endl);

				if( startEdg != originalStartEdge )
				{
					UG_THROW("wir sind nicht am Anfang wieder angekommen" << std::endl);
				}


				if( segmentPart.size() != 0 )
				{
					UG_THROW("die Segmentteile muessen alle verarbeitet sein"  << std::endl);
				}

				UG_LOG("Kreislauf geschlossen" << std::endl);


				// test if the segments and their partition produce sumething useful, for debug purposes

				// als nächstes muss man die Klassen von durch Klüften abgetrennten ordered Faces durchgehen, und die Verschiebevertizes erzeugen
				// als nächstes die verschiedenen Sektionen durch gehen, eventuell nochmal extra Objekte dafür erzeugen
				// oder gleich beim Durchgehen die neuen Vertizes erzeugen, Startsignal durch ein Face mit erster Edge KLuft, und dann die nächste
				// Kluftedge finden, egal ob vom gleihen face oder von einem späteren face im kreis

				// now figure out to which face this next edge belongs, and if this is a fracture edge, then we have the triple and the normal info
				// else we let the normal zero

				// figure out if second edge is also frac edge, i.e. if it belongs to an edge  of the remaining vVFT elements
				// first easily asking if it is marked as frac edge, simplifies research

				// in principle from here on need to loop through all triples and through all faces, find out some way to construct next edge and to
				// build one element after the other of the ordered faces vector, still even the first element is not completed

				// TODO FIXME kreuzende Fractures im Innenraum -> Arte in Reinform implementieren
				// later assign somehow next edge to start edge, or use new variable, when we have figured out next face
				// at end, chech if we have arrived again at original first edge


				int totalSubsNum = sh.num_subsets();

//				int newSubsToAdd = totalSubsNum;

//				for( VertexOfFaceInfo const & vertFracInfo : orderedFaces )
//				{
////					Face * fa = vertFracInfo.getFace();
////
////					sh.assign_subset(fa,newSubsToAdd++);
//				}

				number totAnglsEdg = 0;
				number totAnglsNrm = 0;

				for( VecVertexOfFaceInfo const & segPart : segments )
				{
//					newSubsToAdd++;

					IndexType numbTriangs = segPart.size();

					if( numbAttTripl == 0 )
					{
						UG_THROW("zu wenig Dreiecke " << std::endl);
					}

					VertexOfFaceInfo const & vFISBegin = segPart[0];
					VertexOfFaceInfo const & vFISEnd = segPart[numbTriangs-1];

					std::pair<Edge*, Edge* >  edgesBegin = vFISBegin.getEdge();
					std::pair<Edge*, Edge* >  edgesEnd = vFISEnd.getEdge();

					std::pair<vector3, vector3> normalBegin = vFISBegin.getNormal();
					std::pair<vector3, vector3> normalEnd = vFISEnd.getNormal();

					Edge* edgeFracOne = edgesBegin.first;
					Edge* edgeFracTwo = edgesEnd.second;

					auto subsIndFracOne = sh.get_subset_index(edgeFracOne);
					auto subsIndFracTwo = sh.get_subset_index(edgeFracTwo);

					vector3 normalFracOne = normalBegin.first;
					vector3 normalFracTwo = normalEnd.second;

//					sh.assign_subset(edgeFracOne, newSubsToAdd);
//
//					if( edgeFracTwo != originalStartEdge )
//						sh.assign_subset(edgeFracTwo, newSubsToAdd);

					// neue Punkte erzeugen

					number cosBetweenNormals = VecDot( normalFracOne, normalFracTwo );

					if( subsIndFracOne == subsIndFracTwo )
					{
						if( numFracsCrossAtVrt != 3 )
						{
							UG_THROW("Fracture Segment an beiden Seiten gleiche sudo, aber keine T Kreuzung?" << std::endl);
						}

						// dieselben Methoden wie im Fall von einer durchgehenden Kluft an einem Vertex, dort kopieren
						// bzw Funktion schreiben, die beides macht

						// hier wird der Fall abgezweigt, dass wir auf der durchgehenden Seite einer Kluft sind
						// wenn wir eine T-Kreuzung haben

						std::vector<Vertex *> nextFracVrt;

						IndexType foundThisVrtOne = 0;

						for( size_t i = 0; i < 2; ++i )
						{
							Vertex * vrtEdgEnd = edgeFracOne->vertex(i);

							if( vrtEdgEnd == *iterV )
							{
								foundThisVrtOne++;
							}
							else
							{
								nextFracVrt.push_back( vrtEdgEnd );
							}

						}

						if( foundThisVrtOne != 1 )
						{
							UG_THROW("zu viel zu wenig vertizex one " << std::endl);
						}


						IndexType foundThisVrtTwo = 0;

						for( size_t i = 0; i < 2; ++i )
						{
							Vertex * vrtEdgEnd = edgeFracTwo->vertex(i);

							if( vrtEdgEnd == *iterV )
							{
								foundThisVrtTwo++;
							}
							else
							{
								nextFracVrt.push_back( vrtEdgEnd );
							}

						}

						if( foundThisVrtTwo != 1 )
						{
							UG_THROW("zu viel zu wenig vertizex two " << std::endl);
						}

						Face *  faceBegin = vFISBegin.getFace();
						Face *  faceEnd = vFISEnd.getFace();



						expandSingleFractureAtGivenSide( normalFracTwo, normalFracTwo,
														 edgeFracOne, edgeFracTwo,
														 faceBegin, faceEnd,
														 fracInfosBySubset,
														 posOldVrt,
														 aaPos,
														 grid, sh,
														 assoFaces
														 ,
														 nextFracVrt,
														 aaVrtVecFace,
														 dbg_flachen_passiert,
														 *iterV,
														 crossVrtInf
														);



					}
					else
					{

						// create normal vectors into direction of relevant edges

						vector3 alongEdgeOne;
						vector3 alongEdgeTwo;

						Vertex * vrtEdgeOneBegin = nullptr;
						Vertex * vrtEdgeTwoBegin = nullptr;
						Vertex * vrtEdgeOneEnd = nullptr;
						Vertex * vrtEdgeTwoEnd = nullptr;


						for( size_t i = 0; i < 2; ++i )
						{
							Vertex * vrtFromEdgeOne = edgeFracOne->vertex(i);
							Vertex * vrtFromEdgeTwo = edgeFracTwo->vertex(i);

							if( vrtFromEdgeOne == *iterV )
							{
								vrtEdgeOneBegin = vrtFromEdgeOne;
							}
							else
							{
								vrtEdgeOneEnd = vrtFromEdgeOne;
							}

							if( vrtFromEdgeTwo == *iterV )
							{
								vrtEdgeTwoBegin = vrtFromEdgeTwo;
							}
							else
							{
								vrtEdgeTwoEnd = vrtFromEdgeTwo;
							}

						}

						if( vrtEdgeOneBegin == NULL || vrtEdgeTwoBegin == NULL || vrtEdgeOneEnd == NULL || vrtEdgeTwoEnd == NULL )
						{
							UG_THROW("lauter Nullen vertizes" << std::endl);
						}

						vector3 fracVrtPos = aaPos[*iterV];

						vector3 fracEdgOneEndPos = aaPos[ vrtEdgeOneEnd ];
						vector3 fracEdgTwoEndPos = aaPos[ vrtEdgeTwoEnd ];

						vector3 directionEdgOne;
						VecSubtract(directionEdgOne, fracEdgOneEndPos, fracVrtPos);

						vector3 directionEdgTwo;
						VecSubtract(directionEdgTwo, fracEdgTwoEndPos, fracVrtPos);

						vector3 nrmdVecEdgOne;
						VecNormalize(nrmdVecEdgOne, directionEdgOne);

						vector3 nrmdVecEdgTwo;
						VecNormalize(nrmdVecEdgTwo, directionEdgTwo);

						number cosBetweenEdges = VecDot(nrmdVecEdgOne,nrmdVecEdgTwo);

						number angleEdges = std::acos( cosBetweenEdges );
						number angleNormls = std::acos( cosBetweenNormals );

						totAnglsEdg += angleEdges;
						totAnglsNrm += angleNormls;

						UG_LOG("cosinus Edges Normals " << cosBetweenEdges << "  " << cosBetweenNormals << std::endl);
						UG_LOG("angles edges normals " << angleEdges << "  " << angleNormls << std::endl);

						// prject normal 1 onto edge 2 and normal 2 on edge 1, scale with width one half resp with width two half


						number cosBetweenNrmFraOneEdgTwo = VecDot(normalFracOne,nrmdVecEdgTwo);
						number cosBetweenNrmFraTwoEdgOne = VecDot(normalFracTwo,nrmdVecEdgOne);

						vector3 projectNrmFraOneToEdgTwoDirection;
						VecScale(projectNrmFraOneToEdgTwoDirection, nrmdVecEdgTwo, 1./cosBetweenNrmFraOneEdgTwo);

						vector3 projectNrmFraTwoToEdgOneDirection;
						VecScale(projectNrmFraTwoToEdgOneDirection, nrmdVecEdgOne, 1./cosBetweenNrmFraTwoEdgOne);

	//					auto subsIndFracOne = sh.get_subset_index(edgeFracOne);
	//					auto subsIndFracTwo = sh.get_subset_index(edgeFracTwo);

						number shiftOne = fracInfosBySubset.at( subsIndFracOne ).width / 2. ;
						number shiftTwo = fracInfosBySubset.at( subsIndFracTwo ).width / 2. ;

						vector3 shiftAlongEdgeTwo;
						VecScale(shiftAlongEdgeTwo, projectNrmFraOneToEdgTwoDirection, shiftOne);

						vector3 shiftAlongEdgeOne;
						VecScale(shiftAlongEdgeOne, projectNrmFraTwoToEdgOneDirection, shiftTwo);

						vector3 shiftPart;
						VecAdd(shiftPart, fracVrtPos, shiftAlongEdgeTwo);

						vector3 posNewVrt;
						VecAdd( posNewVrt, shiftPart, shiftAlongEdgeOne);

						UG_LOG("neuer Vertex Kreuzung " << posNewVrt << std::endl );

						Vertex * newShiftVrtx = *grid.create<RegularVertex>();
						aaPos[newShiftVrtx] = posNewVrt;

						//					sh.assign_subset(newShiftVrtx,  newSubsToAdd );
						sh.assign_subset(newShiftVrtx, subsIndFracOne );
						// could also be two, but have to select one, no Kompromiss possible......

						crossVrtInf.addShiftVrtx(newShiftVrtx);

//						UG_LOG("ADDED SHIFT VECTOR " << aaPos[newShiftVrtx] << std::endl);


						for( VertexOfFaceInfo const & vertFracInfoSeg : segPart )
						{
							Face * fac = vertFracInfoSeg.getFace();

	//						sh.assign_subset(fa,newSubsToAdd);


							// ACHTUNG neue Variable Face klein geschrieben im Gegensatz zu Prof. Reiter! nicht später falsche verwenden!
							vector<Vertex*>& newVrts4Fac = aaVrtVecFace[ fac ];

							IndexType vrtxFnd = 0;

							for(size_t indVrt = 0; indVrt < (fac)->num_vertices();  indVrt++ )
							{
								Vertex* facVrt = (fac)->vertex(indVrt);

								if(  facVrt == *iterV )
								{
									newVrts4Fac[ indVrt ] = newShiftVrtx;
									vrtxFnd++;

//									crossVrtInf.addShiftVrtx(newShiftVrtx);
//
//									UG_LOG("ADDED SHIFT VECTOR " << aaPos[newShiftVrtx] << std::endl);

								}
							}

//							crossVrtInf.setShiftVrtx(newVrts4Fac);

							if( vrtxFnd <= 0 )
							{
								UG_THROW("vertex not found kreuzung!" << std::endl);
																}
							else if( vrtxFnd > 1 )
							{
								UG_THROW("vertex zu oft gefunden kreuzung " << vrtxFnd << std::endl );
							}
							else if ( vrtxFnd == 1 )
							{
							}
							else
							{
								UG_THROW("vertex finden komisch kreuzung " << std::endl);
							}

						}
					}
				}

				UG_LOG("sum angles edges normals " << totAnglsEdg << "  " << totAnglsNrm << std::endl);

//				return true;

#if 0
				// das folgende ist vermutlich Unsinn TODO FIXME, waren wohl Versuche am Anfang..... nochmal prüfen!!!!

				// get starting point of the "rotation" around the vertex where fractures are crossing
//				for( auto & attVFT : vVFT ) // not const, as we need to erase found elements!
				for( VecVertFracTrip::iterator itAttVFT = vVFT.begin(); itAttVFT !=  vVFT.end(); )
				{
					Face * facAtVrtWithFrac = itAttVFT->getFace();

					bool facFound = false;

//					for( auto const & ifac : assoFaces ) // not const, as we need to erase found elements!
					for( std::vector<Face*>::iterator itFac = aF.begin(); itFac != aF.end(); )
					{
						if( *itFac == facAtVrtWithFrac )
						{
							// found a starting face
							// copy first the found info, then delete the original one
							auto attVFTCop = *itAttVFT; // copy, not reference!

							vVFT.erase(itAttVFT);
							aF.erase(itFac);

							// TODO FIXME erase ifac and attVFT, how to do?

							Face * startFace = facAtVrtWithFrac;

							// now determine the common edge(s), the first edge of the vector must be a frac edge, the second one might be one

							Edge * startEdg = attVFTCop.getEdge();

							// unnecessary check, but for test purposes at beginning, later should be removed
							if( !FaceContains(facAtVrtWithFrac, startEdg ))
							{
								UG_THROW("face hat ecke nicht, die es haben soll" << std::endl);
							}

							// loop around the edges of the ifac face attached to the vertex

							// determin second edge of the startFace, not stored in the vecVertFracTrip information
							// check number of common edges containing the same vertex

							IndexType fndCommEdg = 0;
							std::vector<Edge*> assoEdg2Fac;

							assoEdg2Fac.push_back( startEdg );

							std::vector<vector3> assoNorm;

							vector3 norm2Frac = attVFTCop.getNormal();

							Edge * secondEdge;
							secondEdge = NULL;

							for( auto const & iE : allAssoEdges ) // werden nicht gelöscht, deswegen Zugriff auf attachment direkt
							{
								if( FaceContains(facAtVrtWithFrac, iE) )
								{
									fndCommEdg++;
									if( iE != startEdg )
									{
										secondEdge = iE;
									}
								}

								if( fndCommEdg != 2 )
								{
									UG_THROW("komische Anzahl gemeinsamer Ecke " << fndCommEdg << std::endl);
								}

								if( secondEdge == NULL )
								{
									UG_THROW("wieso keine zweite Ecke gefunden???? " << std::endl);
								}

								assoEdg2Fac.push_back(secondEdge);

								// check, if second edge belongs to anothter fracture fac, i.e. if it is also a fracture edge

								// check if second edge is edge of a fracture face, i.e. either this edge, or another one is from a fracture

								bool scndEdgIsFracEdg = aaMarkEdgeB[secondEdge];

								if( scndEdgIsFracEdg )
								{
									// TODO FIXME figure out second vertex fracture info, in this case, we have found the next part of the chain!

									for( VecVertFracTrip::iterator itAttVFTTwo = vVFT.begin(); itAttVFTTwo !=  vVFT.end(); )
									{
										// need to erase also this element soon, also in the list of all attached elements

										Face* vFTTwoFac = itAttVFTTwo->getFace();

										if( FaceContains( vFTTwoFac, secondEdge))
										{
											auto vVFT2 = *itAttVFTTwo;
											vVFT.erase( itAttVFTTwo );

											Face * nextFace = vFTTwoFac;

											if( secondEdge != vVFT2.getEdge() )
											{
												UG_THROW("Ecke nicht Ecke " << std::endl);
											}


										}
									}
								}
								else // find the next attached face, either from the
								{

								}
							}

							break;
						}

						if( ! facFound )
						{
							UG_THROW("Gesicht unauffindbar" << std::endl);
						}

						break;
					}
				}

				while( vVFT.size() != 0 )
				{
					while( aF.size() != 0 )
					{

					}
				}

#endif

//				if( numFracsCrossAtVrt == 3 )
//				{
//					teeVrts.push_back(*iterV);
//				}
//				else if( numFracsCrossAtVrt == 4 )
//				{
//					crossVrts.push_back(*iterV);
//				}

				IndexType groesseShiVe = crossVrtInf.getVecShiftedVrts().size();

				auto gro = groesseShiVe;

				UG_LOG("GROESSE SHIFT VECTORS " << gro << std::endl);

				vecCrossVrtInf.push_back(crossVrtInf);

			}

		}
//		// // different treatment for boundary vertizes
		else
		{

			// es muss wohl noch ein Problem mit den Verschiebungen bei boundary Vertizes geben.....


			if( numFracsCrossAtVrt < 1 )
			{
				UG_THROW("no fracs crossing but marked vertex at boundary? << std::endl");
			}
			else if( numFracsCrossAtVrt == 1 ) // no crossing point  at boundary
			{
				// in this case, we have ONE attached edges, the edges has two attached faces
				// the faces have a normal, and based on the normal, we can decide which faces belong to the same side of the edges

				if( numbAttTripl != 2 )
				{
					UG_THROW("Anzahl der angehaengten Triples kann nicht stimmen, Vertex einer Kluft ohne Schnittpunkte  am Rand " << std::endl);
				}

				// Zuordnung der Edges und Faces, die auf der gleichen Seite der fracture sind

				// und gleich auch Erzeugung der neuen Knoten, die dann
				// in einem Doublett zusammen mit ihren Normalen an die alten Vertizes
				// angehängt werden; der Winkel zur Normalen hilft später, die Seite
				// heraus zu finden, Seite von den Edges

				// get  edges adjacent to this vertex which lie on the boundary of the domain

				std::vector<Edge* > adjBndEdgs;

//				std::vector<Edge*> & allAssoEdges = aaVrtInfoAssoEdges[*iterV];

//				for( std::vector<Edge*>::iterator iterBVEdg = grid.associated_edges_begin(*iterV); iterBVEdg != grid.associated_edges_end(*iterV); iterBVEdg++  )
//				{
//					if( IsBoundaryEdge2D(grid,*iterBVEdg) )
//					{
//						adjBndEdgs.push_back( *iterBVEdg );
//					}
//				}
				for( auto const & iBVE : allAssoEdges )
				{
					if( IsBoundaryEdge2D(grid, iBVE ) )
					{
						adjBndEdgs.push_back( iBVE );
					}
				}

				// to compute the normals, compute the vector of the edge and normalize it
				std::vector<vector3> bndEdgeDirection;

				for( auto const & bE : adjBndEdgs )
				{

					// get vertices, i.e. get seocnd vertex, first one must be known

//					std::vector<Vertex* > verticesEdg;

					static_assert( std::is_same< Edge* const &, decltype( bE ) >::value );

					static_assert( std::is_same< Vertex*, decltype( bE->vertex(0) ) >::value );

					IndexType fndIV = 0;

					Vertex * vrtOtherEdg;
					vrtOtherEdg = NULL;

					for( size_t i = 0; i < 2; ++i )
					{
//						verticesEdg.push_back( adjBndEdgs.vertex(i) );

						Vertex * vrtOfEdg = bE->vertex(i);

						if( vrtOfEdg == *iterV )
						{
							fndIV++;
						}
						else
						{
							vrtOtherEdg = vrtOfEdg;
						}
					}

					vector3 posOtherVrt = aaPos[vrtOtherEdg];

					UG_LOG("BOUNDARY EDGE VERTIZES " << posOldVrt << ", " << posOtherVrt << std::endl);

					vector3 fromIterV2Other;

					VecSubtract(fromIterV2Other, posOtherVrt, posOldVrt);

					vector3 nV;

					VecNormalize(nV, fromIterV2Other);

					bndEdgeDirection.push_back(nV);
				}



#if NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX

				IndexType dbg_lim = vecVertFracTrip.size();

				int dbg_cnt = 0;
#endif

				for( VvftIterator vvftAtBnd = vecVertFracTrip.begin();
						vvftAtBnd != vecVertFracTrip.end();
						vvftAtBnd++
				)
				{
#if NOTLOESUNG_EINSCHALTEN_SEGFAULT_CREATE_VERTEX

					if( dbg_lim == dbg_cnt )
					{
						UG_LOG("DARF NICHT SEIN" << std::endl);
						break;
					}

					dbg_cnt++;
#endif

					// Ziel: den parallelen Anteil der Normalen auf die jeweilige Randkante projizieren

					vector3 nrmEdg = vvftAtBnd->getNormal();

					Edge * edgeOfFrac = vvftAtBnd->getEdge();

					// figure out the adjoint boundary edge into the same direction

					// the normal in both directions have to be compared with the vectors in direction of boundary edges
					for( auto bED : bndEdgeDirection )
					{
						// check orientation of boundary edges wrt the normals

						number cosinus = VecDot( nrmEdg, bED );

						UG_LOG("BOUNDARY COSINUS between " << nrmEdg << " and " << bED << " -> " << cosinus << std::endl);

						if( cosinus > 0 )
						{
							// gleiche Seite vermutet

							// muessen wissen, wie lange das gestreckt werden soll

							vector3 alongEdgV;

							auto subsIndEdgOF = sh.get_subset_index(edgeOfFrac);

							number width = fracInfosBySubset.at(subsIndEdgOF).width;

							number scal = width / 2. / cosinus;

							VecScale( alongEdgV, bED, scal );

							vector3 posNewVrtOnBnd;

							VecAdd(posNewVrtOnBnd, posOldVrt, alongEdgV );

							UG_LOG("neuer Vertex Edge " << posNewVrtOnBnd << std::endl );

							Vertex * newShiftEdgVrtx = *grid.create<RegularVertex>();
							aaPos[newShiftEdgVrtx] = posNewVrtOnBnd;

							sh.assign_subset(newShiftEdgVrtx, subsIndEdgOF );

							std::vector<Edge * > attEdg;
							std::vector<Face * > attFac;

							attEdg.push_back(edgeOfFrac);

							Face * facFrac = vvftAtBnd->getFace();

							attFac.push_back( facFrac );

							// we store the boundary edge direction for boundary verzizes rather than the normal, in contrast to inner vertizes, where we store the averaged normal
							ExpandVertexMultiplett vrtMtpl( attEdg, attFac, bED );

							aaVrtExpMP[ *iterV ].push_back( vrtMtpl );

#if 0
							// the attached faces need to know that they need a new vertex to be shifted
							for( std::vector<Face *>::iterator iterFac = grid.associated_faces_begin(*iterV); iterFac != grid.associated_faces_end(*iterV); iterFac++ )
							{
								bool isFromFrac = false;

								for( auto const & facFrac : attFac )
								{
									if( *iterFac == facFrac )
									{
										isFromFrac = true;
									}
								}

								bool atRightSide = false;

								if( isFromFrac )
									atRightSide = true;

								if( !isFromFrac )
								{

									// check if on same side of edge where the normal points to: compute cosinus between vector of face center
									//  perpendicular to the edge
									// TODO FIXME

									vector3 facCenter = CalculateCenter( *iterFac, aaPos );

									vector3 perpendicu;

									DropAPerpendicular(perpendicu, facCenter, aaPos[edgeOfFrac->vertex(0)], aaPos[edgeOfFrac->vertex(1)]);

									vector3 tmpN;

									VecSubtract(tmpN, facCenter, perpendicu );

									VecNormalize(tmpN, tmpN);

									UG_LOG("Normale Boundary zum Face ist " << tmpN << std::endl);

									number cosBetwFracEdgAndDirection2Face = VecDot(tmpN, nrmEdg );

									UG_LOG("Cosinus Boundary zur Normalen ist " << cosBetwFracEdgAndDirection2Face << std::endl);

									if( cosBetwFracEdgAndDirection2Face > 0 )
									{
										UG_LOG("assuming boundary face to be on richt side" << std::endl);

										atRightSide = true;

#if ANSCHAULICH_ERZEUGE_SUDOS_ANHANG
												Vertex * otherFacCent = *grid.create<RegularVertex>();
												aaPos[otherFacCent] = facCenter;
												sh.assign_subset(otherFacCent, 5 );

												Vertex * pp = *grid.create<RegularVertex>();
												aaPos[pp] = perpendicu;
												sh.assign_subset(pp, 6 );

												sh.assign_subset(*iterFac,7);
#endif


									}
									else
									{
										UG_LOG("assuming boundary face to be on wrong side" << std::endl);
									}

								}

								if( atRightSide ) // atRightSide ) NOCH FALSCH TODO FIXME muss nur auf richtiger Seite sein
								{


									vector<Vertex*>& newVrts4Fac = aaVrtVecFace[ * iterFac ];

									IndexType vrtxFnd = 0;

									for(size_t indVrt = 0; indVrt < (*iterFac)->num_vertices();  indVrt++ )
									{
										Vertex* facVrt = (*iterFac)->vertex(indVrt);

										if(  facVrt == *iterV )
										{
											newVrts4Fac[ indVrt ] = newShiftEdgVrtx;
											vrtxFnd++;
										}
									}



									if( vrtxFnd <= 0 )
									{
										UG_THROW("vertex not found bnd!" << std::endl);
									}
									else if( vrtxFnd > 1 )
									{
										UG_THROW("vertex zu oft gefunden bnd " << vrtxFnd << std::endl );
									}
									else if ( vrtxFnd == 1 )
									{
//														UG_LOG("vertex found abgeschlossen" << std::endl);
									}
									else
									{
										UG_THROW("vertex finden bnd komisch " << std::endl);
									}
								}

							}
#else
							for( auto const & ifac : assoFaces )
							{
								bool isFromFrac = false;

								for( auto const & facFrac : attFac )
								{
									if( ifac == facFrac )
									{
										isFromFrac = true;
									}
								}

								bool atRightSide = false;

								if( isFromFrac )
									atRightSide = true;

								if( !isFromFrac )
								{

									// check if on same side of edge where the normal points to: compute cosinus between vector of face center
									//  perpendicular to the edge
									// TODO FIXME

									vector3 facCenter = CalculateCenter( ifac, aaPos );

									vector3 perpendicu;

									DropAPerpendicular(perpendicu, facCenter, aaPos[edgeOfFrac->vertex(0)], aaPos[edgeOfFrac->vertex(1)]);

									vector3 tmpN;

									VecSubtract(tmpN, facCenter, perpendicu );

									VecNormalize(tmpN, tmpN);

									UG_LOG("Normale Boundary zum Face ist " << tmpN << std::endl);

									number cosBetwFracEdgAndDirection2Face = VecDot(tmpN, nrmEdg );

									UG_LOG("Cosinus Boundary zur Normalen ist " << cosBetwFracEdgAndDirection2Face << std::endl);

									if( cosBetwFracEdgAndDirection2Face > 0 )
									{
										UG_LOG("assuming boundary face to be on richt side" << std::endl);

										atRightSide = true;

#if ANSCHAULICH_ERZEUGE_SUDOS_ANHANG
												Vertex * otherFacCent = *grid.create<RegularVertex>();
												aaPos[otherFacCent] = facCenter;
												sh.assign_subset(otherFacCent, 5 );

												Vertex * pp = *grid.create<RegularVertex>();
												aaPos[pp] = perpendicu;
												sh.assign_subset(pp, 6 );

												sh.assign_subset(*iterFac,7);
#endif


									}
									else
									{
										UG_LOG("assuming boundary face to be on wrong side" << std::endl);
									}

								}

								if( atRightSide ) // atRightSide ) NOCH FALSCH TODO FIXME muss nur auf richtiger Seite sein
								{


									vector<Vertex*>& newVrts4Fac = aaVrtVecFace[ ifac ];

									IndexType vrtxFnd = 0;

									for(size_t indVrt = 0; indVrt < (ifac)->num_vertices();  indVrt++ )
									{
										Vertex* facVrt = (ifac)->vertex(indVrt);

										if(  facVrt == *iterV )
										{
											newVrts4Fac[ indVrt ] = newShiftEdgVrtx;
											vrtxFnd++;
										}
									}



									if( vrtxFnd <= 0 )
									{
										UG_THROW("vertex not found bnd!" << std::endl);
									}
									else if( vrtxFnd > 1 )
									{
										UG_THROW("vertex zu oft gefunden bnd " << vrtxFnd << std::endl );
									}
									else if ( vrtxFnd == 1 )
									{
//														UG_LOG("vertex found abgeschlossen" << std::endl);
									}
									else
									{
										UG_THROW("vertex finden bnd komisch " << std::endl);
									}
								}
							}

#endif



						}
					}
				}

			}
			else // fractures are crossing at boundary even
			{

			}


			UG_LOG("END THIS BOUNDARY VERTEX" << std::endl);
		}

		dbg_vertizesPassiert++;

	}


//		// neue Vertizes in der Entfernung der Klüfte von den Klüften weg erzeugen,
//		// basierend auf den Normalen multipliziert mit der halben Kluftdicke
//		//für eine Kluft erstmal nur
//		// die neuen Kanten und Faces erzeugen, die alten falschen Kanten löschen und ebenso die alten Faces
//		// später auf mehr Klüfte ausdehnen, mit Problemstelle Kreuzung, aber erst, wenn eine Kluft funktioniert
//

//	return true;

	// jetzt Seb Sachen beinahe unverändert

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
		for(size_t i_edge = 0; i_edge < sf->num_edges(); ++i_edge)
		{
			Edge* e = grid.get_edge(sf, i_edge);

			if(sel.is_selected(e))
			{
				//	check the associated vertices through the volumes aaVrtVecVol attachment.
				//	If at least one has an associated new vertex and if no edge between the
				//	new vertices already exists, we'll create the new edge.
				size_t ind0 = i_edge;
				size_t ind1 = (i_edge + 1) % sf->num_edges();

				Vertex* nv0 = (aaVrtVecFace[sf])[ind0];
				Vertex* nv1 = (aaVrtVecFace[sf])[ind1];

				if(nv0 || nv1)
				{
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


	std::vector<Face * > newFaces;
	std::vector<int> subsOfNewFaces;

	//	iterate over all surrounding faces and create new vertices.
	//	Since faces are replaced on the fly, we have to take care with the iterator.
	for(FaceIterator iter_sf = sel.faces_begin(); iter_sf != sel.faces_end();)
	{
		Face* sf = *iter_sf;
		++iter_sf;

		std::vector<Vertex*> newVrts = aaVrtVecFace[sf];

		//	all new vertices have been assigned to newVrts.
		//	Note that if newVrts[i] == NULL, then we have to take the
		//	old vertex sf->vertex(i).
		//	now expand the fracture edges of sf to faces.
		for(size_t i_vrt = 0; i_vrt < sf->num_vertices(); ++i_vrt)
		{
			size_t iv1 = i_vrt;
			size_t iv2 = (i_vrt + 1) % sf->num_vertices();

			Edge* tEdge = grid.get_edge(sf->vertex(iv1), sf->vertex(iv2));

			if(tEdge)
			{
				if( aaMarkEdgeB[tEdge] )
				{
					Face* expFace = NULL;
					if(newVrts[iv1] && newVrts[iv2])
					{
						//	create a new quadrilateral
						expFace = *grid.create<Quadrilateral>(
										QuadrilateralDescriptor(sf->vertex(iv1), sf->vertex(iv2),
																newVrts[iv2], newVrts[iv1]));
					}
					else if(newVrts[iv1])
					{
							//	create a new triangle
						expFace = *grid.create<Triangle>(
										TriangleDescriptor(sf->vertex(iv1), sf->vertex(iv2), newVrts[iv1]));
					}
					else if(newVrts[iv2])
					{
						//	create a new triangle
						expFace = *grid.create<Triangle>(
										TriangleDescriptor(sf->vertex(iv1), sf->vertex(iv2), newVrts[iv2]));
					}
					else
					{

////						sh.assign_subset(*iter_sf,10);
////						sh.assign_subset(tEdge,11);
//
//						return true;
//						//	this code-block should never be entered. If it is entered then
//						//	we selected the wrong faces. This is would be a BUG!!!
//						//	remove the temporary attachments and throw an error
//
//						//	remove the temporary attachments
//#if OLD_PROFREITER_STUFF
//						grid.detach_from_vertices(aAdjMarker);
//						grid.detach_from_edges(aAdjMarker);
//#endif
//						grid.detach_from_vertices(aAdjMarkerVFP);
//						grid.detach_from_edges(aAdjMarkerB);
//
//						grid.detach_from_vertices( aAdjInfoAVVFT );
//						grid.detach_from_faces(attVrtVec);
//
//						grid.detach_from_vertices( aAdjInfoEdges );
//						grid.detach_from_vertices( aAdjInfoFaces );
//
//
//						// TODO FIXME auch die weiteren Marker und INfos, alle Attachments, detachen!!!!
//
//						throw(UGError("Implementation error in ExpandFractures2d Arte."));
					}

					// TODO FIXME selektion closen irgendwie, damit auch alle Randkanten zum subset gehoeren!!!

					if( expFace )
					{
						sh.assign_subset(expFace, fracInfosBySubset.at(sh.get_subset_index(tEdge)).newSubsetIndex);

						int subs = fracInfosBySubset.at(sh.get_subset_index(tEdge)).newSubsetIndex;

						subsOfNewFaces.push_back( subs );

						newFaces.push_back( expFace );
					}
				}
			}
		}




		//	now set up a new face descriptor and replace the face.
		if(fd.num_vertices() != sf->num_vertices())
			fd.set_num_vertices(sf->num_vertices());

		for(size_t i_vrt = 0; i_vrt < sf->num_vertices(); ++i_vrt)
		{
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

		if(!aaMarkEdgeB[e])
			grid.erase(e);
	}

	if( subsOfNewFaces.size() != newFaces.size() )
	{
		UG_THROW("andere zahl neue faces als subdoms " << std::endl);
	}

	IndexType nfn = 0;

	for( auto const & nf : newFaces )
	{
		for(size_t i_edge = 0; i_edge < nf->num_edges(); ++i_edge)
		{
			Edge* edg = grid.get_edge(nf, i_edge);

			sh.assign_subset( edg, subsOfNewFaces[nfn] );

		}

		for( size_t iVrt = 0; iVrt < nf->num_vertices(); iVrt++ )
		{
			Vertex * vrt = nf->vertex(iVrt);

			sh.assign_subset( vrt, subsOfNewFaces[nfn] );
		}

		nfn++;
	}

	// sollen die Boundary Edges zur boundary gehören, oder zur Kluft?
	// wie ist es mit den Knoten, sind die alle richtig zugewiesen bezüglich subset?

	// TODO FIXME HHHHHHHHHHHHHHHHHH
	// jetzt muss noch der Diamant erzeugt werden
	// Ziel: KluftInnen erzeugen


	//	remove the temporary attachments

#if OLD_PROFREITER_STUFF
	grid.detach_from_vertices(aAdjMarker);
	grid.detach_from_edges(aAdjMarker);
#endif
//	grid.detach_from_vertices(aAdjMarkerVFP);
//	grid.detach_from_edges(aAdjMarkerB);
//
//	grid.detach_from_vertices( aAdjInfoAVVFT );
//	grid.detach_from_faces(attVrtVec);
//
//	grid.detach_from_vertices( aAdjInfoEdges );
//	grid.detach_from_vertices(aAdjInfoFaces );


	grid.detach_from_vertices(aAdjMarkerVFP );
	grid.detach_from_edges(aAdjMarkerB);
	grid.detach_from_vertices( aAdjInfoAVVFT  );
	grid.detach_from_vertices( aAdjInfoEdges );
	grid.detach_from_vertices( aAdjInfoFaces );
	grid.detach_from_faces(attVrtVec);
	grid.detach_from_vertices( aAdjInfoVVEVM );

//	sel.clear();

//	return true;

	//  alles detachen, was noch attached ist, da ist einiges hinzu gekommen!


	// Keilstruktur entfernen und durch den gewünschten Diamanten ersetzen

	// new vertices which divide the fracture edges

//	AVertex aAdjVert; // used to know if an edge is frac edge and in the old faces
//	grid.attach_to_edges_dv( aAdjVert, nullptr );
//	grid.attach_to_faces_dv( aAdjVert, nullptr );
//	Grid::EdgeAttachmentAccessor<AVertex> aaEdgeVert( grid, aAdjVert );
//	Grid::FaceAttachmentAccessor<AVertex> aaFaceVert( grid, aAdjVert );

	for( auto const & cfi : vecCrossVrtInf )
	{

//		IndexType nuCroFra =  cfi.getNumbCrossFracs();
		CrossVertInf::FracTyp fracTyp =  cfi.getFracTyp();

		Vertex * crossPt = cfi.getCrossVertex();

		std::vector<Vertex * > shiftVrtcs = cfi.getVecShiftedVrts();

//		IndexType shiffsOrig = shiftVrtcs.size();
//
//		auto soc = shiffsOrig;

		VecEdge origFracEdg; // = cfi.getVecOrigFracEdges();

		VecEdge allAssoEdgCP;

		for( std::vector<Edge *>::iterator iterEdg = grid.associated_edges_begin(crossPt); iterEdg != grid.associated_edges_end(crossPt); iterEdg++ )
		{

			allAssoEdgCP.push_back(*iterEdg);

			bool hasShiftedVrtx = false;

			for( IndexType i = 0; i < 2; i++ )
			{
				Vertex * side = (*iterEdg)->vertex(i);

				for( auto const & vrt : shiftVrtcs )
				{
					if( side == vrt )
						hasShiftedVrtx = true;
				}

			}

			if( ! hasShiftedVrtx )
				origFracEdg.push_back(*iterEdg);
		}


//		UG_LOG("Shift Vectors Orig " <<  soc << std::endl);

//		// aim: old cross vertex first, shift vertex second, new established vertex third
//		using VrtxPair = std::pair<Vertex*, Vertex*>;
//
//		using DiamondVertexInfo = VertexFractureTriple<Edge*, Face*, Vertex *>;
//
//		using VecDiamondVertexInfo = std::vector<DiamondVertexInfo>;
//
//		VecDiamondVertexInfo vecDiamVrtInfo;

		//		if( nuCroFra == 3 )
		if( fracTyp == CrossVertInf::TEnd )
		{

			// TODO FIXME hier sind wir HHHHHHHHHHHHHHHHHHHH

		}
		//		else if( nuCroFra == 4 )
		else if( fracTyp == CrossVertInf::XCross )
		{
//			for( auto const & edg : origFracEdg )
			for( auto edg : origFracEdg )
			{

				static_assert( std::is_same< decltype( edg ), Edge * >::value  );
//				static_assert( std::is_same< const_cast<Edge*>(decltype( edg )), Edge * >::value  );

//				//Vertex* vrtSpliEd =
//				if( edg != nullptr )
//				SplitEdge<Vertex>( grid, static_cast<Edge*>(edg) );
////				else
////					UG_LOG("Nullptr getroffen " << std::endl);
//
////				size_t numSubs = sh.num_subsets();
////
////				sh.assign_subset(edg, numSubs );
//
//				//SplitEdge
//
//				return true;

			}

			// alle neuen faces entfernen, die an den Schnittknoten anhängen

			std::vector<Face *> assoFacCross;
			std::vector<IndexType> subdomList;

			for( std::vector<Face *>::iterator iterFac = grid.associated_faces_begin(crossPt); iterFac != grid.associated_faces_end(crossPt); iterFac++ )
			{
				assoFacCross.push_back(*iterFac);

				bool sudoAlreadyThere = false;

				IndexType sudo = sh.get_subset_index(*iterFac);

				for( auto sl : subdomList )
				{
					if( sl == sudo )
						sudoAlreadyThere = true;
				}

				if( !sudoAlreadyThere )
					subdomList.push_back(sudo);
			}

			// need to know the subsets of all faces!




			for( auto const & afc : assoFacCross )
			{
				// figure out respective original frac edge, need to know all edges and all vertices of this fracture

				// was sind die Ecken und Kanten dieses Faces? welche Ecke davon ist die originale Fracture edge?
				// was ist die subdomain der originalen Fracture edge? klar das ist die subdomain des faces

				auto subdom = sh.get_subset_index(afc);

				bool subdomFnd = false;

				for( auto const & sd : subdomList )
				{
					if( subdom == sd )
						subdomFnd = true;
				}

				if( ! subdomFnd )
					UG_THROW("SUBDOM NOT found" << std::endl);

//				for( auto const &  )
//				AssociatedEdgeIterator associated_edges_begin(Face* face);///< DO NOT INVOKE! Subject to change.
//						AssociatedEdgeIterator associated_edges_end(Face* face);///

//				std::vector<Edge*> edgesThisFac;
//
//				for( std::vector<Edge *>::iterator iterEdgF = grid.associated_edges_begin(afc);
//						iterEdgF != grid.associated_edges_end(afc); iterEdgF++ )
//				{
//					edgesThisFac.push_back(*iterEdgF);
//				}


			}

			if( subdomList.size() != 2 )
				UG_THROW("wieviele Subdoms um Kreuz herum? 	" << subdomList.size() << std::endl );

			// sort faces and corresponding edges, start with an oririnal fracture edge and select the next edge which
			// has a newly created shift vertex

			auto assoFacCrossCop = assoFacCross;

			// "durchdrehen"

//			Face * assoFacConsider = *(assoFacCrossCop.begin());
			Face * assoFacConsider = assoFacCrossCop[0];

			// soll einmal am Anfang festgelegt werden und dann bleiben
			Edge * assoFacEdgBeg2Fix = nullptr;
			// soll erst ganz am Schluss festgelegt werden
			Edge * assoFacEdgEnd2Fix = nullptr;

			// soll in jedem Lauf aktualisiert werden
			Edge * firstEdgeFac = assoFacEdgBeg2Fix;;
			Edge * secondEdgeFac = nullptr;

			bool atStartSort = true;

			using ExpandCrossFracInfo = VertexFractureTriple< std::pair<Edge*, Edge*>, Face*, std::pair<Vertex*, Vertex*> >;
			// Vertex nullptr wo original fracture, und shift vertex, wo Keilecke, die weg soll

			using VecExpandCrossFracInfo = std::vector<ExpandCrossFracInfo>;

			VecExpandCrossFracInfo vecExpCrossFI;

			bool boundAtShiftVrtEdg = true;

//			auto shiftVrtcsCop = shiftVrtcs;

//			UG_LOG("starting Rundlauf " << std::endl);

			IndexType dbg_rndl = 0;

			while( assoFacCrossCop.size() != 0 )
			{

//				UG_LOG("Debug Rundlauf " << dbg_rndl << std::endl);

				dbg_rndl++;

				secondEdgeFac = nullptr;

				Edge * fracEdge = nullptr;
				Edge * shiftVrtxEdg = nullptr;

//				IndexType fndEdgEnd = 0;

				std::vector<Edge*> edgesThisFac;

				edgesThisFac.clear();

				IndexType eoeo = edgesThisFac.size();

				auto eiei = eoeo;

//				UG_LOG("Edges this fac Orig Orig " << eiei <<  " " << dbg_rndl << std::endl);


//				UG_LOG("Debug Ecken " << std::endl);

				IndexType dbg_itEd = 0;

				for( std::vector<Edge *>::iterator iterEdgF = grid.associated_edges_begin(assoFacConsider);
						iterEdgF != grid.associated_edges_end(assoFacConsider); iterEdgF++ )
				{
					edgesThisFac.push_back(*iterEdgF);

//					UG_LOG("und noch eines dazu " << dbg_itEd << " " << dbg_rndl << std::endl);

					//IndexType sudos = sh.num_subsets();

					//sh.assign_subset(*iterEdgF,sudos);
				}

				IndexType effsOrig = edgesThisFac.size();

//				auto efeu = effsOrig;

//				UG_LOG("Edges this fac Orig " << efeu <<  dbg_rndl << std::endl);


				// figure out original fracture edge

				IndexType fndFracEdg = 0;

				for( auto const & etf : edgesThisFac )
				{
					for( auto const & ofe : origFracEdg )
					{
						if( etf == ofe )
						{
							fndFracEdg++;
							fracEdge = etf;
						}
					}
				}

//				UG_LOG("Debug Ofen	 " << std::endl);


				if( fracEdge == nullptr || fndFracEdg != 1 )
				{
					UG_LOG("Frac Orig Ecke nicht gefunden oder falsche Zahl " << fndFracEdg << std::endl );
//					return false;
					UG_THROW("Frac Orig Ecke nicht gefunden oder falsche Zahl " << fndFracEdg << std::endl );
				}


				// find expanded shift vertex

				Vertex * shiftVrtxFound = nullptr;
				IndexType fndVrt = 0;

//				IndexType suse = sh.num_subsets();

				//sh.assign_subset(crossPt,suse);

//				for( auto & sv : shiftVrtcsCop )
//				{
//					IndexType suseV = sh.num_subsets();
//					//sh.assign_subset(sv,suseV);
//				}

//				return true;

//				UG_LOG("Debug Entfernene	 " << std::endl);

				IndexType dbg_edgnum = 0;

				IndexType helpVarEdges = 0;

				IndexType effs = edgesThisFac.size();
				IndexType shiffs = shiftVrtcs.size();

//				UG_LOG("Edges this fac " <<  effs <<  dbg_rndl << std::endl);
//				UG_LOG("Shift Vectors " <<  shiffs <<  dbg_rndl << std::endl);


				for( auto const & etf : edgesThisFac )
				{

					if( helpVarEdges >= edgesThisFac.size() )
					{
						UG_LOG("Indexueberschreitung Edges" << std::endl);
						break;
					}

					helpVarEdges++;

					dbg_edgnum++;

					IndexType helpShiVaNum = 0;

					IndexType dbg_shiVe = 0;

//					for( Vertex * const & sv : shiftVrtcs )
					for( auto const & sv : shiftVrtcs )
					{

						if( helpShiVaNum >= shiftVrtcs.size() )
						{
							UG_LOG("Shift Vertex Index Verletzung " << std::endl);
							break;
						}

						helpShiVaNum++;

						dbg_shiVe++;

						for( IndexType i = 0; i < 2; i++ )
						{
//							if( ( etf->vertex(i) == crossPt && etf->vertex((i+1)%2) == sv )  || (etf->vertex((i+1)%2) == crossPt && etf->vertex(i) == sv ))
							if( etf->vertex(i) == crossPt && etf->vertex((i+1)%2) == sv )
							{
								shiftVrtxFound = sv;
								fndVrt++;
								shiftVrtxEdg = etf;

//								UG_LOG("Shift Vertex " << aaPos[shiftVrtxFound] << " " << dbg_edgnum << " " << dbg_shiVe << " " << dbg_rndl << std::endl);
//								UG_LOG("Cross Vertex " << aaPos[crossPt] << " " << dbg_edgnum << " " << dbg_shiVe <<  " " << dbg_rndl <<  std::endl );
//
//								UG_LOG("dbg edgenu shive " << dbg_edgnum << " " << dbg_shiVe <<  " " << dbg_rndl <<  std::endl);
							}
						}
					}
				}

//				UG_LOG("Debug Entfert durch	 " << std::endl);


				if( fndVrt != 1 || shiftVrtxFound == nullptr || shiftVrtxEdg == nullptr )
				{
					UG_LOG("shift vertex komische Zahl oder null " << fndVrt << std::endl);
//					return false;
					UG_THROW("shift vertex komische Zahl oder null " << fndVrt << std::endl);
				}

//				UG_LOG("Debug Entfert Text durch	 " << std::endl);


//				for( std::vector<Vertex*>::iterator itV = shiftVrtcsCop.begin(); itV != shiftVrtcsCop.end(); itV++ )
//				{
//					Vertex * vrt = *itV;
//
//					if( vrt == shiftVrtxFound )
//					{
//						shiftVrtcsCop.erase(itV);
//						break;
//					}
//				}
//
//				UG_LOG("Debug Rasieren durch	 " << std::endl);


				if( atStartSort )
				{
					assoFacEdgBeg2Fix = fracEdge;
					atStartSort = false;
				}

//				Edge * firstEdgeFac = fracEdge;
//				Edge * secondEdgeFac = shiftEdge;
				firstEdgeFac = fracEdge;
				secondEdgeFac = shiftVrtxEdg;


				Vertex * firstVrt = nullptr;
				Vertex * secondVrt = shiftVrtxFound;

				if( !boundAtShiftVrtEdg )
				{
					firstEdgeFac = shiftVrtxEdg;
					secondEdgeFac = fracEdge;

					firstVrt = shiftVrtxFound;
					secondVrt = nullptr;
				}

//				UG_LOG("Debug Paarbildung	 " << std::endl);


				std::pair<Edge*, Edge*> edgesFac( firstEdgeFac, secondEdgeFac );

				std::pair<Vertex*,Vertex*> vrtcsSF( firstVrt, secondVrt );

				ExpandCrossFracInfo startFacInf( edgesFac, assoFacConsider, vrtcsSF );

				vecExpCrossFI.push_back(startFacInf);

//				IndexType sui = sh.num_subsets();
//
//				sh.assign_subset(assoFacConsider,sui);

//				UG_LOG("Debug Paarbildung	Rasieren " << std::endl);

				IndexType dbg_it_er = 0;

				for( std::vector<Face*>::iterator itFac = assoFacCrossCop.begin(); itFac != assoFacCrossCop.end(); itFac++ )
				{
					Face * iFa = *itFac;

//					UG_LOG("interieren " << dbg_it_er << std::endl );

					dbg_it_er++;

//					UG_LOG("ifa " << iFa << std::endl );
//					UG_LOG("assoFac " << assoFacConsider << std::endl );

//					bool enthaltung = FaceContains( iFa, firstEdgeFac );
//
////					UG_LOG("Enthaltung " << std::endl);
//
//					bool entzwei = FaceContains(iFa, secondEdgeFac);
//
////					UG_LOG("Entzweiung " << std::endl);


					if( iFa == assoFacConsider && FaceContains( iFa, firstEdgeFac ) && FaceContains(iFa, secondEdgeFac) )
					{
//						UG_LOG("Erasieren " << std::endl);
						assoFacCrossCop.erase(itFac);
						break;
					}
				}

//				UG_LOG("Debug Paarbildung	Rasieren durch " << std::endl);


				if( assoFacCrossCop.size() == 0 )
				{
					if( secondEdgeFac != assoFacEdgBeg2Fix )
					{
						UG_LOG("Gesichter Diamant leer, aber keine Anfangsecke gefunden" << std::endl);
//						return false;
						UG_THROW("Gesichter Diamant leer, aber keine Anfangsecke gefunden" << std::endl);
					}
					else
					{
						assoFacEdgEnd2Fix = secondEdgeFac;

						break; // while loop zu Ende, raus aus dem while loop, den Rest nicht mehr machen, würde schief gehen zwingendermassen
					}

				}

				// figure out the next face

				firstEdgeFac = secondEdgeFac;
				secondEdgeFac = nullptr;

				IndexType nextFaceFound = 0;

				for( std::vector<Face*>::iterator itFac = assoFacCrossCop.begin(); itFac != assoFacCrossCop.end(); itFac++ )
				{
					Face * iFa = *itFac;

					if( FaceContains(iFa, firstEdgeFac ) )
					{
						nextFaceFound++;
					}
				}

				if( nextFaceFound != 1 )
				{
					UG_LOG("folgendes Gesicht in falscher Anztahl gefunden Diamant " << nextFaceFound << std::endl);
//					return false;
					UG_THROW("folgendes Gesicht in falscher Anztahl gefunden Diamant " << nextFaceFound << std::endl);
				}

				for( std::vector<Face*>::iterator itFac = assoFacCrossCop.begin(); itFac != assoFacCrossCop.end(); itFac++ )
				{
					Face * iFa = *itFac;

					if( FaceContains(iFa, firstEdgeFac ) )
					{
						assoFacConsider = iFa;
						break;
					}
				}


				boundAtShiftVrtEdg = ! boundAtShiftVrtEdg;
			}

			if( assoFacEdgEnd2Fix != assoFacEdgBeg2Fix || assoFacEdgEnd2Fix == nullptr || assoFacEdgBeg2Fix == nullptr )
			{
				UG_THROW("Anfang und Ende stimmen nicht ueberein " << std::endl);
//				return false;
				UG_THROW("Anfang und Ende stimmen nicht ueberein " << std::endl);
			}

//			if( shiftVrtcsCop.size() != 0 )
//			{
//				UG_LOG("Shift vertizes nicht alle gefinden " << std::endl);
//				return false;
//				UG_THROW("Shift vertizes nicht alle gefinden " << std::endl);
//			}

			if( assoFacCrossCop.size() != 0 )
			{
				UG_LOG("nicht alle asso facs gefunden " << std::endl);
//				return false;
				UG_THROW("nicht alle asso facs gefunden " << std::endl);
			}

			UG_LOG("Kreis des Diamanten X Fall geschlossen " << std::endl);

			// create new vertices and new edges and new faces, delete the old ones at end not to forget

//			IndexType vecfis = vecExpCrossFI.size();
//
//			UG_LOG("Punkte werden erzeugt " << vecfis << std::endl);
//
//			for( IndexType i = -1; i < 8; i = i + 2 )
//				UG_LOG("iiii " << i << std::endl);
//
//			for( IndexType indC = -1; indC < 8; indC = indC + 2 )
//			{
//
//				UG_LOG("Punkterzeugung Test " << indC << std::endl );
//
//			}
//
//			return true;

			int vecfis = vecExpCrossFI.size();

//			for( int indC = -1; indC < vecfis; indC = indC + 2 )
//				UG_LOG("Punkterzeugung " << indC << std::endl );
//

//			for( int indC = -1; indC < vecExpCrossFI.size(); indC = indC + 2 )
			for( int indC = -1; indC < vecfis-2; indC = indC + 2 )
			{

//				UG_LOG("Punkterzeugung X " << indC << std::endl );


				IndexType indBefore = ( indC + vecExpCrossFI.size() ) % vecExpCrossFI.size();
				IndexType indAfter = ( indC + 1 + vecExpCrossFI.size() ) % vecExpCrossFI.size();

				ExpandCrossFracInfo & expCFIBeforeFracEdg = vecExpCrossFI[ indBefore ];
				ExpandCrossFracInfo & expCFIAfterFracEdg = vecExpCrossFI[ indAfter ];

				// compute cross point

				std::pair<Edge*, Edge*> edgesCrossSegBefore = expCFIBeforeFracEdg.getEdge();
				std::pair<Edge*, Edge*> edgesCrossSegAfter = expCFIAfterFracEdg.getEdge();

				Edge * beforeShiftEdg = edgesCrossSegBefore.first;
				Edge * beforeFracEdg = edgesCrossSegBefore.second;
				Edge * afterFracEdg = edgesCrossSegAfter.first;
				Edge * afterShiftEdg = edgesCrossSegAfter.second;

				std::pair<Vertex*,Vertex*> vrtcsCrossSegBefore = expCFIBeforeFracEdg.getNormal();
				std::pair<Vertex*,Vertex*> vrtcsCrossSegAfter = expCFIAfterFracEdg.getNormal();

				Vertex * shiftBefore = vrtcsCrossSegBefore.first;
				Vertex * shiftAfter = vrtcsCrossSegAfter.second;

				// zur Zielsetzung
//				Vertex * toBeEstablishedCutEdgeVrtBefore = vrtcsCrossSegBefore.second;
//				Vertex * toBeEstablishedCutEdgeVrtAfter = vrtcsCrossSegAfter.first;

				if(    vrtcsCrossSegBefore.second != nullptr ||  vrtcsCrossSegAfter.first != nullptr
					|| 	shiftBefore == nullptr || shiftAfter == nullptr
//				if(    toBeEstablishedCutEdgeVrtBefore != nullptr ||  toBeEstablishedCutEdgeVrtAfter != nullptr
//					|| 	shiftBefore == nullptr || shiftAfter == nullptr
				  )
					UG_THROW("Nullpointer fehlen oder zu viel " << std::endl);

				if(    beforeFracEdg != afterFracEdg || beforeFracEdg == nullptr || afterFracEdg == nullptr
					|| beforeShiftEdg == nullptr || afterShiftEdg == nullptr
				  )
					UG_LOG("Ecken Nullpunkter " << std::endl);

				// determin cut point of line between the shift vertices and the frac edge

				Edge * fracEdge = beforeFracEdg; // muss gleich sein offenkundig afterFracEdge, ist getestet auch

				// Gerade bestimmen, die durch fracEdge bestimmt wird, und Gerade, die durch Verbindungsvektor shift Verzices bestimmt wird

				// figure out frac vertex which is the cross point

				IndexType fracEdgInd = -1;

				for( IndexType fiv = 0; fiv < 2; fiv++ )
				{
					if( fracEdge->vertex(fiv) == crossPt )
						fracEdgInd = fiv;
				}

				if( fracEdgInd < 0 || fracEdgInd > 1 )
					UG_THROW("cross point nicht Teil von fracture edge" << std::endl );

				Vertex * fracEdgEnd = fracEdge->vertex( ( fracEdgInd + 1 ) % 2 );

				vector3 posCrossPt = aaPos[ crossPt ];
				vector3 posFracEnd = aaPos[ fracEdgEnd ];

				vector3 posShiftBefore = aaPos[ shiftBefore ];
				vector3 posShiftAfter = aaPos[ shiftAfter ];

				vector3 directionFrac;

				VecSubtract(directionFrac, posFracEnd, posCrossPt );

				vector3 directionShiftBefore;

				VecSubtract(directionShiftBefore, posShiftBefore, posCrossPt);

				vector3 directionShiftAfter;

				VecSubtract(directionShiftAfter, posShiftAfter, posCrossPt );

				vector3 sumShift;

				VecAdd(sumShift, directionShiftBefore, directionShiftAfter);

				vector3 halfSumShift;

				VecScale(halfSumShift,sumShift,0.5);

				vector3 posNewVrtOnEdg;

				VecAdd( posNewVrtOnEdg, posCrossPt, halfSumShift );

				Vertex * newEdgVrtx = *grid.create<RegularVertex>();
				aaPos[newEdgVrtx] = posNewVrtOnEdg;

				IndexType sudoEdg = sh.get_subset_index(fracEdge);

				Face * faceBefore = expCFIBeforeFracEdg.getFace();
				Face * faceAfter = expCFIAfterFracEdg.getFace();

				IndexType sudoBefore = sh.get_subset_index(faceBefore);
				IndexType sudoAfter = sh.get_subset_index(faceAfter);

				if( sudoEdg != sudoBefore || sudoBefore != sudoAfter )
					UG_THROW("komische sudos vor Diamant " << std::endl);

				sh.assign_subset(newEdgVrtx, sudoEdg);

				UG_LOG("neuer Diamant Vorbereitungs Vertex " << posNewVrtOnEdg << std::endl);

				//				Vertex * toBeEstablishedCutEdgeVrtBefore = vrtcsCrossSegBefore.second;
				//				Vertex * toBeEstablishedCutEdgeVrtAfter = vrtcsCrossSegAfter.first;



				// gleich neue Faces erzeugen?
//				std::pair<Vertex*,Vertex*> vrtcsCrossSegBefore = expCFIBeforeFracEdg.getNormal();
//				std::pair<Vertex*,Vertex*> vrtcsCrossSegAfter = expCFIAfterFracEdg.getNormal();

				// insert the newly established vertex into the vertex info of the ExpandCrossFracInfo of the face
////				vrtcsCrossSegBefore.second = newEdgVrtx;
////				vrtcsCrossSegAfter.first = newEdgVrtx;
////
//				std::pair<Vertex*,Vertex*> vrtcsCrossSegBeforeNew( vrtcsCrossSegBefore.first, newEdgVrtx );
//				std::pair<Vertex*,Vertex*> vrtcsCrossSegAfterNew( newEdgVrtx, vrtcsCrossSegAfter.second );
//
//
//				expCFIBeforeFracEdg.setNewNormal( vrtcsCrossSegBeforeNew );
//				expCFIAfterFracEdg.setNewNormal( vrtcsCrossSegAfterNew );

				vrtcsCrossSegBefore.second = newEdgVrtx;
				vrtcsCrossSegAfter.first = newEdgVrtx;
				expCFIBeforeFracEdg.setNewNormal( vrtcsCrossSegBefore );
				expCFIAfterFracEdg.setNewNormal( vrtcsCrossSegAfter );


//				 vecExpCrossFI[ indBefore ].setNewEdgVertex(newEdgVrtx);
//				 vecExpCrossFI[ indAfter ].setNewEdgVertex(newEdgVrtx);

//				 DiamondVertexInfo dviBefore();
//				 DiamondVertexInfo dviAfter();

			}


			// create new edges and new faces


//			IndexType diamantSubsNum = sh.num_subsets()+1; // +1 notwendig? TODO FIXME
			IndexType diamantSubsNum = sh.num_subsets(); // +1 notwendig? TODO FIXME

			std::vector<Face*> newFracFaceVec = std::vector<Face*>();
			std::vector<Face*> newDiamFaceVec = std::vector<Face*>();


			boundAtShiftVrtEdg = true;
			atStartSort = true;


			for( auto const & ecf : vecExpCrossFI )
			{
				// get new vertex at the original fracture edge

				std::pair<Edge*, Edge*> edgesCrossSeg = ecf.getEdge();

				Face * facSeg = ecf.getFace();

				std::pair<Vertex*, Vertex*> vertcsCrossSeg = ecf.getNormal();

				std::pair<Vertex*, Vertex*> vertcsCrossSegWithNewV = ecf.getNewNormal();


				if( atStartSort || boundAtShiftVrtEdg )
				{
					if( vertcsCrossSeg.first != nullptr || vertcsCrossSeg.second == nullptr )
						UG_THROW("Verwechslung " << vertcsCrossSeg.first << " " << vertcsCrossSeg.second << std::endl);
				}

				atStartSort = false;

				Edge * fracEdge = boundAtShiftVrtEdg ? edgesCrossSeg.first : edgesCrossSeg.second;
				Edge * shiftVrtEdge = boundAtShiftVrtEdg ?  edgesCrossSeg.second : edgesCrossSeg.first;

				Vertex * fracVrtNew =  boundAtShiftVrtEdg ? vertcsCrossSegWithNewV.first : vertcsCrossSegWithNewV.second; // should be nullptr at first segment, to be assigned afterwards / shifted
				Vertex * shiftVrt = boundAtShiftVrtEdg ? vertcsCrossSeg.second : vertcsCrossSeg.first;
				Vertex * shiftVrtTest = boundAtShiftVrtEdg ? vertcsCrossSegWithNewV.second : vertcsCrossSegWithNewV.first;

				if( shiftVrtTest != shiftVrt )
					UG_THROW("Shift Vertex verloren gegangen " << std::endl);

				IndexType sudoFac = sh.get_subset_index(facSeg);

				IndexType sudoFracEdge = sh.get_subset_index(fracEdge);

				if( sudoFac != sudoFracEdge )
				{
					UG_THROW("subdoms frac edge und face nicht gleich " << std::endl);
				}

				// get all vertices of face, check if both known ones are contained, delete the face, create
				// the additional needed edge, and create new face with new vertex


				std::vector<Vertex*> vrtcsFace;
				// assign first old vertices, then replace
//				std::vector<Vertex*> vrtcsNewFaceFrac = vrtcsFace;
//				std::vector<Vertex*> vrtcsNewFaceDiam = vrtcsFace;

				//	all new vertices have been assigned to newVrts.
				//	Note that if newVrts[i] == NULL, then we have to take the
				//	old vertex sf->vertex(i).
				//	now expand the fracture edges of sf to faces.
				for(size_t iVrt = 0; iVrt < facSeg->num_vertices(); iVrt++ )
				{

					Vertex * vrt = facSeg->vertex(iVrt);
					vrtcsFace.push_back( vrt );
//					vrtcsNewFaceFrac.push_back( vrt );
//					vrtcsNewFaceDiam.push_back( vrt );
				}

				std::vector<Vertex*> vrtcsNewFaceFrac = vrtcsFace;
//				std::vector<Vertex*> vrtcsNewFaceDiam = vrtcsFace;


				// check if known vertices are contained

				IndexType fraVeNuInd = -1;
				IndexType shiftVeNuInd = -1;

				IndexType cntr = 0;

				for( auto const & vf : vrtcsFace )
				{
					if( vf == fracVrtNew )
					{
						fraVeNuInd = cntr;
						UG_THROW("wie kann man hierher kommen?" << std::endl);
					}

					if( vf == crossPt )
					{
						fraVeNuInd = cntr;
					}

//
					if( vf == shiftVrt )
						shiftVeNuInd = cntr;

					cntr++;
				}

				if( fraVeNuInd < 0 || shiftVeNuInd < 0 )
					UG_THROW("frac vertex oder shift vertex not contained " << std::endl);

				UG_LOG("neuer frac vertex " << fraVeNuInd << " " << shiftVeNuInd << std::endl );

				// replace vrtcs
				vrtcsNewFaceFrac[fraVeNuInd] = fracVrtNew;

				// compute shift of center vertex along frac edge

				// check subdom of frac vertex and check if in subdom List of X cross

				IndexType sudoOther = -1;

				IndexType foundSudoOther = 0;
				IndexType foundSudoFac = 0;

				for( auto const & sd : subdomList )
				{
					if( sd != sudoFac )
					{
						sudoOther = sd;
						foundSudoOther++;
					}
					else if( sd == sudoFac )
					{
						foundSudoFac++;
					}
					else
					{
						UG_THROW("Sudo not from frac and not from other?" << std::endl);
					}
				}

				if( foundSudoFac != 1 && foundSudoOther != 1 )
					UG_THROW("sudo zu oft oder zu selten gefunden " << std::endl);

				// establish new edges and new faces

				if( vrtcsNewFaceFrac.size() != 4 )
					UG_LOG("komische Groesse Gesicht " << std::endl);

				for( IndexType i = 0; i < vrtcsNewFaceFrac.size(); i++ )
				{
					if( vrtcsNewFaceFrac[i] == nullptr )
					{
						UG_THROW("null auf " << i << std::endl);
					}
//					else
//					{
//						UG_LOG("kein null auf " << i << std::endl );
//					}
				}


				UG_LOG("neue Gesichter ausserhalb Diamant Ziel " << std::endl);

//				int a = 1 + 2;
//
				Face * newFracFace =
					*grid.create<Quadrilateral>( QuadrilateralDescriptor( vrtcsNewFaceFrac[0], vrtcsNewFaceFrac[1],
																		  vrtcsNewFaceFrac[2], vrtcsNewFaceFrac[3]
												) );

				sh.assign_subset(newFracFace, sh.get_subset_index(facSeg) );
//				sh.assign_subset(newFracFace, diamantSubsNum ); testweise

				newFracFaceVec.push_back(newFracFace);

				boundAtShiftVrtEdg = ! boundAtShiftVrtEdg;

			}

			UG_LOG("neue Gesichter ausserhalb Diamant erzeugt " << std::endl);

			for( int indC = 0; indC < vecfis-1; indC = indC + 2 )
			{
				IndexType indBefore = ( indC + vecExpCrossFI.size() ) % vecExpCrossFI.size();
				IndexType indAfter = ( indC + 1 + vecExpCrossFI.size() ) % vecExpCrossFI.size();

				ExpandCrossFracInfo & expCFIBeforeFracEdg = vecExpCrossFI[ indBefore ];
				ExpandCrossFracInfo & expCFIAfterFracEdg = vecExpCrossFI[ indAfter ];

				Face * facBefore = expCFIBeforeFracEdg.getFace();
				Face * facAfter = expCFIAfterFracEdg.getFace();

				std::vector<Vertex*> vrtcsFaceBefore;
				std::vector<Vertex*> vrtcsFaceAfter;


				for(size_t iVrt = 0; iVrt < facBefore->num_vertices(); iVrt++ )
				{
					Vertex * vrt = facBefore->vertex(iVrt);
					vrtcsFaceBefore.push_back( vrt );
				}

				for(size_t iVrt = 0; iVrt < facAfter->num_vertices(); iVrt++ )
				{
					Vertex * vrt = facAfter->vertex(iVrt);
					vrtcsFaceAfter.push_back( vrt );
				}

				Vertex * newVrtBefore = expCFIBeforeFracEdg.getNewNormal().first;
				Vertex * shiftVrt = expCFIBeforeFracEdg.getNewNormal().second;
				Vertex * newVrtAfter = expCFIAfterFracEdg.getNewNormal().second;

				if(    expCFIBeforeFracEdg.getNewNormal().second != expCFIBeforeFracEdg.getNormal().second
					|| expCFIBeforeFracEdg.getNewNormal().second == nullptr
					|| expCFIBeforeFracEdg.getNewNormal().second != expCFIAfterFracEdg.getNewNormal().first
					|| expCFIAfterFracEdg.getNewNormal().first == nullptr
					|| expCFIAfterFracEdg.getNewNormal().first != expCFIAfterFracEdg.getNormal().first
				)
				{
					UG_THROW("Vektorchaos " << std::endl);
				}

				std::vector<Vertex *> vrtxSmallDiam;

				vrtxSmallDiam.push_back( crossPt );
				vrtxSmallDiam.push_back( newVrtBefore );
				vrtxSmallDiam.push_back( shiftVrt );
				vrtxSmallDiam.push_back( newVrtAfter );

				Face * newFracFace =
					*grid.create<Quadrilateral>( QuadrilateralDescriptor( vrtxSmallDiam[0], vrtxSmallDiam[1],
																		  vrtxSmallDiam[2], vrtxSmallDiam[3]
												) );

				sh.assign_subset(newFracFace, diamantSubsNum );

				newDiamFaceVec.push_back(newFracFace);



			}

//			sh.get_subset_name() XXXXX

//			auto sunam = sh.get_subset_name(subdomList[0]);
//
//			for( auto const & su : subdomList )
//			{
//
//			}

//			using SuNaTyp = decltype( sh.get_subset_name(0) );

			static_assert( std::is_same<char const *, decltype( sh.get_subset_name(subdomList[0]) ) >::value );

			std::string diamNam = std::string("diamant_") + std::string(const_cast<char*>( sh.get_subset_name( subdomList[0] ) ))
					              + std::string("_") + std::string(const_cast<char*>( sh.get_subset_name( subdomList[1] ) ));

			sh.set_subset_name(diamNam.c_str(),diamantSubsNum);

			// TODO FIXME in extra Funktion packen, vielfach aufgerufen in der Art!

			for( auto const & nF : newFracFaceVec )
			{
				for(size_t iEdge = 0; iEdge < nF->num_edges(); iEdge++ )
				{
					Edge* edg = grid.get_edge(nF, iEdge);

					sh.assign_subset( edg, sh.get_subset_index(nF) );

				}

				for( size_t iVrt = 0; iVrt < nF->num_vertices(); iVrt++ )
				{
					Vertex * vrt = nF->vertex(iVrt);

					sh.assign_subset( vrt, sh.get_subset_index(nF) );
				}

			}

			for( auto const & nF : newDiamFaceVec )
			{
				for(size_t iEdge = 0; iEdge < nF->num_edges(); iEdge++ )
				{
					Edge* edg = grid.get_edge(nF, iEdge);

					sh.assign_subset( edg, sh.get_subset_index(nF) );

				}

				for( size_t iVrt = 0; iVrt < nF->num_vertices(); iVrt++ )
				{
					Vertex * vrt = nF->vertex(iVrt);

					sh.assign_subset( vrt, sh.get_subset_index(nF) );
				}

			}


			// at end delete all fracture edges which are too long

			for( auto const & fdel : vecExpCrossFI )
			{
				Face * fac2BeDeleted = fdel.getFace();

				if( fac2BeDeleted != nullptr )
					grid.erase(fac2BeDeleted);
				else
					UG_THROW("hier fehlt ein Gesicht " << std::endl);
			}

//			IndexType subsNumNow = sh.num_subsets();
//
////			IndexType susu = subsNumNow;
//
////			UG_LOG("subs num " << susu << std::endl);
//			UG_LOG("subs num " << subsNumNow << std::endl);
//
			for( auto const & edg : allAssoEdgCP )
			{
//				Edge * e2D = oEdg;

				if( edg != nullptr )
				{
					UG_LOG("will erasieren " << edg << std::endl );
					grid.erase(edg);

//					sh.assign_subset( e2D, subsNumNow );
//					sh.assign_subset( e2D, subsNumNow );
				}
				else
				{
					UG_LOG("hier fehlt eine Ecke " << std::endl);
				}
			}

			UG_LOG("ALles erasiert " << std::endl);

//			for( auto & afc : assoFacCross )
//			{
//				grid.erase(afc);
//			}
//
//			// von den Edges, die vom Schnittknoten ausgehen, die anderen Endvertizes merken, dann edges löschen
//
//			std::vector<Edge *> assoEdgCross;
//			std::vector<Vertex *> endVertices;
//
//			for( std::vector<Edge *>::iterator iterEdg = grid.associated_edges_begin(crossPt); iterEdg != grid.associated_edges_end(crossPt); iterEdg++ )
//			{
//				assoEdgCross.push_back(*iterEdg);
//
//				for( size_t i = 0; i < 2; i++ )
//				{
//					Vertex * vrtEdgEnd = (*iterEdg)->vertex(i);
//
//					if( vrtEdgEnd != crossPt )
//					{
//						endVertices.push_back( vrtEdgEnd );
//					}
//				}
//			}
//
//#if 0
//			vector3 shiftPart;
//			VecAdd(shiftPart, fracVrtPos, shiftAlongEdgeTwo);
//
//			vector3 posNewVrt;
//			VecAdd( posNewVrt, shiftPart, shiftAlongEdgeOne);
//
//			UG_LOG("neuer Vertex Kreuzung " << posNewVrt << std::endl );
//
//			Vertex * newShiftVrtx = *grid.create<RegularVertex>();
//			aaPos[newShiftVrtx] = posNewVrt;
//
//#endif
//
//			for( auto & aec : assoEdgCross )
//			{
//				grid.erase(aec);
//			}


		}
	}

//	grid.detach_from_edges( aAdjVert );

	// die frac vertices entfernen noch

//	for( auto const & cfi : vecCrossVrtInf )
//	{
//		IndexType nuCroFra =  cfi.getNumbCrossFracs();
//
//		VecEdge origFracEdg = cfi.getVecOrigFracEdges();
//
//
//		if( nuCroFra == 3 )
//		{
//
//		}
//		else if( nuCroFra == 4 )
//		{
//			IndexType subsNumNow = sh.num_subsets();
//
//	//			IndexType susu = subsNumNow;
//
//	//			UG_LOG("subs num " << susu << std::endl);
//			UG_LOG("subs num " << subsNumNow << std::endl);
//
//			for( auto const & oEdg : origFracEdg )
//			{
//				Edge * e2D = oEdg;
//
//				if( e2D != nullptr )
//				{
//	//					grid.erase(edg2BeDel);
//					UG_LOG("will erasieren " << e2D << std::endl );
//
//					sh.assign_subset( e2D, subsNumNow );
//	//					sh.assign_subset( e2D, subsNumNow );
//				}
//				else
//				{
//						UG_LOG("hier fehlt eine Ecke " << std::endl);
//				}
//			}
//
//		}
//	}

	UG_LOG("zu Ende gekommen mit Arte 2D" << std::endl);

	return true;

	// ENDE NEUES ZEUG SELEKTION







#if OLD_PROFREITER_STUFF
	// TODO FIXME von diesem Loop kann man noch für oben die calculate crease normal lernen, vielleicht minimal abgewandelt, vielleicht exakt gleich

	//	a callback that returns true if the edge is a fracture edge, neues System
	AttachmentUnequal<Edge, Grid::EdgeAttachmentAccessor<ABool> > isFracEdgeB(aaMarkEdgeB, false);

	//	iterate over all surrounding faces and create new vertices.
	for(FaceIterator iter_sf = sel.faces_begin(); iter_sf != sel.faces_end(); ++iter_sf)
	{
		Face* sf = *iter_sf;

	//	check for each vertex whether it lies in the fracture
	//	(aaMarkVRT > 1 in this case)
	//	if so, we have to copy or create a vertex from/in aaVrtVec[vrt] which is
	//	associated with the crease normal on the side of sf.
		for(size_t i_vrt = 0; i_vrt < sf->num_vertices(); ++i_vrt)
		{
			Vertex* vrt = sf->vertex(i_vrt);
			if(aaMarkVRT[vrt] > 1)
			{
			//	calculate the normal on this side of the frac
				// TODO FIXME so eine Funktion brauchen wir vielleicht oben auch zur Vereinfachung des Codes!!!
				vector3 n_v2 = CalculateCreaseNormal(grid, sf, vrt, isFracEdgeB, aaPos);
				// das calculate crease normal scheint mir ein Schwachsinn zu sein
				// aber vielleicht doch nicht?

				UG_LOG("calculated crease normal v2: " << n_v2 << endl);

		}
	}
#endif







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
