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

#ifndef __H__LIB_GRID__SELECTION_UTIL_IMPL__
#define __H__LIB_GRID__SELECTION_UTIL_IMPL__

#include <vector>
#include <queue>

#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "common/util/metaprogramming_util.h"
#include "common/math/misc/math_constants.h"

namespace ug {

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	selection util methods

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
template <typename TAAPosVRT>
bool CalculateCenter(typename TAAPosVRT::ValueType& centerOut, Selector& sel, TAAPosVRT& aaPos)
{
	if(!sel.grid()){
		throw(UGError("No grid assigned to selector"));
	}
	
	Grid& grid = *sel.grid();
	
//	collect all vertices that are adjacent to selected elements
//	we have to make sure that each vertex is only counted once.
//	we do this by using grid::mark.
	grid.begin_marking();

//	std::vector<Vertex*> vrts;
//	vrts.assign(sel.vertices_begin(), sel.vertices_end());
//	grid.mark(sel.vertices_begin(), sel.vertices_end());

	VecSet(centerOut, 0);
	size_t n = 0;
	for(VertexIterator iter = sel.vertices_begin();
		iter != sel.vertices_end(); ++iter)
	{
		VecAdd(centerOut, centerOut, aaPos[*iter]);
		grid.mark(*iter);
		++n;
	}

	for(EdgeIterator iter = sel.edges_begin();
		iter != sel.edges_end(); ++iter)
	{
		Edge::ConstVertexArray vrts = (*iter)->vertices();
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked(vrts[i])){
				grid.mark(vrts[i]);
				VecAdd(centerOut, centerOut, aaPos[vrts[i]]);
				++n;
			}
		}
	}

	for(FaceIterator iter = sel.faces_begin();
		iter != sel.faces_end(); ++iter)
	{
		Face::ConstVertexArray vrts = (*iter)->vertices();
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked(vrts[i])){
				grid.mark(vrts[i]);
				VecAdd(centerOut, centerOut, aaPos[vrts[i]]);
				++n;
			}
		}
	}

	for(VolumeIterator iter = sel.volumes_begin();
		iter != sel.volumes_end(); ++iter)
	{
		Volume::ConstVertexArray vrts = (*iter)->vertices();
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked(vrts[i])){
				grid.mark(vrts[i]);
				VecAdd(centerOut, centerOut, aaPos[vrts[i]]);
				++n;
			}
		}
	}

	grid.end_marking();

	if(n > 0){
		VecScale(centerOut, centerOut, 1. / static_cast<number>(n));
		return true;
	}
	return false;
}

template <typename TAAPosVRT>
void TranslateSelection(Selector& sel, const typename TAAPosVRT::ValueType& offset,
						TAAPosVRT& aaPos)
{
	if(!sel.grid()){
		throw(UGError("No grid assigned to selector"));
	}

	Grid& grid = *sel.grid();

//	collect all vertices that are adjacent to selected elements
//	we have to make sure that each vertex is only counted once.
//	we do this by using grid::mark.
	grid.begin_marking();

	for(VertexIterator iter = sel.vertices_begin();
		iter != sel.vertices_end(); ++iter)
	{
		VecAdd(aaPos[*iter], aaPos[*iter], offset);
		grid.mark(*iter);
	}

	for(EdgeIterator iter = sel.edges_begin();
		iter != sel.edges_end(); ++iter)
	{
		Edge::ConstVertexArray vrts = (*iter)->vertices();
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked(vrts[i])){
				grid.mark(vrts[i]);
				VecAdd(aaPos[vrts[i]], aaPos[vrts[i]], offset);
			}
		}
	}

	for(FaceIterator iter = sel.faces_begin();
		iter != sel.faces_end(); ++iter)
	{
		Face::ConstVertexArray vrts = (*iter)->vertices();
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked(vrts[i])){
				grid.mark(vrts[i]);
				VecAdd(aaPos[vrts[i]], aaPos[vrts[i]], offset);
			}
		}
	}

	for(VolumeIterator iter = sel.volumes_begin();
		iter != sel.volumes_end(); ++iter)
	{
		Volume::ConstVertexArray vrts = (*iter)->vertices();
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked(vrts[i])){
				grid.mark(vrts[i]);
				VecAdd(aaPos[vrts[i]], aaPos[vrts[i]], offset);
			}
		}
	}

	grid.end_marking();
}

////////////////////////////////////////////////////////////////////////
//	InvertSelection
template <typename TSelector, typename TIterator>
void InvertSelection(TSelector& sel, TIterator begin, TIterator end)
{
	for(TIterator iter = begin; iter != end;){
	//	be careful - since iterator could be an iterator of the selector,
	//	we have to make sure, that we will not invalidate it.
		typename TIterator::value_type v = *iter;
		++iter;
		
		if(sel.is_selected(v))
			sel.deselect(v);
		else
			sel.select(v);
	}
}

////////////////////////////////////////////////////////////////////////
template <typename TElem, typename TIterator>
void
SelectAssociated(ISelector& sel, TIterator begin, TIterator end,
				 ISelector::status_t status)
{
	Grid* pGrid = sel.grid();
	if(!pGrid)
		return;

	Grid& grid = *pGrid;
	std::vector<TElem*> elems;
	for(TIterator iter = begin; iter != end; ++iter){
		CollectAssociated(elems, grid, *iter);
		for(size_t i = 0; i < elems.size(); ++i){
			sel.select(elems[i], status);
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedVertices
template <typename TSelector, typename TElemIterator>
void SelectAssociatedVertices(TSelector& sel, TElemIterator elemsBegin,
							  TElemIterator elemsEnd, ISelector::status_t status)
{
	while(elemsBegin != elemsEnd)
	{
		uint numVrts = (*elemsBegin)->num_vertices();
		for(uint i = 0; i < numVrts; ++i)
			sel.select((*elemsBegin)->vertex(i), status);
		elemsBegin++;
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedEdges
template <typename TSelector, typename TElemIterator>
void SelectAssociatedEdges(TSelector& sel, TElemIterator elemsBegin,
						   TElemIterator elemsEnd, ISelector::status_t status)
{
	Grid* pGrid = sel.grid();
	if(pGrid)
	{
		Grid& grid = *pGrid;
		std::vector<Edge*> vEdges;
		while(elemsBegin != elemsEnd)
		{
			CollectEdges(vEdges, grid, *elemsBegin);
			for(uint i = 0; i < vEdges.size(); ++i)
				sel.select(vEdges[i], status);
			elemsBegin++;
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedFaces
template <typename TSelector, typename TElemIterator>
void SelectAssociatedFaces(TSelector& sel, TElemIterator elemsBegin,
						   TElemIterator elemsEnd, ISelector::status_t status)
{
	Grid* pGrid = sel.grid();
	if(pGrid)
	{
		Grid& grid = *pGrid;
		std::vector<Face*> vFaces;
		while(elemsBegin != elemsEnd)
		{
			CollectFaces(vFaces, grid, *elemsBegin);
			for(uint i = 0; i < vFaces.size(); ++i)
				sel.select(vFaces[i], status);
			elemsBegin++;
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedFaces
template <typename TSelector, typename TElemIterator>
void SelectAssociatedVolumes(TSelector& sel, TElemIterator elemsBegin,
						     TElemIterator elemsEnd, ISelector::status_t status)
{
	SelectAssociated<Volume>(sel, elemsBegin, elemsEnd, status);
}


template <typename TElem, typename TSelector>
void AssignSelectionStateToSides(TSelector& sel, bool recursive)
{
	using TIter = typename TSelector::template traits<TElem>::level_iterator;
	using TSide = typename TElem::side;

	UG_ASSERT(sel.grid(), "A selector has to operate on a grid");

	Grid& g = *sel.grid();
	typename Grid::traits<TSide>::secure_container sides;

	for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl){
		for(TIter iter = sel.template begin<TElem>(lvl);
			iter != sel.template end<TElem>(lvl); ++iter)
		{
			TElem* e = *iter;
			ISelector::status_t elemStatus = sel.get_selection_status(e);
			g.associated_elements(sides, e);

			for(size_t i = 0; i < sides.size(); ++i){
				ISelector::status_t sideStatus = sel.get_selection_status(sides[i]);
				sel.select(sides[i], elemStatus | sideStatus);
			}
		}
	}

	if(recursive && TSide::HAS_SIDES){
		AssignSelectionStateToSides<TSide>(sel, recursive);
	}
}


template <typename TElemIterator>
void SelectBoundaryElements(ISelector& sel, TElemIterator elemsBegin,
						 TElemIterator elemsEnd)
{
	UG_ASSERT(sel.grid(), "A grid has to be associated with the selector");
	Grid& grid = *sel.grid();

	for(TElemIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		typename TElemIterator::value_type e = *iter;
		if(LiesOnBoundary(grid, e)){
			sel.select(e);
		}
	}
}

template <typename TElemIterator>
void SelectInnerElements(ISelector& sel, TElemIterator elemsBegin,
						 TElemIterator elemsEnd)
{
	UG_ASSERT(sel.grid(), "A grid has to be associated with the selector");
	Grid& grid = *sel.grid();

	for(TElemIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		typename TElemIterator::value_type e = *iter;
		if(!LiesOnBoundary(grid, e)){
			sel.select(e);
		}
	}
}



template <typename TSelector, typename TAAPos>
void ExtendSelectionInDirection(
		TSelector& sel,
        size_t extSize,
        const typename TAAPos::ValueType& dir,
        number minAngle,
        number maxAngle,
      	const TAAPos& aaPos,
		ISelector::status_t status)
{
	if(!sel.grid()){
		UG_LOG("ERROR in ExtendSelection: Selector has to be assigned to a grid.\n");
		return;
	}
	
	Grid& grid = *sel.grid();
	
//	first select associated elements of volumes, faces and edges.
//	then select associated elements of selected vertices.
//	do this extSize times.
//	elements that have already been processed are marked.
	
	grid.begin_marking();
	
//	perform iteration
	for(size_t extIters = 0; extIters < extSize; ++extIters)
	{
//TODO: speed-up by only calling SelectAssociatedGridObjects once before the loop.
//		During the loop only newly selected elements should be checked for associated elements.

	//	select associated elements
		SelectAssociatedGridObjects(sel, status);

	//	iterate over all selected vertices.
		for(size_t lvl = 0; lvl < sel.num_levels(); ++lvl){
			for(VertexIterator iter = sel.template begin<Vertex>(lvl);
				iter != sel.template end<Vertex>(lvl); ++iter)
			{
				Vertex* vrt = *iter;
			//	all marked vertices have already been processed.
				if(!grid.is_marked(vrt)){
					grid.mark(vrt);

				//	select associated volumes, faces and edges.
					for(auto asIter = grid.associated_edges_begin(vrt);
						asIter != grid.associated_edges_end(vrt); ++asIter)
					{
						if(CheckDirection(vrt, *asIter, aaPos, dir, minAngle, maxAngle))
							sel.select(*asIter, status);
					}

					for(auto asIter = grid.associated_faces_begin(vrt);
						asIter != grid.associated_faces_end(vrt); ++asIter)
					{
						if(CheckDirection(vrt, *asIter, aaPos, dir, minAngle, maxAngle))
							sel.select(*asIter, status);
					}

					for(auto asIter = grid.associated_volumes_begin(vrt);
						asIter != grid.associated_volumes_end(vrt); ++asIter)
					{
						if(CheckDirection(vrt, *asIter, aaPos, dir, minAngle, maxAngle))
							sel.select(*asIter, status);
					}
				}
			}
		}
	}
	
	grid.end_marking();
}

////////////////////////////////////////////////////////////////////////
template <typename TAAPos>
void SelectEdgesByDirection(
				Selector& sel,
				TAAPos& aaPos,
				const vector3& dir,
				number minDeviationAngle,
				number maxDeviationAngle,
				bool selectFlipped)
{

	UG_COND_THROW(!sel.grid(), "A grid has to be assigned to the given selector");

	Grid& g = *sel.grid();
	vector3 n;
	VecNormalize(n, dir);

	number maxDot = cos(deg_to_rad(minDeviationAngle));
	number minDot = cos(deg_to_rad(maxDeviationAngle));

	for(auto _feI = g.begin<Edge>(); _feI != g.end<Edge>(); ++_feI){
		Edge* e = *_feI;
		vector3 dir;
		VecSubtract(dir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		VecNormalize(dir, dir);
		number d = VecDot(dir, n);
		if((d >= minDot - SMALL && d <= maxDot + SMALL) ||
			(selectFlipped && (-d >= minDot - SMALL && -d <= maxDot + SMALL)))
		{
			sel.select(e);
		}
	}
}

////////////////////////////////////////////////////////////////////////
template <typename TAAPos>
void SelectSubsetEdgesByDirection(
				Selector& sel,
				SubsetHandler& sh,
				int subsetIndex,
				TAAPos& aaPos,
				const vector3& dir,
				number minDeviationAngle,
				number maxDeviationAngle,
				bool selectFlipped)
{

	UG_COND_THROW(!sel.grid(), "A grid has to be assigned to the given selector");

	vector3 n;
	VecNormalize(n, dir);

	number maxDot = cos(deg_to_rad(minDeviationAngle));
	number minDot = cos(deg_to_rad(maxDeviationAngle));

	for(auto _feI = sh.begin<Edge>(subsetIndex); _feI != sh.end<Edge>(subsetIndex); ++_feI){
		Edge* e = *_feI;
		vector3 dir;
		VecSubtract(dir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);
		VecNormalize(dir, dir);
		number d = VecDot(dir, n);
		if((d >= minDot - SMALL && d <= maxDot + SMALL) ||
			(selectFlipped && (-d >= minDot - SMALL && -d <= maxDot + SMALL)))
		{
			sel.select(e);
		}
	}
}


////////////////////////////////////////////////////////////////////////
template <typename TEdgeIterator>
void SelectCreaseEdges(ISelector& sel, TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
						number minAngle, APosition aPos,
						bool ignoreBoundaryEdges, ISelector::status_t state)
{
	if(!sel.grid())
		return;

	Grid& grid = *sel.grid();

//	get the position accessor
	if(!grid.has_vertex_attachment(aPos))
		return;

	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

//	we'll store face normals in those vectors:
	vector3 n[2];

//	associated faces are stored in this array
	Face* f[2];

//	all dot-products between normals lower than minDot mark a crease.
	number minDot = cos(minAngle * M_PI / 180.f);

//	iterate through the edges
	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
	{
		Edge* e = *iter;
		if(!(ignoreBoundaryEdges && IsBoundaryEdge2D(grid, e))){
		//	get the associated faces
		//	all edges that do not have exactly 2 associated edges
		//	are regarded as crease-edges
			if(GetAssociatedFaces(f, grid, e, 2) == 2){
			//	get the normals of the associated faces
				CalculateNormal(n[0], f[0], aaPos);
				CalculateNormal(n[1], f[1], aaPos);
			//	if the dot-product is lower than minDot, then the edge is a crease edge.
				if(VecDot(n[0], n[1]) < minDot)
					sel.select(e, state);
			}
			else{
				sel.select(e, state);
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////
template <typename TIter>
void SelectAreaBoundary(ISelector& sel, const TIter begin, const TIter end)
{
	using TElem = typename Pointer2Value<typename TIter::value_type>::type;
	using TSide = typename TElem::side;

	if(!sel.grid())
		return;

	Grid& grid = *sel.grid();

	grid.begin_marking();

	std::vector<TSide*> sides;
	TIter iter = begin;
	while(iter != end){
		TElem* elem = *iter;
		++iter;
		CollectAssociated(sides, grid, elem);
		for(size_t i = 0; i < sides.size(); ++i){
			TSide* side = sides[i];
			if(!grid.is_marked(side)){
			//	if the side was initially selected, it should stay that way
				if(!sel.is_selected(side)){
					grid.mark(side);
					sel.select(side);
				}
			}
			else{
			//	if the side is marked, then it is an inner side
				sel.deselect(side);
			}
		}
	}

	grid.end_marking();
}

////////////////////////////////////////////////////////////////////////
template <typename TIter>
void SelectInterfaceElements(ISelector& sel, ISubsetHandler& sh,
							 const TIter begin, const TIter end,
							 bool regardSelectedNbrsOnly)
{
	using TElem = typename Pointer2Value<typename TIter::value_type>::type;
	using TNbr = typename TElem::sideof;

	if(!TElem::CAN_BE_SIDE)
		return;

	if(!sel.grid())
		return;

	Grid& grid = *sel.grid();

	std::vector<TNbr*> nbrs;

	for(TIter iter = begin; iter != end;){
		TElem* elem = *iter;
		++iter;

		CollectAssociated(nbrs, grid, elem);

		int si = -2;
		for(size_t i = 0; i < nbrs.size(); ++i){
			if(!regardSelectedNbrsOnly || sel.is_selected(nbrs[i])){
				if(sh.get_subset_index(nbrs[i]) != si){
					if(si == -2)
						si = sh.get_subset_index(nbrs[i]);
					else{
					//	elem is an interface element
						sel.select(elem);
						break;
					}
				}
			}
		}
	}
}

template <typename TElem>
void SelectSubsetElements(ISelector& sel, ISubsetHandler& sh, int subsetIndex,
						  ISelector::status_t status)
{
	using TIter = typename GridObjectCollection::traits<TElem>::iterator;
	GridObjectCollection goc = sh.get_grid_objects_in_subset(subsetIndex);

	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
		for(TIter iter = goc.begin<TElem>(lvl); iter != goc.end<TElem>(lvl); ++iter)
			sel.select(*iter, status);
	}
}


template <typename TGeomObj, typename TAAPos>
bool SelectRegion(Selector& sel, const typename TAAPos::ValueType& p, TAAPos& aaPos,
			   	  typename Grid::traits<typename TGeomObj::side>::callback cbRegionBoundary)
{
	using TIter = typename Grid::traits<TGeomObj>::iterator;

	if(!sel.grid())
		return false;

	Grid& g = *sel.grid();

//	first try to find the element which contains p
	TGeomObj* startElem = nullptr;
	for(TIter iter = g.begin<TGeomObj>(); iter != g.end<TGeomObj>(); ++iter){
		if(ContainsPoint(*iter, p, aaPos)){
			startElem = *iter;
			break;
		}
	}

	if(!startElem)
		return false;

	sel.clear<TGeomObj>();
	sel.select(startElem);
	SelectionFill<TGeomObj>(sel, cbRegionBoundary);

	return true;
}

template <typename TAAPos>
void SelectShortPolychains(ISelector& sel, number maxLength, bool closedChainsOnly,
						   TAAPos aaPos)
{
	if(!sel.grid())
		return;
	Grid& grid = *sel.grid();

//	we'll collect all contained short polychains in this vector before deciding
//	to select them or not. If a polychain is already longer than maxLength
//	we won't add its edges to the vector
	std::vector<Edge*> curChain;
	std::queue<Edge*> nextEdges;
	Grid::edge_traits::secure_container edges;

	std::vector<Vertex*> junctionPoints;

	grid.begin_marking();

	for(EdgeIterator eiter = grid.begin<Edge>(); eiter != grid.end<Edge>(); ++eiter){
		if(grid.is_marked(*eiter))
			continue;

		bool curChainIsClosed = true;
		number curChainLength = 0;

		curChain.clear();
		junctionPoints.clear();
		
		nextEdges.push(*eiter);
		grid.mark(*eiter);

		while(!nextEdges.empty()){
			Edge* e = nextEdges.front();
			nextEdges.pop();
			
			curChainLength += EdgeLength(e, aaPos);
			if(curChainLength <= maxLength)
				curChain.push_back(e);

			for(size_t ivrt = 0; ivrt < 2; ++ivrt){
				Vertex* vrt = e->vertex(ivrt);
				grid.associated_elements(edges, vrt);
				if(edges.size() < 2)
					curChainIsClosed = false;
				else if(edges.size() == 2){
					for(size_t iedge = 0; iedge < 2; ++iedge){
						Edge* nextEdge = edges[iedge];
						if(!grid.is_marked(nextEdge)){
							grid.mark(nextEdge);
							nextEdges.push(nextEdge);
						}
					}
				}
				else{
					junctionPoints.push_back(vrt);
				}
			}
		}

		if((curChainLength <= maxLength)){
			if(closedChainsOnly && curChainIsClosed && !junctionPoints.empty()){
			//	count the number of associated edges of each junction-point
			//	in curChain. If one is associated with != 2 vertices the chain
			//	is considered as not closed 
				for(size_t ivrt = 0; ivrt < junctionPoints.size(); ++ivrt){
					Vertex* vrt = junctionPoints[ivrt];
					size_t numConnectedEdges = 0;
					for(size_t iedge = 0; iedge < curChain.size(); ++iedge){
						if(EdgeContains(curChain[iedge], vrt))
							++numConnectedEdges;
					}

					if(numConnectedEdges != 2){
						curChainIsClosed = false;
						break;
					}
				}

		 	}

		 	if(curChainIsClosed || !closedChainsOnly)
				sel.select(curChain.begin(), curChain.end());
		}
	}

	grid.end_marking();
}

template <typename TElem>
void SelectLinkedElements(ISelector& sel,
		  typename Grid::traits<TElem>::callback cbIsSelectable,
		  typename Grid::traits<typename TElem::side>::callback cbIsTraversable)
{
	using namespace std;
	using ElemIter = typename Grid::traits<TElem>::iterator;
	using Side = typename TElem::side;

	if(!sel.grid())
		return;
	Grid& grid = *sel.grid();
	
	queue<TElem*> qElems;
	
//	add all currently selected elements to the qElems queue
	GridObjectCollection goc = sel.get_grid_objects();
	for(size_t i = 0; i < goc.num_levels(); ++i){
		for(ElemIter iter = goc.begin<TElem>(i); iter != goc.end<TElem>(i); ++iter)
			qElems.push(*iter);
	}
	
	typename Grid::traits<TElem>::secure_container	nbrs;
	typename Grid::traits<Side>::secure_container	sides;
	
	while(!qElems.empty()){
		TElem* e = qElems.front();
		qElems.pop();
		
		grid.associated_elements(sides, e);
		
		for(size_t i_side = 0; i_side < sides.size(); ++i_side){
			Side* side = sides[i_side];
		//	if stopAtSelectedSides is active and if the side is selected,
		//	we won't traverse it.
			if(!cbIsTraversable(side))
				continue;
			
		//	get all neighboring elements of side. Check for each unselected,
		//	whether it lies on a boundary. If it does, select it and push it
		//	to the queue.
			grid.associated_elements(nbrs, side);
			
			for(size_t i_nbr = 0; i_nbr < nbrs.size(); ++i_nbr){
				TElem* nbr = nbrs[i_nbr];
				if(sel.is_selected(nbr))
					continue;
				if(!cbIsSelectable(nbr))
					continue;
				
			//	we found a new linked boundary element
				sel.select(nbr);
				qElems.push(nbr);
			}
		}
	}
}

template <typename TAAPosVRT>
UG_API
number FaceArea(ISelector& sel, TAAPosVRT& aaPos)
{
	number sum = 0.;

	if(!sel.grid()) {
		UG_WARNING("A grid has to be associated with the selector!");
		return sum;
	}

	using FaceIter = Grid::traits<Face>::const_iterator;
	GridObjectCollection goc = sel.get_grid_objects();

	for(size_t i = 0; i < goc.num_levels(); ++i)
		for(FaceIter iter = goc.begin<Face>(i); iter != goc.end<Face>(i); ++iter)
			sum += FaceArea(*iter, aaPos);

	return sum;
}

}// end of namespace

#endif
