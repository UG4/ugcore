// created by Sebastian Reiter
// y10 m05 d18
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SELECTION_UTIL_IMPL__
#define __H__LIB_GRID__SELECTION_UTIL_IMPL__

#include <vector>
#include <queue>
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "common/util/metaprogramming_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	selection util methods

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
template <class TAAPosVRT>
bool CalculateCenter(typename TAAPosVRT::ValueType& centerOut,
					 Selector& sel, TAAPosVRT& aaPos)
{
	if(!sel.grid()){
		throw(UGError("No grid assigned to selector"));
	}
	
	Grid& grid = *sel.grid();
	
//	collect all vertices that are adjacent to selected elements
//	we have to make sure that each vertex is only counted once.
//	we do this by using grid::mark.
	grid.begin_marking();

//	std::vector<VertexBase*> vrts;
//	vrts.assign(sel.vertices_begin(), sel.vertices_end());
//	grid.mark(sel.vertices_begin(), sel.vertices_end());

	VecSet(centerOut, 0);
	size_t n = 0;
	for(VertexBaseIterator iter = sel.vertices_begin();
		iter != sel.vertices_end(); ++iter)
	{
		VecAdd(centerOut, centerOut, aaPos[*iter]);
		grid.mark(*iter);
		++n;
	}

	for(EdgeBaseIterator iter = sel.edges_begin();
		iter != sel.edges_end(); ++iter)
	{
		EdgeBase::ConstVertexArray vrts = (*iter)->vertices();
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
		VecScale(centerOut, centerOut, 1. / (number)n);
		return true;
	}
	return false;
}

template <class TAAPosVRT>
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

	for(VertexBaseIterator iter = sel.vertices_begin();
		iter != sel.vertices_end(); ++iter)
	{
		VecAdd(aaPos[*iter], aaPos[*iter], offset);
		grid.mark(*iter);
	}

	for(EdgeBaseIterator iter = sel.edges_begin();
		iter != sel.edges_end(); ++iter)
	{
		EdgeBase::ConstVertexArray vrts = (*iter)->vertices();
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
template <class TSelector, class TIterator>
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
template <class TElem, class TIterator>
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
template <class TSelector, class TElemIterator>
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
template <class TSelector, class TElemIterator>
void SelectAssociatedEdges(TSelector& sel, TElemIterator elemsBegin,
						   TElemIterator elemsEnd, ISelector::status_t status)
{
	Grid* pGrid = sel.grid();
	if(pGrid)
	{
		Grid& grid = *pGrid;
		std::vector<EdgeBase*> vEdges;
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
template <class TSelector, class TElemIterator>
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
template <class TSelector, class TElemIterator>
void SelectAssociatedVolumes(TSelector& sel, TElemIterator elemsBegin,
						     TElemIterator elemsEnd, ISelector::status_t status)
{
	SelectAssociated<Volume>(sel, elemsBegin, elemsEnd, status);
}


template <class TElem, class TSelector>
void AssignSelectionStateToSides(TSelector& sel, bool recursive)
{
	typedef typename TSelector::template traits<TElem>::level_iterator TIter;
	typedef typename TElem::side TSide;

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


template <class TElemIterator>
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

template <class TElemIterator>
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

////////////////////////////////////////////////////////////////////////
template <class TIter>
void SelectAreaBoundary(ISelector& sel, const TIter begin, const TIter end)
{
	typedef typename Pointer2Value<typename TIter::value_type>::type	TElem;
	typedef typename TElem::side										TSide;

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
template <class TIter>
void SelectInterfaceElements(ISelector& sel, ISubsetHandler& sh,
							 const TIter begin, const TIter end,
							 bool regardSelectedNbrsOnly)
{
	typedef typename Pointer2Value<typename TIter::value_type>::type	TElem;
	typedef typename TElem::sideof										TNbr;

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

template <class TElem>
void SelectSubsetElements(ISelector& sel, ISubsetHandler& sh, int subsetIndex,
						  ISelector::status_t status)
{
	typedef typename GeometricObjectCollection::traits<TElem>::iterator	TIter;
	GeometricObjectCollection goc = sh.get_geometric_objects_in_subset(subsetIndex);

	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
		for(TIter iter = goc.begin<TElem>(lvl); iter != goc.end<TElem>(lvl); ++iter)
			sel.select(*iter, status);
	}
}


template <class TGeomObj, class TAAPos>
bool SelectRegion(Selector& sel, const typename TAAPos::ValueType& p, TAAPos& aaPos,
			   	  typename Grid::traits<typename TGeomObj::side>::callback cbRegionBoundary)
{
	typedef typename Grid::traits<TGeomObj>::iterator	TIter;

	if(!sel.grid())
		return false;

	Grid& g = *sel.grid();

//	first try to find the element which contains p
	TGeomObj* startElem = NULL;
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


template <class TElem>
void SelectLinkedElements(ISelector& sel,
		  typename Grid::traits<TElem>::callback cbIsSelectable,
		  typename Grid::traits<typename TElem::side>::callback cbIsTraversable)
{
	using namespace std;
	typedef typename Grid::traits<TElem>::iterator	ElemIter;
	typedef typename TElem::side Side;

	if(!sel.grid())
		return;
	Grid& grid = *sel.grid();
	
	queue<TElem*> qElems;
	
//	add all currently selected elements to the qElems queue
	GeometricObjectCollection goc = sel.get_geometric_objects();
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

template <class TAAPosVRT>
UG_API
number FaceArea(ISelector& sel, TAAPosVRT& aaPos)
{
	number sum = 0.;

	if(!sel.grid()) {
		UG_WARNING("A grid has to be associated with the selector!");
		return sum;
	}

	typedef Grid::traits<Face>::const_iterator FaceIter;
	GeometricObjectCollection goc = sel.get_geometric_objects();

	for(size_t i = 0; i < goc.num_levels(); ++i)
		for(FaceIter iter = goc.begin<Face>(i); iter != goc.end<Face>(i); ++iter)
			sum += FaceArea(*iter, aaPos);

	return sum;
}

}// end of namespace

#endif
