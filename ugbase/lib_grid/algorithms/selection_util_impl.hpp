// created by Sebastian Reiter
// y10 m05 d18
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SELECTION_UTIL_IMPL__
#define __H__LIB_GRID__SELECTION_UTIL_IMPL__

#include <vector>
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

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
		throw(UGFatalError("No grid assigned to selector"));
	}
	
	Grid& grid = *sel.grid();
	
//	collect all vertices that are adjacent to selected elements
//	we have to make sure that each vertex is only counted once.
//	we do this by using grid::mark.
	grid.begin_marking();

	std::vector<VertexBase*> vrts;
	vrts.assign(sel.vertices_begin(), sel.vertices_end());
	grid.mark(sel.vertices_begin(), sel.vertices_end());

	for(EdgeBaseIterator iter = sel.edges_begin();
		iter != sel.edges_end(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked((*iter)->vertex(i))){
				grid.mark((*iter)->vertex(i));
				vrts.push_back((*iter)->vertex(i));
			}
		}
	}

	for(FaceIterator iter = sel.faces_begin();
		iter != sel.faces_end(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked((*iter)->vertex(i))){
				grid.mark((*iter)->vertex(i));
				vrts.push_back((*iter)->vertex(i));
			}
		}
	}

	for(VolumeIterator iter = sel.volumes_begin();
		iter != sel.volumes_end(); ++iter)
	{
		for(size_t i = 0; i < (*iter)->num_vertices(); ++i){
			if(!grid.is_marked((*iter)->vertex(i))){
				grid.mark((*iter)->vertex(i));
				vrts.push_back((*iter)->vertex(i));
			}
		}
	}

	grid.end_marking();

	if(vrts.size() > 0){
		centerOut = CalculateCenter(vrts.begin(), vrts.end(), aaPos);
		return true;
	}
	return false;
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

}// end of namespace

#endif
