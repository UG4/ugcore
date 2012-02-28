/*
 * mg_solver_util.cpp
 *
 *  Created on: 02.08.2011
 *      Author: andreasvogel
 */

#include "mg_solver_util.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// SurfaceToToplevelMap
////////////////////////////////////////////////////////////////////////////////

/// creates a mapping levIndex = vMap[surfIndex];
template <typename TElem>
static void CreateSurfaceToToplevelMap(std::vector<size_t>& vMap,
                                       ConstSmartPtr<SurfaceDoFDistribution> surfDD,
                                       ConstSmartPtr<LevelDoFDistribution> topDD)
{
//	type of element iterator
	typedef typename SurfaceDoFDistribution::traits<TElem>::const_iterator iter_type;

//	vector of indices
	std::vector<size_t> surfaceInd, levelInd;

	for(int si = 0; si < surfDD->num_subsets(); ++si)
	{
	//	iterators for subset
		iter_type iter = surfDD->begin<TElem>(si);
		iter_type iterEnd = surfDD->end<TElem>(si);

	//	loop all elements of type
		for( ; iter != iterEnd; ++iter)
		{
		//	get elem
			TElem* elem = *iter;

		//	extract all algebra indices for the element on surface
			surfDD->inner_algebra_indices(elem, surfaceInd);

		//	extract all algebra indices for the element on level
			topDD->inner_algebra_indices(elem, levelInd);

		//	check that index sets have same cardinality
			UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

		//	copy all elements of the vector
			for(size_t i = 0; i < surfaceInd.size(); ++i)
			{
			//	copy entries into level vector
				vMap[surfaceInd[i]] = levelInd[i];
			}
		}
	}
}

void CreateSurfaceToToplevelMap(std::vector<size_t>& vMap,
                                ConstSmartPtr<SurfaceDoFDistribution> surfDD,
                                ConstSmartPtr<LevelDoFDistribution> topDD)
{
//	check full refinement
	if(surfDD->num_indices() != topDD->num_indices())
		UG_THROW_FATAL("CreateSurfaceToToplevelMap: This function can only"
				" be applied to a full refined grid, where the surface is the "
				" top level.");

//	resize mapping
	vMap.resize(surfDD->num_indices(), 10000555);

// 	add dofs on elements
	if(surfDD->has_indices_on(VERTEX))
		CreateSurfaceToToplevelMap<VertexBase>(vMap, surfDD, topDD);
	if(surfDD->has_indices_on(EDGE))
		CreateSurfaceToToplevelMap<EdgeBase>(vMap, surfDD, topDD);
	if(surfDD->has_indices_on(FACE))
		CreateSurfaceToToplevelMap<Face>(vMap, surfDD, topDD);
	if(surfDD->has_indices_on(VOLUME))
		CreateSurfaceToToplevelMap<Volume>(vMap, surfDD, topDD);
}


template <typename TElemBase>
void SelectNonShadowsAdjacentToShadowsOnLevel(ISelector& sel,
                                              const SurfaceView& surfView,
                                              int level)
{
//	vectors for associated elements
	std::vector<VertexBase*> vAssVertex;
	std::vector<EdgeBase*> vAssEdge;
	std::vector<Face*> vAssFace;
	std::vector<Volume*> vAssVolume;

//	get grid
	Grid& grid = *sel.grid();

//	get multigrid
	MultiGrid& mg = *dynamic_cast<MultiGrid*>(&grid);

//	check multigrid
	if(&mg == NULL)
		UG_THROW_FATAL("SelectNonShadowsAdjacentToShadowsOnLevel: No "
					"Multigrid given, selection ob level not possible.");

//	check level
	if(level >= (int) mg.num_levels() || level < 0)
		UG_THROW_FATAL("SelectNonShadowsAdjacentToShadowsOnLevel: Requested "
						"level "<<level<<" does not exist in Multigrid.");

//	iterator type
	typename geometry_traits<TElemBase>::const_iterator iter, iterEnd;

	iterEnd = mg.end<TElemBase>(level);

//	loop all base elems
	for(iter = mg.begin<TElemBase>(level); iter != iterEnd; ++iter)
	{
	//	get element
		TElemBase* shadow = *iter;

	//	check if element is a shadow
		if(!surfView.is_shadowed(shadow)) continue;

	//	get adjacent elements
		CollectAssociated(vAssVertex, grid, shadow);
		CollectAssociated(vAssEdge, grid, shadow);
		CollectAssociated(vAssFace, grid, shadow);
		CollectAssociated(vAssVolume, grid, shadow);

	//	select associated elements
		for(size_t i = 0; i < vAssVertex.size(); ++i)
			if(surfView.is_contained(vAssVertex[i]))
				sel.select(vAssVertex[i]);
		for(size_t i = 0; i < vAssEdge.size(); ++i)
			if(surfView.is_contained(vAssEdge[i]))
				sel.select(vAssEdge[i]);
		for(size_t i = 0; i < vAssFace.size(); ++i)
			if(surfView.is_contained(vAssFace[i]))
				sel.select(vAssFace[i]);
		for(size_t i = 0; i < vAssVolume.size(); ++i)
			if(surfView.is_contained(vAssVolume[i]))
				sel.select(vAssVolume[i]);
	}
}



void SelectNonShadowsAdjacentToShadowsOnLevel(ISelector& sel,
                                              const SurfaceView& surfView,
                                              int level)
{
//	clear all marks
	sel.clear();

//	get grid
	Grid& grid = *sel.grid();

//	get multigrid
	MultiGrid& mg = *dynamic_cast<MultiGrid*>(&grid);

//	check multigrid
	if(&mg == NULL)
		UG_THROW_FATAL("SelectNonShadowsAdjacentToShadowsOnLevel: No "
				"Multigrid given, selection ob level not possible.");

//	if level does not exist in multigrid, no element can be selected and we're done
	if(level >= (int) mg.num_levels()) return;

//	check level
	if(level < 0)
		UG_THROW_FATAL("SelectNonShadowsAdjacentToShadowsOnLevel: Requested "
				"level "<<level<<" does not exist in Multigrid.");

//	No loops on Edges/Faces/Volumes are needed, since if an edge/face is in the
//	surface view, then also its vertices. Thus, the adjacend element is marked
//	already by the loop over VertexBase.

//	select elements
	SelectNonShadowsAdjacentToShadowsOnLevel<VertexBase>(sel, surfView, level);
}

template <typename TElemBase>
void SelectNonShadowsAdjacentToShadows(ISelector& sel, const SurfaceView& surfView)
{
//	vectors for associated elements
	std::vector<VertexBase*> vAssVertex;
	std::vector<EdgeBase*> vAssEdge;
	std::vector<Face*> vAssFace;
	std::vector<Volume*> vAssVolume;

//	get grid
	Grid& grid = *sel.grid();

	typename SurfaceView::template traits<TElemBase>::const_iterator iter, iterEnd;

//	loop all base elems
	iterEnd = surfView.end<TElemBase>();
	for(iter = surfView.begin<TElemBase>(); iter != iterEnd; ++iter)
	{
	//	get element
		TElemBase* shadow = *iter;

	//	check if element is a shadow
		if(!surfView.is_shadowed(shadow)) continue;

	//	get adjacent elemens
		CollectAssociated(vAssVertex, grid, shadow);
		CollectAssociated(vAssEdge, grid, shadow);
		CollectAssociated(vAssFace, grid, shadow);
		CollectAssociated(vAssVolume, grid, shadow);

	//	select associated elements
		for(size_t i = 0; i < vAssVertex.size(); ++i)
			if(surfView.is_contained(vAssVertex[i]))
				sel.select(vAssVertex[i]);
		for(size_t i = 0; i < vAssEdge.size(); ++i)
			if(surfView.is_contained(vAssEdge[i]))
				sel.select(vAssEdge[i]);
		for(size_t i = 0; i < vAssFace.size(); ++i)
			if(surfView.is_contained(vAssFace[i]))
				sel.select(vAssFace[i]);
		for(size_t i = 0; i < vAssVolume.size(); ++i)
			if(surfView.is_contained(vAssVolume[i]))
				sel.select(vAssVolume[i]);
	}
}

void SelectNonShadowsAdjacentToShadows(ISelector& sel, const SurfaceView& surfView)
{
//	clear all marks
	sel.clear();

//	get grid
	Grid& grid = *sel.grid();

//	note: the highest dimension of elements must not be loop, since there are
//		  no slaves of the highest dimension
	SelectNonShadowsAdjacentToShadows<VertexBase>(sel, surfView);

	if(grid.num<Face>() > 0 || grid.num<Volume>() > 0)
		SelectNonShadowsAdjacentToShadows<EdgeBase>(sel, surfView);

	if(grid.num<Volume>() > 0)
		SelectNonShadowsAdjacentToShadows<Face>(sel, surfView);
}


} // end namespace ug
