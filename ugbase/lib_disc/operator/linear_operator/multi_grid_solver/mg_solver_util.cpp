/*
 * mg_solver_util.cpp
 *
 *  Created on: 02.08.2011
 *      Author: andreasvogel
 */

#include "mg_solver_util.h"

namespace ug{

template <typename TElemBase>
bool SelectNonShadowsAdjacentToShadowsOnLevel(ISelector& sel,
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
	{
		UG_LOG("ERROR in SelectNonShadowsAdjacentToShadowsOnLevel: No "
				"Multigrid given, selection ob level not possible.\n");
		return false;
	}

//	check level
	if(level >= (int) mg.num_levels() || level < 0)
	{
		UG_LOG("ERROR in SelectNonShadowsAdjacentToShadowsOnLevel: Requested "
				"level "<<level<<" does not exist in Multigrid.\n");
		return false;
	}

//	loop all subsets
	for(int si = 0; si < surfView.num_subsets(); ++si)
	{
	//	iterator type
		typename geometry_traits<TElemBase>::const_iterator iter, iterEnd;
		iterEnd = surfView.end<TElemBase>(si);

	//	loop all base elems
		for(iter = surfView.begin<TElemBase>(si); iter != iterEnd; ++iter)
		{
		//	get element
			TElemBase* elem = *iter;

		//	check if element is a shadow
			if(surfView.shadows(elem))
			{
			//	get shadow
				GeometricObject* shadow = surfView.get_parent(elem);

			//	check if this is the correct level
				if(mg.get_level(shadow) != level) continue;

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

	}

//	we're done
	return true;
}


bool SelectNonShadowsAdjacentToShadowsOnLevel(ISelector& sel,
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
	{
		UG_LOG("ERROR in SelectNonShadowsAdjacentToShadowsOnLevel: No "
				"Multigrid given, selection ob level not possible.\n");
		return false;
	}

//	if level does not exist in multigrid, no element can be selected and we're done
	if(level >= (int) mg.num_levels())
		return true;

//	check level
	if(level < 0)
	{
		UG_LOG("ERROR in SelectNonShadowsAdjacentToShadowsOnLevel: Requested "
				"level "<<level<<" does not exist in Multigrid.\n");
		return false;
	}

//	No loops on Edges/Faces/Volumes are needed, since if an edge/face is in the
//	surface view, then also its vertices. Thus, the adjacend element is marked
//	already by the loop over VertexBase.

//	select elements
	return SelectNonShadowsAdjacentToShadowsOnLevel<VertexBase>(sel, surfView, level);
}

template <typename TElemBase>
bool SelectNonShadowsAdjacentToShadows(ISelector& sel, const SurfaceView& surfView)
{
//	vectors for associated elements
	std::vector<VertexBase*> vAssVertex;
	std::vector<EdgeBase*> vAssEdge;
	std::vector<Face*> vAssFace;
	std::vector<Volume*> vAssVolume;

//	get grid
	Grid& grid = *sel.grid();

//	loop all subsets
	for(int si = 0; si < surfView.num_subsets(); ++si)
	{
	//	iterator type
		typename geometry_traits<TElemBase>::const_iterator iter, iterEnd;
		iterEnd = surfView.end<TElemBase>(si);

	//	loop all base elems
		for(iter = surfView.begin<TElemBase>(si); iter != iterEnd; ++iter)
		{
		//	get element
			TElemBase* elem = *iter;

		//	check if element is a shadow
			if(surfView.shadows(elem))
			{
			//	get shadow
				GeometricObject* shadow = surfView.get_parent(elem);

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

	}

//	we're done
	return true;
}

bool SelectNonShadowsAdjacentToShadows(ISelector& sel, const SurfaceView& surfView)
{
//	clear all marks
	sel.clear();

//	get grid
	Grid& grid = *sel.grid();

//	select elements
	bool bRes = true;

//	note: the highest dimension of elements must not be loop, since there are
//		  no slaves of the highest dimension
	bRes &= SelectNonShadowsAdjacentToShadows<VertexBase>(sel, surfView);

	if(grid.num<Face>() > 0 || grid.num<Volume>() > 0)
		bRes &= SelectNonShadowsAdjacentToShadows<EdgeBase>(sel, surfView);

	if(grid.num<Volume>() > 0)
		bRes &= SelectNonShadowsAdjacentToShadows<Face>(sel, surfView);

//	we're done
	return bRes;
}


} // end namespace ug
