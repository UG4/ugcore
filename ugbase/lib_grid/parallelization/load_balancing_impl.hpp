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

#ifndef __H__UG__load_balancing_impl__
#define __H__UG__load_balancing_impl__

#include <vector>
#include "load_balancing.h"
#include "common/static_assert.h"
#include "lib_grid/algorithms/trees/kd_tree_static.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_distribution
///	@{

////////////////////////////////////////////////////////////////////////////////
template <class TElem, class TIterator, class TAAPos>
bool PartitionElements_RegularGrid(SubsetHandler& shOut,
								TIterator begin, TIterator end,
								int numCellsX, int numCellsY, int numCellsZ,
								TAAPos& aaPos,
								typename Grid::traits<TElem>::callback cbConsiderElem,
								int bucketSubset)
{
	using namespace ug;
	using namespace std;

	typedef typename TAAPos::ValueType vector_t;

//	make sure that the dimension is right
	UG_STATIC_ASSERT(TAAPos::ValueType::Size >= 2,
					TAPosition_has_to_be_at_least_two_dimensional);

	UG_ASSERT(shOut.grid(), "A grid has to be associated with the "
											"specified subset handler.");

	if((numCellsX < 1) || (numCellsY < 1) || (numCellsZ < 1)){
		UG_THROW("At least one cell per direction should be specified.");
	}

	Grid& grid = *shOut.grid();

//	collect all elements which shall be considered for partitioning.
//	All others are assigned to bucketSubset
	vector<TElem*> elems;
	for(TIterator iter = begin; iter != end; ++iter){
		if(cbConsiderElem(*iter))
			elems.push_back(*iter);
		else
			shOut.assign_subset(*iter, bucketSubset);
	}

//	calculate the bounding box
	vector_t tmin(0), tmax(0);
	{
		grid.begin_marking();
		vector<Vertex*>	vrts, associatedVrts;
		for(TIterator iter = begin; iter != end; ++iter){
			CollectAssociated(vrts, grid, *iter);
			for(size_t i = 0; i < vrts.size(); ++i){
				if(!grid.is_marked(vrts[i])){
					associatedVrts.push_back(vrts[i]);
					grid.mark(vrts[i]);
				}
			}
		}
		grid.end_marking();

		CalculateBoundingBox(tmin, tmax, associatedVrts.begin(),
							 associatedVrts.end(), aaPos);
	}

	vector3 min, max;
	VecCopy(min, tmin, 0);
	VecCopy(max, tmax, 0);

	number width = max.x() - min.x();
	number height = max.y() - min.y();
	number depth = max.z() - min.z();

	if((width < SMALL) && numCellsX > 1){
		UG_LOG("Can't execute PartitionElements_Rectangle: Geometry has no width.\n");
		return false;
	}
	if((height < SMALL) && numCellsY > 1){
		UG_LOG("Can't execute PartitionElements_Rectangle: Geometry has no height.\n");
		return false;
	}
	if((depth < SMALL) && numCellsZ > 1){
		UG_LOG("Can't execute PartitionElements_Rectangle: Geometry has no height.\n");
		return false;
	}


//	iterate over all elements and calculate the index at which they shall be
//	inserted into the subset handler
	for(typename vector<TElem*>::iterator iter = elems.begin();
		iter != elems.end(); ++iter)
	{
		TElem* elem = *iter;
		vector3 center;
		VecCopy(center, CalculateCenter(elem, aaPos), 0);

	//	get the cell index
		int xInd = 0;
		if(numCellsX > 1)
				xInd = (int)((number)numCellsX * (center.x() - min.x()) / width);
		int yInd = 0;
		if(numCellsY > 1)
			yInd = (int)((number)numCellsY * (center.y() - min.y()) / height);
		int zInd = 0;
		if(numCellsZ > 1)
			zInd = (int)((number)numCellsZ * (center.z() - min.z()) / depth);

	//	calculate the subset index (one could think of several alternatives here)
		int si = zInd * numCellsX * numCellsY + yInd * numCellsX + xInd;

	//	assign the subset
		shOut.assign_subset(elem, si);
	}
	return true;
}
}//	end of namespace

#endif







//	This is some old and unmaintained code.
//	Maybe it will be useful somewhen...
/**
 * this is a very simple and unoptimized algorithm!
 * Only intended for testing-purposes.
 *
 * It should work for faces
 * (call PartitionElementsByRepeatedIntersection<Face, 2>(...))
 *
 * and for volumes
 * (call PartitionElementsByRepeatedIntersection<Volume, 3>(...))
 *
 * probably even for edges. It won't however work for vertices.
 *
 * The resulting parts are not necessarily connected.
 * The size of the cut between to parts can be quite large.
 * The volume to boundary relation is not optimized.
 *
 * The algorithm simply divides the grid by repeated orthogonal cuts along
 * the main axes. The leafs of the resulting kd-tree are then assigned to
 * the different subsets.
 * \param shOut 	Subset Handler
 * \param grid		Grid
 * \param numSubsets has to be chosen as a power of 2.
 * \param aVrtPos	position attachment
 *
 * This algorithm may possibly fail. However, if the grid is somehow aligned
 * to the main axes and there is a similar amount of elements in each direction
 * (no anisotropy), it should work well.
 */
/*
template <class TElem, int IDimension>
bool PartitionElementsByRepeatedIntersectionKD(ug::SubsetHandler& shOut,
										ug::Grid& grid,
										int numSubsets,
										ug::APosition& aVrtPos)
{
//	TODO: move implementation to separate ..._impl.hpp file.
	using namespace ug;
	using namespace std;

	typedef KDTreeStatic<APosition, IDimension> KDTree;
	typedef typename geometry_traits<TElem>::iterator ElemIterator;

//	get the max-size from the number of subsets that shall be constructed
	int treeDepth = 0;
	int numLeafs = 1;
	while(numLeafs < numSubsets)
	{
		treeDepth++;
		numLeafs *= 2;
	}

	if(numLeafs != numSubsets)
	{
		LOG("PROBLEM in PartitionElementsByRepeatedIntersection: ");
		LOG("numSubsets has to be chosen as a power of 2)\n");
		return false;
	}

//	fill a kd-tree with the vertices. This will help to find the elem-partition later on.
	KDTree kdTree;
	kdTree.create_from_grid(grid, grid.vertices_begin(), grid.vertices_end(),
						aVrtPos, treeDepth, 1, KDSD_CIRCULAR);

//	get the leafs
	vector<typename KDTree::Node*> vLeafs;
	kdTree.get_leafs(vLeafs);

	if(vLeafs.size() != numLeafs)
	{
	//	the algorithm failed, since it is not suited to deal with the
	//	given grid and the given parameters.
	//	You could try to distribute your grid onto less processes.
	//	The problem is, that the amount of vertices on both sides of a cut
	//	is not the same. This could be improved, by not choosing the cut-plane
	//	in the kd-tree in a geometric way, but by choosing it so, that
	//	the amount of vertices on both sides equals.
	//	To avoid this problem you could try to rotate your grid. Try to align it
	//	with the x, y and z axis - or write a better algorithm! This one is only
	//	suited for simple test-problems!!!
		LOG("PROBLEM in PartitionElementsByRepeatedIntersection: ");
		LOG("could not find the correct intersections.\n");
		LOG("This can happen depending on the properties of the grid.\n");
		LOG("Search the source-code for more informations.\n");
		return false;
	}

	shOut.clear();

//	we need a vertex-selector
	VertexSelector sel(grid);

//	iterate through the leafs. They are sorted in a special way.
	for(uint i = 0; i < vLeafs.size(); ++i)
	{
	//	select all vertices of the leaf.
		sel.select(vLeafs[i]->m_pvVertices->begin(), vLeafs[i]->m_pvVertices->end());

	//	now iterate through the elements of the grid.
		for(ElemIterator iter = grid.begin<TElem>(); iter != grid.end<TElem>(); ++iter)
		{
		//	push all elements that do not reference unselected vertices into the i-th subset.
		//	all elements will be assigned to a subset this way, since all vertices are
		//	selected in the last iteration.
			TElem* ele = *iter;
		//	make sure, that the element is not already assigned to a subset.
			if(shOut.get_subset_index(ele) == -1)
			{
				bool bAllSelected = true;
				for(int j = 0; j < ele->num_vertices(); ++j)
				{
					if(!sel.is_selected(ele->vertex(j)))
					{
						bAllSelected = false;
						break;
					}
				}

			//	if all vertices were selected, we'll push the element to the subset.
				if(bAllSelected)
					shOut.assign_subset(ele, i);
			}
		}
	}

//TODO: exchange elements between subsets until all subsets meet a certain criterion.

	return true;
}
*/
