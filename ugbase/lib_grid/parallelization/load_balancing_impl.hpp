// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 21.04.2011 (m,d,y)

#ifndef __H__UG__load_balancing_impl__
#define __H__UG__load_balancing_impl__

#include <vector>
#include "load_balancing.h"
#include "common/static_assert.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/trees/kd_tree_static.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_distribution
///	@{

////////////////////////////////////////////////////////////////////////
//	some classes and methods that help sorting certain elements by position
template <class TElem, class TPos>
class ElemWithPosition
{
	public:
		TElem*	pElem;
		TPos	pos;
};

/// if two values are identical, we'll return the cmp of the next dimension.
template <class TElem, int ICurDim, class TPos>
bool cmp(const ElemWithPosition<TElem, TPos>& ele1,
		 const ElemWithPosition<TElem, TPos>& ele2)
{
	for(size_t i = 0; i < TPos::Size; ++i)
	{
		int ind = (ICurDim + i) % TPos::Size;
		if(ele1.pos[ind] != ele2.pos[ind])
			return ele1.pos[ind] < ele2.pos[ind];
	}
	return false;
}

///	This method is normally not directly used.
/**
 * assigns elements to subsets.
 * Tries to construct subsets that all have the same size.
 * Runs with something like O(n*log(n)).
 * all elements between firstSubset and numSubsets are sorted by their
 * IDim-coordinate. The method then divides the elements into two parts.
 * The size of each part depends on the number of subsets that have to
 * be constructed from it.
 * The method the calls itself recursivly, until it finally assigns the
 * elements to subsets.
 */
template <class TElem, int IDim, class TPos>
bool PartitionElementsByRepeatedIntersection(ug::SubsetHandler& shOut,
										int firstSubset, int numSubsets,
										int firstElemInd, int numElems,
										std::vector<ElemWithPosition<TElem, TPos> >& vElems,
										int actDim)
{
	typename std::vector<ElemWithPosition<TElem, TPos> >::iterator iterBegin =
												vElems.begin() + firstElemInd;
	typename std::vector<ElemWithPosition<TElem, TPos> >::iterator iterEnd =
												iterBegin + numElems;

//	if we have to split vElems into multiple subsets
	if(numSubsets > 1)
	{
	//	sort the elements in the range between firstElemInd and firstElemInd + numElems
	//	by their coordinate-value for actDim.
		switch(actDim)
		{
			case 0:	sort(iterBegin, iterEnd, cmp<TElem, 0, TPos>); break;
			case 1:	sort(iterBegin, iterEnd, cmp<TElem, 1, TPos>); break;
			case 2:	sort(iterBegin, iterEnd, cmp<TElem, 2, TPos>); break;
		}

	//	get the number of elements that go into the first and second subtree
		int numSubsInFirst = numSubsets / 2;
		int numElemsInFirst = (int)(((float)numSubsInFirst / (float)numSubsets) * (float)numElems);

	//	perform recursive call to this method.
		bool retVal;
	//	first subtree
		retVal = PartitionElementsByRepeatedIntersection<TElem, IDim>(shOut, firstSubset,
															numSubsInFirst,
															firstElemInd, numElemsInFirst,
															vElems, (actDim + 1) % IDim);
	//	second subtree
		retVal |= PartitionElementsByRepeatedIntersection<TElem, IDim>(shOut,
															firstSubset + numSubsInFirst,
															numSubsets - numSubsInFirst,
															firstElemInd + numElemsInFirst,
															numElems - numElemsInFirst,
															vElems, (actDim + 1) % IDim);
	//	done
		return retVal;
	}
	else if(numSubsets == 1)
	{
	//	assign all elements in vElems to the subset-handler
		for(; iterBegin != iterEnd; ++iterBegin)
			shOut.assign_subset(iterBegin->pElem, firstSubset);

		return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////
template <class TElem, int IDimension, class TAPosition>
bool PartitionElementsByRepeatedIntersection(ug::SubsetHandler& shOut,
										ug::Grid& grid,
										int numSubsets,
										TAPosition& aVrtPos,
										int startDim)
{
//	TODO: move implementation to separate ..._impl.hpp file.
	using namespace ug;
	using namespace std;

	typedef typename TAPosition::ValueType	TPos;

//	position attachment accessor
	typename Grid::VertexAttachmentAccessor<TAPosition> aaPos(grid, aVrtPos);

//	fill all elements of type TElem of the grid into an ElemWithPosition-vector
	typedef ElemWithPosition<TElem, TPos>	Element;
	vector<Element> vElems(grid.template num<TElem>());

	typename geometry_traits<TElem>::iterator elemIter = grid.template begin<TElem>();

	for(typename vector<Element>::iterator iter = vElems.begin();
		iter != vElems.end(); ++iter, ++elemIter)
	{
		iter->pElem = *elemIter;
		iter->pos = CalculateCenter(*elemIter, aaPos);
	}

//	create the subsets from this vector
	return PartitionElementsByRepeatedIntersection<TElem, IDimension, TPos>(
															shOut, 0, numSubsets,
															0, vElems.size(),
															vElems,
															startDim % IDimension);
}

////////////////////////////////////////////////////////////////////////
template <class TElem, int IDimension, class TAPosition>
bool PartitionElementsByRepeatedIntersection(ug::SubsetHandler& shOut,
										ug::MultiGrid& mg,
										int level,
										int numSubsets,
										TAPosition& aVrtPos,
										int startDim)
{
//	TODO: move implementation to separate ..._impl.hpp file.
	using namespace ug;
	using namespace std;

	typedef typename TAPosition::ValueType	TPos;

//	position attachment accessor
	typename Grid::VertexAttachmentAccessor<TAPosition> aaPos(mg, aVrtPos);

//	fill all elements of type TElem of the grid into an ElemWithPosition-vector
	typedef ElemWithPosition<TElem, TPos>	Element;
	vector<Element> vElems(mg.template num<TElem>(level));

	typename geometry_traits<TElem>::iterator elemIter = mg.template begin<TElem>(level);

	for(typename vector<Element>::iterator iter = vElems.begin();
		iter != vElems.end(); ++iter, ++elemIter)
	{
		iter->pElem = *elemIter;
		iter->pos = CalculateCenter(*elemIter, aaPos);
	}

//	create the subsets from this vector
	return PartitionElementsByRepeatedIntersection<TElem, IDimension, TPos>(
															shOut, 0, numSubsets,
															0, vElems.size(),
															vElems,
															startDim % IDimension);
}

////////////////////////////////////////////////////////////////////////////////
template <class TElem, class TIterator, class TAAPos>
bool PartitionElements_RegularGrid(SubsetHandler& shOut,
								TIterator begin, TIterator end,
								int numCellsX, int numCellsY,
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
	vector_t min, max;
	{
		grid.begin_marking();
		vector<VertexBase*>	vrts, associatedVrts;
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

		CalculateBoundingBox(min, max, associatedVrts.begin(),
							 associatedVrts.end(), aaPos);
	}

	number width = max.x - min.x;
	number height = max.y - min.y;

	if(width < SMALL){
		UG_LOG("Can't execute PartitionElements_Rectangle: Geometry has no width.\n");
		return false;
	}
	if(height < SMALL){
		UG_LOG("Can't execute PartitionElements_Rectangle: Geometry has no height.\n");
		return false;
	}

//	iterate over all elements and calculate the index at which they shall be
//	inserted into the subset handler
	for(typename vector<TElem*>::iterator iter = elems.begin();
		iter != elems.end(); ++iter)
	{
		TElem* elem = *iter;
		vector_t center = CalculateCenter(elem, aaPos);

	//	get the cell index
		int xInd = (int)((number)numCellsX * (center.x - min.x) / width);
		int yInd = (int)((number)numCellsY * (center.y - min.y) / height);

	//	calculate the subset index (one could think of several alternatives here)
		int si = yInd * numCellsX + xInd;

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
