//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d14

#ifndef __H__LIB_GRID__PARALLELL_SUBSET_UTIL_IMPL__
#define __H__LIB_GRID__PARALLELL_SUBSET_UTIL_IMPL__

#include "../distributed_grid.h"
#include "lib_grid/algorithms/subset_util.h"

namespace ug
{
template <class TElem>
void CollectSurfaceViewElements(ISubsetHandler& surfaceViewOut,
                                DistributedGridManager& distGridMgr,
								MultiGridSubsetHandler& mgsh,
								bool clearContainer)
{
//	get multigrid
	MultiGrid* pMG = dynamic_cast<MultiGrid*>(distGridMgr.get_assigned_grid());
	if(!pMG){
		throw(UGFatalError("  Can't create surface-view. A Multigrid is required.\n"));
	}

//	some typedefs
	typedef typename geometry_traits<TElem>::iterator ElemIter;

//	clear the target surfaceView
	if(clearContainer)
		surfaceViewOut.clear();

//	iterate through all levels of the mgsh
	for(size_t level = 0; level < mgsh.num_levels(); ++level){
//	iterate through all subsets on that level
		for(int si = 0; si < mgsh.num_subsets(); ++si){
		//	iterate through all elements in the subset on that level
			for(ElemIter iter = mgsh.begin<TElem>(si, level);
				iter != mgsh.end<TElem>(si, level); ++iter)
			{
				TElem* elem = *iter;
			//	check whether the element has children
				if(!distGridMgr.is_ghost(elem)){
					if(!pMG->has_children(elem)){
					//	the element is a part of the surface-view. add it to the handler
						surfaceViewOut.assign_subset(elem, si);
					}
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////
//	CreateSurfaceView
template <class TSurfaceView>
void CreateSurfaceView(TSurfaceView& surfaceViewOut,
                       DistributedGridManager& distGridMgr,
						MultiGridSubsetHandler& mgsh)
{
//	get multigrid
	MultiGrid* pMG = dynamic_cast<MultiGrid*>(distGridMgr.get_assigned_grid());
	if(!pMG){
		throw(UGFatalError("  Can't create surface-view. A Multigrid is required.\n"));
	}

//	This method clears the surfaceViewOut and assigns all objects of
//	which lie on the surface of the mg to the surface view.
	CollectSurfaceViewElements<VertexBase>(surfaceViewOut, distGridMgr, mgsh, true);
	CollectSurfaceViewElements<EdgeBase>(surfaceViewOut, distGridMgr, mgsh, false);
	CollectSurfaceViewElements<Face>(surfaceViewOut, distGridMgr, mgsh, false);
	CollectSurfaceViewElements<Volume>(surfaceViewOut, distGridMgr, mgsh, false);

//	assign associated elements of lower dimension to the surface view
	bool assignSidesOnly = true;
	if(mgsh.num<Volume>() > 0 && !pMG->option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
		assignSidesOnly = false;
	else if(mgsh.num<Volume>() > 0 && !pMG->option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
		assignSidesOnly = false;
	else if(mgsh.num<Face>() > 0 && !pMG->option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
		assignSidesOnly = false;

	if(assignSidesOnly){
		AssignAssociatedSidesToSubsets<Volume>(surfaceViewOut, mgsh);
		AssignAssociatedSidesToSubsets<Face>(surfaceViewOut, mgsh);
		AssignAssociatedSidesToSubsets<EdgeBase>(surfaceViewOut, mgsh);
	}
	else
	{
		UG_LOG("INFO in CreateSurfaceView: Performing AssignAssociatedLowerDimElemsToSubsets ");
		UG_LOG("for all elements (Small Performance drawback).\n");

		AssignAssociatedLowerDimElemsToSubsets<Volume>(surfaceViewOut, mgsh);
		AssignAssociatedLowerDimElemsToSubsets<Face>(surfaceViewOut, mgsh);
		AssignAssociatedLowerDimElemsToSubsets<EdgeBase>(surfaceViewOut, mgsh);
	}

//	set num subsets for surface view
//	this needed, when some elements subsets do not appear in this part of the
//	surface grid. But still, num_subsets() has to return the correct number.
	if(mgsh.num_subsets() > 0)
		surfaceViewOut.subset_required(mgsh.num_subsets() - 1);
}


}//	end of namespace

#endif
