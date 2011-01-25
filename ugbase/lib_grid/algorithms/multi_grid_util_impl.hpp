//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d18

#ifndef __H__UG__MULTI_GRID_UTIL_IMPL__
#define __H__UG__MULTI_GRID_UTIL_IMPL__

#include "subset_util.h"

namespace ug
{

template <class TElem>
void CollectSurfaceViewElements(ISubsetHandler& surfaceViewOut,
								MultiGrid& mg,
								MultiGridSubsetHandler& mgsh)
{
//	some typedefs
	typedef typename geometry_traits<TElem>::iterator ElemIter;

//	clear the target surfaceView
	surfaceViewOut.clear();
		
//	iterate through all levels of the mgsh
	for(size_t level = 0; level < mgsh.num_levels(); ++level){
//	iterate through all subsets on that level
		for(int si = 0; si < mgsh.num_subsets(); ++si){
		//	iterate through all elelements in the subset on that level
			for(ElemIter iter = mgsh.begin<TElem>(si, level);
				iter != mgsh.end<TElem>(si, level); ++iter)
			{
				TElem* elem = *iter;
			//	check whether the element has children
				if(!mg.has_children(elem)){
				//	the element is a part of the surface-view. add it to the handler
					surfaceViewOut.assign_subset(elem, si);
				}
			}
		}
	}
}

template <class TElem, class TSurfaceView>
void CreateSurfaceView(TSurfaceView& surfaceViewOut,
						MultiGrid& mg,
						MultiGridSubsetHandler& mgsh)
{
//	This method clears the surfaceViewOut and assigns all objects of
//	type TElem which lie on the surface of the mg to the surface view.
	CollectSurfaceViewElements<TElem>(surfaceViewOut, mg, mgsh);
	
//	assign associated elements of lower dimension to the surface view
	AssignAssociatedLowerDimElemsToSubsets<TElem>(surfaceViewOut, mgsh);	
}

/*



*/

}//	end of namespace

#endif
