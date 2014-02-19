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
								MultiGridSubsetHandler& mgsh,
								bool clearContainer)
{
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
				if(!mg.has_children(elem)){
				//	the element is a part of the surface-view. add it to the handler
					surfaceViewOut.assign_subset(elem, si);
				}
			}
		}
	}
}

template <class TSurfaceView>
void CreateSurfaceView(TSurfaceView& surfaceViewOut,
						MultiGrid& mg,
						MultiGridSubsetHandler& mgsh)
{
//	This method clears the surfaceViewOut and assigns all objects of
//	which lie on the surface of the mg to the surface view.
	CollectSurfaceViewElements<Vertex>(surfaceViewOut, mg, mgsh, true);
	CollectSurfaceViewElements<EdgeBase>(surfaceViewOut, mg, mgsh, false);
	CollectSurfaceViewElements<Face>(surfaceViewOut, mg, mgsh, false);
	CollectSurfaceViewElements<Volume>(surfaceViewOut, mg, mgsh, false);
	
//	assign associated elements of lower dimension to the surface view
	bool assignSidesOnly = true;
	if(mgsh.num<Volume>() > 0 && !mg.option_is_enabled(VOLOPT_AUTOGENERATE_FACES))
		assignSidesOnly = false;
	else if(mgsh.num<Volume>() > 0 && !mg.option_is_enabled(VOLOPT_AUTOGENERATE_EDGES))
		assignSidesOnly = false;
	else if(mgsh.num<Face>() > 0 && !mg.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
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
}

template <class TElem>
bool IsSubSurfaceElement(MultiGrid& mg, TElem* e, bool checkSides)
{
	typedef typename TElem::grid_base_object TBaseElem;

	size_t numChildren = mg.num_children<TBaseElem>(e);
	for(size_t i = 0; i < numChildren; ++i){
		TBaseElem* child = mg.get_child<TBaseElem>(e, i);
		assert(child);
		if(mg.num_children<TBaseElem>(child) != 0)
			return false;
	}

//	ok all children of the same base type are surface elements.
	if(checkSides){
	//	Now we have to check whether all sides are surface elements, too.
	//	Since vertices do not have sides, this recursion will terminate.
		for(size_t i = 0; i < e->num_sides(); ++i){
			if(!IsSubSurfaceElement(mg, mg.get_side(e, i)))
				return false;
		}
	}

	return true;
}

}//	end of namespace

#endif
