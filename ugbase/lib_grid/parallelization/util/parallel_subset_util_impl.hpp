#ifndef __H__LIB_GRID__PARALLELL_SUBSET_UTIL_IMPL__
#define __H__LIB_GRID__PARALLELL_SUBSET_UTIL_IMPL__

#include "../distributed_grid.h"
#include "compol_subset.h"
#include "lib_grid/algorithms/subset_util.h"
#include "pcl/pcl_interface_communicator.h"

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
		throw(UGError("  Can't create surface-view. A Multigrid is required.\n"));
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
		throw(UGError("  Can't create surface-view. A Multigrid is required.\n"));
	}

//	This method clears the surfaceViewOut and assigns all objects of
//	which lie on the surface of the mg to the surface view.
	CollectSurfaceViewElements<Vertex>(surfaceViewOut, distGridMgr, mgsh, true);
	CollectSurfaceViewElements<Edge>(surfaceViewOut, distGridMgr, mgsh, false);
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
		AssignAssociatedSidesToSubsets<Edge>(surfaceViewOut, mgsh);
	}
	else
	{
		UG_LOG("INFO in CreateSurfaceView: Performing AssignAssociatedLowerDimElemsToSubsets ");
		UG_LOG("for all elements (Small Performance drawback).\n");

		AssignAssociatedLowerDimElemsToSubsets<Volume>(surfaceViewOut, mgsh);
		AssignAssociatedLowerDimElemsToSubsets<Face>(surfaceViewOut, mgsh);
		AssignAssociatedLowerDimElemsToSubsets<Edge>(surfaceViewOut, mgsh);
	}

//	set num subsets for surface view
//	this needed, when some elements subsets do not appear in this part of the
//	surface grid. But still, num_subsets() has to return the correct number.
	if(mgsh.num_subsets() > 0)
		surfaceViewOut.subset_required(mgsh.num_subsets() - 1);

//	for a parallel subset handler we have to copy the subset-indices to
//	avoid problems at interfaces.
//	This happens e.g. in 2d, where a triangle is an element of the surface view only
//	on one process. Associated vertices on other processes wouldn't know that
//	they are surface view element, too. This has to be communicated.
	ComPol_Subset<VertexLayout> cpSubsetVRT(surfaceViewOut);
	ComPol_Subset<EdgeLayout> cpSubsetEDGE(surfaceViewOut);
	ComPol_Subset<FaceLayout> cpSubsetFACE(surfaceViewOut);
	ComPol_Subset<VolumeLayout> cpSubsetVOL(surfaceViewOut);
	pcl::InterfaceCommunicator<VertexLayout> comVRT;
	pcl::InterfaceCommunicator<EdgeLayout> comEDGE;
	pcl::InterfaceCommunicator<FaceLayout> comFACE;
	pcl::InterfaceCommunicator<VolumeLayout> comVOL;

	comVRT.send_data(distGridMgr.grid_layout_map().template get_layout<Vertex>(INT_H_SLAVE), cpSubsetVRT);
	comEDGE.send_data(distGridMgr.grid_layout_map().template get_layout<Edge>(INT_H_SLAVE), cpSubsetEDGE);
	comFACE.send_data(distGridMgr.grid_layout_map().template get_layout<Face>(INT_H_SLAVE), cpSubsetFACE);
	comVOL.send_data(distGridMgr.grid_layout_map().template get_layout<Volume>(INT_H_SLAVE), cpSubsetVOL);

	comVRT.receive_data(distGridMgr.grid_layout_map().template get_layout<Vertex>(INT_H_MASTER), cpSubsetVRT);
	comEDGE.receive_data(distGridMgr.grid_layout_map().template get_layout<Edge>(INT_H_MASTER), cpSubsetEDGE);
	comFACE.receive_data(distGridMgr.grid_layout_map().template get_layout<Face>(INT_H_MASTER), cpSubsetFACE);
	comVOL.receive_data(distGridMgr.grid_layout_map().template get_layout<Volume>(INT_H_MASTER), cpSubsetVOL);

	comVRT.communicate();
	comEDGE.communicate();
	comFACE.communicate();
	comVOL.communicate();

	comVRT.send_data(distGridMgr.grid_layout_map().template get_layout<Vertex>(INT_H_MASTER), cpSubsetVRT);
	comEDGE.send_data(distGridMgr.grid_layout_map().template get_layout<Edge>(INT_H_MASTER), cpSubsetEDGE);
	comFACE.send_data(distGridMgr.grid_layout_map().template get_layout<Face>(INT_H_MASTER), cpSubsetFACE);
	comVOL.send_data(distGridMgr.grid_layout_map().template get_layout<Volume>(INT_H_MASTER), cpSubsetVOL);

	comVRT.receive_data(distGridMgr.grid_layout_map().template get_layout<Vertex>(INT_H_SLAVE), cpSubsetVRT);
	comEDGE.receive_data(distGridMgr.grid_layout_map().template get_layout<Edge>(INT_H_SLAVE), cpSubsetEDGE);
	comFACE.receive_data(distGridMgr.grid_layout_map().template get_layout<Face>(INT_H_SLAVE), cpSubsetFACE);
	comVOL.receive_data(distGridMgr.grid_layout_map().template get_layout<Volume>(INT_H_SLAVE), cpSubsetVOL);

	comVRT.communicate();
	comEDGE.communicate();
	comFACE.communicate();
	comVOL.communicate();
}


}//	end of namespace

#endif
