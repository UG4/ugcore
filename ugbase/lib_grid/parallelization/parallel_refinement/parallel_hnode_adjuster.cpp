// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 16.11.2012 (d,m,y)

#include "parallel_hnode_adjuster.h"
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug{

template <class TStdVector>
static bool ContainsInterfaceElem(const TStdVector& elems,
								  DistributedGridManager& distGridMgr)
{
	for(size_t i = 0; i < elems.size(); ++i){
		if(distGridMgr->is_interface_element(elems[i]))
			return true;
	}
	return false;
}


ParallelHNodeAdjuster::AdjustRetVal ParallelHNodeAdjuster::
ref_marks_changed(IRefiner& ref,
			   	  const std::vector<VertexBase*>& vrts,
			   	  const std::vector<EdgeBase*>& edges,
			   	  const std::vector<Face*>& faces,
			   	  const std::vector<Volume*>& vols)
{
	UG_ASSERT(ref.grid(), "A refiner has to operate on a grid, before marks can be adjusted!");
	if(!ref.grid()){
		return CONTINUE_IF_MARKED_NEW;
	}
	
	if(!grid.is_parallel())
		return CONTINUE_IF_MARKED_NEW;

	DistributedGridManager& distGridMgr = *grid.distributed_grid_manager();

//	check whether new interface elements have been selected
	bool newInterfaceVrtsMarked = ContainsInterfaceElem(vrts, distGridMgr);
	bool newInterfaceEdgeMarked = ContainsInterfaceElem(edges, distGridMgr);
	bool newInterfaceFacesMarked = ContainsInterfaceElem(faces, distGridMgr);
	bool newInterfaceVolsMarked = ContainsInterfaceElem(vols, distGridMgr);

	bool newlyMarkedElems = newInterfaceVrtsMarked ||
							newInterfaceEdgeMarked ||
							newInterfaceFacesMarked ||
							newInterfaceVolsMarked;

	bool exchangeFlag = pcl::OneProcTrue(newlyMarkedElems);

	if(exchangeFlag){

	}

	return CONTINUE_IF_MARKED_NEW;
}
}// end of namespace
