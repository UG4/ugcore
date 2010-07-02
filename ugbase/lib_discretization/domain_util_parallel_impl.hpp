//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d30

#ifndef __H__LIBDISCRETIZATION__DOMAIN_UTIL_PARALLEL_IMPL__
#define __H__LIBDISCRETIZATION__DOMAIN_UTIL_PARALLEL_IMPL__

#include <sstream>
#include "lib_grid/lg_base.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/parallelization/parallelization.h"
#include "lib_grid/algorithms/attachment_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
//	CreateTopView
///	Collects all geom-objects between iterBegin and iterEnd that don't have any children.
/**
 * Ghost objects are ignored.
 *
 * In the current implementation all top-view-elements will be assigned
 * to subset 0 - this will change in future versions (as well as the parameters).
 *
 * TIterator has to be an STL compatible iterator, whose value-type is a
 * pointer to an VertexBase, EdgeBase, Face, Volume or derived class.
 *
 * make sure that all elements between iterBegin and iterEnd are members
 * of the given MultiGrid.
 *
 * This method will extend the subsets. The caller is responsible for
 * clearing them before calling this method.
 */

template <class TIterator>
void CreateTopView(DistributedGridManager& distGridMgr,
					SubsetHandler& sh,
					TIterator iterBegin, TIterator iterEnd)
{
	MultiGrid* pMG = dynamic_cast<MultiGrid*>(distGridMgr.get_assigned_grid());
	if(!pMG){
		UG_LOG("  Can't create top-view. A Multigrid is required.\n");
		return;
	}

	while(iterBegin != iterEnd)
	{
		if(!distGridMgr.is_ghost(*iterBegin)){
			if(!pMG->has_children(*iterBegin))
				sh.assign_subset(*iterBegin, 0);
		}
		++iterBegin;
	}
}

////////////////////////////////////////////////////////////////////////
template <class TDomain>
bool PrepareDomain(TDomain& domainOut, SubsetHandler& shTopViewOut,
					const char* filename,
					int numProcs,
					bool keepSrcGrid,
					size_t numPreRefinements,
					size_t numPostRefinements,
					int autoAssignInnerObjectsToSubset,
					int autoAssignBoundaryObjectsToSubset)
{
//	load the base-grid of the multi-grid.
//	we're enabling GRIDOPT_FULL_INTERCONNECTION here since the
//	multi-grid refiner would autoenable it otherwise.
	typedef typename TDomain::grid_type grid_type;
	grid_type& mg = domainOut.get_grid();;

//  Create Subset Handler for all elements of Multigrid
	typedef typename TDomain::subset_handler_type subset_handler_type;
	subset_handler_type& sh = domainOut.get_subset_handler();

	typename TDomain::distributed_grid_manager_type* pDistGridMgr =
						domainOut.get_distributed_grid_manager();

	if(!pDistGridMgr){
		assert(!"  Distributed Grid Manager required in PrepareDomain!\n");
		UG_LOG("  Distributed Grid Manager required in PrepareDomain!\n");
		return false;
	}
	
	typename TDomain::distributed_grid_manager_type& distGridMgr = *pDistGridMgr;
	
//TODO: add partition-method
	if(!LoadAndDistributeGrid(*pDistGridMgr, sh, numProcs, filename,
							keepSrcGrid,
							AdjustGrid_AutoAssignSubsetsAndRefine(
								autoAssignInnerObjectsToSubset,
								autoAssignBoundaryObjectsToSubset,
								numPreRefinements)))
	{
		return false;
	}

//	if post-refinement is required then do it now
//	initialise the multi-grid-refiner
	if(numPostRefinements > 0){
		ParallelGlobalMultiGridRefiner refiner(distGridMgr);

	//	perform the refinement
		PROFILE_BEGIN(postRefinement);
		for(size_t i = 0; i < numPostRefinements; ++i)
		{
			size_t level = mg.num_levels() - 1;
			LOG("refining level " << level << " ...");
			
		//	perform refinement.
			refiner.refine();
			LOG(" done\n");
		}
		PROFILE_END();
	}

//	select the top-view
	LOG("creating top view ...");
	CreateTopView(distGridMgr, shTopViewOut, mg.faces_begin(), mg.faces_end());
	LOG(" done\n");

//	save the top-view for debug purposes
	if(pcl::GetProcRank() == pcl::GetOutputProcRank())
		SaveGridToFile(mg, "mg_top_view.obj", shTopViewOut);

//	save the local geometry to a file
	std::stringstream ss;
	ss << "gridOnProc_" << pcl::GetProcRank() << ".obj";
	SaveGridToFile(mg, ss.str().c_str(), mg.get_hierarchy_handler());
	
/*
//	debug - check subsets
//	make sure that all elements are in the same subset as their parents
	for(size_t l = 1; l <= mg.num_levels() - 1; ++l)
	{
		for(geometry_traits<VertexBase>::iterator iter = mg.template begin<VertexBase>(l);
			iter != mg.template end<VertexBase>(l); ++iter)
		{
			if(sh.get_subset_index(*iter) !=
				sh.get_subset_index(mg.get_parent(*iter)))
			{
				UG_LOG("VERTEX-CHILD AND PARENT SUBSET-INDEX MISMATCH!\n");
			}
		}

		for(geometry_traits<EdgeBase>::iterator iter = mg.template begin<EdgeBase>(l);
			iter != mg.template end<EdgeBase>(l); ++iter)
		{
			if(sh.get_subset_index(*iter) !=
				sh.get_subset_index(mg.get_parent(*iter)))
			{
				UG_LOG("EDGE-CHILD AND PARENT SUBSET-INDEX MISMATCH!\n");
			}
		}

		for(geometry_traits<Face>::iterator iter = mg.template begin<Face>(l);
			iter != mg.template end<Face>(l); ++iter)
		{
			if(sh.get_subset_index(*iter) !=
				sh.get_subset_index(mg.get_parent(*iter)))
			{
				UG_LOG("FACE-CHILD AND PARENT SUBSET-INDEX MISMATCH!\n");
			}
		}
	}
*/

/*
//debug: check interfaces
	GridLayoutMap& gridLayoutMap = distGridMgr.grid_layout_map();
	if(gridLayoutMap.has_layout<VertexBase>(INT_VERTICAL_MASTER)){
		UG_LOG(">>> has vertical masters\n");
	//	get the vertex layout (this is a multi-level-layout!)
		VertexLayout& vl = gridLayoutMap.get_layout<VertexBase>(INT_VERTICAL_MASTER);
		for(size_t l = 0; l < vl.num_levels(); ++l)
		{
			for(VertexLayout::iterator iter = vl.begin(l);
				iter != vl.end(l); ++iter)
			{
				VertexLayout::Interface& interface = vl.interface(iter);
				UG_LOG("interface size on level " << l << ": " << interface.size() << std::endl);
			}
		}
	}
	if(gridLayoutMap.has_layout<VertexBase>(INT_VERTICAL_SLAVE)){
		UG_LOG(">>> has vertical slaves\n");
	//	get the vertex layout (this is a multi-level-layout!)
		VertexLayout& vl = gridLayoutMap.get_layout<VertexBase>(INT_VERTICAL_SLAVE);
		for(size_t l = 0; l < vl.num_levels(); ++l)
		{
			for(VertexLayout::iterator iter = vl.begin(l);
				iter != vl.end(l); ++iter)
			{
				VertexLayout::Interface& interface = vl.interface(iter);
				UG_LOG("interface size on level " << l << ": " << interface.size() << std::endl);
			}
		}
	}
	if(gridLayoutMap.has_layout<VertexBase>(INT_MASTER)){
		UG_LOG(">>> has horizontal masters\n");
	//	get the vertex layout (this is a multi-level-layout!)
		VertexLayout& vl = gridLayoutMap.get_layout<VertexBase>(INT_MASTER);
		for(size_t l = 0; l < vl.num_levels(); ++l)
		{
			for(VertexLayout::iterator iter = vl.begin(l);
				iter != vl.end(l); ++iter)
			{
				VertexLayout::Interface& interface = vl.interface(iter);
				UG_LOG("interface size on level " << l << ": " << interface.size() << std::endl);
			}
		}
	}
	if(gridLayoutMap.has_layout<VertexBase>(INT_SLAVE)){
		UG_LOG(">>> has horizontal slaves\n");
	//	get the vertex layout (this is a multi-level-layout!)
		VertexLayout& vl = gridLayoutMap.get_layout<VertexBase>(INT_SLAVE);
		for(size_t l = 0; l < vl.num_levels(); ++l)
		{
			for(VertexLayout::iterator iter = vl.begin(l);
				iter != vl.end(l); ++iter)
			{
				VertexLayout::Interface& interface = vl.interface(iter);
				UG_LOG("interface size on level " << l << ": " << interface.size() << std::endl);
			}
		}
	}
*/

//	we have to copy or convert the grids position attachment.
//	if we're in 2d, the third dimension will be skipped.
//	if we're in 3d, the values are simply copied.
	typename TDomain::position_attachment_type& aPos =
				domainOut.get_position_attachment();
	ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition, aPos);
	mg.detach_from_vertices(aPosition);

	return true;
}


} // end namespace ug


#endif /* __H__LIBDISCRETIZATION__DOMAIN_UTIL_PARALLEL_IMPL__ */
