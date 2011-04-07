//	created by Andreas Vogel

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_UTIL_GENERAL_IMPL__
#define __H__LIB_DISCRETIZATION__DOMAIN_UTIL_GENERAL_IMPL__

#include <string>
#include <sstream>
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/tools/subset_handler_multi_grid.h"
#include "lib_grid/tools/subset_handler_grid.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/attachment_util.h"

#include "domain_util.h"
#include "lib_discretization/reference_element/reference_element.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/parallelization/parallelization_util.h"
#include "lib_grid/parallelization/parallel_refinement/parallel_global_multi_grid_refiner.h"
#endif

namespace ug{

//\todo: remove this function
#ifdef UG_PARALLEL
////////////////////////////////////////////////////////////////////////
template <class TDomain>
bool PrepareDomain(TDomain& domainOut, SubsetHandler& shTopViewOut,
					const char* filename,
					int numProcs,
					bool keepSrcGrid,
					size_t numPreRefinements,
					size_t numPostRefinements,
					bool writeProcessGrids,
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
	if(!LoadAndDistributeGrid(distGridMgr, sh, numProcs, filename,
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
/*
	LOG("creating top view ...");
	CreateSurfaceView(shTopViewOut, distGridMgr, sh,
					  mg.faces_begin(), mg.faces_end());
	LOG(" done\n");
*/

//	save the top-view for debug purposes
	if(writeProcessGrids)
	{
		if(pcl::GetProcRank() == pcl::GetOutputProcRank())
			SaveGridToFile(mg, "mg_top_view.obj", shTopViewOut);

		//	save the local geometry to a file
		std::stringstream ss;
		ss << "gridOnProc_" << pcl::GetProcRank() << ".obj";
		SaveGridToFile(mg, ss.str().c_str(), mg.get_hierarchy_handler());
	}

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
	if(gridLayoutMap.has_layout<VertexBase>(INT_V_MASTER)){
		UG_LOG(">>> has vertical masters\n");
	//	get the vertex layout (this is a multi-level-layout!)
		VertexLayout& vl = gridLayoutMap.get_layout<VertexBase>(INT_V_MASTER);
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
	if(gridLayoutMap.has_layout<VertexBase>(INT_V_SLAVE)){
		UG_LOG(">>> has vertical slaves\n");
	//	get the vertex layout (this is a multi-level-layout!)
		VertexLayout& vl = gridLayoutMap.get_layout<VertexBase>(INT_V_SLAVE);
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
	if(gridLayoutMap.has_layout<VertexBase>(INT_H_MASTER)){
		UG_LOG(">>> has horizontal masters\n");
	//	get the vertex layout (this is a multi-level-layout!)
		VertexLayout& vl = gridLayoutMap.get_layout<VertexBase>(INT_H_MASTER);
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
	if(gridLayoutMap.has_layout<VertexBase>(INT_H_SLAVE)){
		UG_LOG(">>> has horizontal slaves\n");
	//	get the vertex layout (this is a multi-level-layout!)
		VertexLayout& vl = gridLayoutMap.get_layout<VertexBase>(INT_H_SLAVE);
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
	if(TDomain::dim != 3)
	{
		typename TDomain::position_attachment_type& aPos =
			domainOut.get_position_attachment();
		ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition, aPos);
		mg.detach_from_vertices(aPosition);
	}

	return true;
}

#else
////////////////////////////////////////////////////////////////////////
template <class TDomain>
bool PrepareDomain(TDomain& domainOut, SubsetHandler& shTopViewOut,
					const char* filename,
					int numProcs,
					bool keepSrcGrid,
					size_t numPreRefinements,
					size_t numPostRefinements,
					bool writeProcessGrids,
					int autoAssignInnerObjectsToSubset,
					int autoAssignBoundaryObjectsToSubset)
{
	UG_LOG("PrepareDomain is currently not implemented for the serial case.\n");
	UG_LOG("  This has to be done in the near future!\n");
	return false;
}
#endif

/// returns if a subset is a regular grid
inline bool SubsetIsRegularGrid(const SubsetHandler& sh, int si)
{
//	check for constraining/constrained elements
	if(sh.num<HangingVertex>(si) > 0) return false;
	if(sh.num<ConstrainedEdge>(si) > 0) return false;
	if(sh.num<ConstrainingEdge>(si) > 0) return false;
	if(sh.num<ConstrainedTriangle>(si) > 0) return false;
	if(sh.num<ConstrainingTriangle>(si) > 0) return false;
	if(sh.num<ConstrainedQuadrilateral>(si) > 0) return false;
	if(sh.num<ConstrainingQuadrilateral>(si) > 0) return false;

//	if not found, subset describes a regular grid
	return true;
}

/// returns if a subset is a regular grid
inline bool SubsetIsRegularGrid(const MGSubsetHandler& sh, int si)
{
//	check for constraining/constrained elements
	if(sh.num<HangingVertex>(si) > 0) return false;
	if(sh.num<ConstrainedEdge>(si) > 0) return false;
	if(sh.num<ConstrainingEdge>(si) > 0) return false;
	if(sh.num<ConstrainedTriangle>(si) > 0) return false;
	if(sh.num<ConstrainingTriangle>(si) > 0) return false;
	if(sh.num<ConstrainedQuadrilateral>(si) > 0) return false;
	if(sh.num<ConstrainingQuadrilateral>(si) > 0) return false;

//	if not found, subset describes a regular grid
	return true;
}

/// returns if a subset is a regular grid
inline bool SubsetIsRegularGrid(const ISubsetHandler& ish, int si)
{
//	test SubsetHandler
	const SubsetHandler* sh = dynamic_cast<const SubsetHandler*>(&ish);
	if(sh != NULL)
		return SubsetIsRegularGrid(*sh, si);

//	test MGSubsetHandler
	const MGSubsetHandler* mgsh = dynamic_cast<const MGSubsetHandler*>(&ish);
	if(mgsh != NULL)
		return SubsetIsRegularGrid(*mgsh, si);

//	unknown type of subset handler
	throw(UGFatalError("Unknown SubsetHandler type."));
	return false;
}

///	returns the current dimension of the subset
inline int DimensionOfSubset(const SubsetHandler& sh, int si)
{
	// choose dimension
	if(sh.num<Volume>(si) > 0) return 3;
	if(sh.num<Face>(si) > 0) return 2;
	if(sh.num<EdgeBase>(si) > 0) return 1;
	if(sh.num<VertexBase>(si) > 0) return 0;
	else return -1;
}

///	returns the current dimension of the subset
inline int DimensionOfSubset(const MGSubsetHandler& sh, int si)
{
	// choose dimension
	if(sh.num<Volume>(si) > 0) return 3;
	if(sh.num<Face>(si) > 0) return 2;
	if(sh.num<EdgeBase>(si) > 0) return 1;
	if(sh.num<VertexBase>(si) > 0) return 0;
	else return -1;
}

///	returns the current dimension of the subset
inline int DimensionOfSubset(const ISubsetHandler& ish, int si)
{
//	test SubsetHandler
	const SubsetHandler* sh = dynamic_cast<const SubsetHandler*>(&ish);
	if(sh != NULL)
		return DimensionOfSubset(*sh, si);

//	test MGSubsetHandler
	const MGSubsetHandler* mgsh = dynamic_cast<const MGSubsetHandler*>(&ish);
	if(mgsh != NULL)
		return DimensionOfSubset(*mgsh, si);

//	unknown type of subset handler
	return -1;
}

inline int DimensionOfSubsets(const ISubsetHandler& sh)
{
//	dimension to be computed
	int dim = -1;

//	loop subsets
	for(int si = 0; si < sh.num_subsets(); ++si)
	{
	//	get dimension of subset
		int siDim = DimensionOfSubset(sh, si);

	//	if no dimension available, return -1
		if(siDim == -1) return -1;

	//	check if dimension is higher than already checked subsets
		if(dim < siDim)
			dim = siDim;
	}
	return dim;
}

///	returns the current dimension of the subset
template <typename TDomain>
inline int DimensionOfSubset(const TDomain& domain, int si)
{
	// extract subset handler
	const typename TDomain::subset_handler_type& sh = domain.get_subset_handler();

	return DimensionOfSubset(sh, si);
}

//	returns the corner coordinates of a geometric object
template <typename TElem, typename TAAPos>
void CollectCornerCoordinates(	std::vector<typename TAAPos::ValueType>& vCornerCoordsOut,
								const TElem& elem, const TAAPos& aaPos, bool clearContainer)
{
	if(clearContainer)
		vCornerCoordsOut.clear();

	// number of vertices of element
	const size_t numVertices = elem.num_vertices();

	// loop vertices
	for(size_t i = 0; i < numVertices; ++i)
	{
		// get element
		VertexBase* vert = elem.vertex(i);

		// write corner coordinates
		vCornerCoordsOut.push_back(aaPos[vert]);
	}
}

///	returns the corner coordinates of a geometric object
template <typename TElem, typename TDomain>
void CollectCornerCoordinates(	std::vector<typename TDomain::position_type>& vCornerCoordsOut,
								const TElem& elem, const TDomain& domain, bool clearContainer)
{
	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.get_position_accessor();

	CollectCornerCoordinates(vCornerCoordsOut, elem, aaPos, clearContainer);
}

////////////////////////////////////////////////////////////////////////
///	returns the size of a geometric object
template <typename TElem, typename TPosition>
number ElementSize(const TElem& elem, const Grid::VertexAttachmentAccessor<Attachment<TPosition> >& aaPos)
{
	// corner coords
	std::vector<TPosition> vCornerCoords;

	// load corner coords
	CollectCornerCoordinates(vCornerCoords, elem, aaPos);

	// get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type TRefElem;

	// dimension of Positions
	static const int dim = TPosition::Size;

	// return Element Size
	return ElementSize<TRefElem, dim>(&vCornerCoords[0]);
}

///	returns the size of a geometric object
template <typename TElem, typename TDomain>
number ElementSize(const TElem& elem, const TDomain& domain)
{
	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.get_position_accessor();

	return ElementSize(elem, aaPos);
}

// writes domain to *.ugx file
template <typename TDomain>
bool WriteDomainToUGX(const char* filename, const TDomain& domain)
{
	// filename
	std::string strName = filename;

	// check filename
	if(strName.find(" ") != std::string::npos)
		{UG_LOG("Filename must not include spaces. Cannot write domain."); return false;}

	// check if filename has already ending (if not add it)
	if(strName.find(".ugx") == std::string::npos)
	{
		if(strName.find(".") != std::string::npos)
		{
			UG_LOG("Filename must not include dots. Cannot write domain.");
			return false;
		}
		else
		{
			strName = strName + ".ugx";
		}
	}

	// types
	typedef typename TDomain::grid_type GridType;
	typedef typename TDomain::subset_handler_type SubsetHandlerType ;

	// extract grid and subset handler
	GridType& grid = *const_cast<GridType*>(&domain.get_grid());
	SubsetHandlerType& sh = *const_cast<SubsetHandlerType*>(&domain.get_subset_handler());

	// save grid
	if(!SaveGridToUGX(grid, sh, strName.c_str()))
		{UG_LOG("WriteDomainToUGX: Cannot save grid.\n"); return false;}

	return true;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOMAIN_UTIL_GENERAL_IMPL__ */
