//	created by Andreas Vogel

#ifndef __H__UG__LIB_DISC__DOMAIN_UTIL_GENERAL_IMPL__
#define __H__UG__LIB_DISC__DOMAIN_UTIL_GENERAL_IMPL__

#include <string>
#include <sstream>
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/tools/subset_handler_multi_grid.h"
#include "lib_grid/tools/subset_handler_grid.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/attachment_util.h"

#include "domain_util.h"
#include "lib_disc/reference_element/reference_element.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/parallelization/parallelization_util.h"
#include "lib_grid/parallelization/parallelization.h"
#endif

namespace ug{

template <typename TDomain>
bool LoadDomain(TDomain& domain, const char* filename)
{
	return LoadDomain(domain, filename, 0);
}


template <typename TDomain>
bool LoadDomain(TDomain& domain, const char* filename, int procId)
{
#ifdef UG_PARALLEL
	if((procId >= 0 ) && (pcl::GetProcRank() != procId))
		return true;
#endif

	if(!LoadGridFromFile(domain.grid(), domain.subset_handler(),
						 filename, domain.position_attachment()))
	{
		return false;
	}

	domain.update_local_subset_dim_property();

	return true;
}


template <typename TDomain>
bool SaveDomain(TDomain& domain, const char* filename)
{
	return SaveGridToFile(domain.grid(), domain.subset_handler(),
						  filename, domain.position_attachment());
}



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
inline int DimensionOfSubset(const ISubsetHandler& ish, int si
#ifdef UG_PARALLEL
                             , const pcl::ProcessCommunicator* pProcCom
#endif
							)
{
	int dim = DIM_SUBSET_EMPTY_GRID;
	if(ish.contains_volumes(si))
		dim = 3;
	else if(ish.contains_faces(si))
		dim = 2;
	else if(ish.contains_edges(si))
		dim = 1;
	else if(ish.contains_vertices(si))
		dim = 0;

//	in parallel, we have to check if another proc has a higher dimension
	#ifdef UG_PARALLEL
		if(pProcCom) dim = pProcCom->allreduce(dim, PCL_RO_MAX);
	#endif

	return dim;
}

inline int DimensionOfSubsets(const ISubsetHandler& sh
#ifdef UG_PARALLEL
                             , const pcl::ProcessCommunicator* pProcCom
#endif
							)
{
//	dimension to be computed
	int dim = DIM_SUBSET_EMPTY_GRID;

//	loop subsets
	for(int si = 0; si < sh.num_subsets(); ++si)
	{
	//	get dimension of subset
		int siDim = DimensionOfSubset(sh, si);

	//	if empty grid given, skip
		if(siDim == DIM_SUBSET_EMPTY_GRID) continue;

	//	check if dimension is higher than already checked subsets
		if(dim < siDim)
			dim = siDim;
	}

//	get globally highest subset
#ifdef UG_PARALLEL
	if(pProcCom) dim = pProcCom->allreduce(dim, PCL_RO_MAX);
#endif

//	return computed domain
	return dim;
}

///	returns the current dimension of the subset
template <typename TDomain>
inline int DimensionOfDomainSubset(const TDomain& domain, int si
#ifdef UG_PARALLEL
                             , const pcl::ProcessCommunicator* pProcCom
#endif
)
{
	// extract subset handler
	const typename TDomain::subset_handler_type& sh = domain.subset_handler();

#ifdef UG_PARALLEL
	return DimensionOfSubset(sh, si, pProcCom);
#else
	return DimensionOfSubset(sh, si);
#endif
}

///	returns the current dimension of the domain
template <typename TDomain>
inline int DimensionOfDomain(const TDomain& domain
#ifdef UG_PARALLEL
                             , const pcl::ProcessCommunicator* pProcCom
#endif
							)
{
//	get grid
	const typename TDomain::grid_type& grid = domain.grid();

// 	get local dimension of subset
	int locDim = DIM_SUBSET_EMPTY_GRID;
	if(grid.template num<VertexBase>() > 0) locDim = 0;
	if(grid.template num<EdgeBase>() > 0) locDim = 1;
	if(grid.template num<Face>() > 0) locDim = 2;
	if(grid.template num<Volume>() > 0) locDim = 3;

//	in parallel, we have to check if another proc has a higher dimension
#ifdef UG_PARALLEL
	if(pProcCom) locDim = pProcCom->allreduce(locDim, PCL_RO_MAX);
#endif

//	return result
	return locDim;
}

//	returns the corner coordinates of a geometric object
template <typename TElem, typename TAAPos>
void CollectCornerCoordinates(	std::vector<typename TAAPos::ValueType>& vCornerCoordsOut,
								TElem* elem, const TAAPos& aaPos, bool clearContainer)
{
	if(clearContainer)
		vCornerCoordsOut.clear();

	// number of vertices of element
	const size_t numVertices = NumVertices(elem);

	// loop vertices
	for(size_t i = 0; i < numVertices; ++i)
	{
		// get element
		VertexBase* vert = GetVertex(elem, i);

		// write corner coordinates
		vCornerCoordsOut.push_back(aaPos[vert]);
	}
}

//	returns the corner coordinates of a geometric object
template <typename TElem, typename TAAPos>
void CollectCornerCoordinates(	std::vector<typename TAAPos::ValueType>& vCornerCoordsOut,
								const TElem& elem, const TAAPos& aaPos, bool clearContainer)
{
//	cast constness away
	TElem* pElem = const_cast<TElem*>(&elem);

//	forward
	return CollectCornerCoordinates(vCornerCoordsOut, pElem, aaPos, clearContainer);
}

///	returns the corner coordinates of a geometric object
template <typename TElem, typename TDomain>
void CollectCornerCoordinates(	std::vector<typename TDomain::position_type>& vCornerCoordsOut,
								const TElem& elem, const TDomain& domain, bool clearContainer)
{
	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.position_accessor();

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
	const typename TDomain::position_accessor_type& aaPos = domain.position_accessor();

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
	GridType& grid = *const_cast<GridType*>(&domain.grid());
	SubsetHandlerType& sh = *const_cast<SubsetHandlerType*>(&domain.subset_handler());

	// save grid
	if(!SaveGridToUGX(grid, sh, strName.c_str()))
		{UG_LOG("WriteDomainToUGX: Cannot save grid.\n"); return false;}

	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOMAIN_UTIL_GENERAL_IMPL__ */
