//	created by Andreas Vogel

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_UTIL_GENERAL_IMPL__
#define __H__LIB_DISCRETIZATION__DOMAIN_UTIL_GENERAL_IMPL__

#include <string>
#include "./domain_util.h"

namespace ug{

///	returns the current dimension of the subset
template <typename TDomain>
int DimensionOfSubset(const TDomain& domain, int si)
{
	// extract subset handler
	const typename TDomain::subset_handler_type& sh = domain.get_subset_handler();

	// choose dimension
	if(sh.template num<Volume>(si) > 0) return 3;
	if(sh.template num<Face>(si) > 0) return 2;
	if(sh.template num<EdgeBase>(si) > 0) return 1;
	if(sh.template num<VertexBase>(si) > 0) return 0;
	else return -1;
}

////////////////////////////////////////////////////////////////////////
///	returns the corner coordinates of a geometric object
template <typename TElem, typename TDomain>
void CollectCornerCoordinates(	std::vector<typename TDomain::position_type>& vCornerCoordsOut,
								const TElem& elem, const TDomain& domain, bool clearContainer)
{
	if(clearContainer)
		vCornerCoordsOut.clear();

	// get position accessor
	const typename TDomain::position_accessor_type& aaPos = domain.get_position_accessor();

	// number of vertices of element
	const size_t numVertices = elem.num_vertices();

	// loop vertices
	for(int i = 0; i < numVertices; ++i)
	{
		// get element
		VertexBase* vert = elem.vertex(i);

		// write corner coordinates
		vCornerCoordsOut.push_back(aaPos[vert]);
	}
}


////////////////////////////////////////////////////////////////////////
/// writes domain to *.ugx file
template <typename TDomain>
bool WriteDomainToUGX(const char* filename, const TDomain& domain)
{
	// extract grid
	typename TDomain::grid_type& grid = domain.get_grid();

	// extract subset handler
	const typename TDomain::subset_handler_type& sh = domain.get_subset_handler();

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

	// save grid
	if(!SaveGridToUGX(grid, sh, strName.c_str()))
		{UG_LOG("WriteDomainToUGX: Cannot save grid.\n"); return false;}

	return true;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOMAIN_UTIL_GENERAL_IMPL__ */
