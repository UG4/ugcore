#ifndef __H__UG__MULTI_GRID_UTIL__
#define __H__UG__MULTI_GRID_UTIL__

#include "lib_grid/lg_base.h"

namespace ug
{

/**	iterates over the multi-grid and assigns all surface-elements to surfaceViewOut.*/
template <class TElem>
void CollectSurfaceViewElements(ISubsetHandler& surfaceViewOut,
								MultiGrid& mg,
								MGSubsetHandler& mgsh,
								bool clearContainer = true);

/**	calls CollectSurfaceViewElements and then assigns all associated elements
 *	of lower dimension to the surface-view, too.
 *
 *	TSurfaceView has to be a SubsetHandler or MGSubsetHandler compatible type.*/
template <class TElem, class TSurfaceView>
void CreateSurfaceView(TSurfaceView& surfaceViewOut,
						MultiGrid& mg,
						MGSubsetHandler& mgsh);

///	returns true, if the element lies one level below the surface
/**	If checkSides == false, then only children of the same base type as TElem
 * will be regarded. If checkSides == true, also children of the elements sides
 * will be regarded, too.
 * If all regarded children lie on the surface (i.e. do not have children them selfs),
 * then the element is regarded as a surface element.
 */
template <class TElem>
bool IsSubSurfaceElement(MultiGrid& mg, TElem* e, bool checkSides = false);


}//	end of namespace

////////////////////////////////
//	include implementation
#include "multi_grid_util_impl.hpp"

#endif
