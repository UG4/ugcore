//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d18

#ifndef __H__UG__MULTI_GRID_UTIL__
#define __H__UG__MULTI_GRID_UTIL__

#include "lib_grid/lg_base.h"

namespace ug
{

/**	iterates over the multi-grid and assigns all surface-elements to surfaceViewOut.*/
template <class TElem>
void CollectSurfaceViewElements(ISubsetHandler& surfaceViewOut,
								MultiGrid& mg,
								MGSubsetHandler& mgsh);

/**	calls CollectSurfaceViewElements and then assigns all associated elements
 *	of lower dimension to the surface-view, too.
 *
 *	TSurfaceView has to be a SubsetHandler or MGSubsetHandler compatible type.*/
template <class TElem, class TSurfaceView>
void CreateSurfaceView(TSubsetHandler& surfaceViewOut,
						MultiGrid& mg,
						MGSubsetHandler& mgsh);
						
}//	end of namespace

#endif
