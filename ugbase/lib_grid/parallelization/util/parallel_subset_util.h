//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d14

#ifndef __H__LIB_GRID__PARALLELL_SUBSET_UTIL__
#define __H__LIB_GRID__PARALLELL_SUBSET_UTIL__

#include "../distributed_grid.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CreateSurfaceView
///	Collects all elements between iterBegin and iterEnd that don't have any children.
/**
 * Elements which are on the surface of the multi-grid-hierarchy
 * (elements that don't have children) are assigned to a subset of the
 * shSurfaceViewOut. The subset-index is taken from sh.
 *
 * Special care is taken for so called ghost-elements. Ghost elements
 * are elements that lie on the surface of the multi-grid but shall
 * not be considered for the suface-view (e.g. vertical-masters).
 *
 * TIterator has to be an STL compatible iterator, whose value-type is a
 * pointer to a VertexBase, EdgeBase, Face, Volume or derived class.
 *
 * make sure that all elements between iterBegin and iterEnd are members
 * of the MultiGrid of the given distGridMgr.
 *
 * If the distGridMgr is not operating on a multi-grid, then the method won't add
 * elements to the surface-view.
 *
 * This method will extend the surface-view. The caller is responsible for
 * clearing it before calling this method.
 */
template <class TIterator>
void CreateSurfaceView(SubsetHandler& shSurfaceViewOut,
						DistributedGridManager& distGridMgr,
						ISubsetHandler& sh, TIterator iterBegin,
						TIterator iterEnd);

}//	end of namespace

////////////////////////////////
//	include implementation
#include "parallel_subset_util_impl.hpp"

#endif
