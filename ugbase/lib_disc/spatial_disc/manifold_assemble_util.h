

#ifndef __H_MANIFOLD_ASSEMBLE_UTIL__
#define __H_MANIFOLD_ASSEMBLE_UTIL__

#include "bindings/lua/lua_user_data.h"
#include "common/util/smart_pointer.h"
#include "lib_disc/spatial_disc/domain_disc.h"

// other ug4 modules
#include "common/common.h"
#include "lib_grid/tools/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"


namespace ug {


/**
 * This function marks all grid elements of type TElem for assembly BUT H-slaves in a BoolMarker.
 *
 * \param[out]		grid		Grid
 * \param[out]		bMarker		BoolMarker containing all marked elements of type TElem
 * 								BUT H-slaves
 */
template <typename TElem>
static void MarkAllElemsForAssemblyButHSlaves(Grid& grid, BoolMarker& bMarker)
{
	typedef typename Grid::traits<TElem>::iterator ElemIter;

//	Loop all elements and unmark H slaves
	for(ElemIter eIter = grid.begin<TElem>(); eIter != grid.end<TElem>(); ++eIter)
	{
		TElem* e = *eIter;
		if(grid.distributed_grid_manager()->contains_status(e, ES_H_SLAVE))
			bMarker.unmark(e);
		else
			bMarker.mark(e);
	}
}

/**
 * This function marks all elements for assembly BUT H-slaves in a BoolMarker and sets it
 * in the assemble tuner of the given assemble object. Necessary to avoid multiple assembly
 * of manifold elements when using inner_boundary discretization in parallel mode.
 *
 * \param[out]		grid		Smart pointer to assemble object with assemble tuner as member
 * \param[out]		ass			Grid
 */
template <typename TAlgebra>
static void MarkAllElemsForAssemblyButHSlaves(SmartPtr<IAssemble<TAlgebra> > ass, Grid& grid)
{
//	TODO: encapsulate BoolMarker in SmartPtr in ass_tuner.h

#ifdef UG_PARALLEL
	BoolMarker* p_bMarker = new BoolMarker(grid);

//	Unmark H slaves in all element types
	MarkAllElemsForAssemblyButHSlaves<Vertex>(grid, *p_bMarker);
	MarkAllElemsForAssemblyButHSlaves<Edge>(grid, *p_bMarker);
	MarkAllElemsForAssemblyButHSlaves<Face>(grid, *p_bMarker);
	MarkAllElemsForAssemblyButHSlaves<Volume>(grid, *p_bMarker);

//	Set marker in assembly tuner
	ass->ass_tuner()->set_marker(p_bMarker);
#endif
}

} // end namespace ug

#endif /*__H_MANIFOLD_ASSEMBLE_UTIL__*/
