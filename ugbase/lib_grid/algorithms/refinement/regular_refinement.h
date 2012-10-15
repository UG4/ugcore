// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m02 d09

#ifndef __H__LIB_GRID__REGULAR_REFINER__
#define __H__LIB_GRID__REGULAR_REFINER__

#include "lib_grid/lg_base.h"
#include "refinement_callbacks.h"
#include "common/ug_config.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

////////////////////////////////////////////////////////////////////////
//	Refine
///	refines selected faces and edges regularily and builds a closure on adjacent unselected faces.
/**
 * Selected faces will be refined with regular refinement.
 * Faces which are adjacent to selected faces or selected edges and
 * which are not selected, will be refined too to build a closure.
 * In this case however, the refinement is not regular - faces are
 * refined in a way to avoid hanging nodes.
 *
 * Pass a grid and a selector (which is working on the grid). The
 * provided aInt is required by the algorithm to store temporary
 * values. You have to pass it to the algorithm to allow
 * maximal speed for repeated small refinements.
 * 
 * All involved geometric objects will be selected after the method
 * terminated (This includes vertices adjacent to selected edges and faces).
 *
 * aInt has to be attached to the edges of the grid.
 *
 * If you are interested in rare, big refinements, you may also use
 * the overloaded version of Refine, which only takes a Grid and a 
 * Selector.
 *
 * \sa ug::RegularRefiner, ug::HangingNodeRefiner
 */
UG_API
bool Refine(Grid& grid, Selector& sel, AInt& aInt,
			IRefinementCallback* refCallback = NULL);

///	refines selected faces and edges regularily and builds a closure on adjacent unselected faces.
/**
 * This method overloads Refine(Grid& grid, Selector& sel, AInt& aInt).
 * It is slower that the full version of refine, since it tempiroraily
 * attaches the required aInt attachment to the edges of aInt before
 * it calls the original version.
 *
 * This method should only be used if only very few refinement steps
 * are performed and if speed is not cruical.
 *
 * \sa ug::RegularRefiner, ug::HangingNodeRefiner
 */
UG_API
bool Refine(Grid& grid, Selector& sel,
			IRefinementCallback* refCallback = NULL);

/// @}
}// end of namespace

#endif
