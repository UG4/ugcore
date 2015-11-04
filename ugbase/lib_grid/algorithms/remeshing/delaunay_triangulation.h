// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.09.2011 (m,d,y)

#ifndef __H__UG__delaunay_triangulation__
#define __H__UG__delaunay_triangulation__

#include <queue>
#include <vector>
#include <sstream>
#include "delaunay_info.h"
#include "common/ug_config.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"
#include "lib_grid/callbacks/basic_callbacks.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
template <class TAAPos>
bool MakeDelaunay(DelaunayInfo<TAAPos>& info);

////////////////////////////////////////////////////////////////////////////////
///	Transforms the given triangle-set into a delaunay set
/** Creates a delaunay triangulation. If a minAngle greater than 0 is specified,
 * then additional vertices are introduced, if required to generate fulfill
 * the min-angle-condition.
 *
 * Make sure that init_marks was performed on the specified deaunay-info object.
 * If the triangulation is not already delaunay, all interior edges should also
 * be pushed to the candidates pool (can be done automatically in DelaunayInfo::init_marks)
 */
template <class TAAPos>
bool QualityGridGeneration(Grid& grid, DelaunayInfo<TAAPos>& info,
						   number minAngle = 0,
				  	  	   int maxSteps = -1/*remove this*/);

template <class TriIter, class TAAPos>
bool QualityGridGeneration(Grid& grid, TriIter trisBegin, TriIter trisEnd,
						   TAAPos& aaPos, number minAngle = 0,
				  	  	   Grid::edge_traits::callback cbConstrainedEdge = ConsiderNone(),
				  	  	   int maxSteps = -1/*remove this*/)
{
//	set up a delaunay-info structure
	DelaunayInfo<TAAPos> info(grid, aaPos, cbConstrainedEdge);
	info.init_marks(trisBegin, trisEnd, true);
	return QualityGridGeneration(grid, info, minAngle, maxSteps);	
}

}//	end of namespace

#endif
