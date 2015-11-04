#ifndef __H__UG__check_associated_elements__
#define __H__UG__check_associated_elements__

#include "lib_grid/lg_base.h"

namespace ug{
namespace grid_unit_tests{
/** THIS METHOD USES Grid::mark!
 * iterates over all volumes in g and checks for each whether the edges stored
 * in the associated container are the same as the edges which can be found
 * in the grid and which are contained by the volume...
 *
 * If something is wrong, the method throws an instance of UGError.
 */
void CheckAssociatedEdgesOfVolumes(Grid& g);

/** THIS METHOD USES Grid::mark!
 * iterates over all edges in g and checks for each whether the volumes stored
 * in the associated container are the same as the edges which can be found
 * in the grid and which are contained by the volume...
 *
 * If something is wrong, the method throws an instance of UGError.
 */
void CheckAssociatedVolumesOfEdges(Grid& g);
}//	end of namespace
}//	end of namespace

#endif
