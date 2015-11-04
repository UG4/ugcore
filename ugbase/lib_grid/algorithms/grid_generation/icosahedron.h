#ifndef __H__UG__icosahedron__
#define __H__UG__icosahedron__

#include "lib_grid/lg_base.h"

namespace ug
{

///	Creates an Icosahedron with the given radius. (outer circle)
void GenerateIcosahedron(Grid& grid, const vector3& center,
						 number radius, AVector3& aPos);

/// Creates a ico-sphere by repeatedly refining an icosahedron
/**	Make sure not to choose a too high number of refinements. Number of triangles
 * produced equals 20 * 4^numRefinements.
 *
 * You may optionally specify a selector. If you do so, the selector will be
 * used for internal calculations. The whole sphere will be selected when the
 * algorithm is done (all vertices, edges and faces).
 *
 * If you won't specify a selector, an internal selector has to be created,
 * which introduces a runtime overhead of O(n). This could be avoided by a
 * more elaborate implementation.
 */
void GenerateIcosphere(Grid& grid, const vector3& center, number radius,
						int numRefinements, AVector3& aPos, Selector* psel = NULL);

}//	end of namespace

#endif
