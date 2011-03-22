//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d17

#ifndef __H__LIB_GRID__TETRAHEDRALIZATION__
#define __H__LIB_GRID__TETRAHEDRALIZATION__

#include "lib_grid/lg_base.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	fills a closed surface-grid with tetrahedrons.
/**	You may specify a quality parameter. If this parameter is <= 0, no
 *	inner vertices will be created. The default value for quality is 5.
 *	The quality resembles the minimal valid dihedral angle.
 *	The algorithm should always terminate for this quality. If you choose
 *	a lower quality parameter (careful with quality < 1), the algotithm may
 *	not terminate. The lower the quality parameter (but > 0), the
 *	better the tetrahedron quality.
 *	\{
 */
bool Tetrahedralize(Grid& grid, number quality = 5,
					bool preserveBnds = false,
					bool preserveAll = false,
					APosition& aPos = aPosition);

bool Tetrahedralize(Grid& grid, SubsetHandler& sh,
					number quality = 5,
					bool preserveBnds = false,
					bool preserveAll = false,
					APosition& aPos = aPosition);
///	\}

}//	end of namespace

#endif
