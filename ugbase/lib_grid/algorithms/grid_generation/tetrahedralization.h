#ifndef __H__LIB_GRID__TETRAHEDRALIZATION__
#define __H__LIB_GRID__TETRAHEDRALIZATION__

#include "lib_grid/lg_base.h"

namespace ug
{

/// \addtogroup lib_grid_algorithms_grid_generation
///	@{

////////////////////////////////////////////////////////////////////////
///	fills a closed surface-grid with tetrahedrons.
/**	You may specify a quality parameter. If this parameter is <= 0, no
 *	inner vertices will be created. The default value for quality is 5.
 *	The quality resembles the minimal valid dihedral angle.
 *	The algorithm should always terminate for this quality. If you choose
 *	a lower quality parameter (careful with quality < 1), the algotithm may
 *	not terminate. The lower the quality parameter (but > 0), the
 *	better the tetrahedron quality.
 *
 *	Using tetgen by Hang Si.
 *
 *  \param verbosity	number between 0 and 3 indicating how detailed the
 *						verbosity should be
 *	\{
 */
bool Tetrahedralize(Grid& grid, number quality = 5,
					bool preserveBnds = false,
					bool preserveAll = false,
					APosition& aPos = aPosition,
					int verbosity = 0);

bool Tetrahedralize(Grid& grid, ISubsetHandler& sh,
					number quality = 5,
					bool preserveBnds = false,
					bool preserveAll = false,
					APosition& aPos = aPosition,
					int verbosity = 0);
///	\}

///	If tetrahedrons are already present, this method refines them based on the given volume constraints.
/**	A negative volume constraint implies no constraint for that element.*/
bool Retetrahedralize(Grid& grid, SubsetHandler& sh,
					ANumber& aVolumeConstraint,
					number quality = 5,
					bool preserveBnds = false,
					bool preserveAll = false,
					APosition& aPos = aPosition,
					bool applyVolumeConstraint = true,
					int verbosity = 0);
/**@}*/ // end of doxygen defgroup command

}//	end of namespace

#endif
