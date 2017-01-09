/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
