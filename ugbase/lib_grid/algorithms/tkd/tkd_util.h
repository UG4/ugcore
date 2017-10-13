/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_tkd_util
#define __H__UG_tkd_util

#include "tkd_info.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/common_attachments.h"

namespace ug {

/** Creates a set of TKD vertices (inner and outer). This method is
 * typically only used internally.
 * \param vrtsOut	(optional) Array of size TKDInfo::NUM_COORDS*/
void CreateTKDVertices (const TKDInfo& tkdInfo,
                        Grid& g,
                        APosition3& aPos,
                        Vertex** vrtsOut = NULL);

///	Creates a tkd mesh based on the given tkdInfo object
/** \param volsOut	(optional) Array of size TKDInfo::NUM_INNER_ELEMENTS
 *					Created volumes are written to this array, if specified.*/
void CreateTKD (const TKDInfo& tkdInfo,
                Grid& g,
				APosition3& aPos,
				Volume** volsOut = NULL);

///	Creates a tkd mesh with a surrounding layer based on the given tkdInfo object
/** \param volsOut	(optional) Array of size TKDInfo::NUM_ELEMENTS.
 *					Created volumes are written to this array, if specified.
 *					The first TKDInfo::NUM_INNER_VERTICES volumes are inner volumes,
 *					the rest ar outer volumes.*/
void CreateTKDWithOuterLayer (const TKDInfo& tkdInfo,
	                          Grid& g,
	                          APosition3& aPos,
							  Volume** volsOut = NULL);

}//	end of namespace

#endif	//__H__UG_tkd_util
