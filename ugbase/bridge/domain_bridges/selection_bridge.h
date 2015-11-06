/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Stephan Grein
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

#include <vector>
#include <string>
#include <sstream>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"

using namespace std;

namespace ug{

/**
 * \defgroup selection_bridge Selection Bridge
 * \ingroup domain_bridge
 * \{
 */

////////////////////////////////////////////////////////////////////////////////
///	Selects / Deselects all elements
void SelectDomainElements(ISelector& sel, bool bSelect, bool selectVrts,
				  bool selectEdges, bool selectFaces, bool selectVolumes);


////////////////////////////////////////////////////////////////////////////////
///	Selects / Deselects associated (low dimensional) elements.
void SelectAssociatedElements(ISelector& sel, bool bSelect, bool selectVrts,
				  	  	  	  bool selectEdges, bool selectFaces);


////////////////////////////////////////////////////////////////////////////////
///	Selects / Deselects elements in the given subset of the given domain
template <class TDomain>
void SelectDomainSubset(ISelector& sel, TDomain& dom, int subsetIndex,
						bool bSelect, bool selectVrts, bool selectEdges,
						bool selectFaces, bool selectVolumes);
}
