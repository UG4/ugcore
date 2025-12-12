/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__registry_util__
#define __H__UG__registry_util__

#include "common/util/string_util.h"

namespace ug {

/// \addtogroup registry
/// \{

/**
 * Checks whether the specified name is a valid registry identifier name.
 * <p>
 * <b>Note:</b> identifiers starting with <code>F_</code>, <code>C_</code>,
 * <code>I_</code> or containing <code>__</code>
 *  (double underscore) or being equal to 'constructor' are invalid.
 * </p>
 * @param name name to check
 * @return  <code>true</code> if the specified name is valid;
 *          <code>false</code> otherwise
 */
UG_API bool IsValidRegistryIdentifier(const std::string& name);

/**
 * Returns a message describing which registry identifiers are valid and which are not.
 * @return message string
 */
UG_API std::string GetRegistryIdentifierMessage();

// end group registry
/// \}

}//	end of namespace

#endif
