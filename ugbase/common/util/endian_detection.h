/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__COMMON__UTIL__ENDIAN__
#define __H__UG__COMMON__UTIL__ENDIAN__

namespace ug{

/// \addtogroup ugbase_common_util
/// \{

/**
 * Returns true if the system endianess is little endian, i.e. the least
 * significant byte is stored first. If returned false the system endianess is
 * big endian, i.e. most significant byte first. Please note, that this is a
 * runtime check.
 *
 * @return	true if little endian, false if big endian
 */
inline bool IsLittleEndian()
{
    short int number = 0x1;
    char *numPtr = (char*)&number;
    return (numPtr[0] == 1);
}

/**
 * Returns false if the system endianess is little endian, i.e. the least
 * significant byte is stored first. If returned true the system endianess is
 * big endian, i.e. most significant byte first. Please note, that this is a
 * runtime check.
 *
 * @return	false if little endian, true if big endian
 */
inline bool IsBigEndian()
{
	return !IsLittleEndian();
}

// end group ugbase_common_util
/// \}

} // end namespace ug

#endif
