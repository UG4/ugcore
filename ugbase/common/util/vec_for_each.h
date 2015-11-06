/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_for_each_in_vec
#define __H__UG_for_each_in_vec


///	Allows iteration over all members of an std::vector compatible type
/**	Use e.g. like this:
 * \code
 * std::vector<T>	vec;
 * //...
 * for_each_in_vec(T& t, vec){
 *    t.memberFunc();
 * }end_for
 * \endcode
 * The '{' and '}' brackets are hereby optional, since a new block is automatically
 * defined between 'for_each_in_vec' and 'end_for'.
 *
 * The specified vector has to feature methods 'T& operator[](size_t i)' and
 * 'size_t size()'.
 * \{ */
#define for_each_in_vec(_vfeDecl, _vfeVec) \
			for(size_t _vfeI = 0; _vfeI < _vfeVec.size(); ++_vfeI){\
				_vfeDecl = _vfeVec[_vfeI];

#define end_for	}
/** \} */

#endif	//__H__UG_for_each_in_vec
