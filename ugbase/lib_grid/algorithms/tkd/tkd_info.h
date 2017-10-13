/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Martin Scherer, Rebecca Wittum, Sebastian Reiter
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

#ifndef __H__UG_tkd_info
#define __H__UG_tkd_info

#include <vector>
#include "common/math/ugmath_types.h"

namespace ug {

class TKDInfo
{
public:
	/**
	 * \param a	...?
	 * \param w ...?
	 * \param h		height of the (inner) TKD
	 * \param d		thickness of the outer layer
	 */
	TKDInfo (number a, number w, number h, number d);

	static const int NUM_COORDS			= 52;
	static const int NUM_INNER_COORDS	= 26;
	static const int NUM_OUTER_COORDS	= 26;

	static const int NUM_ELEMENTS		= 54;
	static const int NUM_INNER_ELEMENTS	= 18;
	static const int NUM_OUTER_ELEMENTS	= 36;

	size_t num_coords () const				{return (size_t)NUM_COORDS;}
	size_t num_inner_coords () const		{return (size_t)NUM_INNER_COORDS;}
	size_t num_outer_coords () const		{return (size_t)NUM_OUTER_COORDS;}

	const vector3* coords () const			{return &m_coords[0];}
	const vector3* inner_coords () const	{return &m_coords[0];}
	const vector3* outer_coords () const	{return &m_coords[26];}

	size_t num_elements () const {return NUM_ELEMENTS;}
	size_t num_inner_elements () const {return NUM_INNER_ELEMENTS;}
	size_t num_outer_elements () const {return NUM_OUTER_ELEMENTS;}

	/** returns a list of vertex numbers and vertex indices descriping elements
	 * of the tkd as follows:
	 * \code
	 * numVrts_0, vrt_0_0, vrt_0_1, ..., numVrts_1, vrt_1_0, vrt_1_1, ...
	 * \endcode
	 * numVrts can be mapped to element types since only tetrahedra, prisms, and
	 * hexahedra are used.
	 * \{ */
	const int* inner_element_indices () const;
	const int* outer_element_indices () const;
	/** \} */


private:
	inline void init_coords (	vector3* coordsOut,
								number a,
								number w,
								number h);

	std::vector<vector3> m_coords;
};

}//	end of namespace

#endif	//__H__UG_tkd_info
