/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

///	Generates a list of triangle-corners for the given Ico-Sphere
/**
 * \param trisOut	std::vector of points. Contents are not cleared!
 *					Computed vertices are written to this vector.
 *					Three consecutive points form a triangle.
 * \param center	Center of the sphere.
 * \param radius	Radius of the sphere.
 * \param refinements	number of refinements used to generate the Icosphere.
 */
void GenerateIcosphere(	std::vector<vector3>& trisOut,
						const vector3& center,
						number radius,
						size_t refinements);
}//	end of namespace

#endif
