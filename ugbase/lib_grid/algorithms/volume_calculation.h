/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Martin Stepniewski, Martin Scherer
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

#ifndef __H__UG__volume_calculation__
#define __H__UG__volume_calculation__

#include "lib_grid/grid_objects/grid_objects_3d.h"

namespace ug {

/**
 * \brief	contains methods for the calculation of the volumes of
 * 			edge, face and volume elements
 *
 * \defgroup lib_grid_algorithms_volume_calculation volume_calculation
 * \ingroup lib_grid_algorithms
 * \{
 */

///	Calculates the volume of the given element.
/**	Make sure that aaPos has at least the same dimensionality as the specified element!
 * \note	If Volume-elements feature bended sides, the returned value is only an
 * 			approximation to the actual volume.
 * \{ */
//	VOLUMES
template <typename TAAPos>
inline
number CalculateVolume(Volume* elem, TAAPos aaPos);

template <typename TAAPos>
inline
number CalculateVolume(Tetrahedron* elem, TAAPos aaPos);

template <typename TAAPos>
inline
number CalculateVolume(Pyramid* elem, TAAPos aaPos);

template <typename TAAPos>
inline
number CalculateVolume(Prism* elem, TAAPos aaPos);

template <typename TAAPos>
inline
number CalculateVolume(Hexahedron* elem, TAAPos aaPos);

template <typename TAAPos>
inline
number CalculateVolume(Octahedron* elem, TAAPos aaPos);

//	FACES
template <typename TAAPos>
inline
number CalculateVolume(FaceVertices* elem, TAAPos aaPos);

//	EDGES
template <typename TAAPos>
inline
number CalculateVolume(EdgeVertices* elem, TAAPos aaPos);

template <typename TAAPos>
inline
number CalculateVolume(Vertex* elem, TAAPos aaPos);

template <typename TIterator, typename TAAPos>
inline
number CalculateVolume(TIterator begin, TIterator end, TAAPos aaPos);
/** \} */

template<int dim>
void CalculateBoundingBox(size_t npoints, const MathVector<dim> points[], MathVector<dim> &vMinBB, MathVector<dim> &vMaxBB);

/** \} */	//end of ingroup

}// end of namespace


////////////////////////////////////////
//	include implementation
#include "volume_calculation_impl.hpp"

#endif
