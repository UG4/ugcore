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

#ifndef __H__UG__volume_calculation_impl__
#define __H__UG__volume_calculation_impl__

#include "common/math/misc/math_util.h"
#include "geom_obj_util/face_util.h"

namespace ug{


template <typename TAAPos>
number CalculateVolume(Volume* elem, TAAPos aaPos)
{
	switch (elem->reference_object_id()) {
	case ROID_TETRAHEDRON:
		return CalculateVolume(static_cast<Tetrahedron*>(elem), aaPos);
	case ROID_PRISM:
		return CalculateVolume(static_cast<Prism*>(elem), aaPos);
	case ROID_PYRAMID:
		return CalculateVolume(static_cast<Pyramid*>(elem), aaPos);
	case ROID_HEXAHEDRON:
		return CalculateVolume(static_cast<Hexahedron*>(elem), aaPos);
	case ROID_OCTAHEDRON:
		return CalculateVolume(static_cast<Octahedron*>(elem), aaPos);
	default:
		UG_THROW("Unknown volume type");
		break;
	}

	return NAN;
}

template <typename TAAPos>
number CalculateVolume(Tetrahedron* elem, TAAPos aaPos)
{
	return CalculateTetrahedronVolume(aaPos[elem->vertex(0)],
									aaPos[elem->vertex(1)],
									aaPos[elem->vertex(2)],
									aaPos[elem->vertex(3)]);
}

template <typename TAAPos>
number CalculateVolume(Pyramid* elem, TAAPos aaPos)
{
	return CalculatePyramidVolume(aaPos[elem->vertex(0)],
								aaPos[elem->vertex(1)],
								aaPos[elem->vertex(2)],
								aaPos[elem->vertex(3)],
								aaPos[elem->vertex(4)]);
}

template <typename TAAPos>
number CalculateVolume(Prism* elem, TAAPos aaPos)
{
	return CalculatePrismVolume(aaPos[elem->vertex(0)],
								aaPos[elem->vertex(1)],
								aaPos[elem->vertex(2)],
								aaPos[elem->vertex(3)],
								aaPos[elem->vertex(4)],
								aaPos[elem->vertex(5)]);
}

template <typename TAAPos>
number CalculateVolume(Hexahedron* elem, TAAPos aaPos)
{
	return CalculateHexahedronVolume(aaPos[elem->vertex(0)],
									aaPos[elem->vertex(1)],
									aaPos[elem->vertex(2)],
									aaPos[elem->vertex(3)],
									aaPos[elem->vertex(4)],
									aaPos[elem->vertex(5)],
									aaPos[elem->vertex(6)],
									aaPos[elem->vertex(7)]);
}

template <typename TAAPos>
number CalculateVolume(Octahedron* elem, TAAPos aaPos)
{
	return CalculateOctahedronVolume(aaPos[elem->vertex(0)],
								aaPos[elem->vertex(1)],
								aaPos[elem->vertex(2)],
								aaPos[elem->vertex(3)],
								aaPos[elem->vertex(4)],
								aaPos[elem->vertex(5)]);
}

template <typename TAAPos>
number CalculateVolume(FaceVertices* elem, TAAPos aaPos)
{
	return FaceArea(elem, aaPos);
}


template <typename TAAPos>
number CalculateVolume(EdgeVertices* elem, TAAPos aaPos)
{
	return EdgeLength(elem, aaPos);
}

template <typename TAAPos>
number CalculateVolume(Vertex*, TAAPos)
{
	return 0;
}

template <typename TIterator, typename TAAPos>
number CalculateVolume(TIterator begin, TIterator end, TAAPos aaPos)
{
	number totalVolume = 0;
	for(TIterator iter = begin; iter != end; ++iter){
		totalVolume += CalculateVolume(*iter, aaPos);
	}
	return totalVolume;
}

//! Determine the bounding box for a set of points.
template<int dim>
void CalculateBoundingBox(size_t npoints, const MathVector<dim> points[], MathVector<dim> &vMinBB, MathVector<dim> &vMaxBB)
{
	// determine bounding box
	vMinBB= points[0];
	vMaxBB = points[0];

	for(size_t ii = 1; ii < npoints; ++ii)
	{
		for(int i = 0; i < dim; ++i)
		{
			const MathVector<dim>& v = points[ii];
			if(v[i] < vMinBB[i]) vMinBB[i] = v[i];
			else if(v[i] > vMaxBB[i]) vMaxBB[i] = v[i];
		}
	}
}

}// end of namespace

#endif
