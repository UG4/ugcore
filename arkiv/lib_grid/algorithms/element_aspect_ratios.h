øunused
/*
 * Copyright (c) 2012-2018:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

#ifndef __H__UG_ELEMENT_ASPECT_RATIOS
#define __H__UG_ELEMENT_ASPECT_RATIOS

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#endif

#include <lib_grid/algorithms/bounding_box_util.h>

/* system includes */
#include <stddef.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "lib_grid/lib_grid.h"


using namespace std;


namespace ug {


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinTriangleHeight
////////////////////////////////////////////////////////////////////////////////////////////

template <typename TAAPosVRT>
number CalculateMinTriangleHeight(Face* face, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
	if(face->num_vertices() == 3)
	{
	//	Get type of vertex attachment in aaPos and define it as ValueType
		using ValueType = typename TAAPosVRT::ValueType;

		number minHeight, tmpMinHeight;
		ValueType v = aaPos[face->vertex(2)];
		ValueType dir;

	//	Calculate start height and set to minHeight
		VecSubtract(dir, aaPos[face->vertex(1)], aaPos[face->vertex(0)]);
		minHeight = DistancePointToRay(	v, aaPos[face->vertex(0)], dir);

		for(uint i = 1; i < 3; ++i)
		{
			v = aaPos[face->vertex((i+2)%3)];
			VecSubtract(dir, aaPos[face->vertex((i+1)%3)], aaPos[face->vertex((i))]);
			tmpMinHeight = DistancePointToRay(v, aaPos[face->vertex(i )], dir);

			if(tmpMinHeight < minHeight)
			{
				minHeight = tmpMinHeight;
			}
		}

		return minHeight;
	}
	else
		UG_ASSERT(false, "Error. Face is not a triangle.");

	return NAN;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateAspectRatio
////////////////////////////////////////////////////////////////////////////////////////////

///	An unimplemented version, so that a compile error occurs if no overload exists.
template <typename TElem, typename TAAPosVRT>
number CalculateAspectRatio(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

/// Quadrilaterals
template <typename TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Quadrilateral* quad, TAAPosVRT& aaPos) {
	using box_t = AABox<typename TAAPosVRT::ValueType>;
	Vertex* const* vrts = quad->vertices();
	box_t box(aaPos[vrts[0]], aaPos[vrts[0]]);
	for (size_t i = 1; i < quad->num_vertices(); ++i){
		box = box_t(box, aaPos[vrts[i]]);
	}
	const number xDist = box.max.x() - box.min.x();
	const number yDist = box.max.y() - box.min.y();
	return xDist > yDist ? xDist/yDist : yDist/xDist;
}

///	Triangles and Constrained Triangles
template <typename TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Triangle* tri, TAAPosVRT& aaPos)
{
	number maxEdgeLength;

	//	Collect element edges, find longest edge and calculate its length
	vector<Edge*> edges;
	CollectAssociated(edges, grid, tri);
	Edge* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
	maxEdgeLength = EdgeLength(longestEdge, aaPos);

	/*
	 * optimal aspect ratio of a regular tetrahedron with edge lengths a:
	 * Q = hmin/lmax = sqrt(3)/2*a / a = 0.866...
	 *
	 * Info: return value is normalized by factor 2/sqrt(3)
	 * (s. Shewchuk, "What is a Good Linear Element? Interpolation, Conditioning, and Quality Measures?", 2002)
	 */

	//	Calculate minimal triangle height
	number minTriangleHeight = CalculateMinTriangleHeight(tri, aaPos);

	//	Calculate the aspect ratio
	return  2/std::sqrt(3.0) * minTriangleHeight / maxEdgeLength;
}

/// CalculateAspectRatio for faces
template <typename TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Face* face, TAAPosVRT& aaPos)
{
	switch (face->reference_object_id())
	{
		case ROID_TRIANGLE:
		{
			return CalculateAspectRatio(grid, static_cast<Triangle*>(face), aaPos);
		}

		case ROID_QUADRILATERAL:
		{
			return CalculateAspectRatio(grid, static_cast<Quadrilateral*>(face), aaPos);
		}

		default:
		{
		 	UG_THROW("Note: Currently only faces of type triangle and quadrilateral"
		 			" supported in aspect ratio calculation.");
		 	break;
		}
	}
	return NAN;
}

///	Tetrahedron
template <typename TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	/*
	 * optimal Aspect Ratio of a regular tetrahedron with edge lengths a:
	 * Q = hmin/lmax = sqrt(2/3)*a / a = 0.8164...
	 *
	 * Info: return value is normalized by factor sqrt(3/2)
	 * (s. Shewchuk, "What is a Good Linear Element? Interpolation, Conditioning, and Quality Measures?", 2002)
	 */
	return CalculateTetrahedronAspectRatio(grid, tet, aaPos);
}

/// Hexahedron
template <typename TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Hexahedron* hex, TAAPosVRT& aaPos)
{
	return CalculateHexahedronAspectRatio(grid, hex, aaPos);
}

/// Pyramid
template <typename TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Pyramid* pyr, TAAPosVRT& aaPos) {
	return CalculatePyramidAspectRatio(grid, pyr, aaPos);
}

///	Volume
template <typename TAAPosVRT>
number CalculateAspectRatio(Grid& grid, Volume* vol, TAAPosVRT& aaPos)
{
	switch (vol->reference_object_id())
	{
		case ROID_TETRAHEDRON:
		{
			return CalculateAspectRatio(grid, static_cast<Tetrahedron*>(vol), aaPos);
		}

		case ROID_HEXAHEDRON:
		{
			return CalculateAspectRatio(grid, static_cast<Hexahedron*>(vol), aaPos);
		}

		case ROID_PYRAMID:
		{
			return CalculateAspectRatio(grid, static_cast<Pyramid*>(vol), aaPos);
		}

		default:
		{
		 	UG_THROW("Note: Currently only volumes of type tetrahedron, hexahedron"
		 			"and pyramids supported in aspect ratio calculation.");
		 	break;
		}
	}

	return NAN;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateVolToRMSFaceAreaRatio
////////////////////////////////////////////////////////////////////////////////////////////

///	An unimplemented version, so that a compile error occurs if no overload exists.
template <typename TElem, typename TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, TElem* elem, TAAPosVRT& aaPos);

///	Face (Triangles and Constrained Triangles supported)
template <typename TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, Face* face, TAAPosVRT& aaPos)
{
	UG_THROW("CalculateVolToRMSFaceAreaRatio: Currently only volumes of type tetrahedron"
			" supported in volume to root-mean-square face area ratio calculation.");
	return NAN;
}

///	Tetrahedron
template <typename TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, Tetrahedron* tet, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
	number ratio;

	/*
	 * optimal volume to root-mean-square face area ratio of a
	 * regular tetrahedron with edge lengths a:
	 * Q = V/A_rms^(3/2)
	 *
	 * Info: return value is normalized by factor pow(3, 7/4.0) / 2.0 / sqrt(2);
	 * (s. Shewchuk, "What is a Good Linear Element? Interpolation, Conditioning, and Quality Measures?", 2002)
	 */

	//	Calculate the ratio
	ratio = CalculateTetrahedronVolToRMSFaceAreaRatio(grid, tet, aaPos);

	return ratio;
}

/// Hexahedron
template <typename TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, Hexahedron* hex, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
	number ratio;

	//	Calculate the ratio
	ratio = CalculateHexahedronVolToRMSFaceAreaRatio(grid, hex, aaPos);

	return ratio;
}

///	Volume
template <typename TAAPosVRT>
number CalculateVolToRMSFaceAreaRatio(Grid& grid, Volume* vol, TAAPosVRT& aaPos)
{
	switch (vol->reference_object_id())
	{
		case ROID_TETRAHEDRON:
		{
			return CalculateVolToRMSFaceAreaRatio(grid, static_cast<Tetrahedron*>(vol), aaPos);
		}

		case ROID_HEXAHEDRON:
		{
			return CalculateVolToRMSFaceAreaRatio(grid, static_cast<Hexahedron*>(vol), aaPos);
		}

		default:
		{
			UG_THROW("CalculateVolToRMSFaceAreaRatio: Currently only volumes of type tetrahedron"
					" supported in volume to root-mean-square face area ratio calculation.");
		 	break;
		}
	}

	return NAN;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindLargestFace
template <typename TIterator, typename TAAPosVRT>
Face* FindLargestFace(TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
	//	if facesBegin equals facesEnd, then the list is empty and we can
	//	immediately return nullptr
	if(facesBegin == facesEnd)
		return nullptr;

	//	the first face is the first candidate for the smallest face.
		Face* largestFace = *facesBegin;
		number largestArea = FaceArea(largestFace, aaPos);
		++facesBegin;

		for(; facesBegin != facesEnd; ++facesBegin){
			Face* curFace = *facesBegin;
			number curArea = FaceArea(curFace, aaPos);
			if(curArea > largestArea){
				largestFace = curFace;
				largestArea = curArea;
			}
		}

		return largestFace;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindSmallestVolumeElement
template <typename TIterator, typename TAAPosVRT>
typename TIterator::value_type
FindSmallestVolume(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return nullptr
	if(volumesBegin == volumesEnd)
		return nullptr;

//	Initializations
	typename TIterator::value_type smallestVolume = *volumesBegin;
	number smallestVolumeVolume = CalculateVolume(smallestVolume, aaPos);
	++volumesBegin;

//	compare all tetrahedrons and find minimal volume
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curVolumeVolume = CalculateVolume(curVolume, aaPos);

		if(curVolumeVolume < smallestVolumeVolume)
		{
			smallestVolume = curVolume;
			smallestVolumeVolume = curVolumeVolume;
		}
	}

	return smallestVolume;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindLargestVolumeElement
template <typename TIterator, typename TAAPosVRT>
typename TIterator::value_type
FindLargestVolume(TIterator volumesBegin, TIterator volumesEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return nullptr
	if(volumesBegin == volumesEnd)
		return nullptr;

//	Initializations
	typename TIterator::value_type largestVolume = *volumesBegin;
	number largestVolumeVolume = CalculateVolume(largestVolume, aaPos);
	++volumesBegin;

//	compare all tetrahedrons and find minimal volume
	for(; volumesBegin != volumesEnd; ++volumesBegin)
	{
		Volume* curVolume = *volumesBegin;
		number curVolumeVolume = CalculateVolume(curVolume, aaPos);

		if(curVolumeVolume > largestVolumeVolume)
		{
			largestVolume = curVolume;
			largestVolumeVolume = curVolumeVolume;
		}
	}

	return largestVolume;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithSmallestAspectRatio
template <typename TIterator, typename TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestAspectRatio(Grid& grid, 	TIterator elemsBegin,
												TIterator elemsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return nullptr
	if(elemsBegin == elemsEnd)
		return nullptr;

//	Initializations
	typename TIterator::value_type elementWithSmallestAspectRatio = *elemsBegin;
	number smallestAspectRatio = CalculateAspectRatio(grid, elementWithSmallestAspectRatio, aaPos);
	++elemsBegin;

//	compare all tetrahedrons and find that one with minimal aspect ratio
	for(; elemsBegin != elemsEnd; ++elemsBegin)
	{
		typename TIterator::value_type curElement = *elemsBegin;
		//TElem* curElement = *elemsBegin;
		number curSmallestAspectRatio = CalculateAspectRatio(grid, curElement, aaPos);

		if(curSmallestAspectRatio < smallestAspectRatio)
		{
			elementWithSmallestAspectRatio = curElement;
			smallestAspectRatio = curSmallestAspectRatio;
		}
	}

	return elementWithSmallestAspectRatio;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithLargestAspectRatio
template <typename TIterator, typename TAAPosVRT>
typename TIterator::value_type
FindElementWithLargestAspectRatio(Grid& grid,  	TIterator elemsBegin,
												TIterator elemsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return nullptr
	if(elemsBegin == elemsEnd)
		return nullptr;

//	Initializations
	typename TIterator::value_type elementWithLargestAspectRatio = *elemsBegin;
	number largestAspectRatio = CalculateAspectRatio(grid, elementWithLargestAspectRatio, aaPos);
	++elemsBegin;

//	compare all tetrahedrons and find that one with maximal aspect ratio
	for(; elemsBegin != elemsEnd; ++elemsBegin)
	{
		typename TIterator::value_type curElement = *elemsBegin;
		//TElem* curElement = *elemsBegin;
		number curSmallestAspectRatio = CalculateAspectRatio(grid, curElement, aaPos);

		if(curSmallestAspectRatio > largestAspectRatio)
		{
			elementWithLargestAspectRatio = curElement;
			largestAspectRatio = curSmallestAspectRatio;
		}
	}

	return elementWithLargestAspectRatio;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithSmallestVolToRMSFaceAreaRatio
template <typename TIterator, typename TAAPosVRT>
typename TIterator::value_type
FindElementWithSmallestVolToRMSFaceAreaRatio(Grid& grid, TIterator elemsBegin,
											 TIterator elemsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return nullptr
	if(elemsBegin == elemsEnd)
		return nullptr;

//	Initializations
	typename TIterator::value_type elementWithSmallestRatio = *elemsBegin;
	number smallestRatio = CalculateVolToRMSFaceAreaRatio(grid, elementWithSmallestRatio, aaPos);
	++elemsBegin;

//	compare all tetrahedrons and find that one with minimal aspect ratio
	for(; elemsBegin != elemsEnd; ++elemsBegin)
	{
		typename TIterator::value_type curElement = *elemsBegin;
		//TElem* curElement = *elemsBegin;
		number curSmallestRatio = CalculateVolToRMSFaceAreaRatio(grid, curElement, aaPos);

		if(curSmallestRatio < smallestRatio)
		{
			elementWithSmallestRatio = curElement;
			smallestRatio = curSmallestRatio;
		}
	}

	return elementWithSmallestRatio;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	FindElementWithLargestVolToRMSFaceAreaRatio
template <typename TIterator, typename TAAPosVRT>
typename TIterator::value_type
FindElementWithLargestVolToRMSFaceAreaRatio(Grid& grid, TIterator elemsBegin,
											TIterator elemsEnd, TAAPosVRT& aaPos)
{
	//PROFILE_FUNC();
//	if volumesBegin equals volumesBegin, then the list is empty and we can
//	immediately return nullptr
	if(elemsBegin == elemsEnd)
		return nullptr;

//	Initializations
	typename TIterator::value_type elementWithLargestRatio = *elemsBegin;
	number largestRatio = CalculateVolToRMSFaceAreaRatio(grid, elementWithLargestRatio, aaPos);
	++elemsBegin;

//	compare all tetrahedrons and find that one with maximal aspect ratio
	for(; elemsBegin != elemsEnd; ++elemsBegin)
	{
		typename TIterator::value_type curElement = *elemsBegin;
		//TElem* curElement = *elemsBegin;
		number curSmallestRatio = CalculateVolToRMSFaceAreaRatio(grid, curElement, aaPos);

		if(curSmallestRatio > largestRatio)
		{
			elementWithLargestRatio = curElement;
			largestRatio = curSmallestRatio;
		}
	}

	return elementWithLargestRatio;
}
}	 
#endif
