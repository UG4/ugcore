/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__MISC_UTIL__IMPL__
#define __H__LIB_GRID__MISC_UTIL__IMPL__

#include "misc_util.h"
#include "lib_grid/grid/grid_util.h"
#include "vertex_util.h"
#include "edge_util.h"
#include "face_util.h"
#include "volume_util.h"

#include <limits>

namespace ug
{

////////////////////////////////////////////////////////////////////////////////////////////
template <typename TElem>
UG_API
typename TElem::side*
GetSharedSide(Grid& grid, TElem* e1, TElem* e2)
{
	if(!TElem::HAS_SIDES)
		return nullptr;

	using side_t = typename TElem::side;
	typename Grid::traits<side_t>::secure_container	sides1;
	typename Grid::traits<side_t>::secure_container	sides2;

	grid.associated_elements(sides1, e1);
	grid.associated_elements(sides2, e2);

	for(size_t i1 = 0; i1 < sides1.size(); ++i1){
		for(size_t i2 = 0; i2 < sides2.size(); ++i2){
			if(sides1[i1] == sides2[i2])
				return sides1[i1];
		}
	}

	return nullptr;
}


////////////////////////////////////////////////////////////////////////
//	CalculateGridObjectCenter
template <typename TAAPosVRT>
UG_API
inline
typename TAAPosVRT::ValueType
CalculateGridObjectCenter(const GridObject* o, TAAPosVRT& aaPosVRT)
{
	switch(o->base_object_id()){
		case VERTEX:	return CalculateCenter(static_cast<const Vertex*>(o), aaPosVRT);
		case EDGE:		return CalculateCenter(static_cast<const Edge*>(o), aaPosVRT);
		case FACE:		return CalculateCenter(static_cast<const Face*>(o), aaPosVRT);
		case VOLUME:	return CalculateCenter(static_cast<const Volume*>(o), aaPosVRT);
		default:
			UG_THROW("Unknown geometric-object type.");
	}
}

template<typename TAAPosVRT, typename TAAWeightVRT>
UG_API
inline
typename TAAPosVRT::ValueType
CalculateGridObjectCenter(const GridObject* o, TAAPosVRT& aaPosVRT,
							   TAAWeightVRT& aaWeight)
{
	switch(o->base_object_id()){
		case VERTEX:
			return CalculateCenter(static_cast<const Vertex*>(o), aaPosVRT, aaWeight);
		case EDGE:
			return CalculateCenter(static_cast<const Edge*>(o), aaPosVRT, aaWeight);
		case FACE:
			return CalculateCenter(static_cast<const Face*>(o), aaPosVRT, aaWeight);
		case VOLUME:
			return CalculateCenter(static_cast<const Volume*>(o), aaPosVRT, aaWeight);
		default:
			UG_THROW("Unknown geometric-object type.");
	}
}


////////////////////////////////////////////////////////////////////////
//	CalculateCenter
template <typename TIterator, typename TAAPosVRT>
typename TAAPosVRT::ValueType
CalculateCenter(TIterator begin, TIterator end, TAAPosVRT& aaPos)
{
	int counter = 0;
	typename TAAPosVRT::ValueType center;
	VecSet(center, 0);
	for(TIterator iter = begin; iter != end; ++iter, ++counter)
		VecAdd(center, center, CalculateCenter(*iter, aaPos));
		
	if(counter > 0)
		VecScale(center, center, 1./(number)counter);
		
	return center;
}

////////////////////////////////////////////////////////////////////////
//	FindClosestByCoordinate
template<typename TElem, typename TVertexPositionAttachmentAccessor>
TElem* FindClosestByCoordinate(const typename TVertexPositionAttachmentAccessor::ValueType& coord,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						TVertexPositionAttachmentAccessor& aaPosVRT)
{
	if(iterBegin == iterEnd)
		return nullptr;

	typename geometry_traits<TElem>::iterator iter = iterBegin;
	TElem* bestElem = *iter;
	number bestDistSq = VecDistanceSq(coord, CalculateCenter(bestElem, aaPosVRT));
	iter++;

	while(iter != iterEnd)
	{
		number distSq = VecDistanceSq(coord, CalculateCenter(*iter, aaPosVRT));
		if(distSq < bestDistSq)
		{
			bestDistSq = distSq;
			bestElem = *iter;
		}

		++iter;
	}

	return bestElem;
}

////////////////////////////////////////////////////////////////////////
template<typename vector_t, typename TIterator, typename TAAPos>
void CalculateBoundingBox(vector_t& vMinOut, vector_t& vMaxOut,
						  TIterator begin, TIterator end,
						  TAAPos& aaPos)
{
	if(begin == end){
		VecSet(vMinOut, 0);
		VecSet(vMaxOut, 0);
		return;
	}

	vMinOut = vMaxOut = aaPos[GetVertex(*begin, 0)];

	const int dim = vector_t::Size;

//	iterate over all elements and find min and max values
	vector_t tmin, tmax;
	for(TIterator iter = begin; iter != end; ++iter){
		typename TIterator::value_type elem = *iter;

		for(size_t i_vrt = 0; i_vrt < NumVertices(elem); ++i_vrt){
			for(int i = 0; i < dim; ++i){
				vector_t& v = aaPos[GetVertex(elem, i_vrt)];
				if(v[i] < vMinOut[i])
					vMinOut[i] = v[i];
				else if(v[i] > vMaxOut[i])
					vMaxOut[i] = v[i];
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	NumSharedVertices
template <typename TElemPtr1, typename TElemPtr2>
size_t NumSharedVertices(Grid& grid, TElemPtr1 elem1, TElemPtr2 elem2)
{
	grid.begin_marking();
//	first mark all vertices of elem1
	for(size_t i = 0; i < elem1->num_vertices(); ++i)
		grid.mark(elem1->vertex(i));

//	now count how many of vertex 2 are marked.
	size_t counter = 0;
	for(size_t i = 0; i < elem2->num_vertices(); ++i){
		if(grid.is_marked(elem2->vertex(i)))
			++counter;
	}
	
	grid.end_marking();
	
	return counter;
}

////////////////////////////////////////////////////////////////////////
//	EraseElements
template <typename TElem>
void EraseElements(Grid& grid, typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd)
{
//	be careful to not invalidate the iterators.
	while(iterBegin != iterEnd)
	{
		TElem* e = *iterBegin;
		iterBegin++;
		grid.erase(e);
	}
}

////////////////////////////////////////////////////////////////////////
//	AssignIndices
template <typename TElem>
void AssignIndices(typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					Grid::AttachmentAccessor<TElem, AInt>& aaInt)
{
	int index = 0;
	while(iterBegin != iterEnd)
	{
		aaInt[*iterBegin] = index++;
		iterBegin++;
	}
}

////////////////////////////////////////////////////////////////////////
//	ElementDiameter
template <typename TElem, typename TAAPos>
number ElementDiameterSq(Grid& grid,
                         TAAPos& aaPos,
					     TElem* elem)
{
	PointerConstArray<Vertex*> vVert;
	grid.associated_elements(vVert, elem);

	number max = 0.0;
	for(size_t i = 0; i < vVert.size(); ++i)
		for(size_t j = i+1; j < vVert.size(); ++j)
			max = std::max(max, VecDistanceSq(aaPos[vVert[i]], aaPos[vVert[j]]));

	return max;
}

template <typename TAAPos>
number ElementDiameterSq(Grid& grid,
						 TAAPos& aaPos,
						 GridObject* elem)
{
	switch(elem->base_object_id()){
		case VERTEX: return ElementDiameterSq(grid, aaPos, static_cast<Vertex*>(elem));
		case EDGE: return ElementDiameterSq(grid, aaPos, static_cast<Edge*>(elem));
		case FACE: return ElementDiameterSq(grid, aaPos, static_cast<Face*>(elem));
		case VOLUME: return ElementDiameterSq(grid, aaPos, static_cast<Volume*>(elem));
		default: UG_THROW("ElementDiameterSq: Element type not found.")
	}
}

template <typename TElem, typename TAAPos>
number ElementDiameter(Grid& grid,
                       TAAPos& aaPos,
					   TElem* elem)
{
	return std::sqrt(ElementDiameterSq(grid, aaPos, elem));
}

template <typename TAAPos, typename TIterator>
number MaxElementDiameter(Grid& grid, TAAPos& aaPos,
                          TIterator iter, TIterator iterEnd)
{
	number max = 0.0;
	for(; iter != iterEnd; ++iter)
		max = std::max(max, ElementDiameterSq(grid, aaPos, *iter));

#ifdef UG_PARALLEL
	// share value between all procs
	pcl::ProcessCommunicator com;
	max = com.allreduce(max, PCL_RO_MAX);
#endif

	return std::sqrt(max);
}

template <typename TAAPos, typename TIterator>
number MinElementDiameter(Grid& grid, TAAPos& aaPos,
                          TIterator iter, TIterator iterEnd)
{
	number min = std::numeric_limits<number>::max();
	for(; iter != iterEnd; ++iter)
		min = std::min(min, ElementDiameterSq(grid, aaPos, *iter));

#ifdef UG_PARALLEL
	// share value between all procs
	pcl::ProcessCommunicator com;
	min = com.allreduce(min, PCL_RO_MIN);
#endif

	return std::sqrt(min);
}


template <typename TElem1, typename TElem2, typename TAAPos>
typename TAAPos::ValueType
GetDirection (TElem1* e1, TElem2* e2, const TAAPos& aaPos)
{
	using vector_t = typename TAAPos::ValueType;

	vector_t c1 = CalculateCenter (e1, aaPos);
	vector_t c2 = CalculateCenter (e2, aaPos);

	c2 -= c1;
	return c2;
}

template <typename TElem1, typename TElem2, typename TAAPos>
bool CheckDirection (TElem1* e1,
                     TElem2* e2,
                     const TAAPos& aaPos,
                     const typename TAAPos::ValueType& dir,
                     number minAngle,
                     number maxAngle)
{
	using vector_t = typename TAAPos::ValueType;

	const vector_t v = GetDirection (e1, e2, aaPos);
	const number angle = rad_to_deg(VecAngle(v, dir));
	return (minAngle <= angle) && (maxAngle >= angle);
}

}//	end of namespace

#endif
