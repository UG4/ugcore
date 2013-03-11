//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d03

#ifndef __H__LIB_GRID__MISC_UTIL__IMPL__
#define __H__LIB_GRID__MISC_UTIL__IMPL__

#include "misc_util.h"
#include "lib_grid/grid/grid_util.h"
#include "vertex_util.h"
#include "edge_util.h"
#include "face_util.h"
#include "volume_util.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	CalculateGeometricObjectCenter
template<class TAAPosVRT>
UG_API
inline
typename TAAPosVRT::ValueType
CalculateGeometricObjectCenter(const GeometricObject* o, TAAPosVRT& aaPosVRT)
{
	switch(o->base_object_id()){
		case VERTEX:	return CalculateCenter(static_cast<const VertexBase*>(o), aaPosVRT);
		case EDGE:		return CalculateCenter(static_cast<const EdgeBase*>(o), aaPosVRT);
		case FACE:		return CalculateCenter(static_cast<const Face*>(o), aaPosVRT);
		case VOLUME:	return CalculateCenter(static_cast<const Volume*>(o), aaPosVRT);
		default:
			UG_THROW("Unknown geometric-object type.");
	}
}

template<class TAAPosVRT, class TAAWeightVRT>
UG_API
inline
typename TAAPosVRT::ValueType
CalculateGeometricObjectCenter(const GeometricObject* o, TAAPosVRT& aaPosVRT,
							   TAAWeightVRT& aaWeight)
{
	switch(o->base_object_id()){
		case VERTEX:
			return CalculateCenter(static_cast<const VertexBase*>(o), aaPosVRT, aaWeight);
		case EDGE:
			return CalculateCenter(static_cast<const EdgeBase*>(o), aaPosVRT, aaWeight);
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
template <class TIterator, class TAAPosVRT>
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
template<class TElem, class TVertexPositionAttachmentAccessor>
TElem* FindClosestByCoordinate(const typename TVertexPositionAttachmentAccessor::ValueType& coord,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						TVertexPositionAttachmentAccessor& aaPosVRT)
{
	if(iterBegin == iterEnd)
		return NULL;

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
template<class vector_t, class TIterator, class TAAPos>
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
template <class TElemPtr1, class TElemPtr2>
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
template <class TElem>
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
template <class TElem>
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
template <class TElem, class TAAPos>
number ElementDiameterSq(Grid& grid,
                         TAAPos& aaPos,
					     TElem* elem)
{
	PointerConstArray<VertexBase*> vVert;
	grid.associated_elements(vVert, elem);

	number max = 0.0;
	for(size_t i = 0; i < vVert.size(); ++i)
		for(size_t j = i+1; j < vVert.size(); ++j)
			max = std::max(max, VecDistanceSq(aaPos[vVert[i]], aaPos[vVert[j]]));

	return max;
}

template <class TAAPos>
number ElementDiameterSq(Grid& grid,
						 TAAPos& aaPos,
						 GeometricObject* elem)
{
	switch(elem->base_object_id()){
		case VERTEX: return ElementDiameterSq(grid, aaPos, static_cast<VertexBase*>(elem));
		case EDGE: return ElementDiameterSq(grid, aaPos, static_cast<EdgeBase*>(elem));
		case FACE: return ElementDiameterSq(grid, aaPos, static_cast<Face*>(elem));
		case VOLUME: return ElementDiameterSq(grid, aaPos, static_cast<Volume*>(elem));
		default: UG_THROW("ElementDiameterSq: Element type not found.")
	}
}

template <class TElem, class TAAPos>
number ElementDiameter(Grid& grid,
                       TAAPos& aaPos,
					   TElem* elem)
{
	return std::sqrt(ElementDiameterSq(grid, aaPos, elem));
}

template <class TAAPos, class TIterator>
number MaxElementDiameter(Grid& grid, TAAPos& aaPos,
                          TIterator iter, TIterator iterEnd)
{
	number max = 0.0;
	for(; iter != iterEnd; ++iter)
		max = std::max(max, ElementDiameterSq(grid, aaPos, *iter));
	return std::sqrt(max);
}

}//	end of namespace

#endif
