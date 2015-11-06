/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__FACE_UTIL_IMPL__
#define __H__LIB_GRID__FACE_UTIL_IMPL__

#include "face_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TIterator>
number AreaFaceQuality(TIterator facesBegin, TIterator facesEnd,
					   Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
//	if the area is empty return 0 (bad)
	if(facesBegin == facesEnd)
		return 0;

//	get the first
	number q = FaceQuality(*facesBegin, aaPos);
	++facesBegin;

//	iterate over the others and find a worse one
	for(; facesBegin != facesEnd; ++facesBegin){
		number tq = FaceQuality(*facesBegin, aaPos);
		if(tq < q)
			q = tq;
	}

//	return the quality
	return q;
}

////////////////////////////////////////////////////////////////////////
inline void Triangulate(Grid& grid,
						Grid::VertexAttachmentAccessor<APosition>* paaPos)
{
	Triangulate(grid, grid.begin<Quadrilateral>(),
				grid.end<Quadrilateral>(), paaPos);
}

////////////////////////////////////////////////////////////////////////
//	FixFaceOrientation
template <class TFaceIterator>
void FixFaceOrientation(Grid& grid, TFaceIterator facesBegin,
						TFaceIterator facesEnd)
{
	using namespace std;
//	we use marks to differentiate between processed and unprocessed faces
	grid.begin_marking();

//	we have to mark all faces between facesBegin and facesEnd initially,
//	to differentiate them from the other faces in the grid.
	for(TFaceIterator iter = facesBegin; iter != facesEnd; ++iter)
		grid.mark(*iter);

//	this edge descriptor will be used multiple times
	EdgeDescriptor ed;
	
//	containers to store neighbour elements.
	vector<Face*> vNeighbours;
	
//	stack that stores candidates
	stack<Face*> stkFaces;
	
//	we'll iterate through all faces
	while(facesBegin != facesEnd)
	{
	//	candidates are empty at this point.
	//	if the face is unprocessed it is a new candidate
		if(grid.is_marked(*facesBegin))
		{
		//	mark it as candidate (by removing the mark)
			grid.unmark(*facesBegin);
			stkFaces.push(*facesBegin);
			
		//	while the stack is not empty
			while(!stkFaces.empty())
			{
			//	get the candidate
				Face* f = stkFaces.top();
				stkFaces.pop();

			//	get the neighbours for each side
				for(size_t i = 0; i < f->num_edges(); ++i)
				{
					f->edge_desc(i, ed);
					GetNeighbours(vNeighbours, grid, f, i);

				//	fix orientation of unprocessed neighbours.
					for(size_t j = 0; j < vNeighbours.size(); ++j)
					{
						Face* fn = vNeighbours[j];
						if(grid.is_marked(fn))
						{
						//	check whether the orientation of f and fn differs.
							if(EdgeOrientationMatches(&ed, fn))
							{
							//	the orientation of ed is the same as the orientation
							//	of an edge in fn.
							//	the faces thus have different orientation.
								grid.flip_orientation(fn);
							}

						//	mark the face as processed and add it to the stack
							grid.unmark(fn);
							stkFaces.push(fn);
						}
					}
				}
			}
		}
		
	//	check the next face
		++facesBegin;
	}

	grid.end_marking();
}

////////////////////////////////////////////////////////////////////////
//	InvertOrientation
template <class TFaceIterator>
void InvertOrientation(Grid& grid, TFaceIterator facesBegin,
					   TFaceIterator facesEnd)
{
	for(TFaceIterator iter = facesBegin; iter != facesEnd; ++iter)
		grid.flip_orientation(*iter);
}


////////////////////////////////////////////////////////////////////////
template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const FaceVertices* f, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	uint numVrts = f->num_vertices();
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

	FaceVertices::ConstVertexArray vrts = f->vertices();

//	sum up
	for(uint i = 0; i < numVrts; ++i)
	{
		VecAdd(v, v, aaPosVRT[vrts[i]]);
	}

//	average
	if(numVrts > 0)
		VecScale(v, v, 1./(number)numVrts);

	return v;
}


////////////////////////////////////////////////////////////////////////
template<class TAAPosVRT, class TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const FaceVertices* f, TAAPosVRT& aaPos, TAAWeightVRT& aaWeight)
{
	uint numVrts = f->num_vertices();
	typename TAAPosVRT::ValueType v;
	typedef typename TAAWeightVRT::ValueType weight_t;
//	init v with 0.
	VecSet(v, 0);

	FaceVertices::ConstVertexArray vrts = f->vertices();

//	sum up
	weight_t totalWeight = 0;
	for(uint i = 0; i < numVrts; ++i)
	{
		weight_t w = aaWeight[vrts[i]];
		VecScaleAppend(v, w, aaPos[vrts[i]]);
		totalWeight += w;
	}

//	average
	if(totalWeight != 0)
		VecScale(v, v, 1./(number)totalWeight);

	return v;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t, class TAAPos>
bool
ContainsPoint(const FaceVertices* f, const vector_t& p, TAAPos aaPos)
{
	switch(f->num_vertices()){
		case 3: return PointIsInsideTriangle(p, aaPos[f->vertex(0)],
											 aaPos[f->vertex(1)],
											 aaPos[f->vertex(2)]);
		case 4: return PointIsInsideQuadrilateral(p, aaPos[f->vertex(0)],
												  aaPos[f->vertex(1)],
												  aaPos[f->vertex(2)],
												  aaPos[f->vertex(3)]);
		default:
			UG_THROW("Unknown face type with " << f->num_vertices()
					<< " vertices encountered in ContainsPoint(...).");
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
//	project points to surface 
template <class TTriangleIterator, class TAAPosVRT>
bool ProjectPointToSurface(vector3& vOut, const vector3& v, const vector3& n,
						   TTriangleIterator trisBegin, TTriangleIterator trisEnd,
						   TAAPosVRT& aaPos, bool compareNormals)
{
	bool gotOne = false;
	number bestDist = 0;	// value doesn't matter - will be overwritten later on.

//	iterate through all triangles and find the closest intersection
	for(TTriangleIterator iter = trisBegin; iter != trisEnd; ++iter)
	{
		Triangle* tri = *iter;
		
		const vector3& p1 = aaPos[tri->vertex(0)];
		const vector3& p2 = aaPos[tri->vertex(1)];
		const vector3& p3 = aaPos[tri->vertex(2)];
		
		vector3 tn;
		CalculateTriangleNormalNoNormalize(tn, p1, p2, p3);
		
	//	if normal-check is enabled, we have to make sure, that the points
	//	normal points into the same direction as the triangles normal.
		if(compareNormals){
			if(VecDot(tn, n) <= 0)
				continue;
		}
		
		number bc1, bc2;
		vector3 vTmp;
		number distance = DistancePointToTriangle(vTmp, bc1, bc2,
													v, p1, p2, p3, tn);
	
		if(gotOne){
			if(distance < bestDist){
				bestDist = distance;
				vOut = vTmp;
			}
		}
		else{
			gotOne = true;
			vOut = vTmp;
			bestDist = distance;
		}
	}

	return gotOne;
}

template <class TAAPosVRT>
int PointFaceTest(vector3& v, Face* f, TAAPosVRT& aaPos)
{
	vector3 n;
	vector3 dir;
	
	CalculateNormal(n, f, aaPos);
	VecSubtract(dir, v, aaPos[f->vertex(0)]);
	
	number d = VecDot(dir, n);
	if(d > 0)
		return 1;
	if(d < 0)
		return -1;
	return 0;
}

template <class TAAPosVRT>
bool IsDegenerated(Face* f, TAAPosVRT& aaPos, number threshold)
{
	number threshSQ = threshold * threshold;
	size_t numVrts = f->num_vertices();

	for(size_t i = 0; i < numVrts; ++i){
		if(VecDistanceSq(aaPos[f->vertex(i)], aaPos[f->vertex((i+1) % numVrts)])
			< threshSQ)
		{
			return true;
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
number FaceArea(FaceVertices* f, TAAPosVRT& aaPos)
{
	number area = 0;
	for(size_t i = 2; i < f->num_vertices(); ++i){
		area += TriangleArea(aaPos[f->vertex(0)],
							 aaPos[f->vertex(i - 1)],
							 aaPos[f->vertex(i)]);
	}
	return area;
}

////////////////////////////////////////////////////////////////////////
template <class TIterator, class TAAPosVRT>
number FaceArea(TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos)
{
	number sum = 0.;

	for (; facesBegin != facesEnd; ++facesBegin)
		sum += FaceArea(*facesBegin, aaPos);

	return sum;
}

////////////////////////////////////////////////////////////////////////
template <class TIterator, class TAAPosVRT>
Face* FindSmallestFace(TIterator facesBegin, TIterator facesEnd, TAAPosVRT& aaPos)
{
	//	if facesBegin equals facesEnd, then the list is empty and we can
	//	immediately return NULL
		if(facesBegin == facesEnd)
			return NULL;

	//	the first face is the first candidate for the smallest face.
		Face* smallestFace = *facesBegin;
		number smallestArea = FaceArea(smallestFace, aaPos);
		++facesBegin;

		for(; facesBegin != facesEnd; ++facesBegin){
			Face* curFace = *facesBegin;
			number curArea = FaceArea(curFace, aaPos);
			if(curArea < smallestArea){
				smallestFace = curFace;
				smallestArea = curArea;
			}
		}

		return smallestFace;
}


}//	end of namespace

#endif
