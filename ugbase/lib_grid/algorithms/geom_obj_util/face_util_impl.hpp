//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d05

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
//	FixOrientation
template <class TFaceIterator>
void FixOrientation(Grid& grid, TFaceIterator facesBegin,
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
				for(uint i = 0; i < f->num_edges(); ++i)
				{
					GetNeighbours(vNeighbours, grid, f, i);

				//	fix orientation of unprocessed neighbours.
					for(size_t j = 0; j < vNeighbours.size(); ++j)
					{
						Face* fn = vNeighbours[j];
						if(grid.is_marked(fn))
						{
						//	check whether the orientation of f and fn differs.
							f->edge(i, ed);
							if(EdgeOrientationMatches(ed, fn))
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
//	CalculateFaceCenter
template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateFaceCenter(Face* f, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	uint numVrts = f->num_vertices();
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

//	sum up
	for(uint i = 0; i < numVrts; ++i)
	{
		VecAdd(v, v, aaPosVRT[f->vertex(i)]);
	}

//	average
	if(numVrts > 0)
		VecScale(v, v, 1./(number)numVrts);

	return v;
}

////////////////////////////////////////////////////////////////////////
template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(Face* f, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	uint numVrts = f->num_vertices();
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

//	sum up
	for(uint i = 0; i < numVrts; ++i)
	{
		VecAdd(v, v, aaPosVRT[f->vertex(i)]);
	}

//	average
	if(numVrts > 0)
		VecScale(v, v, 1./(number)numVrts);

	return v;
}


////////////////////////////////////////////////////////////////////////
//	FindFaceByCoordinate
template<class TVertexPositionAttachmentAccessor>
Face* FindFaceByCoordinate(const typename TVertexPositionAttachmentAccessor::ValueType& coord,
							FaceIterator iterBegin, FaceIterator iterEnd,
							TVertexPositionAttachmentAccessor& aaPosVRT)
{
	if(iterBegin == iterEnd)
		return NULL;

	FaceIterator iter = iterBegin;
	Face* bestFace = *iter;
	number bestDistSq = VecDistanceSq(coord, CalculateFaceCenter(bestFace, aaPosVRT));
	iter++;

	while(iter != iterEnd)
	{
		number distSq = VecDistanceSq(coord, CalculateFaceCenter(*iter, aaPosVRT));
		if(distSq < bestDistSq)
		{
			bestDistSq = distSq;
			bestFace = *iter;
		}

		++iter;
	}

	return bestFace;
}

////////////////////////////////////////////////////////////////////////
//	project points to surface 
template <class TTriangleIterator, class TAAPosVRT>
bool ProjectPointToSurface(vector3& vOut, const vector3& v, const vector3& n,
						   TTriangleIterator trisBegin, TTriangleIterator trisEnd,
						   TAAPosVRT& aaPos, bool compareNormals)
{
	vector3 vInter;
	bool gotOne = false;
	number b1, b2, t;
	number tBest = 0;	// value doesn't matter - will be overwritten later on.

//	iterate through all triangles and find the closest intersection
	for(TTriangleIterator iter = trisBegin; iter != trisEnd; ++iter)
	{
		Triangle* tri = *iter;
		if(RayTriangleIntersection(vInter, b1, b2, t, aaPos[tri->vertex(0)],
								aaPos[tri->vertex(1)], aaPos[tri->vertex(2)], v, n))
		{
		//	check the normal
			vector3 tn = n;
			if(compareNormals){
				CalculateTriangleNormal(tn, aaPos[tri->vertex(0)], aaPos[tri->vertex(1)],
										aaPos[tri->vertex(2)]);
			}
		//	the triangle normal and the point - normal should match
			if(VecDot(tn, n) > 0){
				if(gotOne){
					if(fabs(t) < tBest){
						vOut = vInter;
						tBest = fabs(t);
					}
				}
				else{
					gotOne = true;
					vOut = vInter;
					tBest = fabs(t);
				}
			}
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

}//	end of namespace

#endif
