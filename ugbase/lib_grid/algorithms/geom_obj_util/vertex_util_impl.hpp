//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d05

#ifndef __H__LIB_GRID__VERTEX_UTIL_IMPL__
#define __H__LIB_GRID__VERTEX_UTIL_IMPL__

#include "vertex_util.h"
#include "face_util.h"
#include "../trees/kd_tree_static.h"
#include "misc_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TAAPos>
number VertexDistanceSq(VertexBase* v0, VertexBase* v1, TAAPos& aaPos)
{
	return VecDistanceSq(aaPos[v0], aaPos[v1]);
}

////////////////////////////////////////////////////////////////////////
template <class TAAPos>
number VertexDistance(VertexBase* v0, VertexBase* v1, TAAPos& aaPos)
{
	return VecDistance(aaPos[v0], aaPos[v1]);
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
void CalculateVertexNormal(vector3& nOut, Grid& grid, VertexBase* vrt, TAAPosVRT& aaPos)
{
//	set all normal to zero
	nOut = vector3(0, 0, 0);

//	loop through all associated faces, calculate their normal and add them to thee normal
	Grid::AssociatedFaceIterator iterEnd = grid.associated_faces_end(vrt);
	for(Grid::AssociatedFaceIterator iter = grid.associated_faces_begin(vrt);
		iter != iterEnd; iter++)
	{
		vector3 vN;
		CalculateNormal(vN, *iter, aaPos);
		VecAdd(nOut, nOut, vN);
	}

	VecNormalize(nOut, nOut);
}

////////////////////////////////////////////////////////////////////////
template <class TIterator, class AAPosVRT>
void LaplacianSmooth(Grid& grid, TIterator vrtsBegin,
					TIterator vrtsEnd, AAPosVRT& aaPos,
					number alpha, int numIterations)
{
	for(int iteration = 0; iteration < numIterations; ++iteration){
	//	iterate through all vertices
		for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		//	smooth each one
			VertexBase* vrt = *iter;
			vector3 v(0, 0, 0);
			int num = 0;

			Grid::AssociatedEdgeIterator edgesEnd = grid.associated_edges_end(vrt);
			for(Grid::AssociatedEdgeIterator eIter = grid.associated_edges_begin(vrt);
				eIter != edgesEnd; ++eIter)
			{
				VecAdd(v, v, aaPos[GetConnectedVertex(*eIter, vrt)]);
				++num;
			}

			if(num > 0){
				VecScale(v, v, 1. / (number)num);
				VecSubtract(v, v, aaPos[vrt]);
				VecScale(v, v, alpha);
				VecAdd(aaPos[vrt], aaPos[vrt], v);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
template<class TVertexPositionAttachmentAccessor>
inline
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(VertexBase* v, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	return aaPosVRT[v];
}

template <class TVrtIterator>
VertexBase* MergeMultipleVertices(Grid& grid, TVrtIterator vrtsBegin,
						  	  	  TVrtIterator vrtsEnd)
{
	if(vrtsBegin == vrtsEnd)
		return NULL;

	VertexBase* v = *vrtsBegin;
	++vrtsBegin;
	while(vrtsBegin != vrtsEnd){
		VertexBase* v2 = *vrtsBegin;
		++vrtsBegin;
		MergeVertices(grid, v, v2);
	}
	return v;
}

////////////////////////////////////////////////////////////////////////
//TODO:	replace KDTreeStatic by a dynamic kd-tree.
//TODO: Better support for various iterators.
template <int dim, class TVrtIterator>
void RemoveDoubles(Grid& grid, const TVrtIterator& iterBegin,
					const TVrtIterator& iterEnd, Attachment<MathVector<dim> >& aPos,
					number threshold)
{
	if(!grid.has_vertex_attachment(aPos))
		return;

	typedef Attachment<MathVector<dim> > attachment_type;

	Grid::VertexAttachmentAccessor<attachment_type> aaPos(grid, aPos);

	KDTreeStatic<attachment_type, dim, MathVector<dim> > kdTree;
	kdTree.create_from_grid(grid, iterBegin, iterEnd, aPos, 20, 10, KDSD_LARGEST);

//	we need temporary attachments:
//	a vector<VertexBase*> attachment, that stores for each vertex all other vertices
//	closer than threshold, which have higher attachment data index.
	typedef Attachment<std::list<VertexBase*> >	AVertexList;
	AVertexList aVertexList;
	grid.attach_to_vertices(aVertexList);
	Grid::VertexAttachmentAccessor<AVertexList> aaVL(grid, aVertexList);

//	we'll store in this attachment whether a vertex will be merged or not.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);
	{
		for(TVrtIterator iter = iterBegin; iter != iterEnd; ++iter)
			aaInt[*iter] = 0;
	}

//	compare squares.
	threshold *= threshold;
//	iterate over all vertices and collect all that have aInt == 0 and are within range.
	for(TVrtIterator iter = iterBegin; iter != iterEnd; ++iter)
	{
		VertexBase* v = *iter;
		if(aaInt[v] == 0)
		{//	the vertex will not be removed during merge
		//	find all vertices closer than threshold
			std::list<VertexBase*> neighbours;
			uint numClosest = 3;
			while(numClosest < grid.num_vertices())
			{
				neighbours.clear();
				kdTree.get_neighbourhood(neighbours, aaPos[v], numClosest);

				if(VecDistanceSq(aaPos[neighbours.back()], aaPos[v]) < threshold)
					numClosest *= 2;
				else
					break;
			}

		//	store them in the vertexVec attachment
			if(!neighbours.empty())
			{
				for(std::list<VertexBase*>::iterator nIter = neighbours.begin();
					nIter != neighbours.end(); ++nIter)
				{
					VertexBase* nv = *nIter;
					if(aaInt[nv] == 0)
					{
						if(nv != v)
						{
							if(VecDistanceSq(aaPos[v], aaPos[nv]) < threshold)
							{
								aaVL[v].push_back(nv);
								aaInt[nv] = 1;
							}
							else
								break;
						}
					}
				}
			}
		}
	}

//	iterate over all vertices again and merge collected ones
//	This iteration only works, if the iterators stem from a class
//	like Selector or SubsetHandler or Grid etc, which automatically
//	handle removed elements.
//TODO:	This should be improved somehow!
	{
		TVrtIterator iter = iterBegin;
		while(iter != iterEnd)
		{
			VertexBase* v = *iter;
			if(!aaVL[v].empty())
			{
				std::list<VertexBase*>::iterator nIter = aaVL[v].begin();
				while(nIter != aaVL[v].end())
				{
					VertexBase* delVrt = *nIter;
					nIter++;
					MergeVertices(grid, v, delVrt);
				}
			}

			++iter;
		}
	}

	grid.detach_from_vertices(aVertexList);
	grid.detach_from_vertices(aInt);
}


////////////////////////////////////////////////////////////////////////
template<class TAAPos> inline
void TransformVertex(VertexBase* vrt, matrix33& m, TAAPos& aaPos)
{
//	todo: avoid unnecessary copy
	vector3 oldPos = aaPos[vrt];
	MatVecMult(aaPos[vrt], m, oldPos);
}

////////////////////////////////////////////////////////////////////////
template<class TIterator, class TAAPos>
void TransformVertices(TIterator vrtsBegin, TIterator vrtsEnd,
					   matrix33& m, TAAPos& aaPos)
{
	for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter)
		TransformVertex(*iter, m, aaPos);
}

}//	end of namespace

#endif
