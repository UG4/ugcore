//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d05

#ifndef __H__LIB_GRID__VERTEX_UTIL_IMPL__
#define __H__LIB_GRID__VERTEX_UTIL_IMPL__

#include "vertex_util.h"
#include "edge_util.h"
#include "face_util.h"
#include "../trees/kd_tree_static.h"
#include "misc_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TAAPos>
number VertexDistanceSq(Vertex* v0, Vertex* v1, TAAPos& aaPos)
{
	return VecDistanceSq(aaPos[v0], aaPos[v1]);
}

////////////////////////////////////////////////////////////////////////
template <class TAAPos>
number VertexDistance(Vertex* v0, Vertex* v1, TAAPos& aaPos)
{
	return VecDistance(aaPos[v0], aaPos[v1]);
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
void CalculateVertexNormal(vector3& nOut, Grid& grid, Vertex* vrt, TAAPosVRT& aaPos)
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
template <class TAAPosVRT>
void CalculateBoundaryVertexNormal2D(typename TAAPosVRT::ValueType& nOut,
									 Grid& grid, Vertex* vrt,
						   	   	     TAAPosVRT& aaPos)
{
//	The algorithm is a little cumbersome. However, through this setup, we
//	make sure that the orientation of the normal indeed points outwards,
//	based only on the topology.

//	set nOut to 0
	VecSet(nOut, 0);

//	iterate over associated faces
	std::vector<Face*> faces;
	CollectAssociated(faces, grid, vrt);

	EdgeDescriptor ed;
	for(size_t i_face = 0; i_face < faces.size(); ++i_face){
		Face* f = faces[i_face];
	//	check for each side of f whether it is a boundary edge
		for(size_t i_side = 0; i_side < f->num_sides(); ++i_side){
			if(IsBoundaryEdge2D(grid, grid.get_edge(f, i_side))){
				f->edge_desc(i_side, ed);
			//	make sure that e contains the specified vertex
				if(!EdgeContains(&ed, vrt))
					continue;

			//	the normal pointing outwards is clearly defined from the
			//	orientation of the edge descriptor
				nOut.x() += aaPos[ed.vertex(1)].y() - aaPos[ed.vertex(0)].y();
				nOut.y() += aaPos[ed.vertex(0)].x() - aaPos[ed.vertex(1)].x();
			}
		}
	}

	VecNormalize(nOut, nOut);
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
void CalculateBoundaryVertexNormal3D(vector3& nOut, Grid& grid, Vertex* vrt,
						   	   	     TAAPosVRT& aaPos)
{
//	The algorithm is a little cumbersome. However, through this setup, we
//	make sure that the orientation of the normal indeed points outwards,
//	based only on the topology.

//	set nOut to 0
	VecSet(nOut, 0);

//	iterate over associated volumes
	std::vector<Volume*> vols;
	CollectAssociated(vols, grid, vrt);

	FaceDescriptor fd;
	for(size_t i_vol = 0; i_vol < vols.size(); ++i_vol){
		Volume* v = vols[i_vol];
	//	check for each side of f whether it is a boundary edge
		for(size_t i_side = 0; i_side < v->num_sides(); ++i_side){
			if(IsBoundaryFace3D(grid, grid.get_face(v, i_side))){
				v->face_desc(i_side, fd);

			//	make sure that fd contains the given vertex
				if(!FaceContains(&fd, vrt))
					continue;

			//	the normal pointing outwards is clearly defined from the
			//	orientation of the face descriptor
				vector3 n;
				CalculateNormal(n, &fd, aaPos);
				VecAdd(nOut, nOut, n);
			}
		}
	}

	VecNormalize(nOut, nOut);
}

////////////////////////////////////////////////////////////////////////
template <class TVrtIter, class TAPosition>
void
CalculateBoundingBox(typename TAPosition::ValueType& vMinOut,
					 typename TAPosition::ValueType& vMaxOut,
					 TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					 Grid::AttachmentAccessor<Vertex, TAPosition>& aaPos)
{
	size_t dim = TAPosition::ValueType::Size;

    if(vrtsBegin != vrtsEnd)
    {
		vMinOut = aaPos[*vrtsBegin];
		vMaxOut = vMinOut;

    	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter)
    	{
    		for(size_t i = 0; i < dim; ++i){
				vMinOut[i] = std::min(vMinOut[i], aaPos[*iter][i]);
				vMaxOut[i] = std::max(vMaxOut[i], aaPos[*iter][i]);
    		}
    	}
    }
}

////////////////////////////////////////////////////////////////////////
template<class TVertexPositionAttachmentAccessor>
inline
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const Vertex* v, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	return aaPosVRT[v];
}

template<class TAAPosVRT, class TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const Vertex* v, TAAPosVRT& aaPosVRT, TAAWeightVRT&)
{
	return aaPosVRT[v];
}

////////////////////////////////////////////////////////////////////////
template <class TVrtIter, class TAPosition>
typename TAPosition::ValueType
CalculateCenter(TVrtIter vrtsBegin, TVrtIter vrtsEnd,
				Grid::AttachmentAccessor<Vertex, TAPosition>& aaPos)
{
	typename TAPosition::ValueType vMin, vMax;
	CalculateBoundingBox(vMin, vMax, vrtsBegin, vrtsEnd, aaPos);
	typename TAPosition::ValueType vRet;
	VecScaleAdd(vRet, 0.5, vMin, 0.5, vMax);
	return vRet;
}

////////////////////////////////////////////////////////////////////////
template <class TVrtIter, class TAPosition>
typename TAPosition::ValueType
CalculateBarycenter(TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					Grid::VertexAttachmentAccessor<TAPosition>& aaPos)
{
	typename TAPosition::ValueType v;
	VecSet(v, 0);
	int counter = 0;
	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter)
	{
		VecAdd(v,v,aaPos[*iter]);
		counter++;
	}

	if(counter>0)
		VecScale(v,v,1.f/counter);
	return v;
}

////////////////////////////////////////////////////////////////////////
template <class TIterator, class AAPosVRT>
void LaplacianSmooth(Grid& grid, TIterator vrtsBegin,
					TIterator vrtsEnd, AAPosVRT& aaPos,
					number alpha, int numIterations)
{
	std::vector<Face*> faces;
	std::vector<Volume*> vols;

	for(int iteration = 0; iteration < numIterations; ++iteration){
	//	iterate through all vertices
		for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		//	smooth each one
			Vertex* vrt = *iter;
			vector3 v;
			VecSet(v, 0);
			int num = 0;

		//	calculate smoothing vector relative to neighbors
			CollectAssociated(faces, grid, vrt);
			for(size_t i = 0; i < faces.size(); ++i)
			{
				VecAdd(v, v, CalculateGridObjectCenter(grid.get_opposing_object(vrt, faces[i]), aaPos));
				++num;
			}

			if(grid.num<Volume>() > 0){
				CollectAssociated(vols, grid, vrt);
				for(size_t i = 0; i < vols.size(); ++i)
				{
					VecAdd(v, v, CalculateGridObjectCenter(grid.get_opposing_object(vrt, vols[i]), aaPos));
					++num;
				}
			}

			if(num > 0){
				VecScale(v, v, 1. / (number)num);
				VecSubtract(v, v, aaPos[vrt]);
				VecScale(v, v, alpha);
				VecAdd(aaPos[vrt], aaPos[vrt], v);
			}
			/*
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
			}*/
		}
	}
}

////////////////////////////////////////////////////////////////////////
template <class TVrtIterator>
Vertex* MergeMultipleVertices(Grid& grid, TVrtIterator vrtsBegin,
						  	  	  TVrtIterator vrtsEnd)
{
	if(vrtsBegin == vrtsEnd)
		return NULL;

	Vertex* v = *vrtsBegin;
	++vrtsBegin;
	while(vrtsBegin != vrtsEnd){
		Vertex* v2 = *vrtsBegin;
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
//	a vector<Vertex*> attachment, that stores for each vertex all other vertices
//	closer than threshold, which have higher attachment data index.
	typedef Attachment<std::list<Vertex*> >	AVertexList;
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
		Vertex* v = *iter;
		if(aaInt[v] == 0)
		{//	the vertex will not be removed during merge
		//	find all vertices closer than threshold
			std::list<Vertex*> neighbours;
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
				for(std::list<Vertex*>::iterator nIter = neighbours.begin();
					nIter != neighbours.end(); ++nIter)
				{
					Vertex* nv = *nIter;
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
			Vertex* v = *iter;
			if(!aaVL[v].empty())
			{
				std::list<Vertex*>::iterator nIter = aaVL[v].begin();
				while(nIter != aaVL[v].end())
				{
					Vertex* delVrt = *nIter;
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
void TransformVertex(Vertex* vrt, matrix33& m, TAAPos& aaPos)
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

////////////////////////////////////////////////////////////////////////
template<class TIterator, class TAAPos> inline
void MoveVertices(TIterator vrtsBegin, TIterator vrtsEnd, TAAPos aaPos,
				  const typename TAAPos::ValueType& offset)
{
	for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter)
		aaPos[*iter] += offset;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t, class TAAPos>
UG_API bool
ContainsPoint(const Vertex* v, const vector_t& p, TAAPos aaPos)
{
	const vector_t& pv = aaPos[v];
	for(size_t i = 0; i < vector_t::Size; ++i){
		if(pv[i] != p[i])
			return false;
	}
	return true;
}

}//	end of namespace

#endif
