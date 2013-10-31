// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Oct 31, 2013

#ifndef __H__UG__manifold_smoothing__
#define __H__UG__manifold_smoothing__

#include "common/types.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"

namespace ug{

///	Smoothes vertices in a 2d-manifold in 3d-space be moving vertices in the
///	tangential plane, only.
/**	USES Grid::mark*/
template <class TVrtIter, class TAAPos3>
void TangentialSmoothSimple(Grid& g, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					  TAAPos3 aaPos, number alpha, size_t numIterations)
{
	Grid::traits<Face>::secure_container	faces;
	for(size_t iteration = 0; iteration < numIterations; ++iteration){
		for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter){
			VertexBase* vrt = *iter;
			g.associated_elements(faces, vrt);
			if(faces.empty())
				continue;

		//	calculate the normal at the current vertex
			vector3 n(0, 0, 0);
			for(size_t i = 0; i < faces.size(); ++i){
				vector3 tn;
				CalculateNormal(tn, faces[i], aaPos);
				VecAdd(n, n, tn);
			}

		//	calculate the center of connected vertices
			vector3 c(0, 0, 0);
			size_t numConnected = 0;
			g.begin_marking();
			for(size_t i_face = 0; i_face < faces.size(); ++i_face){
				Face::ConstVertexArray vrts = faces[i_face]->vertices();
				const size_t numVrts = faces[i_face]->num_vertices();
				for(size_t i_vrt = 0; i_vrt < numVrts; ++i_vrt){
					VertexBase* v = vrts[i_vrt];
					if((v != vrt) && (!g.is_marked(v))){
						g.mark(v);
						VecAdd(c, c, aaPos[v]);
						++numConnected;
					}
				}
			}
			g.end_marking();
			if(numConnected == 0)
				continue;

			VecScale(c, c, 1. / (number)numConnected);

		//	project the center of connected vertices to the plane through vrt
			vector3 cp;
			ProjectPointToPlane(cp, c, aaPos[vrt], n);

		//	move the vertex
			VecScaleAdd(aaPos[vrt], (1. - alpha), aaPos[vrt], alpha, cp);
		}
	}
}

///	Smoothes vertices in a 2d-manifold in 3d-space be moving vertices in the
///	tangential plane, only.
/**	USES Grid::mark
 * Computes an offset vector for each vertex and smoothes this offset by averaging
 * with offsets of adjacent vertices. The offset is then projected back into the
 * tangential plane and each vertex is relocated using those smoothed offsets.*/
template <class TVrtIter, class TAAPos3>
void TangentialSmooth(Grid& g, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					  TAAPos3 aaPos, number alpha, size_t numIterations)
{
//	attach an integer index to all vertices
	AInt aInt;
	g.attach_to_vertices_dv(aInt, -1);
	Grid::VertexAttachmentAccessor<AInt> aaInt(g, aInt);

	int counter = 0;
	std::vector<vector3*> vrtPos;
	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter){
		aaInt[*iter] = counter++;
		vrtPos.push_back(&aaPos[*iter]);
	}

	std::vector<int> adjVrtInds;
	std::vector<int> adjVrtOffsets;
	std::vector<Face*> adjFaces;
	std::vector<int> adjFaceOffsets;

	adjVrtInds.reserve(counter * 7);//just a guess
	adjVrtOffsets.reserve(counter + 1);
	vrtPos.reserve(counter);
	adjFaces.reserve(counter * 7);//just a guess
	adjFaceOffsets.reserve(counter + 1);

//	build adjacency graph
	Grid::traits<Face>::secure_container	faces;
	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter){
		VertexBase* vrt = *iter;
		adjVrtOffsets.push_back(adjVrtInds.size());
		adjFaceOffsets.push_back(adjFaces.size());

		g.associated_elements(faces, vrt);
		if(faces.empty())
			continue;

	//	create adjacencies
		g.begin_marking();
		for(size_t i_face = 0; i_face < faces.size(); ++i_face){
			adjFaces.push_back(faces[i_face]);
			Face::ConstVertexArray vrts = faces[i_face]->vertices();
			const size_t numVrts = faces[i_face]->num_vertices();
			for(size_t i_vrt = 0; i_vrt < numVrts; ++i_vrt){
				VertexBase* v = vrts[i_vrt];
				if((v != vrt) && (!g.is_marked(v))){
					g.mark(v);
					if(aaInt[v] == -1){
					//	positions of vertices which are not to be smoothed but which are
					//	adjacent to smoothed ones, are inserted on the fly.
						aaInt[v] = (int)vrtPos.size();
						vrtPos.push_back(&aaPos[v]);
					}
					adjVrtInds.push_back(aaInt[v]);
				}
			}
		}
		g.end_marking();
	}

	adjVrtOffsets.push_back(adjVrtInds.size());
	adjFaceOffsets.push_back(adjFaces.size());
	g.detach_from_vertices(aInt);

	std::vector<vector3> offsets(vrtPos.size());
	std::vector<vector3> newOffsets(vrtPos.size());

	for(size_t iteration = 0; iteration < numIterations; ++iteration){
	//	calculate the initial offset vectors for each vertex
		for(size_t i_vrt = 0; i_vrt + 1 < adjVrtOffsets.size(); ++i_vrt){
		//	calculate plane normal
			int adjFacesBegin = adjFaceOffsets[i_vrt];
			int adjFacesEnd = adjFaceOffsets[i_vrt+1];

			if(adjFacesBegin == adjFacesEnd){
				offsets[i_vrt] = vector3(0, 0, 0);
				continue;
			}

			vector3 n(0, 0, 0);
			for(int i_adj = adjFacesBegin; i_adj < adjFacesEnd; ++i_adj){
				vector3 tn;
				CalculateNormal(tn, adjFaces[i_adj], aaPos);
				VecAdd(n, n, tn);
			}

		//	calculate center of adjacent vertices
			int adjVrtsBegin = adjVrtOffsets[i_vrt];
			int adjVrtsEnd = adjVrtOffsets[i_vrt+1];

			if(adjVrtsBegin == adjVrtsEnd){
				offsets[i_vrt] = vector3(0, 0, 0);
				continue;
			}

			vector3 c(0, 0, 0);
			for(int i_adj = adjVrtsBegin; i_adj < adjVrtsEnd; ++i_adj){
				int adjVrt = adjVrtInds[i_adj];
				VecAdd(c, c, *vrtPos[adjVrt]);
			}
			VecScale(c, c, 1. / number(adjVrtsEnd - adjVrtsBegin));

		//	project center to plane and calculate the offset
			vector3 cp;
			ProjectPointToPlane(cp, c, *vrtPos[i_vrt], n);
			VecSubtract(offsets[i_vrt], cp, *vrtPos[i_vrt]);
		}

	//	now smooth the offset vectors by repeatedly taking the average of
	//	neighbored vectors
		const size_t numOffsetSmoothSteps = 1;
		for(size_t i_step = 0; i_step < numOffsetSmoothSteps; ++i_step){
			for(size_t i_vrt = 0; i_vrt + 1 < adjVrtOffsets.size(); ++i_vrt){
				int adjVrtsBegin = adjVrtOffsets[i_vrt];
				int adjVrtsEnd = adjVrtOffsets[i_vrt+1];

				if(adjVrtsBegin == adjVrtsEnd)
					continue;

				int numAdj = adjVrtsEnd - adjVrtsBegin;
				//vector3 o = offsets[i_vrt];
				vector3 o;
				VecScale(o, offsets[i_vrt], numAdj);
				for(int i_adj = adjVrtsBegin; i_adj < adjVrtsEnd; ++i_adj){
					int adjVrt = adjVrtInds[i_adj];
					VecAdd(o, o, offsets[adjVrt]);
				}
				VecScale(newOffsets[i_vrt], o, 1. / number(2 * numAdj));

			//	restrict offset to plane
				int adjFacesBegin = adjFaceOffsets[i_vrt];
				int adjFacesEnd = adjFaceOffsets[i_vrt+1];

				if(adjFacesBegin == adjFacesEnd){
					offsets[i_vrt] = vector3(0, 0, 0);
					continue;
				}

				vector3 n(0, 0, 0);
				for(int i_adj = adjFacesBegin; i_adj < adjFacesEnd; ++i_adj){
					vector3 tn;
					CalculateNormal(tn, adjFaces[i_adj], aaPos);
					VecAdd(n, n, tn);
				}

				VecNormalize(n, n);
				number dot = VecDot(n, newOffsets[i_vrt]);
				VecScale(n, n, dot);
				VecSubtract(newOffsets[i_vrt], newOffsets[i_vrt], n);
			}
			offsets.swap(newOffsets);
		}

	//	and now move the vertices by their smoothed offset vectors
		for(size_t i_vrt = 0; i_vrt + 1 < adjVrtOffsets.size(); ++i_vrt){
			VecScaleAdd(*vrtPos[i_vrt], 1, *vrtPos[i_vrt], alpha, offsets[i_vrt]);
		}
	}
}

}// end of namespace

#endif
