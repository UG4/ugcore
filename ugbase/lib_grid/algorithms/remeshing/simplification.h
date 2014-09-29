// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_simplification
#define __H__UG_simplification

#include "common/math/ugmath.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/algorithms/polychain_util.h"
#include "lib_grid/grid/grid.h"

namespace ug{

inline bool Valence3VertexIsObsolete(const vector3& p, const vector3& c0,
							  const vector3& c1, const vector3& c2,
							  number maxSquaredHeightToBaseAreaRatio)
{
	vector3 proj, n;
	CalculateTriangleNormalNoNormalize(n, c0, c1, c2);
	if(RayTriangleIntersection(proj, c0, c1, c2, p, n)){
		number hsq = VecDistanceSq(proj, p);
		number area = TriangleArea(c0, c1, c2);

		if(area > 0)
			return (hsq / area) < maxSquaredHeightToBaseAreaRatio;
	}
	return false;
}

/** diagonal is always assumed to be between c0 and c2*/
inline bool Valence4VertexIsObsolete(const vector3& p, const vector3& c0,
									  const vector3& c1, const vector3& c2,
									  const vector3& c3,
									  number maxSquaredHeightToBaseAreaRatio)
{
	vector3 proj, n, n1, n2;
	CalculateTriangleNormalNoNormalize(n1, c0, c1, c2);
	CalculateTriangleNormalNoNormalize(n2, c0, c2, c3);

	VecAdd(n, n1, n2);

	if(RayTriangleIntersection(proj, c0, c1, c2, p, n)){
		number hsq = VecDistanceSq(proj, p);
		number area = TriangleArea(c0, c1, c2) + TriangleArea(c0, c2, c3);

		if(area > 0)
			return (hsq / area) < maxSquaredHeightToBaseAreaRatio;
	}

	if(RayTriangleIntersection(proj, c0, c2, c3, p, n)){
		number hsq = VecDistanceSq(proj, p);
		number area = TriangleArea(c0, c1, c2) + TriangleArea(c0, c2, c3);

		if(area > 0)
			return (hsq / area) < maxSquaredHeightToBaseAreaRatio;
	}

	return false;
}

///	returns true if the normals of the resulting triangles point into the same direction
/**	the quad is always assumed to be split along the diagonal c0, c2*/
inline bool QuadDiagonalIsValid(const vector3& c0, const vector3& c1,
								const vector3& c2, const vector3& c3)
{
	vector3 n0, n1;
	CalculateTriangleNormalNoNormalize(n0, c0, c1, c2);
	CalculateTriangleNormalNoNormalize(n1, c0, c2, c3);
	return VecDot(n0, n1) >= 0;
}


template <class TVrtIter, class TAAPos>
void ReplaceLowValenceVertices(Grid& g, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
							   number maxSquaredHeightToBaseAreaRatio, TAAPos aaPos)
{
	using namespace std;

	Grid::face_traits::secure_container	faces;
	Vertex* conVrts[4];

//	used for valence-4 only
	vector<Vertex*> poly;
	vector<Edge*> conEdges;

//	when we remove a vertex, the valencies of adjacent vertices change.
//	We thus operate on two candidate arrays. One which holds all vertices
//	which are current candidates, and one which records all vertices, which have
//	not been removed. After each run over alll candidates we compare both arrays.
//	If they have the same length, we're done. If not, we'll run again over all
//	new candidates.
	std::vector<Vertex*> candidates;
	std::vector<Vertex*> newCandidates;

	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter)
		newCandidates.push_back(*iter);

	while(newCandidates.size() != candidates.size())
	{
		candidates.swap(newCandidates);
		newCandidates.clear();

		for(size_t icand = 0; icand < candidates.size(); ++icand)
		{
			Vertex* vrt = candidates[icand];

		//	if vrt will be removed, we'll call 'pop_back'
			newCandidates.push_back(vrt);

			g.associated_elements(faces, vrt);

			if(faces.size() < 3 || faces.size() > 4)
				continue;

		//	get the connected vertices. If there are more or less than 3 connected
		//	vertices, then the connected triangles do not form a valid closed surface patch
			size_t numConVrts = 0;
			for(size_t i = 0; i < faces.size(); ++i){
				Face* f = faces[i];
				if(f->num_vertices() != 3){
					numConVrts = 0;
					break;
				}

				int baseInd = GetVertexIndex(f, vrt);

				for(size_t j = 1; j < 3; ++j){
					Vertex* cv = f->vertex((baseInd + j) % 3);
					if(cv != vrt){
						bool gotOne = false;
						for(size_t k = 0; k < numConVrts; ++k){
							if(conVrts[k] == cv){
								gotOne = true;
								break;
							}
						}
						if(!gotOne){
							if(numConVrts < 4)
								conVrts[numConVrts] = cv;
							++numConVrts;
						}
					}
				}
			}

			if(numConVrts != faces.size())
				continue;

			if(numConVrts == 3){
				if(Valence3VertexIsObsolete(aaPos[vrt], aaPos[conVrts[0]], 
											aaPos[conVrts[1]], aaPos[conVrts[2]],
											maxSquaredHeightToBaseAreaRatio))
				{
				//	vrt is currently the last entry in newCandidates and has to be removed
					newCandidates.pop_back();
					g.create<Triangle>(TriangleDescriptor(conVrts[0], conVrts[1], conVrts[2]), faces[0]);
					g.erase(vrt);
				}
			}
			else if(numConVrts == 4){
				conEdges.clear();
				poly.clear();
				for(size_t i = 0; i < faces.size(); ++i){
					std::pair<GridBaseObjectId, int> opObj = faces[i]->get_opposing_object(vrt);
					if(opObj.first == EDGE){
						Edge* e = g.get_edge(faces[i], opObj.second);
						if(e)
							conEdges.push_back(e);
					}
				}
				
				if(conEdges.size() != 4)
					continue;

				CreatePolyChain(poly, g, conEdges.begin(), conEdges.end());

			//	we first try triangulate along the shorter diagonal.
			//	If this doesn't work, we'll try the other one.
			//	diag 0: create diagonal between vrts 0,2
			//	diag 1: between vrts 1,3
				int firstDiag = 0;
				if(VecDistanceSq(aaPos[poly[0]], aaPos[poly[2]])
					> VecDistanceSq(aaPos[poly[1]], aaPos[poly[3]]))
				{
					firstDiag = 1;
				}

			//	we iterate
				for(int idiag = 0; idiag < 2; ++idiag){
					int diag = (firstDiag + idiag) % 2;
					int i0 = diag;
					int i1 = (diag + 1);
					int i2 = (diag + 2);
					int i3 = (diag + 3) % 4;

					if(!QuadDiagonalIsValid(aaPos[poly[i0]], aaPos[poly[i1]],
											aaPos[poly[i2]], aaPos[poly[i3]]))
					{
						continue;
					}

					if(Valence4VertexIsObsolete(aaPos[vrt],
												aaPos[poly[i0]], 
												aaPos[poly[i1]],
												aaPos[poly[i2]],
												aaPos[poly[i3]],
												maxSquaredHeightToBaseAreaRatio))
					{
					//	vrt is currently the last entry in newCandidates and has to be removed
						newCandidates.pop_back();

						// FaceDescriptor fd(vrt, poly[i0], poly[i1]);
						// Face* f1 = g.get_face(fd);
						// if(!f1){
						// 	UG_LOG("num vertices in fd: " << fd.num_vertices() << endl);
						// 	UG_LOG("Couldn't find face with corners: "
						// 			<< aaPos[fd.vertex(0)]
						// 			<< ", " << aaPos[fd.vertex(1)]
						// 			<< ", " << aaPos[fd.vertex(2)] << endl);

						// 	UG_LOG("Existing faces:\n");
						// 	for(FaceIterator iter = g.begin<Face>(); iter != g.end<Face>(); ++iter){
						// 		UG_LOG(ElementDebugInfo(g, *iter) << endl);
						// 	}
						// }

						FaceDescriptor fd0(vrt, poly[i0], poly[i1]);
						FaceDescriptor fd1(vrt, poly[i2], poly[i3]);

						g.create<Triangle>(TriangleDescriptor(poly[i0],
															  poly[i1],
															  poly[i2]),
											g.get_face(fd0));
						g.create<Triangle>(TriangleDescriptor(poly[i0],
															  poly[i2],
															  poly[i3]),
											g.get_face(fd1));
						g.erase(vrt);
						break;
					}
				}
			}
		}
	}
}

}//	end of namespace

#endif	//__H__simplification
