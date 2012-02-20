// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.09.2011 (m,d,y)

#ifndef __H__UG__delaunay_triangulation__
#define __H__UG__delaunay_triangulation__

#include <queue>
#include "lib_grid/lg_base.h"

namespace ug
{

///	Transforms the given triangle-set into a delaunay set
/**	THIS METHOD USES Grid::mark
 *
 */

template <class TriIter, class TAAPos>
bool MakeDelaunay(Grid& grid, TriIter trisBegin, TriIter trisEnd, TAAPos& aaPos,
				  CB_ConsiderEdge cbConstrainedEdge = ConsiderNoEdge)
{
	using namespace std;
	typedef typename TAAPos::ValueType vector_t;

//	helper to collect neighbors
	Face* nbrFaces[2];
	vector<EdgeBase*> edges;

//	flip candidates
	queue<EdgeBase*> candidates;
//	all flip candidates are marked in this attachment
	AByte aIsCandidate;
	grid.attach_to_edges_dv(aIsCandidate, 0, false);
	Grid::AttachmentAccessor<EdgeBase, AByte> aaIsCandidate(grid, aIsCandidate);

	grid.attach_to_faces_dv(aIsCandidate, 0, false);
	Grid::AttachmentAccessor<Face, AByte> aaIsMarkedFACE(grid, aIsCandidate);

//	first mark all triangles
	for(TriIter triIter = trisBegin; triIter != trisEnd; ++triIter)
		aaIsMarkedFACE[*triIter] = 1;

//	Collect all candidates for flips (only edges with two neighbors, both marked).
	for(TriIter triIter = trisBegin; triIter != trisEnd; ++triIter){
		Face* t = *triIter;
		CollectEdges(edges, grid, t);
		for(size_t i = 0; i < edges.size(); ++i){
			EdgeBase* e = edges[i];
			if(!(cbConstrainedEdge(e) || aaIsCandidate[e])){
				int numNbrs = GetAssociatedFaces(nbrFaces, grid, e, 2);
			//	two neighbors, both marked
				if(numNbrs == 2 && aaIsMarkedFACE[nbrFaces[0]]
					&& aaIsMarkedFACE[nbrFaces[1]])
				{
				//	the edge is a flip candidate
					candidates.push(e);
					aaIsCandidate[e] = 1;
				}
			}
		}
	}

	grid.end_marking();

	while(!candidates.empty()){
		EdgeBase* e = candidates.front();
		candidates.pop();
		aaIsCandidate[e] = 0;

	//	we only perform swaps on regular manifolds.
		if(GetAssociatedFaces(nbrFaces, grid, e, 2) == 2){
		//	make sure that both neighbors are triangles
			if(nbrFaces[0]->num_vertices() != 3 || nbrFaces[1]->num_vertices() != 3)
				continue;

		//	This section is just temporary...
			VertexBase* conVrt0 = GetConnectedVertex(e, nbrFaces[0]);
			VertexBase* conVrt1 = GetConnectedVertex(e, nbrFaces[1]);

			vector3& v0 = aaPos[e->vertex(0)];
			vector3& v1 = aaPos[e->vertex(1)];
			vector3& v2 = aaPos[conVrt0];
			vector3& v3 = aaPos[conVrt1];

			vector3 cc;
			vector3 testPoint = v3;
			if(!TriangleCircumcenter(cc, v0, v1, v2)){
			//	check the other triangle
				testPoint = v2;
				if(!TriangleCircumcenter(cc, v0, v1, v3)){
					UG_LOG("TriangleCircumcenter failed! Excpect non-delaunay output!\n");
					UG_LOG("  This is most likely caused by two degenerated triangles which "
							"share an edge.\n");
					continue;
				}
			}

			//if(VecDistanceSq(cc, v0) < VecDistanceSq(cc, v3) + SMALL*SMALL){
			if(VecDistanceSq(cc, v0) <= VecDistanceSq(cc, testPoint)){
			//	the edge is fine - it doesn't have to be swapped.
				continue;
			}

		//	ok - everything is fine. Now swap the edge
			e = SwapEdge(grid,  e);

			if(!e){
				UG_LOG("An edge-swap failed. Expect degenerated or flipped triangles "
						"and a non-delaunay output!\n");
				continue;
			}

		//	all edges of associated triangles are candidates again (except e)
			GetAssociatedFaces(nbrFaces, grid, e, 2);
			for(size_t i = 0; i < 2; ++i){
				CollectAssociated(edges, grid, nbrFaces[i]);
				for(size_t j = 0; j < edges.size(); ++j){
					if(edges[j] != e && !(cbConstrainedEdge(edges[j])
										  || aaIsCandidate[edges[j]]))
					{
						candidates.push(edges[j]);
						aaIsCandidate[edges[j]] = 1;
					}
				}
			}
		}
	}

	grid.detach_from_faces(aIsCandidate);
	grid.detach_from_edges(aIsCandidate);

	return true;
}

}//	end of namespace

#endif
