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
/*
template <class TriIter, class TAAPos>
bool MakeDelaunay(Grid& grid, TriIter trisBegin, TriIter trisEnd, TAAPos& aaPos)
{
	using namespace std;

	bool isDelaunay = false;

//	helper to collect neighbors
	Face* nbrFaces[2];
	vector<EdgeBase*> edges;

//	flipCandidates
	queue<EdgeBase*> flipCandidates;

	grid.begin_marking();

//	first mark all triangles
	for(triIter = trisBegin; triIter != trisEnd; ++triIter)
		grid.mark(*triIter);

//	Collect all candidates for flips (only edges with two neighbors, both marked).
	for(triIter = trisBegin; triIter != trisEnd; ++triIter){
		Triangle* t = *triIter;
		CollectEdges(edges, grid, t);
		for(size_t i = 0; i < edges.size(); ++i){
			EdgeBase* e = edges[i];
			int numNbrs = GetAssociatedFaces(nbrFaces, grid, e, 2);
		//	two neighbors, both marked
			if(numNbrs == 2 && grid.is_marked(nbrFaces[0])
				&& grid.is_marked(nbrFaces[1]))
			{
			//	the edge is a flip candidate
				flipCandidates.push(e);
				grid.mark(e);
			}
		}
	}


	grid.end_marking();
	return isDelaunay;
}
*/
}//	end of namespace

#endif
