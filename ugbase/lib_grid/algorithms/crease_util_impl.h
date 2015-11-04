#ifndef __H__UG_crease_util_impl
#define __H__UG_crease_util_impl

#include <algorithm>
#include "lib_grid/iterators/associated_elements_iterator.h"

namespace ug{

template <class TVrtIter, class TAAPos>
void SelectKinkVertices(Selector& sel, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
						number thresholdAngle, bool selectDarts, TAAPos aaPos,
						Grid::edge_traits::callback cbConsiderEdge)
{
	using std::max;
	using std::endl;
	typedef typename TAAPos::ValueType vector_t;

	if(!sel.grid())
		return;

	Grid& grid = *sel.grid();

	number thresholdRad = deg_to_rad(thresholdAngle);
	AssocElemIter<Vertex, Edge> eiter(cbConsiderEdge);

	for(TVrtIter viter = vrtsBegin; viter != vrtsEnd; ++viter){
		Vertex* v = *viter;
		Vertex* cv[2];
		int numCons = 0;
		for(eiter.reinit(grid, v); eiter.valid(); ++eiter){
			Edge* e = *eiter;
			if(numCons < 2)
				cv[numCons] = GetConnectedVertex(e, v);
			++numCons;
		}
		if(numCons == 2){
			vector_t dir0, dir1;
			VecSubtract(dir0, aaPos[v], aaPos[cv[0]]);
			VecSubtract(dir1, aaPos[cv[1]], aaPos[v]);
			if(VecAngle(dir0, dir1) >= thresholdRad)
				sel.select(v);
		}
		else if(selectDarts){
			sel.select(v);
		}
	}
}

}//	end of namespace

#endif	//__H__UG_crease_util_impl
