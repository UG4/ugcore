// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_delaunay_info_impl
#define __H__UG_delaunay_info_impl

#include "lib_grid/algorithms/geom_obj_util/edge_util.h"

namespace ug{

template <class TAAPos>	
template <class TIter>
void DelaunayInfo<TAAPos>::
init_marks(TIter trisBegin, TIter trisEnd, bool pushFlipCandidates)
{
	using namespace std;

//	first mark all triangles and associated vertices with mark INNER
	for(TIter iter = trisBegin; iter != trisEnd; ++iter){
		Face* f = *iter;
		if(f->num_vertices() != 3)
			continue;

		set_mark(f, INNER);
		for(size_t i = 0; i < 3; ++i){
			set_mark(f->vertex(i), INNER);
		}
	}

//	Collect all candidates for flips (only edges with two neighbors, both marked).
	Grid::edge_traits::secure_container	edges;
	for(TIter triIter = trisBegin; triIter != trisEnd; ++triIter){
		Face* t = *triIter;
		m_grid.associated_elements(edges, t);
		for(size_t i = 0; i < edges.size(); ++i){
			Edge* e = edges[i];

		//	treat edges only once
			if(mark(e) != NONE)
				continue;

		//	mark inner and outer segments
			if(m_cbConstrainedEdge(e)){
				set_mark(e, SEGMENT);
				set_mark(e->vertex(0), SEGMENT);
				set_mark(e->vertex(1), SEGMENT);
			}
			else{
				Face* nbrFaces[2];
				int numNbrs = GetAssociatedFaces(nbrFaces, m_grid, e, 2);
			//	two neighbors, both marked
				if(numNbrs == 2 && is_inner(nbrFaces[0])
					&& is_inner(nbrFaces[1]))
				{
				//	the edge is a flip candidate
					if(pushFlipCandidates)
						push_candidate(e);
				}
				else{
				//	the edge lies on the rim and has to be marked as a segment
					set_mark(e, SEGMENT);
					set_mark(e->vertex(0), SEGMENT);
					set_mark(e->vertex(1), SEGMENT);
				}
			}
		}
	}

//todo:	iterate over triangles again, this time examining their corners.
//		(use grid::mark to only visit each once).
//		Check for each whether it is connected to two segments with a smaller
//		angle than Pi/3. If this is the case, it will be marked as DART
}


template <class TAAPos>	
template <class TElem>
bool DelaunayInfo<TAAPos>::
is_inner(TElem* e)
{
	Mark m = mark(e);
	return (m == INNER) || (m == NEW_INNER);
}


template <class TAAPos>	
template <class TElem>
bool DelaunayInfo<TAAPos>::
is_segment(TElem* e)
{
	Mark m = mark(e);
	return (m == SEGMENT) || (m == NEW_SEGMENT) || (m == DART) || (m == SHELL);
}


}//	end of namespace

#endif	//__H__UG_delaunay_info_impl
