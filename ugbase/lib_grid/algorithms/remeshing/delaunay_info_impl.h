/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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
	AAPos aaPos = position_accessor();

	m_grid.begin_marking();
	for(TIter triIter = trisBegin; triIter != trisEnd; ++triIter){
		Face* f = *triIter;
		Face::ConstVertexArray vrts = f->vertices();
		const size_t numVrts = f->num_vertices();
		for(size_t ivrt = 0; ivrt < numVrts; ++ivrt){
			Vertex* vrt = vrts[ivrt];
			if(m_grid.is_marked(vrt))
				continue;
			m_grid.mark(vrt);
			
			if(mark(vrt) != SEGMENT)
				continue;

			vector_t vrtPos = aaPos[vrt];

			bool searching = true;
			m_grid.associated_elements(edges, vrt);
			for(size_t iedge = 0; (iedge < edges.size()) && searching; ++iedge){
				Edge* edge = edges[iedge];
				if(mark(edge) != SEGMENT)
					continue;

				vector_t dir1;
				VecSubtract(dir1, aaPos[GetConnectedVertex(edge, vrt)], vrtPos);

			//	check angles between all segment-edges. If one is found which
			//	is smaller than PI/3, then the vertex will be marked as DART vertex.
				for(size_t iotherEdge = iedge + 1; iotherEdge < edges.size(); ++iotherEdge){
					Edge* otherEdge = edges[iotherEdge];
					if(mark(otherEdge) != SEGMENT)
						continue;

					vector_t dir2;
					VecSubtract(dir2, aaPos[GetConnectedVertex(otherEdge, vrt)], vrtPos);

					if(VecAngle(dir1, dir2) < PI / 3. + SMALL){
						searching = false;
						set_mark(vrt, DART);
						break;
					}
				}
			}
		}
	}
	m_grid.end_marking();
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
