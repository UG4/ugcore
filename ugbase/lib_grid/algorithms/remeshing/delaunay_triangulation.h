// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 08.09.2011 (m,d,y)

#ifndef __H__UG__delaunay_triangulation__
#define __H__UG__delaunay_triangulation__

#include <queue>
#include <vector>
#include <sstream>
#include "delaunay_info.h"
#include "common/ug_config.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
template <class TAAPos>
bool MakeDelaunay(DelaunayInfo<TAAPos>& info);

////////////////////////////////////////////////////////////////////////////////
///	Transforms the given triangle-set into a delaunay set
/** Creates a delaunay triangulation. If a minAngle greater than 0 is specified,
 * then additional vertices are introduced, if required to generate fulfill
 * the min-angle-condition.
 */
template <class TAAPos>
bool QualityGridGeneration(Grid& grid, DelaunayInfo<TAAPos>& info,
						   number minAngle = 0,
				  	  	   int maxSteps = -1/*remove this*/);

template <class TriIter, class TAAPos>
bool QualityGridGeneration(Grid& grid, TriIter trisBegin, TriIter trisEnd,
						   TAAPos& aaPos, number minAngle = 0,
				  	  	   Grid::edge_traits::callback cbConstrainedEdge = Grid::edge_traits::cb_consider_none,
				  	  	   int maxSteps = -1/*remove this*/)
{
	using namespace std;

//	helper to collect neighbors
	Face* nbrFaces[2];
	vector<Edge*> edges;

//	set up a delaunay-info structure
	DelaunayInfo<TAAPos> info(grid, aaPos, cbConstrainedEdge);

//	first mark all triangles
//	mark all vertices of those tris
	for(TriIter iter = trisBegin; iter != trisEnd; ++iter){
		Face* f = *iter;
		if(f->num_vertices() != 3)
			continue;

		info.mark(f);
		for(size_t i = 0; i < 3; ++i){
			info.mark(f->vertex(i));
		}
	}

//	Collect all candidates for flips (only edges with two neighbors, both marked).
	for(TriIter triIter = trisBegin; triIter != trisEnd; ++triIter){
		Face* t = *triIter;
		CollectEdges(edges, grid, t);
		for(size_t i = 0; i < edges.size(); ++i){
			Edge* e = edges[i];
		//	unmark associated vertices of constrained edges
		//...
			if(info.is_constrained(e)){
			//	unmark associated vertices
				info.mark(e->vertex(0), false);
				info.mark(e->vertex(1), false);
			}
			else{
				if(!info.is_candidate(e)){
					int numNbrs = GetAssociatedFaces(nbrFaces, grid, e, 2);
				//	two neighbors, both marked
					if(numNbrs == 2 && info.is_marked(nbrFaces[0])
						&& info.is_marked(nbrFaces[1]))
					{
					//	the edge is a flip candidate
						info.push_candidate(e);
					}
					else{
					//	unmark associated vertices
						info.mark(e->vertex(0), false);
						info.mark(e->vertex(1), false);
						info.mark_as_constrained(e);
					}
				}
			}
		}
	}


	return QualityGridGeneration(grid, info, minAngle, maxSteps);	
}

}//	end of namespace

#endif
