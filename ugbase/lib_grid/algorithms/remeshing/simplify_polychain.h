#ifndef __H__UG_simplify_polychain
#define __H__UG_simplify_polychain

#include <algorithm>
#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/smoothing/manifold_smoothing.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "lib_grid/iterators/associated_elements_iterator.h"

namespace ug{

struct IndCmp{
	IndCmp(const std::vector<number>& vals) : m_vals(vals) {}
	bool operator() (size_t i1, size_t i2){
		return m_vals[i1] < m_vals[i2];
	}
	const std::vector<number>& m_vals;
};

template <class TEdgeIter, class TAAPos>
void SimplifyPolylines(Grid& grid, TEdgeIter edgesBegin, TEdgeIter edgesEnd,
					   number curvatureThreshold, TAAPos aaPos)
{
	using namespace std;
	typedef typename TAAPos::ValueType	vector_t;

	grid.begin_marking();

	number curvatureThresholdRad = deg_to_rad(curvatureThreshold);

	vector<Vertex*>	candidates;

	for(TEdgeIter eiter = edgesBegin; eiter != edgesEnd; ++eiter){
		Edge* e = *eiter;
		grid.mark(e);
		for(int i = 0; i < 2; ++i){
			if(!grid.is_marked(e->vertex(i))){
				grid.mark(e->vertex(i));
				candidates.push_back(e->vertex(i));
			}
		}
	}

//	some predeclarations
	vector<Vertex*>	 vrts;
	vrts.reserve(candidates.size());
	vector<number> curvatures;
	curvatures.reserve(candidates.size());
	vector<size_t> indices;
	indices.reserve(candidates.size());
	Grid::edge_traits::callback cb = IsMarked(grid);
	AssocElemIter<Vertex, Edge> assIter(cb);
	
	
	while(!candidates.empty()){
	//	find the subset of candidates which actually could be removed and add them to
	//	the vrts array
		vrts.clear();
		curvatures.clear();
		for(size_t icandidate = 0; icandidate < candidates.size(); ++icandidate){
			Vertex* vrt = candidates[icandidate];
			grid.unmark(vrt);
			Vertex* cvrt[2];
			int numConVrts = 0;
			for(assIter.reinit(grid, vrt); assIter.valid(); ++assIter){
				if(numConVrts < 2)
					cvrt[numConVrts] = GetConnectedVertex(*assIter, vrt);
				++numConVrts;
			}

			if(numConVrts == 2){
				vector_t dir[2];
				VecSubtract(dir[0], aaPos[vrt], aaPos[cvrt[0]]);
				VecSubtract(dir[1], aaPos[cvrt[1]], aaPos[vrt]);
				number angle = VecAngle(dir[0], dir[1]);
				if(angle <= curvatureThresholdRad){
					vrts.push_back(vrt);
					curvatures.push_back(angle);
				}
			}
		}

	//	if we havn't found any vertices to remove, we're done.
		if(vrts.empty())
			break;

	//	set up index-permutation array
		indices.clear();
		for(size_t i = 0; i < vrts.size(); ++i)
			indices.push_back(i);

	//	this class helps in creating a permutation of indices sorted by curvature
		IndCmp cmp(curvatures);

		sort(indices.begin(), indices.end(), cmp);

	//	perform the actual removal.
	//	At this point all vertices in vrts are unmarked.
	//	candidates are collected anew
		candidates.clear();
		for(size_t ivrt = 0; ivrt < vrts.size(); ++ivrt){
			Vertex* vrt = vrts[indices[ivrt]];
		//	marked vertices are candidates for the next iteration
			if(grid.is_marked(vrt))
				continue;

		//	the vertex has to be removed. Create a new edge between
		//	the two connected vertices and mark it. Unmark connected vertices.
			Edge* cedge[2] = {NULL, NULL};
			Vertex* cvrt[2] = {NULL, NULL};
			for(assIter.reinit(grid, vrt); assIter.valid(); ++assIter){
				if(!cvrt[0]){
					cedge[0] = *assIter;
					cvrt[0] = GetConnectedVertex(*assIter, vrt);
				}
				else{
					cedge[1] = *assIter;
					cvrt[1] = GetConnectedVertex(*assIter, vrt);
				}
			}

			if(!cvrt[1]){
				UG_LOG("ALGORITHM ERROR in SimplifyPolylines: two connected vertices "
					   "should always be found for a marked candidate... "
					   "ignoring vertex.\n");
				continue;
			}
			if(grid.get_edge(cvrt[0], cvrt[1]))
				continue;

		//todo: one could check for separate subsets and prohibit removal of vrt in this case.
			Edge* e = *grid.create<RegularEdge>(EdgeDescriptor(cvrt[0], cvrt[1]), cedge[0]);
			grid.mark(e);
			grid.erase(vrt);
			for(size_t i = 0; i < 2; ++i){
				if(!grid.is_marked(cvrt[i])){
					grid.mark(cvrt[i]);
					candidates.push_back(cvrt[i]);
				}
			}
		}
	}

	grid.end_marking();
}


template <class TEdgeIter, class TAAPos>
void SimplifySmoothedPolylines(Grid& grid, TEdgeIter edgesBegin, TEdgeIter edgesEnd,
					   		   number curvatureThreshold,
					   		   TAAPos aaPos,
					   		   number smoothingAlpha,
					   		   int smoothingIterations)
{
	Attachment<typename TAAPos::ValueType> aSmoothPos;
	grid.attach_to_vertices(aSmoothPos);
	TAAPos aaSmoothPos(grid, aSmoothPos);

	for(VertexIterator iter = grid.begin<Vertex>(); iter != grid.end<Vertex>(); ++iter)
		aaSmoothPos[*iter] = aaPos[*iter];

	LaplacianSmooth(grid, grid.begin<Vertex>(), grid.end<Vertex>(), aaSmoothPos,
					smoothingAlpha, smoothingIterations);

	SimplifyPolylines(grid, edgesBegin, edgesEnd, curvatureThreshold, aaSmoothPos);

	grid.detach_from_vertices(aSmoothPos);
}

}//	end of namespace

#endif	//__H__UG_simplify_polychain
