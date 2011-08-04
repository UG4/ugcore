// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.05.2011 (m,d,y)

#include <vector>
#include <queue>
#include "fractured_media_refiner.h"
#include "common/assert.h"

using namespace std;

namespace ug{

template <class TGrid, class TAPosition>
FracturedMediaRefiner<TGrid, TAPosition>::
FracturedMediaRefiner(IRefinementCallback* refCallback) :
	BaseClass(refCallback),
	m_degeneratedEdgeThreshold(SMALL)
{
}

template <class TGrid, class TAPosition>
FracturedMediaRefiner<TGrid, TAPosition>::
FracturedMediaRefiner(TGrid& g, IRefinementCallback* refCallback) :
	BaseClass(g, refCallback),
	m_degeneratedEdgeThreshold(SMALL)
{
}

template <class TGrid, class TAPosition>
FracturedMediaRefiner<TGrid, TAPosition>::
~FracturedMediaRefiner()
{
}

template <class TGrid, class TAPosition>
void
FracturedMediaRefiner<TGrid, TAPosition>::
set_degenerated_edge_threshold(number threshold)
{
	m_degeneratedEdgeThreshold = threshold;
}

template <class TGrid, class TAPosition>
void
FracturedMediaRefiner<TGrid, TAPosition>::
set_position_attachment(TAPosition& aPos)
{
	UG_ASSERT(BaseClass::get_associated_grid(),
			  "The refiner has to be registered at a grid");
	m_aaPos.access(*BaseClass::get_associated_grid(), aPos);
}

template <class TGrid, class TAPosition>
bool
FracturedMediaRefiner<TGrid, TAPosition>::
mark(Face* f, RefinementMark refMark)
{
//	make sure that the position accessor is valid
	UG_ASSERT(m_aaPos.valid(),
			  "Set a position attachment before refining!");

	bool wasMarked = BaseClass::is_marked(f);
	if(!BaseClass::mark(f, refMark))
		return false;

	if(!wasMarked){
		if(IsDegenerated(f, m_aaPos, m_degeneratedEdgeThreshold)){
			m_queDegeneratedFaces.push(f);
		}
	}
	return true;
}

template <class TGrid, class TAPosition>
void
FracturedMediaRefiner<TGrid, TAPosition>::
collect_objects_for_refine()
{
//	if possible, we'll work with the squared threshold.
	number degThresholdSq = m_degeneratedEdgeThreshold * m_degeneratedEdgeThreshold;

//	get the grid on which we'll operate
	UG_ASSERT(BaseClass::get_associated_grid(),
			  "The refiner has to be registered at a grid");

	Grid& grid = *BaseClass::get_associated_grid();

//	make sure that the position accessor is valid
	UG_ASSERT(m_aaPos.valid(),
			  "Set a position attachment before refining!");

//	push all marked degenerated faces to a queue.
//	pop elements from that queue, mark them anisotropic and unmark associated
//	degenerated edges.
//	Furthermore we'll push degenerated faces, which are connected to the current
//	face through a regular edge to the queue (only unprocessed ones).

	Selector& sel = BaseClass::get_refmark_selector();

	vector<EdgeBase*> edges;
	vector<Face*> faces;

//	we need two while-loops. The outer is required to process changes which
//	stem from the base-class implementation.
//todo:	This is a lot of processing due to repeated calls to collect_objects_for_refine.
	do{
		while(!m_queDegeneratedFaces.empty())
		{
			Face* f = m_queDegeneratedFaces.front();
			m_queDegeneratedFaces.pop();

		//	mark as anisotropic
			if(BaseClass::get_mark(f) != RM_ANISOTROPIC)
				BaseClass::mark(f, RM_ANISOTROPIC);

		//	check edges
			CollectAssociated(edges, grid, f);

		//	degenerated neighbors of non-degenerated edges have to be selected.
		//	degenerated edges may not be selected
			size_t numDeg = 0;
			for(size_t i_edge = 0; i_edge< edges.size(); ++i_edge){
				EdgeBase* e = edges[i_edge];
				if(EdgeLengthSq(e, m_aaPos) > degThresholdSq){
				//	non-degenerated edge
				//	make sure it is selected
					if(BaseClass::get_mark(e) != RM_REGULAR)
						mark(e, RM_REGULAR);

				//	this edge possibly connects to an unselected degenerated neighbor.
				//	If this is the case, we'll have to mark it and push it to the queue.
					CollectAssociated(faces, grid, e);
					for(size_t i_face = 0; i_face < faces.size(); ++i_face){
						Face* nbr = faces[i_face];
						if(!sel.is_selected(nbr)){
							if(IsDegenerated(nbr, m_aaPos, m_degeneratedEdgeThreshold)){
							//	push it to the queue.
								m_queDegeneratedFaces.push(nbr);
							}
						}
					}
				}
				else{
				//	degenerated edge. unmark it
					mark(e, RM_NONE);
					++numDeg;
				}
			}

		//	if all edges are degenerate, we will have to perform regular refinement
			if(numDeg == edges.size()){
				BaseClass::mark(f, RM_REGULAR);
				for(size_t i = 0; i < edges.size(); ++i)
					mark(edges[i], RM_REGULAR);
			}
		}

	//	now call the base implementation. If degenerated faces are selected during
	//	that step, then we have to process them too.
		BaseClass::collect_objects_for_refine();
	}while(!m_queDegeneratedFaces.empty());
}


////////////////////////////////////
//	explicit instantiation
template class FracturedMediaRefiner<Grid, APosition>;
template class FracturedMediaRefiner<Grid, APosition2>;
template class FracturedMediaRefiner<Grid, APosition1>;

template class FracturedMediaRefiner<MultiGrid, APosition>;
template class FracturedMediaRefiner<MultiGrid, APosition2>;
template class FracturedMediaRefiner<MultiGrid, APosition1>;

}// end of namespace
