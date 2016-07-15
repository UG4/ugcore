/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include <vector>
#include <queue>
#include "fractured_media_refiner.h"
#include "common/assert.h"

using namespace std;

namespace ug{

template <class TGrid, class TAPosition>
FracturedMediaRefiner<TGrid, TAPosition>::
FracturedMediaRefiner(SPRefinementProjector projector) :
	BaseClass(projector),
	m_aspectRatioThreshold(SMALL)
{
}

template <class TGrid, class TAPosition>
FracturedMediaRefiner<TGrid, TAPosition>::
FracturedMediaRefiner(TGrid& g, SPRefinementProjector projector) :
	BaseClass(g, projector),
	m_aspectRatioThreshold(SMALL)
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
set_aspect_ratio_threshold(number threshold)
{
	m_aspectRatioThreshold = threshold;
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
		if(aspect_ratio(f) < m_aspectRatioThreshold)
			m_queDegeneratedFaces.push(f);
	}
	return true;
}

template <class TGrid, class TAPosition>
number FracturedMediaRefiner<TGrid, TAPosition>::
aspect_ratio(Face* f)
{
	if(!m_aaPos.valid())
		UG_THROW("A position attachment has to be specified before this method is called.");

	EdgeDescriptor ed;
	f->edge_desc(0, ed);

	number eMin = EdgeLength(&ed, m_aaPos);
	number eMax = eMin;

	for(size_t i = 1; i < f->num_edges(); ++i){
		f->edge_desc(i, ed);
		number len = EdgeLength(&ed, m_aaPos);
		if(len < eMin)
			eMin = len;
		else if(len > eMax)
			eMax = len;
	}

	if(eMax <= 0)
		return 0;

	return eMin / eMax;
}

template <class TGrid, class TAPosition>
void
FracturedMediaRefiner<TGrid, TAPosition>::
collect_objects_for_refine()
{
//	get the grid on which we'll operate
	if(!BaseClass::get_associated_grid())
		UG_THROW("No grid has been set for the refiner.");

	Grid& grid = *BaseClass::get_associated_grid();

//	make sure that the position accessor is valid
	if(!m_aaPos.valid())
		UG_THROW("A position attachment has to be specified before this method is called.");

//	push all marked degenerated faces to a queue.
//	pop elements from that queue, mark them anisotropic and unmark associated
//	degenerated edges.
//	Furthermore we'll push degenerated faces, which are connected to the current
//	face through a regular edge to the queue (only unprocessed ones).

	typename BaseClass::selector_t& sel = BaseClass::get_refmark_selector();

//	some helpers
	vector<Edge*> edges;
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

		//	get the edge with the maximal length
			number eMax = 0;
			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
				number len = EdgeLength(edges[i_edge], m_aaPos);
				if(len > eMax)
					eMax = len;
			}

			if(eMax <= 0)
				eMax = SMALL;

		//	degenerated neighbors of non-degenerated edges have to be selected.
		//	degenerated edges may not be selected
			size_t numDeg = 0;
			for(size_t i_edge = 0; i_edge< edges.size(); ++i_edge){
				Edge* e = edges[i_edge];
				if(EdgeLength(e, m_aaPos) / eMax >= m_aspectRatioThreshold){
				//	non-degenerated edge
				//	make sure it is selected
					if(BaseClass::get_mark(e) != RM_REFINE)
						BaseClass::mark(e, RM_REFINE);

				//	this edge possibly connects to an unselected degenerated neighbor.
				//	If this is the case, we'll have to mark it and push it to the queue.
					CollectAssociated(faces, grid, e);
					for(size_t i_face = 0; i_face < faces.size(); ++i_face){
						Face* nbr = faces[i_face];
						if(!sel.is_selected(nbr)){
							if(aspect_ratio(f) < m_aspectRatioThreshold){
							//	push it to the queue.
								m_queDegeneratedFaces.push(nbr);
							}
						}
					}
				}
				else{
				//	degenerated edge. unmark it
					BaseClass::mark(e, RM_NONE);
					++numDeg;
				}
			}

		//	if all edges are degenerate, we will have to perform regular refinement
			if(numDeg == edges.size()){
				BaseClass::mark(f, RM_REFINE);
				for(size_t i = 0; i < edges.size(); ++i)
					BaseClass::mark(edges[i], RM_REFINE);
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
