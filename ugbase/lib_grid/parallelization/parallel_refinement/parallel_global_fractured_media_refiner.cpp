/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "parallel_global_fractured_media_refiner.h"

//ø #include <vector>

#include "lib_grid/parallelization/util/compol_boolmarker.h"

namespace ug {

ParallelGlobalFracturedMediaRefiner::
ParallelGlobalFracturedMediaRefiner(DistributedGridManager& distGridMgr,
									SPRefinementProjector projector) :
	GlobalFracturedMediaRefiner(*distGridMgr.get_assigned_grid(), projector),
	m_distGridMgr(distGridMgr)
{
}

bool
ParallelGlobalFracturedMediaRefiner::
refinement_is_allowed(Vertex* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}


bool ParallelGlobalFracturedMediaRefiner::
refinement_is_allowed(Edge* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}


bool ParallelGlobalFracturedMediaRefiner::
refinement_is_allowed(Face* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}


bool ParallelGlobalFracturedMediaRefiner::
refinement_is_allowed(Volume* elem)
{
	return !m_distGridMgr.is_ghost(elem);
}


void ParallelGlobalFracturedMediaRefiner::
refinement_step_begins()
{
	m_distGridMgr.begin_ordered_element_insertion();
}


void ParallelGlobalFracturedMediaRefiner::
refinement_step_ends()
{
	m_distGridMgr.end_ordered_element_insertion();
}


void ParallelGlobalFracturedMediaRefiner::
communicate_marks(BoolMarker& marker)
{
	GridLayoutMap& layoutMap = m_distGridMgr.grid_layout_map();

//	we have to communicate side marks. In 3d we also have to communicate edge marks.
//	we'll simply communicate edge and face marks. This is no overhead,
//	since no face interfaces don't exist in 2d anyway. The 1d case is ignored.
	ComPol_BoolMarker_AddMarks<EdgeLayout> compolMarkerEDGE(marker);
	ComPol_BoolMarker_AddMarks<FaceLayout> compolMarkerFACE(marker);

//	SLAVE->MASTER
	m_intfComEDGE.exchange_data(layoutMap, InterfaceNodeTypes::INT_H_SLAVE, InterfaceNodeTypes::INT_H_MASTER,
								compolMarkerEDGE);
	m_intfComFACE.exchange_data(layoutMap, InterfaceNodeTypes::INT_H_SLAVE, InterfaceNodeTypes::INT_H_MASTER,
								compolMarkerFACE);

	m_intfComEDGE.communicate();
	m_intfComFACE.communicate();

//	MASTER->SLAVE
	m_intfComEDGE.exchange_data(layoutMap, InterfaceNodeTypes::INT_H_MASTER, InterfaceNodeTypes::INT_H_SLAVE,
								compolMarkerEDGE);
	m_intfComFACE.exchange_data(layoutMap, InterfaceNodeTypes::INT_H_MASTER, InterfaceNodeTypes::INT_H_SLAVE,
								compolMarkerFACE);

	m_intfComEDGE.communicate();
	m_intfComFACE.communicate();
}

}//	end of namespace
