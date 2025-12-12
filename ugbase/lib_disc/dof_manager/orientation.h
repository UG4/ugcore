/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__ORIENTATION__
#define __H__UG__LIB_DISC__DOF_MANAGER__ORIENTATION__

#include <vector>

#include "lib_grid/grid/grid_base_objects.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"

namespace ug {

/// returns the orientation offsets if required
/**
 * If more than on DoF is located on a geometric sub-element for some Local
 * finite element space, the DoF order must be orientated sometimes in order to
 * ensure continuity of the spaces (only spaces with some continuity requirement
 * place dofs on sub elements). This can only be the case when the sub-element
 * is of lower dimension than the LFE-Space.
 *
 * This function computes the orientation offsets if required. I.e. if orientation
 * is needed, the size of the output-vector will be equal to the number of dofs
 * on the subelement and contain the numbers {0, ..., numDoFsOnSub-1} in the
 * order the dofs must be sorted. If no orientation is required, the vector is
 * returned empty
 *
 * @param vOrientOffset		orientation offset (or empty if no orientation required)
 * @param Elem				the element
 * @param SubElem			the sub-element where dofs are to be oriented
 * @param nrSub				the number of the sub-element in reference element counting
 * @param lfeid				the local finite element space
 */
/// \{
void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              GridObject* Elem, GridObject* SubElem, size_t nrSub,
                              const LFEID& lfeid);

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Volume* volume, Face* face, size_t nrFace,
                              const LFEID& lfeid);

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Volume* vol, Edge* edge, size_t nrEdge,
                              const LFEID& lfeid);

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Face* face, Edge* edge, size_t nrEdge,
                              const LFEID& lfeid);
/// \}

void ComputeOrientationOffsetLagrange(std::vector<size_t>& vOrientOffset,
                                      EdgeDescriptor& ed, EdgeVertices* edge, size_t p);

void MapLagrangeMultiIndexQuad(std::vector<size_t>& vOrientOffset,
                               int id0, bool sameOrientation,
                               size_t pOuter);

void MapLagrangeMultiIndexTriangle(std::vector<size_t>& vOrientOffset,
                                   int id0, bool sameOrientation,
                                   size_t pOuter);


} // end namespace ug

#endif