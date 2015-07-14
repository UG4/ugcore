/*
 * orientation.h
 *
 *  Created on: 25.07.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__ORIENTATION__
#define __H__UG__LIB_DISC__DOF_MANAGER__ORIENTATION__

#include <vector>
#include "lib_grid/grid/grid_base_objects.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"

namespace ug{

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
                                      EdgeDescriptor& ed, EdgeVertices* edge, const size_t p);

void MapLagrangeMultiIndexQuad(std::vector<size_t>& vOrientOffset,
							   const int id0, bool sameOrientation,
							   const size_t pOuter);

void MapLagrangeMultiIndexTriangle(std::vector<size_t>& vOrientOffset,
								   const int id0, bool sameOrientation,
								   const size_t pOuter);


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__ORIENTATION__ */
