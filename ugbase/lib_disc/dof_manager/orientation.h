/*
 * orientation.h
 *
 *  Created on: 25.07.2013
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__ORIENTATION__
#define __H__UG__LIB_DISC__DOF_MANAGER__ORIENTATION__

#include <vector>
#include "lib_grid/grid/geometric_base_objects.h"

namespace ug{

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              GeometricObject* Elem, GeometricObject* SubElem, size_t nrSub,
                              const size_t p);

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Volume* volume, Face* face, size_t nrFace,
                              const size_t p);

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Volume* vol, EdgeBase* edge, size_t nrEdge,
                              const size_t p);

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Face* face, EdgeBase* edge, size_t nrEdge,
                              const size_t p);

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__ORIENTATION__ */
