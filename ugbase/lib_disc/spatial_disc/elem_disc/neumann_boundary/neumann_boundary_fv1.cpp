/*
 * neumann_boundary_fv1.cpp
 *
 *  Created on: 14.10.2010
 *      Author: andreasvogel
 */

#include "neumann_boundary.h"
#include "neumann_boundary_common.h"
#include "neumann_boundary_fv1.h"

namespace ug{

template class NeumannBoundary<Domain1d>;
template class NeumannBoundary<Domain2d>;
template class NeumannBoundary<Domain3d>;

} // namespace ug

