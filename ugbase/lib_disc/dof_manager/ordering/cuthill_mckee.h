/*
 * cuthill_mckee.h
 *
 *  Created on: 21.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__CUTHILL_MCKEE__
#define __H__UG__LIB_DISC__DOF_MANAGER__CUTHILL_MCKEE__

#include <vector>
#include "lib_disc/function_spaces/approximation_space.h"
#include "common/cuthill_mckee.h"

namespace ug{

/// orders the dof distribution using Cuthill-McKee
inline void OrderCuthillMcKee(DoFDistribution& dofDistr, bool bReverse);

/// orders the all DofDistributions of the ApproximationSpace using Cuthill-McKee
template <typename TDomain>
void OrderCuthillMcKee(ApproximationSpace<TDomain>& approxSpace, bool bReverse);

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__CUTHILL_MCKEE__ */
