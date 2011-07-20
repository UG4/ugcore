/*
 * convection_diffusion.cpp
 *
 *  Created on: 20.07.2011
 *      Author: andreasvogel
 */

#include "convection_diffusion_common.h"
#include "convection_diffusion_fv1.h"
#include "convection_diffusion_fe1.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class ConvectionDiffusionElemDisc<Domain1d>;
template class ConvectionDiffusionElemDisc<Domain2d>;
template class ConvectionDiffusionElemDisc<Domain3d>;

} // end namespace ug
