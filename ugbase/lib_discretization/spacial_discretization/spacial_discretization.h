/*
 * domain_discretization.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SPACIAL_DISCRETIZATION__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SPACIAL_DISCRETIZATION__

// element discs
#include "lib_discretization/spacial_discretization/elem_disc/elem_disc.h"

// plug in discs
#include "lib_discretization/spacial_discretization/plug_in_disc/convection_diffusion_equation/convection_diffusion_assemble.h"
#include "lib_discretization/spacial_discretization/plug_in_disc/density_driven_flow/density_driven_flow_assemble.h"

// domain discretization
#include "lib_discretization/spacial_discretization/dirichlet_bnd_values.h"
#include "lib_discretization/spacial_discretization/plug_in_domain_discretization.h"
#include "lib_discretization/spacial_discretization/domain_discretization.h"

////////////////////////
// coupling data
////////////////////////

#include "lib_discretization/spacial_discretization/system_discretization/coupled_system_domain_discretization.h"
#include "lib_discretization/spacial_discretization/disc_coupling/element_data.h"
#include "lib_discretization/spacial_discretization/coupled_plug_in_disc/convection_diffusion_equation/convection_diffusion_assemble.h"
#include "lib_discretization/spacial_discretization/coupled_plug_in_disc/density_driven_flow/density_driven_flow_assemble.h"



#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__SPACIAL_DISCRETIZATION__ */
