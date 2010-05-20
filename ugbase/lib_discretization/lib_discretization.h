
// include all files from the lib_discretization module

/////////////////////
// basics
/////////////////////

// domain description
#include "lib_discretization/domain.h"

// degree of freedom managers
#include "lib_discretization/dof_manager/dof_manager.h"
#include "lib_discretization/dof_manager/p1conform_dof_manager/p1conform_dof_manager.h"
#include "lib_discretization/dof_manager/general_dof_manager/general_dof_manager.h"

// function spaces
#include "lib_discretization/function_spaces/grid_function_space.h"

// reference elements
#include "lib_discretization/reference_element/reference_elements.h"

// quadratures
#include "lib_discretization/quadrature/quadrature.h"

// local shape functions
#include "lib_discretization/local_shape_function_set/local_shape_function_set_factory.h"

// assembling interface
#include "lib_discretization/assemble.h"

////////////////////////
// function spaces
////////////////////////

#include "lib_discretization/function_spaces/function_spaces.h"

#include "lib_discretization/operator/operator.h"

////////////////////////
// spacial discretizations
////////////////////////

// plug in discs
#include "lib_discretization/domain_discretization/plug_in_disc/convection_diffusion_equation/convection_diffusion_assemble.h"
#include "lib_discretization/domain_discretization/plug_in_disc/density_driven_flow/density_driven_flow_assemble.h"

// domain discretization
#include "lib_discretization/domain_discretization/plug_in_domain_discretization.h"

// coupled system discretization
#include "lib_discretization/domain_discretization/system_discretization/coupled_system_domain_discretization.h"

////////////////////////
// time discretizations
////////////////////////

// time step
#include "lib_discretization/time_discretization/timestep.h"

////////////////////////
// geometric linear solvers
////////////////////////

#include "lib_discretization/multi_grid_solver/mg_solver.h"

////////////////////////
// non-linear solvers
////////////////////////

#include "lib_discretization/non_linear_solver/newton.h"


////////////////////////
// output
////////////////////////

#include "lib_discretization/io/vtkoutput.h"


////////////////////////
// element data
////////////////////////

#include "lib_discretization/domain_discretization/disc_coupling/element_data.h"

