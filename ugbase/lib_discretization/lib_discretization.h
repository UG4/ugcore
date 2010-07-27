
#ifndef __H__LIB_DISCRETIZATION__LIB_DISCRETIZATION__
#define __H__LIB_DISCRETIZATION__LIB_DISCRETIZATION__

/////////////////////
// basics
/////////////////////

// domain description
#include "./domain.h"
#include "./domain_util.h"

// degree of freedom managers
#include "./dof_manager/dof_manager.h"
#include "./dof_manager/p1conform_dof_manager/p1conform_dof_manager.h"
#include "./dof_manager/general_dof_manager/general_dof_manager.h"

// function spaces
#include "./function_spaces/grid_function.h"

#include "./function_spaces/grid_function_space.h"

// reference elements
#include "./reference_element/reference_element.h"

// quadratures
#include "./quadrature/quadrature.h"

// local shape functions
#include "./local_shape_function_set/local_shape_function_set_factory.h"

// assembling interface
#include "./assemble.h"

////////////////////////
// function spaces
////////////////////////

#include "./function_spaces/function_spaces.h"

#include "./operator/operator.h"

////////////////////////
// spacial discretizations
////////////////////////

#include "./spacial_discretization/spacial_discretization.h"

////////////////////////
// time discretizations
////////////////////////

#include "./time_discretization/time_discretization.h"

////////////////////////
// output
////////////////////////

#include "./io/vtkoutput.h"

#endif /* __H__LIB_DISCRETIZATION__LIB_DISCRETIZATION__ */
