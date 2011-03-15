/*
 * lib_discretization.h
 *
 *  Created on: 09.03.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__LIB_DISCRETIZATION__
#define __H__UG__LIB_DISCRETIZATION__LIB_DISCRETIZATION__

//////////////////////////////////////////////////////////////////////////////
//
// This file is intended to include all parts of lib_discretization
//
//////////////////////////////////////////////////////////////////////////////

// Domain
#include "domain.h"
#include "domain_util.h"

// Assemble interface
#include "assemble_interface.h"

// Common
#include "common/function_group.h"
#include "common/geometry_util.h"
#include "common/groups_util.h"
#include "common/local_algebra.h"
#include "common/multi_index.h"
#include "common/subset_group.h"

// DoF Manager
#include "dof_manager/dof_distribution.h"
#include "dof_manager/function_pattern.h"
#include "dof_manager/mg_dof_manager.h"
#include "dof_manager/p1conform/p1conform.h"

// Function Spaces
#include "function_spaces/grid_function_space.h"
#include "function_spaces/grid_function_util.h"
#include "function_spaces/grid_function.h"
#include "function_spaces/interpolate.h"

// IO
#include "io/vtkoutput.h"

// local shape function set
#include "local_shape_function_set/local_dof_pattern.h"
#include "local_shape_function_set/local_shape_function_set.h"
#include "local_shape_function_set/local_shape_function_set_provider.h"
#include "local_shape_function_set/local_shape_function_set_id.h"
#include "local_shape_function_set/LagrangeP1/lagrange.h"
#include "local_shape_function_set/LagrangeP1/lagrangep1.h"
#include "local_shape_function_set/common/lagrange1d.h"
#include "local_shape_function_set/common/polynomial1d.h"

// Operator
#include "operator/linear_operator/assembled_linear_operator.h"
#include "operator/linear_operator/projection_operator.h"
#include "operator/linear_operator/prolongation_operator.h"
#include "operator/linear_operator/multi_grid_solver/mg_solver.h"

#include "operator/non_linear_operator/assembled_non_linear_operator.h"
#include "operator/non_linear_operator/line_search.h"
#include "operator/non_linear_operator/newton_solver/newton.h"

// Parallelization
#ifdef UG_PARALLEL
#include "parallelization/parallel_dof_manager.h"
#include "parallelization/parallel_grid_function.h"
#include "parallelization/parallelization_util.h"
#endif

// Quadrature
#include "quadrature/quadrature.h"

// Reference Elements
#include "reference_element/reference_element.h"

// Spatial Discretization
#include "spatial_discretization/domain_discretization_interface.h"
#include "spatial_discretization/domain_discretization.h"
#include "spatial_discretization/subset_assemble_util.h"

#include "spatial_discretization/disc_helper/finite_element_geometry.h"
#include "spatial_discretization/disc_helper/finite_volume_geometry.h"
#include "spatial_discretization/disc_helper/finite_volume_output.h"
#include "spatial_discretization/disc_helper/finite_volume_util.h"
#include "spatial_discretization/disc_helper/geometry_provider.h"
#include "spatial_discretization/disc_helper/hanging_finite_volume_geometry.h"

#include "spatial_discretization/elem_disc/elem_disc_interface.h"
#include "spatial_discretization/elem_disc/elem_disc_assemble_util.h"
#include "spatial_discretization/elem_disc/convection_diffusion/fv1/convection_diffusion.h"
#include "spatial_discretization/elem_disc/convection_diffusion/fe1/fe1_convection_diffusion.h"
#include "spatial_discretization/elem_disc/density_driven_flow/fv1/density_driven_flow.h"
#include "spatial_discretization/elem_disc/inner_boundary/fv/inner_boundary.h"
#include "spatial_discretization/elem_disc/linear_elasticity/fe1_linear_elasticity.h"
#include "spatial_discretization/elem_disc/navier_stokes/fv/navier_stokes.h"
#include "spatial_discretization/elem_disc/neumann_boundary/fv/neumann_boundary.h"

#include "spatial_discretization/ip_data/user_data_interface.h"
#include "spatial_discretization/ip_data/user_data.h"
#include "spatial_discretization/ip_data/ip_data.h"
#include "spatial_discretization/ip_data/const_user_data.h"
#include "spatial_discretization/ip_data/data_export.h"
#include "spatial_discretization/ip_data/data_import.h"
#include "spatial_discretization/ip_data/data_linker.h"

#include "spatial_discretization/post_process/post_process_interface.h"
#include "spatial_discretization/post_process/post_process_util.h"
#include "spatial_discretization/post_process/constraints/p1_constraints_post_process.h"
#include "spatial_discretization/post_process/dirichlet_boundary/p1_dirichlet_boundary.h"

// Time discretization
#include "time_discretization/time_discretization_interface.h"
#include "time_discretization/theta_time_step.h"

#endif /* __H__UG__LIB_DISCRETIZATION__LIB_DISCRETIZATION__ */
