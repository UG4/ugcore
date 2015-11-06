/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__LIB__DISC__
#define __H__UG__LIB_DISC__LIB__DISC__

//////////////////////////////////////////////////////////////////////////////
//
// This file is intended to include all parts of lib_discretization
//
//////////////////////////////////////////////////////////////////////////////

/**
 * \defgroup lib_discretization lib_discretization
 */

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

// DoF Manager
#include "dof_manager/function_pattern.h"

// Function Spaces
#include "function_spaces/approximation_space.h"
#include "function_spaces/grid_function_util.h"
#include "function_spaces/grid_function.h"
#include "function_spaces/interpolate.h"

// IO
#include "io/vtkoutput.h"

// local shape function set
#include "local_finite_element/local_finite_element_provider.h"
#include "local_finite_element/local_finite_element_provider.h"
#include "local_finite_element/local_finite_element_id.h"
#include "local_finite_element/lagrange/lagrange.h"
#include "local_finite_element/lagrange/lagrangep1.h"
#include "local_finite_element/common/lagrange1d.h"
#include "local_finite_element/common/polynomial1d.h"

// Operator
#include "operator/linear_operator/assembled_linear_operator.h"
#include "operator/linear_operator/std_injection.h"
#include "operator/linear_operator/std_transfer.h"
#include "operator/linear_operator/multi_grid_solver/mg_solver.h"

#include "operator/non_linear_operator/assembled_non_linear_operator.h"
#include "operator/non_linear_operator/line_search.h"
#include "operator/non_linear_operator/newton_solver/newton.h"
#include "operator/non_linear_operator/nl_gauss_seidel/nl_gauss_seidel.h"
#include "operator/non_linear_operator/nl_jacobi/nl_jacobi.h"

// Parallelization
#ifdef UG_PARALLEL
#include "parallelization/parallelization_util.h"
#endif

// Quadrature
#include "quadrature/quadrature.h"

// Reference Elements
#include "reference_element/reference_element.h"

// Spatial Discretization
#include "spatial_disc/domain_disc_interface.h"
#include "spatial_disc/domain_disc.h"
#include "spatial_disc/subset_assemble_util.h"

#include "spatial_disc/disc_util/fe_geom.h"
#include "spatial_disc/disc_util/fv1_geom.h"
#include "spatial_disc/disc_util/fvcr_geom.h"
#include "spatial_disc/disc_util/fv_output.h"
#include "spatial_disc/disc_util/fv_util.h"
#include "spatial_disc/disc_util/hfv1_geom.h"
#include "spatial_disc/disc_util/hfvcr_geom.h"

#include "spatial_disc/elem_disc/elem_disc_interface.h"
#include "spatial_disc/elem_disc/elem_disc_assemble_util.h"
#include "spatial_disc/elem_disc/inner_boundary/inner_boundary.h"

#include "spatial_disc/user_data/const_user_data.h"
#include "spatial_disc/user_data/user_data.h"
#include "spatial_disc/user_data/data_export.h"
#include "spatial_disc/user_data/data_import.h"
#include "spatial_disc/user_data/linker/linker.h"

#include "spatial_disc/constraints/constraint_interface.h"
#include "spatial_disc/constraints/continuity_constraints/p1_continuity_constraints.h"
#include "spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

// Time discretization
#include "time_disc/time_disc_interface.h"
#include "time_disc/theta_time_step.h"

#endif /* __H__UG__LIB_DISC__LIB__DISC__ */
