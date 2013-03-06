cmake_minimum_required(VERSION 2.6)

####
# ug_bridge Library
####

project(P_UGBRIDGE)

include("../../ug_cmake_includes.txt")

set(srcUGBridge    	bridge.cpp
					domain_bridges/domain_bridge.cpp
					domain_bridges/lib_grid_bridge.cpp
					domain_bridges/load_balancing_bridge.cpp
					domain_bridges/periodic_boundary_bridge.cpp
					domain_bridges/refinement_bridge.cpp
					domain_bridges/selection_bridge.cpp
					domain_bridges/transform_bridge.cpp
					misc_bridges/test_bridge.cpp
					misc_bridges/profiler_bridge.cpp
					misc_bridges/util_bridge.cpp
					misc_bridges/misc_bridge.cpp
					misc_bridges/vec_math_bridge.cpp)

set(srcUGBridge     ${srcUGBridge}
					misc_bridges/pcl_bridge.cpp)

# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(buildAlgebra)
	# add source files that depend on the algebra and discretization system here.
	set(srcUGBridge ${srcUGBridge}
					algebra_bridges/common_bridge.cpp
					algebra_bridges/preconditioner_bridge.cpp
					algebra_bridges/solver_bridge.cpp
					algebra_bridges/eigensolver_bridge.cpp
					algebra_bridges/domain_dependent_preconditioner_bridge.cpp
					
					disc_bridges/domain_disc_bridge.cpp
					disc_bridges/interpolate_bridge.cpp
					disc_bridges/grid_function_bridge.cpp
					disc_bridges/ordering_bridge.cpp
					disc_bridges/integrate_bridge.cpp
					disc_bridges/adaptive_tools_bridge.cpp
					disc_bridges/common_bridge.cpp
					disc_bridges/algebra_bridge.cpp
					disc_bridges/elem_discs_bridge.cpp
					disc_bridges/constraints_bridge.cpp
					disc_bridges/multigrid_bridge.cpp
					disc_bridges/output_bridge.cpp
					disc_bridges/finite_volume_bridge.cpp
					disc_bridges/user_data_bridge.cpp)
else(buildAlgebra)

endif(buildAlgebra)

set(Sources ${srcUGBridge})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/bridge  ${Sources})
else(BUILD_ONE_LIB)
	add_library(bridge  ${Sources})
endif(BUILD_ONE_LIB)
