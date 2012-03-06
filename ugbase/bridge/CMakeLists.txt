cmake_minimum_required(VERSION 2.6)

####
# ug_bridge Library
####

project(P_UGBRIDGE)

include("../../ug_cmake_includes.txt")

set(srcUGBridge    	bridge.cpp
					vec_math_bridge.cpp
					lib_grid_bridge.cpp
					test_bridge.cpp
					profiler_bridge.cpp
					util_bridge.cpp
					misc_bridge.cpp
					domain_bridge.cpp
					refinement_bridge)

set(srcUGBridge     ${srcUGBridge}
					pcl_bridge.cpp)

# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(buildAlgebra)
	# add source files that depend on the algebra and discretization system here.
	set(srcUGBridge ${srcUGBridge}
					lib_algebra_bridges/lib_algebra_bridge.cpp
					lib_algebra_bridges/amg_bridge.cpp
					lib_algebra_bridges/eigensolver_bridge.cpp
					
					lib_disc_bridges/adaptive_tools_bridge.cpp
					lib_disc_bridges/common_bridge.cpp
					lib_disc_bridges/algebra_bridge.cpp
					lib_disc_bridges/domain_bridge.cpp
					lib_disc_bridges/elem_discs_bridge.cpp
					lib_disc_bridges/constraints_bridge.cpp
					lib_disc_bridges/multigrid_bridge.cpp
					lib_disc_bridges/output_bridge.cpp
					lib_disc_bridges/user_data_bridge.cpp)
else(buildAlgebra)

endif(buildAlgebra)

set(Sources ${srcUGBridge})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(bridge  ${Sources})
else(BUILD_ONE_LIB)
	add_library(bridge  ${Sources})
endif(BUILD_ONE_LIB)
