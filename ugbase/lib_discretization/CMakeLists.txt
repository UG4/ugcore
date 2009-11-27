################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_discretization
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

find_path (HYPRELIB_LIB_DIR libHYPRE.a 
		PATHS ENV PATH
		PATH_SUFFIXES src/hypre/lib/ )

if(HYPRELIB_LIB_DIR)

# add include path
include_directories (${HYPRELIB_LIB_DIR}/../include/)
include_directories(${P_UGBASE_SOURCE_DIR})

set(srcDiscretization		differentialoperator.cpp
				discretizationscheme.cpp
				dofhandler.cpp
				elementdata.cpp
				equation.cpp
				numericalsolution.cpp
				spacialdiscretization.cpp
		                subsetdiscretization.cpp
				systemdiscretization.cpp
				rhs.cpp
				trialspace.cpp
				vtkoutput.cpp
			     	newton.cpp)

add_library(discretization 	${srcDiscretization})

add_dependencies(discretization algebra grid common)

endif(HYPRELIB_LIB_DIR)
