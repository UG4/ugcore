################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_discretization
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

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

include_directories(/Users/andreasvogel/Software/hypre-2.0.0/src/hypre/include/)
link_directories(/Users/andreasvogel/Software/hypre-2.0.0/src/hypre/lib)

add_library(discretization 	${srcDiscretization})


add_dependencies(discretization grid common)
