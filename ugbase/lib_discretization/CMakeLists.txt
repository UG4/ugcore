################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_discretization
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcDiscretization		dofpattern.cpp
				numericalsolution.cpp
				trialspace.cpp
				#vtkoutput.cpp
			     	newton.cpp)
				#differentialoperator.cpp
				#dofhandler.cpp
				#elementdata.cpp
				#equation.cpp
				#spacialdiscretization.cpp
		                #subsetdiscretization.cpp
				#systemdiscretization.cpp

add_library(discretization 	${srcDiscretization})

add_dependencies(discretization algebra grid common)
