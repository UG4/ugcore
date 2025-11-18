# Copyright (c) 2012:  G-CSC, Goethe University Frankfurt
# Author: Torbjörn Klatt
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# - Run doxygen on source files as a custom target
#
#  include(DoxygenTargets)
#  add_doxygen(<doxyfile>
#    [OUTPUT_DIRECTORY <outputdir>]
#    [INSTALL_DESTINATION <installdir>
#    [INSTALL_COMPONENT <installcomponent>]
#    [INSTALL_PDF_NAME <installpdfname>] ]
#    [DOC_TARGET <targetname>]
#    [WORKING_DIRECTORY] <dirname>]
#    [PROJECT_NUMBER <versionnumber>]
#    [NO_WARNINGS]
#    [NO_PDF]
#    [QUIET])
# Options:
#  OUTPUT_DIRECTORY     name of directory, where the documentation should be 
#                       generated, defaults to 'docs-generated' (relative to
#                       CMAKE_CURRENT_BINARY_DIR)
#  INSTALL_DESTINATION  where the documentation should be coppied afterwars
#  INSTALL_PDF_NAME     name of the final PDF documentation
#  DOC_TARGET           name of the make target to generate the documentation,
#                       defaults to 'doc'
#  WORKING_DIRECTORY    path from which Doxygen should be invoced, defaults to
#                       CMAKE_CURRENT_SOURCE_DIR
#  PROJECT_NUMBER       project version to appear in documentation
#  NO_WARNINGS          don't print warnings while running Doxygen, default TRUE
#  NO_PDF               don't generate PDF version of the documentation
#  QUIET                do not let Doxygen print anything to STDOUT,
#                       defaults to TRUE
#
# Requires these CMake modules:
#  FindDoxygen, FindLATEX, CMakeParseArguments
#
# Requires CMake 2.6 or newer (uses the 'function' command)
#
# Original Author:
#  2009-2010 Ryan Pavlik <rpavlik@iastate.edu> <abiryan@ryand.net>
#  http://academic.cleardefinition.com
#  Iowa State University HCI Graduate Program/VRAC
# 
# Updates by:
#  2012 Torbjörn Klatt <opensource at torbjoern minus klatt dot de>
#  http://torbjoern-klatt.de
#
# Copyright Iowa State University 2009-2010.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)

if (__add_doxygen)
	return()
endif ()
set (__add_doxygen YES)

# We must run the following at "include" time, not at function call time,
# to find the path to this module rather than the path to a calling list file
get_filename_component (_doxygenmoddir ${CMAKE_CURRENT_LIST_FILE} PATH)

if (APPLE)
	list(APPEND CMAKE_PREFIX_PATH "/usr/texbin")
endif ()

if (NOT DOXYGEN_FOUND)
	find_package(Doxygen QUIET)
endif ()

set (DOXYGEN_LATEX "NO")
set (DOXYGEN_PDFLATEX "NO")
set (DOXYGEN_DOT "NO")

if (DOXYGEN_DOT_EXECUTABLE)
	set (DOXYGEN_DOT "YES")
endif ()

find_package (LATEX QUIET)
if (LATEX_COMPILER AND MAKEINDEX_COMPILER)
	set (DOXYGEN_LATEX "YES")
endif ()

if (PDFLATEX_COMPILER)
	set (DOXYGEN_PDFLATEX "YES")
endif ()

# An optional single-file install that supports cmake older than 2.8.0
# For internal use
function (_dt_install_file target filename dest rename)
	if (CMAKE_VER VERSION_LESS 2.8.0)
		set (INSTALL_CODE  "
			if (EXISTS \"${filename}\")
				message (STATUS \"Found: ${filename}\")
				file (INSTALL
					DESTINATION \"\${CMAKE_INSTALL_PREFIX}/${dest}\"
					TYPE FILE
					RENAME \"${rename}\"
					FILES \"${filename}\")
			else ()
				message (STATUS \"Skipping (build '${target}' to create): ${filename}\")
			endif ()
			")
		if (NOT ARGN STREQUAL "")
			set (INSTALL_COMPONENT "${ARGN}")
			set (INSTALL_CODE "
			if (NOT CMAKE_INSTALL_COMPONENT OR \"\${CMAKE_INSTALL_COMPONENT}\" STREQUAL \"${INSTALL_COMPONENT}\")
				${INSTALL_CODE}
			endif ()
			")
		endif ()
		install (CODE "${INSTALL_CODE}")
	else ()
		set (COMPONENT_ARGS)
		if (NOT ARGN STREQUAL "")
			set (COMPONENT_ARGS COMPONENT "${ARGN}")
		endif ()
		install (FILES
			"${filename}"
			DESTINATION
			"${dest}"
			RENAME "${rename}"
			${COMPONENT_ARGS}
			OPTIONAL)
	endif ()

endfunction ()

# An optional single-directory install that supports cmake older than 2.8.0
# For internal use
function (_dt_install_dir target dir dest)
	if (CMAKE_VER VERSION_LESS 2.8.0)
		set (INSTALL_CODE  "
			if (EXISTS \"${dir}\")
				message (STATUS \"Found: ${dir}\")
				file(INSTALL
					DESTINATION \"\${CMAKE_INSTALL_PREFIX}/${dest}\"
					TYPE DIRECTORY
					FILES \"${dir}\")
			else ()
				message (STATUS \"Skipping (build '${target}' to create): ${dir}\")
			endif ()
			")
		if (NOT ARGN STREQUAL "")
			set (INSTALL_COMPONENT "${ARGN}")
			set (INSTALL_CODE "

			if (NOT CMAKE_INSTALL_COMPONENT OR \"\${CMAKE_INSTALL_COMPONENT}\" STREQUAL \"${INSTALL_COMPONENT}\")
				${INSTALL_CODE}
			endif ()
			")
		endif ()
		install (CODE "${INSTALL_CODE}")
	else ()
		set (COMPONENT_ARGS)
		if (NOT ARGN STREQUAL "")
			set (COMPONENT_ARGS COMPONENT "${ARGN}")
		endif ()
		install (DIRECTORY
			"${dir}"
			DESTINATION
			"${dest}"
			${COMPONENT_ARGS}
			OPTIONAL)
	endif ()

endfunction ()

function (add_doxygen _doxyfile)
	# parse arguments
	set (options NO_WARNINGS NO_PDF QUIET)
	set (oneValueArgs OUTPUT_DIRECTORY INSTALL_DESTINATION INSTALL_COMPONENT WORKING_DIRECTORY INSTALL_PDF_NAME DOC_TARGET PROJECT_NUMBER)
	set (multiValueArgs)
	cmake_parse_arguments (_doxy "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

	if (NOT _doxyfile)
		message (FATAL_ERROR "Syntax error in use of add_doxygen!")
	endif ()

	if (NOT _doxy_DOC_TARGET)
		set (_doxy_DOC_TARGET doc)
	endif ()

	if (NOT _doxy_OUTPUT_DIRECTORY)
		set (_doxy_OUTPUT_DIRECTORY "docs-generated")
	endif ()
	file (MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_doxy_OUTPUT_DIRECTORY}")

	if (NOT _doxy_INSTALL_PDF_NAME)
		set (_doxy_INSTALL_PDF_NAME "docs-generated.pdf")
	endif ()

	if (NOT _doxy_WORKING_DIRECTORY)
		set ( _doxy_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	endif ()

	if (NOT _doxy_PROJECT_NUMBER)
		set (_doxy_PROJECT_NUMBER "${CPACK_PACKAGE_VERSION}")
	endif ()
	
	# kind of wierd as we are using the negation of Doxygen's variable here
	if (_doxy_NO_WARNINGS)
		set (_doxy_NO_WARNINGS "NO")
	else ()
		set (_doxy_NO_WARNINGS "YES")
	endif ()
	
	if (_doxy_QUIET)
		set (_doxy_QUIET "YES")
	else ()
		set (_doxy_QUIET "NO")
	endif ()

	if (DOXYGEN_FOUND)
		if (NOT TARGET ${_doxy_DOC_TARGET})

			if (NOT IN_DASHBOARD_SCRIPT)
				add_custom_target (${_doxy_DOC_TARGET})
				set_target_properties (${_doxy_DOC_TARGET}
					PROPERTIES
					EXCLUDE_FROM_ALL
					TRUE)
				set_target_properties (${_doxy_DOC_TARGET}
					PROPERTIES
					EXCLUDE_FROM_DEFAULT_BUILD
					TRUE)
			else ()
				add_custom_target (${_doxy_DOC_TARGET} ALL)
			endif ()

		endif ()

		if (NOT IS_ABSOLUTE "${_doxy_OUTPUT_DIRECTORY}")
			get_filename_component (OUTPUT_DIRECTORY
				"${CMAKE_CURRENT_BINARY_DIR}/${_doxy_OUTPUT_DIRECTORY}"
				ABSOLUTE)
		endif ()

		set_property (DIRECTORY
			APPEND
			PROPERTY
			ADDITIONAL_MAKE_CLEAN_FILES
			"${_doxy_OUTPUT_DIRECTORY}/html"
			"${_doxy_OUTPUT_DIRECTORY}/latex")

		if (NOT TARGET ${_doxy_DOC_TARGET}_open)
			# Create a target to open the generated HTML file.
			if (WIN32)
				set (DOXYGEN_LAUNCHER_COMMAND start "Documentation")
			elseif (NOT APPLE)
				set (DOXYGEN_LAUNCHER_COMMAND xdg-open)
			endif ()
			if (DOXYGEN_LAUNCHER_COMMAND)
				add_custom_target (${_doxy_DOC_TARGET}_open
					COMMAND ${DOXYGEN_LAUNCHER_COMMAND} "${_doxy_OUTPUT_DIRECTORY}/html/index.html")
				set_target_properties (${_doxy_DOC_TARGET}_open
					PROPERTIES
					EXCLUDE_FROM_ALL
					TRUE)
				set_target_properties (${_doxy_DOC_TARGET}_open
					PROPERTIES
					EXCLUDE_FROM_DEFAULT_BUILD
					TRUE)
				add_dependencies (${_doxy_DOC_TARGET}_open ${_doxy_DOC_TARGET})
			endif ()
		endif ()

		get_filename_component (_doxyfileabs "${_doxyfile}" ABSOLUTE)
		get_filename_component (_doxy_INCLUDE_FILE "${_doxyfileabs}" NAME)
		get_filename_component (_doxy_INCLUDE_PATH "${_doxyfileabs}" PATH)

		# Doesn't currently work on Windows, so don't bother
		if (DOXYGEN_LATEX AND NOT _doxy_NO_PDF AND NOT WIN32)
			set (_doxy_MAKE_PDF "YES")
			set (_doxy_GENERATE_LATEX "YES")
		else ()
			set (_doxy_MAKE_PDF NO)
			set (_doxy_GENERATE_LATEX "NO")
		endif ()

		if (DOXYGEN_PDFLATEX AND _doxy_MAKE_PDF)
			set (_doxy_USE_PDFLATEX "YES")
		else ()
			set (_doxy_USE_PDFLATEX "NO")
		endif ()

		if (DOXYGEN_DOT)
			set (_doxy_HAVE_DOT "YES")
			set (_doxy_DOT_PATH ${DOXYGEN_DOT_PATH})
		else ()
			set (_doxy_HAVE_DOT "NO")
			set (_doxy_DOT_PATH)
		endif ()

		# See http://www.cmake.org/pipermail/cmake/2006-August/010786.html
		# for info on this variable
		if ("${CMAKE_BUILD_TOOL}" MATCHES "(msdev|devenv)")
			set (_doxy_WARN_FORMAT "\"$file($line) : $text \"")
		else ()
			set (_doxy_WARN_FORMAT "\"$file:$line: $text \"")
		endif ()

		configure_file ("${_doxygenmoddir}/DoxygenTargets.doxyfile.in"
			"${CMAKE_CURRENT_BINARY_DIR}/${_doxyfile}.additional"
			@ONLY)

		add_custom_command (TARGET
			${_doxy_DOC_TARGET}
			COMMAND
			${DOXYGEN_EXECUTABLE}
			"${CMAKE_CURRENT_BINARY_DIR}/${_doxyfile}.additional"
			WORKING_DIRECTORY
			"${_doxy_WORKING_DIRECTORY}"
			#MAIN_DEPENDENCY ${_doxy_DOC_TARGET}
			COMMENT
			"Running Doxygen with configuration ${_doxyfile}..."
			VERBATIM)

		if (_doxy_MAKE_PDF)
			add_custom_command (TARGET
				${_doxy_DOC_TARGET}
				POST_BUILD
				COMMAND
				${CMAKE_MAKE_PROGRAM}
				WORKING_DIRECTORY
				"${_doxy_OUTPUT_DIRECTORY}/latex"
				COMMENT
				"Generating PDF using PDFLaTeX..."
				VERBATIM)
		endif ()

		if (_doxy_INSTALL_DESTINATION)
			if (_doxy_INSTALL_COMPONENT)
				_dt_install_dir ("${_doxy_DOC_TARGET}" "${_doxy_OUTPUT_DIRECTORY}/html" "${_doxy_INSTALL_DESTINATION}" "${_doxy_INSTALL_COMPONENT}")
				if (_doxy_MAKE_PDF)
					_dt_install_file ("${_doxy_DOC_TARGET}" "${_doxy_OUTPUT_DIRECTORY}/latex/refman.pdf" "${_doxy_INSTALL_DESTINATION}" "${_doxy_INSTALL_PDF_NAME}" "${_doxy_INSTALL_COMPONENT}")
				endif ()

			else ()
				_dt_install_dir ("${_doxy_DOC_TARGET}" "${_doxy_OUTPUT_DIRECTORY}/html" "${_doxy_INSTALL_DESTINATION}")
				if (_doxy_MAKE_PDF)
					_dt_install_file ("${_doxy_DOC_TARGET}" "${_doxy_OUTPUT_DIRECTORY}/latex/refman.pdf" "${_doxy_INSTALL_DESTINATION}" "${_doxy_INSTALL_PDF_NAME}")
				endif ()
			endif ()
		endif ()

	endif ()
endfunction()
