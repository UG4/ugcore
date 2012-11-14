MESSAGE(STATUS "Info: Adding make target 'gen_docu' for generating Doxygen documentation.")

# Add the directory with the additional CMake module.
LIST(APPEND CMAKE_MODULE_PATH "../cmake")
# The FindDoxygen CMake module is shipped with CMake since ages; nothing new and shiny.
FIND_PACKAGE(Doxygen QUIET)
# Use the DoxygenTargets CMake module to provide the add_doxygen command
INCLUDE(DoxygenTargets)

# Add the make target while using cmake/DoxygenTargets.doxyfile.in, substitute
# the CMake variables therein and use it as an additional Doxyfile, which overrides
# all Doxygen variables from other Doxyfiles.
ADD_DOXYGEN(
	## the main Doxygen configuration file to be used
	## TODO: maybe make this configurable through a CMAKE argument
	${CMAKE_CURRENT_SOURCE_DIR}/doxy_config_ug4_mathjax.txt
	## where should the documentation be generated to (/html, /latex ...)
	OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ug4
	## the make target name to generate the documentation
	DOC_TARGET gen_docu
	## the directory from where Doxygen should be called (this is docs/)
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	## TODO: add UG's version number to documentation
# 	PROJECT_NUMBER ${UG_VERSION}
	## Uncomment to suppress all Doxygen warnings
	## Leave it commented to write warnings to stderr (WARN_LOGFILE from Doxyfile will be ignored)
# 	NO_WARNINGS
	NO_PDF
	## uncomment QUIET to suppress the verbose Doxygen outout
	## TODO: probably only enable it when DEBUG=ON
	QUIET
)
