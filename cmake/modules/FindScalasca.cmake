#
# - Find Scalasca
#
# This module finds the Scalasca include directory and library
#
# It sets the following variables:
#  SCALASCA_FOUND       - Set to false, or undefined, if Scalasca isn't found.
#  SCALASCA_INCLUDE_DIR - The Scalasca include directory.
#  SCALASCA_COMMAND     - The full path command for 'scalasca'
#  SCALASCA_KCONFIG     - The full path command for 'kconfig'
#  SCALASCA_USER_CFLAGS - C/C++ flags needed to compile using "scalasca -instrument -user"
	
# 1.a) It may be the case the path hints are given in environment, so check for
#      those hints first

# If SCALASCA_ROOT is defined in the environment, use it.
if ( NOT SCALASCA_ROOT AND NOT $ENV{SCALASCA_ROOT} STREQUAL "" )
  set( SCALASCA_ROOT $ENV{SCALASCA_ROOT} )
endif ( NOT SCALASCA_ROOT AND NOT $ENV{SCALASCA_ROOT} STREQUAL "" )

# If SCALASCA_INCLUDEDIR is defined in the environment, use it.
if ( NOT SCALASCA_INCLUDEDIR AND NOT $ENV{SCALASCA_INCLUDEDIR} STREQUAL "" )
  set( SCALASCA_INCLUDEDIR $ENV{SCALASCA_INCLUDEDIR} )
endif ( NOT SCALASCA_INCLUDEDIR AND NOT $ENV{SCALASCA_INCLUDEDIR} STREQUAL "" )

if( SCALASCA_ROOT )
    set(SCALASCA_ENVVAR_INCLUDE_SEARCH_DIRS ${SCALASCA_ROOT}/include )
endif( SCALASCA_ROOT )

if( SCALASCA_INCLUDEDIR )
  file(TO_CMAKE_PATH ${SCALASCA_INCLUDEDIR} SCALASCA_INCLUDEDIR)
  SET(SCALASCA_ENVVAR_INCLUDE_SEARCH_DIRS ${SCALASCA_INCLUDEDIR} )
endif( SCALASCA_INCLUDEDIR )

# 1.b) Lets try also invoking PkgConfig to get hints

find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    # this will set 
    # PKG_CONFIG_LIBSCALASCA_INCLUDE_DIRS
    pkg_check_modules(PKG_CONFIG_LIBSCALASCA QUIET scalasca)
endif(PKG_CONFIG_FOUND)

# 1.c) Lets try to get command line tool an execute it

# search command line tools
find_program(SCALASCA_COMMAND scalasca)
find_program(SCALASCA_KCONFIG kconfig)

# run scalasca and try to extract include hints
execute_process(COMMAND ${SCALASCA_KCONFIG} --user --cflags RESULT_VARIABLE SCALASCA_KCONFIG_CFLAGS_RES OUTPUT_VARIABLE SCALASCA_USER_CFLAGS)

if(NOT SCALASCA_KCONFIG_CFLAGS_RES)
    string(REGEX MATCHALL "[-][I]([^ ;\n])+" SCALASCA_COMMAND_INCLUDE_SEARCH_DIRS "${SCALASCA_USER_CFLAGS}")
    string(REGEX REPLACE "-I" "" SCALASCA_COMMAND_INCLUDE_SEARCH_DIRS "${SCALASCA_COMMAND_INCLUDE_SEARCH_DIRS}")
endif(NOT SCALASCA_KCONFIG_CFLAGS_RES)


# 2. Now we search for the includes and libraries at the hints

# now find SCALASCA_INCLUDE_DIR
find_path(SCALASCA_INCLUDE_DIR 
          NAMES epik_user.h
          HINTS ${SCALASCA_ENVVAR_INCLUDE_SEARCH_DIRS}
                ${SCALASCA_COMMAND_INCLUDE_SEARCH_DIRS}
                ${PKG_CONFIG_LIBSCALASCA_INCLUDE_DIRS})

message("Scalasca: SCALASCA_INCLUDE_DIR: ${SCALASCA_INCLUDE_DIR}")
message("Scalasca: SCALASCA_COMMAND_INCLUDE_SEARCH_DIRS: ${SCALASCA_COMMAND_INCLUDE_SEARCH_DIRS}")
message("Scalasca: SCALASCA_COMMAND: ${SCALASCA_COMMAND}")
message("Scalasca: SCALASCA_KCONFIG: ${SCALASCA_KCONFIG}")
message("Scalasca: SCALASCA_USER_CFLAGS: ${SCALASCA_USER_CFLAGS}")

# 3. set the variables
# handle the QUIETLY and REQUIRED arguments and set SCALASCA_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCALASCA  DEFAULT_MSG  SCALASCA_INCLUDE_DIR SCALASCA_COMMAND SCALASCA_KCONFIG SCALASCA_USER_CFLAGS)

mark_as_advanced(SCALASCA_INCLUDE_DIR SCALASCA_COMMAND SCALASCA_KCONFIG SCALASCA_USER_CFLAGS)
