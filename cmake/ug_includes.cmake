################################################################################
# created by Sebastian Reiter, Andreas Vogel, Martin Rupp, Michael Hoffer, ...
# s.b.reiter@googlemail.com
#
# This file handles all options, with which ug can be compiled. Depending on
# those options, it then resolves all dependencies by adding required
# include-paths and by linking to required libraries.
# It also sets associated defines which can be queried in C and C++ Code.
#
# If you're creating an executable or a library which mainly depends on ug4,
# including this file in your CMakeLists.txt may be a good idea.
################################################################################

################################################################################
# Make sure code is only executed once.
# (no indent here, since it affects the whole file)
if(UG_CMAKE_INCLUDES_INCLUDED)
	return()
endif()
set(UG_CMAKE_INCLUDES_INCLUDED on)

################################################################################
# include ug header and library path
get_filename_component(UG_ROOT_CMAKE_PATH ${CMAKE_CURRENT_LIST_FILE} PATH)
set(UG_ROOT_PATH ${UG_ROOT_CMAKE_PATH}/../)
include_directories(${UG_ROOT_PATH}/ugbase)
include_directories(${CMAKE_BINARY_DIR})
link_directories(${UG_ROOT_PATH}/lib)
################################################################################
# include cmake functions
set(CMAKE_MODULE_PATH ${UG_ROOT_PATH}/cmake/modules)
include(${UG_ROOT_PATH}/cmake/compiler_flags.cmake)
# reset flags, because of maybe switched DEBUG option
reset_cxx_flags()

################################################################################
# We want the code to be built into one library (except the plugins, of course)
set(BUILD_ONE_LIB ON)


# Those variables control the build process. They are later modified depending
# on the value of TARGET.
set(buildUGShell OFF)
set(buildForVRL OFF)
set(buildForLUA OFF)
set(buildAlgebra OFF)
set(buildPluginSystem OFF)
set(buildEmbeddedPlugins OFF)
set(buildDynamicLib OFF)
set(buildGrid OFF)
set(buildDisc OFF)
set(buildBridge OFF)
set(buildBindings OFF)
set(buildRegistry OFF)
set(buildCompileInfo ON)
set(buildMetis OFF)
set(buildParmetis OFF)

# Here we'll store libs, if we find and need them
set(linkLibraries)

################################################################################
# In this section we'll define default variables
	
# Values for the TARGET option
set(targetOptions "ugshell, vrl, libug4, libgrid, ugplugin, gridshell, amg")
set(targetDefault "ugshell")
set(targetExecutableName ugshell)
set(targetLibraryName ug4)

# Values for the DIM option
set(dimOptions "1, 2, 3, ALL, \"1\;2\", \"1\;3\", \"2\;3"\")
set(dimDefault "ALL")

# Values for the CPU option
set(cpuOptions "1, 2, 3, 4, VAR, ALL, \"2\;4\", \"1\;3\;4\" , ..." )
set(cpuDefault "ALL")

# Values for the blas / lapack option
set(blasDefault ON)
set(lapackDefault ON)

# If enabled, the boost version contained in the externals folder will be used
set(internalBoostDefault ON)

# Option to add svn head revision, compile date and build host into ugshell's initial output
set(svnDefault ON)

# Precision of the number type
set(precisionDefault "double")
set(precisionOptions "single, double")

# Values for the PROFILER option
set(profilerOptions "None, Shiny, Scalasca, Vampir")
set(profilerDefault "None")

# If we run the script the first time, search for MPI to determine the default value
if(LOCAL_OPENMPI)
	message("local openmpi")
	set(MPI_C_INCLUDE_PATH "~/local/openmpi/used/include")
	set(MPI_CXX_INCLUDE_PATH "~/local/openmpi/used/include")
	set(MPI_C_LIBRARIES "mpi")
	set(MPI_CXX_LIBRARIES "mpi;mpi_cxx")
endif(LOCAL_OPENMPI)
if(BUILTIN_MPI)
	set(MPI_FOUND YES)
else(BUILTIN_MPI)
	if(NOT mpiHasBeenSearched)
		find_package(MPI)
		set(mpiHasBeenSearched ON CACHE BOOL "This var is set to true after mpi has been searched the first time.")
	endif(NOT mpiHasBeenSearched)
endif(BUILTIN_MPI)


set(posixDefault OFF)
if(WIN32)
	add_definitions(-DUG_WIN32)
	if(CYGWIN)
		add_definitions(-DUG_CYGWIN)
	endif(CYGWIN)
else(WIN32)
	if(UNIX)
		set(posixDefault ON)
	endif(UNIX)
endif(WIN32)


################################################################################
# All available options should be defined here:
# cmake-options can either be on or off.
# note: none-on/off-variables are set below, and the docstring is set at the end of this file


# the following options are real cmake-options
option(STATIC "Enables static linking. Valid options are: ON, OFF" OFF)
option(DEBUG "Enables debugging. Valid options are: ON, OFF" OFF)
option(DEBUG_LOGS "Enables debug output. Valid options are: ON, OFF" OFF)
option(PARALLEL "Enables parallel compilation. Valid options are: ON, OFF" ${MPI_FOUND})
option(PROFILE_PCL "Enables profiling of the pcl-library. Valid options are ON, OFF" OFF)
option(PROFILE_BRIDGE "Enables profiling of bridge objects. Valid options are ON, OFF" OFF)
option(PCL_DEBUG_BARRIER "Enables debug barriers in the pcl-library. Valid options are ON, OFF" OFF)
option(LAPACK "Lapack won't be used, even if available. Valid options are ON, OFF" ${lapackDefault})
option(BLAS "Blas won't be used, even if available. Valid options are ON, OFF" ${blasDefault})
option(METIS "Embeds the Metis library into UG. Valid options are ON, OFF" OFF)
option(PARMETIS "Embeds the Parmetis library for non commercial use. Valid options are ON, OFF" OFF)
option(INTERNAL_BOOST "If enabled, the boost version found in the externals directory will be used. Valid options are ON, OFF" ${internalBoostDefault})
option(BUILTIN_BLAS "BLAS is built into compiler" OFF)
option(BUILTIN_LAPACK "LAPACK is built into compiler" OFF)
option(BUILTIN_MPI "MPI is built into compiler" OFF)
option(OPENMP "Enables use of OpenMP. Valid options are ON, OFF" OFF)
option(CXX11 "Enables compilation with C++11 standard. Valid options are ON, OFF" OFF)
option(EMBEDDED_PLUGINS "Plugin sources are directly included in libug4. No dynamic loading required. Valid options are ON, OFF " OFF)
option(COMPILE_INFO "Embeds information on compile revision and date. Requires relinking of all involved libraries. Valid options are ON, OFF " ${buildCompileInfo})
option(POSIX "If enabled and available, some additional functionality may be available. Valid options are ON, OFF " ${posixDefault})
option(BUILD_UGDOCU "If enabled, every build builds a new completion file for ugIDE" OFF)
option(CRS_ALGEBRA "Use the CRS Sparse Matrix" OFF)
option(CPU_ALGEBRA "Use the old CPU Sparse Matrix" ON)

if(APPLE)
	option(USE_LUA2C "Use LUA2C" ON)
else(APPLE)
	option(USE_LUA2C "Use LUA2C" OFF)
endif(APPLE)


################################################################################
# set default values for pseudo-options
if(NOT TARGET)
	if(BUILDING_PLUGIN)
		set(TARGET ugplugin)
	else(BUILDING_PLUGIN)
		set(TARGET ${targetDefault})
	endif(BUILDING_PLUGIN)
endif(NOT TARGET)

if(NOT DIM)
	set(DIM ${dimDefault})
endif(NOT DIM)

if(NOT CPU)
	set(CPU ${cpuDefault})
endif(NOT CPU)

if(NOT PRECISION)
	set(PRECISION ${precisionDefault})
endif(NOT PRECISION)

if(NOT PROFILER)
	set(PROFILER ${profilerDefault})
endif(NOT PROFILER)

# convert the DIM and CPU sets to readable sets
# otherwise "2;3" gets "23" .
# todo: make this a function/macro
set(DIMReadable "")
foreach(d ${DIM})
	if("${DIMReadable}" STREQUAL "")
		set(DIMReadable ${d})
	else("${DIMReadable}" STREQUAL "")
		set(DIMReadable "${DIMReadable}" "\;" ${d})
	endif("${DIMReadable}" STREQUAL "")
endforeach(d)
set(CPUReadable "")
foreach(d ${CPU})
	if("${CPUReadable}" STREQUAL "")
		set(CPUReadable ${d})
	else("${CPUReadable}" STREQUAL "")
		set(CPUReadable "${CPUReadable}" "\;" ${d})
	endif("${CPUReadable}" STREQUAL "")
endforeach(d)



################################################################################
# We'll output the current options-setting in this section
message(STATUS "")
message(STATUS "Info: Current options:")
message(STATUS "Info: TARGET:            "${TARGET}" (options are: "${targetOptions}")")
message(STATUS "Info: DIM:               "${DIMReadable}" (options are: "${dimOptions}")")
message(STATUS "Info: CPU:               "${CPUReadable}" (options are: "${cpuOptions}")")
message(STATUS "Info: PRECISION:         "${PRECISION}" (options are: "${precisionOptions}")")
message(STATUS "Info: STATIC:            "${STATIC}" (options are: ON, OFF)")
message(STATUS "Info: DEBUG:             "${DEBUG}" (options are: ON, OFF)")
message(STATUS "Info: DEBUG_LOGS:        "${DEBUG_LOGS}" (options are: ON, OFF)")
message(STATUS "Info: PARALLEL:          "${PARALLEL}" (options are: ON, OFF)")
message(STATUS "Info: PCL_DEBUG_BARRIER: "${PCL_DEBUG_BARRIER}" (options are: ON, OFF)")
message(STATUS "Info: PROFILER:          "${PROFILER}" (options are: "${profilerOptions}")")
message(STATUS "Info: PROFILE_PCL:       "${PROFILE_PCL}" (options are: ON, OFF)")
message(STATUS "Info: PROFILE_BRIDGE:    "${PROFILE_BRIDGE}" (options are: ON, OFF)")
message(STATUS "Info: LAPACK:            "${LAPACK}" (options are: ON, OFF)")
message(STATUS "Info: BLAS:              "${BLAS}" (options are: ON, OFF)")
message(STATUS "Info: METIS:             "${METIS}" (options are: ON, OFF)")
message(STATUS "Info: PARMETIS:          "${PARMETIS}" (options are: ON, OFF)")
message(STATUS "Info: INTERNAL_BOOST:    "${INTERNAL_BOOST}" (options are: ON, OFF)")
message(STATUS "Info: EMBEDDED_PLUGINS   "${EMBEDDED_PLUGINS}" (options are: ON, OFF)")
message(STATUS "Info: COMPILE_INFO       "${COMPILE_INFO}" (options are: ON, OFF)")
message(STATUS "")
message(STATUS "Info: External libraries (path which contains the library or ON if you used uginstall):")
message(STATUS "Info: TETGEN:   "${TETGEN})
message(STATUS "Info: HYPRE:    "${HYPRE})
message(STATUS "Info: HLIBPRO:  "${HLIBPRO})
message(STATUS "")
message(STATUS "Info: C Compiler ID: ${CMAKE_C_COMPILER_ID}, C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}") 
message(STATUS "")

################################################################################
# Options are processed in this section.

########################################
# TARGET
add_definitions(-DUG_TARGET="${TARGET}") # (dummy) preprocessor directive used for displaying build configuration


# The target option enables specific build variables
if("${TARGET}" STREQUAL "ugshell")
	set(buildUGShell ON)
	set(buildForLUA ON)
	set(buildAlgebra ON)
	set(buildPluginSystem ON)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)
	
elseif("${TARGET}" STREQUAL "vrl")
	# The vrl works only, if a dynamic library is built.
	if(STATIC)
		message(FATAL_ERROR "ug4 for vrl can only be build as a dynamic library. Please set STATIC = OFF.")
	endif(STATIC)
	
	set(buildAlgebra ON)
	set(buildForVRL ON)
	set(buildPluginSystem ON)
	#todo: rename targetLibraryName to ug4_vrl
	#set(targetLibraryName ug4_vrl)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)
	
elseif("${TARGET}" STREQUAL "libug4")
	set(buildAlgebra ON)
	set(buildPluginSystem ON)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)	
	
elseif("${TARGET}" STREQUAL "libgrid")
	set(targetLibraryName grid)
	set(buildGrid ON)

elseif("${TARGET}" STREQUAL "ugplugin")
	set(buildAlgebra ON)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildForLUA ON)
	set(buildPluginSystem ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)

elseif("${TARGET}" STREQUAL "gridshell")
	set(targetExecutableName gridshell)
	set(buildUGShell ON)
	set(buildPluginSystem ON)
	set(buildGrid ON)
	set(buildDisc ON)
	set(buildBridge ON)
	set(buildBindings ON)
	set(buildRegistry ON)
	
elseif("${TARGET}" STREQUAL "amg")
	set(targetLibraryName ugamg)
	set(buildAlgebra ON)
	set(buildGrid OFF)
else("${TARGET}" STREQUAL "ugshell")
	message(FATAL_ERROR "Unsupported TARGET: "${TARGET}". Options are: "${targetOptions})
	
endif("${TARGET}" STREQUAL "ugshell")

########################################
# PRECISION
if("${PRECISION}" STREQUAL "single")
	add_definitions(-DUG_SINGLE_PRECISION)
elseif("${PRECISION}" STREQUAL "double")
	# Nothing to do here. double-precision is the default
else("${PRECISION}" STREQUAL "single")
	message(FATAL_ERROR "Unsupported PRECISION: "${PRECISION}". Options are: "${precisionOptions})
endif("${PRECISION}" STREQUAL "single")

########################################
# STATIC
# If STATIC is enabled, then a static ug-lib will be built. This will add a
# _s suffix to the resulting lib name.
# Note that no plugins are built if STATIC is enabled.


set(buildEmbeddedPlugins ${EMBEDDED_PLUGINS})

if(STATIC)
	set(buildDynamicLib OFF)
	set(buildEmbeddedPlugins ON)
	set(targetLibraryName ${targetLibraryName}_s)
	add_definitions(-DUG_STATIC)
else(STATIC)
	set(buildDynamicLib ON)
	set(UG_SHARED ON)
endif(STATIC)

if(buildEmbeddedPlugins)
	add_definitions(-DUG_EMBEDDED_PLUGINS)
endif(buildEmbeddedPlugins)


########################################
# DIM
# The dim option sets defines for C and C++
if("${DIM}" STREQUAL "ALL")
	add_definitions(-DUG_DIM_1 -DUG_DIM_2 -DUG_DIM_3)
	
else("${DIM}" STREQUAL "ALL")
	# DIM is a string of dimensions (e.g. "1;2")
	# loop dims
	foreach(d ${DIM})
		# check if dim is valid
		if(d GREATER 3 OR d LESS 1)
			message(FATAL_ERROR "ERROR: Cannot build world dimension ${d}. "
								"Valid options are: "${dimOptions})
		endif(d GREATER 3 OR d LESS 1)
	
		add_definitions(-DUG_DIM_${d})
	endforeach(d)
endif("${DIM}" STREQUAL "ALL")

########################################
# CPU
# The cpu option sets defines for C and C++
if(CRS_ALGEBRA)
	MESSAGE(STATUS "Info: Using CRS Matrix Algebra.")

    if("${CPU}" STREQUAL "ALL")
        # todo checks for 4, VAR (-DUG_CPU_4 -DUG_CPU_VAR)!
        add_definitions(-DUG_CRS_1 -DUG_CRS_2 -DUG_CRS_3)

    else("${CPU}" STREQUAL "ALL")
        # CPU is a string of numbers (e.g. "1;2")
        # loop dims
        foreach(d ${CPU})
            # check if dim is valid
            if(d GREATER 4 OR d LESS 1)
                message(FATAL_ERROR "ERROR: Cannot build cpu blocksize ${d}. "
                                    "Valid options are: "${cpuOptions})
            endif(d GREATER 4 OR d LESS 1)

            add_definitions(-DUG_CRS_${d})
        endforeach(d)
    endif("${CPU}" STREQUAL "ALL")
endif()
if(CPU_ALGEBRA)
    MESSAGE(STATUS "Info: Using CPU Matrix Algebra.")
    if("${CPU}" STREQUAL "ALL")
        # todo checks for 4, VAR (-DUG_CPU_4 -DUG_CPU_VAR)!
        add_definitions(-DUG_CPU_1 -DUG_CPU_2 -DUG_CPU_3)

    else("${CPU}" STREQUAL "ALL")
        # CPU is a string of numbers (e.g. "1;2")
        # loop dims
        foreach(d ${CPU})
            # check if dim is valid
            if(d GREATER 4 OR d LESS 1)
                message(FATAL_ERROR "ERROR: Cannot build cpu blocksize ${d}. "
                                    "Valid options are: "${cpuOptions})
            endif(d GREATER 4 OR d LESS 1)

            add_definitions(-DUG_CPU_${d})
        endforeach(d)
    endif("${CPU}" STREQUAL "ALL")
endif()

########################################
# COMPILER specific flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
	# Cray Compiler: add support for gnu extensions
	add_cxx_flag("-h gnu")
	# remove warning "The controlling expression is constant"
	add_cxx_flag("-hnomessage=236")
else()
	#add_cxx_flag("-Wall")
	#§set(CMAKE_CPP_FLAGS	"${CMAKE_CPP_FLAGS} -Wno-overloaded-virtual -Wno-autological-compare" CACHE STRING "overriden flags!" FORCE)
endif()

########################################
# DEBUG
########################################
if(DEBUG)
	add_definitions(-DUG_DEBUG)
	# if no build type set, add default debug flags
	if("${CMAKE_BUILD_TYPE}" STREQUAL "None" OR
		"${CMAKE_BUILD_TYPE}" STREQUAL "")
		# use cmake standard flags for debug builds (-g)
		add_cxx_flag(${CMAKE_CXX_FLAGS_DEBUG})
		# if user specified additional debug flags add them
		if(DEBUG_FORMAT)
			add_cxx_flag(${DEBUG_FORMAT})
			message(STATUS "Info: Debug Information is ${DEBUG_FORMAT}") 
		endif(DEBUG_FORMAT)
	elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
		message(WARNING "Build type set to Release, but DEBUG build wanted.")
	elseif("${CMAKE_BUILD_TYPE}" STREQUAL "MinSizeRel")
		message(WARNING "Build type set to MinSizeRel, but DEBUG build wanted.")
	endif()
	
	# This code would enable strict bounds checking for STL objects like in vector::operator[]. 
	# however, GLIBCXX_DEBUG and strstream don't work together on mac (bug: http://bit.ly/cH78bC). 
	# when this bug is fixed, one could set those flags (or similar) depending on the compiler.
	
	# I disabled the following definitions, since they cause several additional problems.
	# First they also cause problems on MinGW builds.
	# Secondly they have to be defined also for all executables, which link against
	# ug. Since they are well hidden, this may cause quite some headache.
	# One should think about a flag DEBUG_STL or something like this, which explicitly
	# enables those flags. (sreiter)
	#if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT APPLE)
	#	add_definitions(-D_GLIBCXX_DEBUG=1 -D_GLIBCXX_DEBUG_PEDANTIC=1)
	#ENDIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT APPLE)
else(DEBUG)
	# add release definitions
	add_definitions(-DBOOST_UBLAS_NDEBUG)
	if(NOT CMAKE_BUILD_TYPE)
		# if DEBUG=OFF also add standard cmake release flags (-O3 -DNDEBUG)
		add_cxx_flag(${CMAKE_CXX_FLAGS_RELEASE})
	elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
		message(WARNING "Build type set to Debug, but DEBUG=OFF. Leads to strange cflags!")
	endif()

	# compiler specific release flags
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		add_cxx_flag("-funroll-loops -ftree-vectorize")
	elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
		add_cxx_flag("-hipa5 -hunroll2")
	endif()
endif(DEBUG)

if(DEBUG_LOGS)
	add_definitions(-DUG_ENABLE_DEBUG_LOGS)
endif(DEBUG_LOGS)

# if build type is set print own cflags and flags from build type 
if(CMAKE_BUILD_TYPE)
	string(TOUPPER ${CMAKE_BUILD_TYPE} bt_upper)
	message(STATUS "Info: compiling with cxx flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${bt_upper}}")
else()
	message(STATUS "Info: compiling with cxx flags: ${CMAKE_CXX_FLAGS}")
endif()

########################################
# PROFILER
if(NOT("${PROFILER}" STREQUAL "None"))
    if("${PROFILER}" STREQUAL "Shiny")
    	add_definitions(-DUG_PROFILER_SHINY)    
    
    # Scalasca
    elseif("${PROFILER}" STREQUAL "Scalasca")
        find_package(Scalasca)
        if(SCALASCA_FOUND)
            message("-- Info: Scalasca: using scalasca command: ${SCALASCA_COMMAND}")
            message("-- Info: Scalasca: using inlcude dir: ${SCALASCA_INCLUDE_DIR}")
        else(SCALASCA_FOUND)
        	message(FATAL_ERROR "PROFILER: ${PROFILER}: Cannot find required "
        	        "binary scalasca/kconfig. Make sure "
        	        "that PATH contains scalasca and kconfig executable.")
        endif(SCALASCA_FOUND)
        
        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "scalasca -instrument -comp=none -user ")
    	add_definitions(-DUG_PROFILER_SCALSACA)    
    
    # Vampir
    elseif("${PROFILER}" STREQUAL "Vampir")
        find_package(VampirTrace)
        if(VAMPIRTRACE_FOUND)
            message("-- Info: Vampir: using compiler wrapper c++: ${VAMPIRTRACE_CXX}")
            message("-- Info: Vampir: using compiler wrapper cc : ${VAMPIRTRACE_CC}")
            message("-- Info: Vampir: using inlcude dir: ${VAMPIRTRACE_INCLUDE_DIR}")
            message("-- Info: Vampir: using library dir: ${VAMPIRTRACE_LIBRARIES}")
        else(VAMPIRTRACE_FOUND)
        	message(FATAL_ERROR "PROFILER: ${PROFILER}: Cannot find required "
        	        "binary vtcxx and/or library or include paths. Make sure "
        	        "that PATH contains vtcxx executable.")
        endif(VAMPIRTRACE_FOUND)
    	add_definitions(-DUG_PROFILER_VAMPIR)    

    # wrong string in compiler
    else("${PROFILER}" STREQUAL "Shiny")
    	message(FATAL_ERROR "Unsupported PROFILER: ${PROFILER}. Options are: "${profilerOptions})
    endif("${PROFILER}" STREQUAL "Shiny")
    
    # if one profiler is used this flag will be set
	add_definitions(-DUG_PROFILER)    # add to c++ flags
	set(UG_PROFILER ON)               # add Cmake variable
	
endif(NOT("${PROFILER}" STREQUAL "None"))


########################################
# PROFILE_PCL
if(PROFILE_PCL)
	add_definitions(-DPROFILE_PCL)
endif(PROFILE_PCL)

# PROFILE_BRIDGE
if(PROFILE_BRIDGE)
	add_definitions(-DPROFILE_BRIDGE)
endif(PROFILE_BRIDGE)


########################################
# PCL_DEBUG_BARRIER
if(PCL_DEBUG_BARRIER)
	add_definitions(-DPCL_DEBUG_BARRIER_ENABLED)
endif(PCL_DEBUG_BARRIER)


########################################
# OPENMP
IF(OPENMP)
	# writing to CMAKE_CXX_FLAGS directly can cause problems on some platforms.
	# Please use add_definitions instead. Hope that works in this case, too. sreiter.
  # The linker option '-lgomp' or '-liomp5' can not be added to add_definitions
  # as these are not passed to the link then. But they have to. tklatt.
	#	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")
  IF(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_cxx_flags("-fopenmp")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgomp")
    ADD_DEFINITIONS(-DUG_OPENMP)
    MESSAGE(STATUS "Info: Using OpenMP (experimental)")
  ELSEIF(CMAKE_C_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_cxx_flags("-fopenmp")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -liomp5")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -liomp5")
    ADD_DEFINITIONS(-DUG_OPENMP)
    MESSAGE(STATUS "Info: Using OpenMP (experimental)")
  ELSEIF(CMAKE_C_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    MESSAGE(WARNING "Clang does not support OpenMP yet.")
  ELSE()
    MESSAGE(WARNING "Don't know compiler type, thus don't know how to enable OpenMP")
  ENDIF()
ENDIF(OPENMP)




########################################
# C++11
# Note 1: C++11 is still experimental with most compilers, though the most
# important features are already supported.
# Note 2: Though Clang supports C++11, adding the compiler flag globally (as done
# for GCC below) will break the build on trying to compile C/ObjC files with this
# flag. GCC only emits warnings, Clang throws an error.
# TODO: Workaround Note 2. Probably requires larger refactoring of all CMake files.
IF(CXX11)
	IF(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
		# Check for the GCC's C++11 capabilities
		INCLUDE(CheckCXXCompilerFlag)
		CHECK_CXX_COMPILER_FLAG(-std=c++0x HAVE_CXX0X)
		# since GCC4.7 (c++0x will be removed in future versions of GCC)
		CHECK_CXX_COMPILER_FLAG(-std=c++11 HAVE_CXX11)
		# Add appropriate compiler flags
		IF(HAVE_CXX11)
			SET(CXX11_FLAG "-std=c++11")
		ELSEIF(HAVE_CXX0X)
			SET(CXX11_FLAG "-std=c++0x")
		ENDIF()

		IF(CXX11_FLAG)
			ADD_DEFINITIONS(-DUG_CXX11)
			add_cxx_flag(${CXX11_FLAG})
			MESSAGE(STATUS "Info: C++11 enabled. (flag: ${CXX11_FLAG})")
		ELSE()
			SET(CXX11 OFF)
			MESSAGE(STATUS "Info: Compiler does not support C++11 standard.")
		ENDIF()
	ELSE()
		MESSAGE(STATUS "Info: Enabling C++11 is only supported with GCC currently.")
	ENDIF()
ENDIF(CXX11)


########################################
# buildAlgebra
if(buildAlgebra)
	add_definitions(-DUG_ALGEBRA)
	
	# we'll check for lapack and blas here
	if(LAPACK OR BLAS)
		# On OSX, we know where lapack and blas are located. To avoid errors
		# with Fortran compilers on OSX, we'll add APPLE as a special case here.
		if(APPLE)
			if(LAPACK)
				set(LAPACK_LIBRARIES "-framework vecLib" CACHE STRING "LAPACK library" FORCE)
				FIND_PATH(LAPACK_INCLUDE_PATH clapack.h /usr/local/include/ /usr/include /include)
				if(LAPACK_INCLUDE_PATH)
		  			set(LAPACK_FOUND YES)
		  		endif(LAPACK_INCLUDE_PATH)
			endif(LAPACK)
			
			if(BLAS)
				set(BLAS_LIBRARIES "-framework vecLib" CACHE STRING "CBLAS library" FORCE)
				find_path(BLAS_INCLUDE_PATH cblas.h /usr/local/include/ /usr/include /include)
				set(BLAS_FOUND YES)
				if(BLAS_INCLUDE_PATH)
					set(BLAS_FOUND YES)
				endif(BLAS_INCLUDE_PATH)
			endif(BLAS)
			
		elseif(NOT BUILTIN_BLAS AND NOT BUILTIN_LAPACK)			
			# a fortran compiler is required to find the packages
			# Sadly there seems to be a cmake-bug, hence the following workaround.
			# ENABLE_LANGUAGE workaround begin (issue 0009220)
			message(STATUS "Info: If problems with the Fortran compiler occur, consider deactivating LAPACK and BLAS")
			if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
			  set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
			endif(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
			# ENABLE_LANGUAGE workaround end (issue 0009220)
			
			ENABLE_LANGUAGE(Fortran OPTIONAL)
			if(CMAKE_Fortran_COMPILER_WORKS)
				if(LAPACK)
					find_package(LAPACK)
				endif(LAPACK)
				
				if(BLAS)
					find_package(BLAS)
				endif(BLAS)
				
			endif(CMAKE_Fortran_COMPILER_WORKS)
		endif(APPLE)
		
	# We'll output whether lapack and blas are used, to avoid misconceptions
		if(LAPACK_FOUND)
			message(STATUS "Info: Using Lapack")
			include_directories (${LAPACK_INCLUDE_PATH})
			set(linkLibraries ${linkLibraries} ${LAPACK_LIBRARIES})
			add_definitions(-DLAPACK_AVAILABLE)
		elseif(BUILTIN_LAPACK)
			message(STATUS "Info: Using Builtin Lapack")
			add_definitions(-DLAPACK_AVAILABLE)
		else(LAPACK_FOUND)	
			message(STATUS "Info: Not using Lapack. No package found.")
		endif(LAPACK_FOUND)
		
		if(BLAS_FOUND)
			message(STATUS "Info: Using Blas")
			include_directories (${BLAS_INCLUDE_PATH})
			set(linkLibraries ${linkLibraries} ${BLAS_LIBRARIES})
			add_definitions(-DBLAS_AVAILABLE)
		elseif(BUILTIN_BLAS)
			message(STATUS "Info: Using Builtin Blas")
			add_definitions(-DBLAS_AVAILABLE)
		else(BLAS_FOUND)	
			message(STATUS "Info: Not using Blas. No package found.")
		endif(BLAS_FOUND)
		
		if(BLAS_goto2_LIBRARY)
			message(STATUS "Info: GotoBLAS2 found: (${BLAS_goto2_LIBRARY}). Adding -lgfortran.")
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgfortran")
		endif()
	endif(LAPACK OR BLAS)
endif(buildAlgebra)


########################################
# METIS
if(METIS OR PARMETIS)
	add_definitions(-DUG_METIS)
	include_directories(${UG_ROOT_PATH}/externals/metis-5.0.2/include)
	set(buildMetis ON)
	set(linkLibraries ${linkLibraries} metis)
endif(METIS OR PARMETIS)

########################################
# PARMETIS
if(PARMETIS)
	add_definitions(-DUG_PARMETIS)
	include_directories(${UG_ROOT_PATH}/externals/parmetis-4.0.2/include)
	set(buildParmetis ON)
	set(linkLibraries ${linkLibraries} parmetis)
endif(PARMETIS)

########################################
# TETGEN
if(TETGEN)
	find_library(TETGEN_LIBS NAMES tet PATHS ${TETGEN})
	if(TETGEN_LIBS-NOTFOUND)
		message(FATAL_ERROR "ERROR: Couldn't find TETGEN in the specified path.")
	else(TETGEN_LIBS-NOTFOUND)
		add_definitions(-DUG_TETGEN -DTETLIBRARY)
		include_directories(${TETGEN})
		set(linkLibraries ${linkLibraries} ${TETGEN_LIBS})
	endif(TETGEN_LIBS-NOTFOUND)
endif(TETGEN)

########################################
# HYPRE
if(HYPRE)
	find_library(HYPRE_LIBS HYPRE PATHS ${HYPRE})

	if(HYPRE_LIBS-NOTFOUND)
		message(FATAL_ERROR "ERROR: Couldn't find HYPRE in the specified path.")
	else(HYPRE_LIBS-NOTFOUND)
		add_definitions(-DUG_HYPRE -DHYPRELIB_DIR) # TODO: Is '-DHYPRELIB_DIR' used?
		include_directories(${HYPRE}/../include/)
		set(linkLibraries ${linkLibraries} ${HYPRE_LIBS})
	endif(HYPRE_LIBS-NOTFOUND)
endif(HYPRE)

########################################
# HLIBPRO
if(HLIBPRO)
#	find_library(HLIBPRO_LIBS NAMES hpro PATHS ${HLIBPRO})
	find_library(HLIBPRO_LIBS NAMES libhpro.dylib PATH_SUFFIXES ../hlibpro-0.13.6/lib/ ${HLIBPRO}) # hlib built as shared lib via 'scons shared=1 static=0' (26092011ih)
	find_path (HLIBPROLIB_DIR libhpro.dylib
		PATHS ENV PATH
		PATH_SUFFIXES ../hlibpro-0.13.6/lib/ )
#	message(STATUS "INFO: Using HLibPro; content of 'HLIBPRO_LIBS' is '${HLIBPRO_LIBS}', content of 'HLIBPROLIB_DIR' is '${HLIBPROLIB_DIR}'.")
	if(HLIBPRO_LIBS-NOTFOUND)
		message(FATAL_ERROR "ERROR: Couldn't find HLIBPRO in the specified path.")
	else(HLIBPRO_LIBS-NOTFOUND)
#		add_definitions(-DHLIBPROLIB_DIR)
		add_definitions(-DUG_HLIBPRO)
		include_directories(${HLIBPROLIB_DIR}/../include/)
		include_directories(${HLIBPROLIB_DIR}/../src/include/)
		set(linkLibraries ${linkLibraries} ${HLIBPRO_LIBS})
	endif(HLIBPRO_LIBS-NOTFOUND)
endif(HLIBPRO)

#########################################
# CUDA
if(CUDA)
	SET(MY_HOSTNAME $ENV{HOSTNAME})
	if(MY_HOSTNAME STREQUAL "tales")
        SET(CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda-5.0)
    endif()    
    find_package(CUDA)
	if(CUDA_FOUND)
		SET(CUDA_NVCC_FLAGS "-arch=sm_13")
		MESSAGE(STATUS "Info: CUDA ${CUDA_VERSION}. Toolkit root dir: ${CUDA_TOOLKIT_ROOT_DIR}, nvcc: ${CUDA_NVCC_EXECUTABLE}")

		file(WRITE ${CMAKE_BINARY_DIR}/cuda_jit.h
			"#define CUDA_TOOLKIT_ROOT_DIR \"${CUDA_TOOLKIT_ROOT_DIR}\"\n"
			"#define CUDA_NVCC_EXECUTABLE \"${CUDA_NVCC_EXECUTABLE}\"\n"
			"#define CUDA_VERSION \"${CUDA_VERSION}\"\n"
			"#define CUDA_NVCC_FLAGS \"${CUDA_NVCC_FLAGS}\"\n")

		include_directories(${CUDA_TOOLKIT_INCLUDE})
		include_directories(${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)
        
		set(linkLibraries ${linkLibraries} ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY})   
		add_definitions(-DCUDA_AVAILABLE)
	else(CUDA_FOUND)
		message(FATAL_ERROR "Error: Couldn't find CUDA")
	endif()
endif()

#########################################
# OpenCL
if(OpenCL)
	IF(APPLE)
		MESSAGE(STATUS "Info: Found OpenCL: -framework OpenCL (Apple)")
		SET(OPENCL_FOUND YES)
		set(linkLibraries "${linkLibraries} -framework OpenCL")   
		add_definitions(-DOPENCL_AVAILABLE)
	ELSE(APPLE)
		INCLUDE ("${UG_ROOT_PATH}/cmake/FindOpenCL.cmake")	
		if(OPENCL_FOUND)
			MESSAGE(STATUS "Info: Found OpenCL. Inc: ${OPENCL_INCLUDE_DIRS}, lib: ${OPENCL_LIB_DIR}")
			include_directories(${OPENCL_INCLUDE_DIRS})
			set(linkLibraries ${linkLibraries} ${OPENCL_LIB_DIR})   
			add_definitions(-DOPENCL_AVAILABLE)
		else(OPENCL_FOUND)
			MESSAGE(FATAL_ERROR " OpenCL not found! Aborting.")
		endif(OPENCL_FOUND)
	ENDIF(APPLE)
endif(OpenCL)

#########################################
if(USE_LUA2C)
	MESSAGE(STATUS "Info: Using LUA2C")
	add_definitions(-DUSE_LUA2C)
endif(USE_LUA2C)


################################################################################
# find and collect required libraries

########################################
# boost (required)
if(INTERNAL_BOOST)
	include_directories(${UG_ROOT_PATH}/externals/boost_1_48_0/)
	message(STATUS "Info: Using externals/boost_1_48_0")
else(INTERNAL_BOOST)
	find_package(Boost REQUIRED)
	include_directories(${Boost_INCLUDE_DIRS})

	if(Boost_FOUND)
		if(Boost_MAJOR_VERSION GREATER 1 OR Boost_MINOR_VERSION GREATER 39)
			#message(STATUS "Info: Using BOOST (version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}, libs at '${Boost_INCLUDE_DIRS}').")
			
		else(Boost_MAJOR_VERSION GREATER 1 OR Boost_MINOR_VERSION GREATER 39)
			# we require a newer version of boost
			message(FATAL_ERROR " BOOST in ${Boost_INCLUDE_DIRS} is not compatible (required 1.40.0, but is ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}) . NO Boost available. Aborting.")	
		endif(Boost_MAJOR_VERSION GREATER 1 OR Boost_MINOR_VERSION GREATER 39)
	else(Boost_FOUND)
		message(FATAL_ERROR " BOOST not found! NO Boost available, but at least boost 1.40.0 is required. Aborting.")
	endif(Boost_FOUND)
	
endif(INTERNAL_BOOST)


########################################
# dynamic linking
if(UNIX)
    if(STATIC)
        set(linkLibraries ${linkLibraries})
    else(STATIC)
        set(linkLibraries ${linkLibraries} dl)
    endif(STATIC)
elseif(WIN32)
	set(linkLibraries ${linkLibraries} Kernel32)
endif(UNIX)


########################################
# MPI
if(PARALLEL)
	if(BUILTIN_MPI)
		add_definitions(-DUG_PARALLEL)
	else(BUILTIN_MPI)
		# search mpi
		find_package(MPI)
		# MPI is required for parallel builds
		if(MPI_FOUND)
			add_definitions(-DUG_PARALLEL)
			include_directories(${MPI_INCLUDE_PATH})
			# Add mpi libraries:
			# Standard case: add cxx libraries
			if(MPI_CXX_LIBRARIES)
				set(linkLibraries ${linkLibraries} ${MPI_CXX_LIBRARIES})
			# Depreciated case: In order to support cmake versions < 2.8.6, 
			#                   where MPI_CXX_LIBRARIES cannot be used
			else(MPI_CXX_LIBRARIES)				
				set(linkLibraries ${linkLibraries} ${MPI_LIBRARY})
				if(MPI_EXTRA_LIBRARY)
				set(linkLibraries ${linkLibraries} ${MPI_EXTRA_LIBRARY})
				endif(MPI_EXTRA_LIBRARY)
			endif(MPI_CXX_LIBRARIES)				
		else(MPI_FOUND)
			message(FATAL_ERROR "MPI not found. Please set PARALLEL to OFF (run cmake -DPARALLEL=OFF ...)")
		endif(MPI_FOUND)
	endif(BUILTIN_MPI)
endif(PARALLEL)


########################################
set(buildCompileInfo ${COMPILE_INFO})
if(buildCompileInfo)
	message(STATUS "Info: COMPILE_INFO enabled. Causes relinking on each run of make.")
endif(buildCompileInfo)


########################################
if(POSIX)
	add_definitions(-DUG_POSIX)
endif(POSIX)

########################################
if(HERMIT_EXPERIMENTAL)
	add_definitions(-DUG_PARALLEL)
	add_definitions(-DBLAS_AVAILABLE)
	add_definitions(-DLAPACK_AVAILABLE)
	message(STATUS "Info: Using experimental stuff on Hermit.")
endif(HERMIT_EXPERIMENTAL)

########################################
# JNI
if(buildForVRL)
	find_package(JNI REQUIRED)
    include_directories(${JNI_INCLUDE_DIRS})
	add_definitions(-DUG_FOR_VRL)
endif(buildForVRL)

########################################
# LUA
if(buildForLUA)
	add_definitions(-DUG_FOR_LUA)
endif(buildForLUA)


################################################################################
# This is a hidden option, which is currently only required for builds on jugene.
# Note that the associated option is enabled in toolchain_file_jugene.cmake
if(enableDynamicOption)
	add_definitions(-dynamic)
endif(enableDynamicOption)


################################################################################
# Those options are temporary and should be removed in future builds.
# They are left from the old build script.
if(buildPluginSystem)
	add_definitions(-DUG_PLUGINS)
endif(buildPluginSystem)
if(buildBridge)
	add_definitions(-DUG_BRIDGE)
endif(buildBridge)
set(UG_DEBUG ${DEBUG})
set(USE_NEW_CMAKE_INCLUDES ON)


################################################################################
# link against required libraries
link_libraries(${linkLibraries})

################################################################################
# Finally declare a method that allows all sub-cmake-files to add their sources
# to a common library.

#######################
# Export sources to global variable.
# PURPOSE: use this function to add local package sources to global
#          ugSources property which is used to build libug4
# @param prefix current directory prefix
# @param sources sources list to export
function(ExportSources prefix sources)
    # iterate over all arguments and insert given prefix
    foreach(l ${ARGV})
        # FIXME: this is a hack to omit the first argument
        #        which is the prefix. Shall we use boolean or index variable?
        if(NOT "${l}" STREQUAL "${prefix}")
            # retrieve the global property ugSources and store it
            # in tmp variable
            # NOTE: properties must be assigned to variables before being used
            get_property(tmp GLOBAL PROPERTY ugSources)
            # append tmp to the global ugSources property using the correct prefix
            if("${prefix}" STREQUAL "")
				set_property(GLOBAL PROPERTY ugSources ${tmp} "${l}")
			else("${prefix}" STREQUAL "")
				set_property(GLOBAL PROPERTY ugSources ${tmp} "${prefix}/${l}")
			endif("${prefix}" STREQUAL "")
        endif(NOT "${l}" STREQUAL "${prefix}")
    endforeach(l)
endfunction(ExportSources)


######################################################################################################################
# the following options are pseudo cmake-options (normal options only support ON and OFF).
# Their default values are defined in the section above.
# we need to put this here so we only have to set the doc string and the type once.
# note: do not change the variables after this

set(TARGET ${TARGET} CACHE STRING "Set the target of this build process. Valid options are: ${targetOptions}")
set(DIM ${DIM} CACHE STRING "Set the dimension for which ug4 is build. Valid options are: ${dimOptions}")
set(CPU ${CPU} CACHE STRING "Set block sizes for which ug4 is build. Valid options are: ${cpuOptions}")
set(PRECISION ${PRECISION} CACHE STRING "Set the precision of the number type. Valid options are: ${precisionOptions}")
set(PROFILER ${PROFILER} CACHE STRING "Set the a profiler. Valid options are: ${profilerOptions}")

# the following options too are pseudo cmake-options. However, they should
# contains pathes, if set.
set(TETGEN ${TETGEN} CACHE PATH "Sets the path in which tetgen shall be searched.")
set(HYPRE ${HYPRE} CACHE PATH "Sets the path in which hypre shall be searched.")
set(HLIBPRO ${HLIBPRO} CACHE PATH "Sets the path in which hlibpro shall be searched.")

set(DEBUG_FORMAT ${DEBUG_FORMAT} CACHE STRING "Debug format options like -g, -gstabs, -ggbd. If not set, debug format is -g.")


# mark some stuff as advanced, so it won't show up in ccmake / CMakeSetup.exe
mark_as_advanced(CMAKE_INSTALL_PREFIX)
IF(APPLE)
	mark_as_advanced(CMAKE_OSX_ARCHITECTURES)
	mark_as_advanced(CMAKE_OSX_DEPLOYMENT_TARGET)
	mark_as_advanced(CMAKE_OSX_SYSROOT)
ENDIF(APPLE)
mark_as_advanced(LAPACK_INCLUDE_PATH)
mark_as_advanced(LAPACK_LIBRARIES)
mark_as_advanced(libReadline)
mark_as_advanced(mpiHasBeenSearched)
mark_as_advanced(BLAS_INCLUDE_PATH)
mark_as_advanced(BLAS_LIBRARIES)
#mark_as_advanced(CMAKE_BUILD_TYPE)
mark_as_advanced(CMAKE_CXX_FLAGS)
mark_as_advanced(CMAKE_C_FLAGS)

mark_as_advanced(BUILTIN_LAPACK)
mark_as_advanced(BUILTIN_BLAS)
mark_as_advanced(BUILTIN_MPI)

