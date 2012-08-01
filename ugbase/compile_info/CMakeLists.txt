cmake_minimum_required(VERSION 2.6)

################################################
# ug4 - compile-info
################################################

project(P_UG4_COMPILE_INFO)

# set global output paths for UG4
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include("../../ug_cmake_includes.txt")


IF(EMBED_SVN_REV)
	include(FindSubversion)

	IF(Subversion_FOUND)
		#MESSAGE(STATUS "Info: Embedding svn-revision number.")
	        EXECUTE_PROCESS(
	                COMMAND "svnversion" "--no-newline" ${CMAKE_SOURCE_DIR}
	                OUTPUT_VARIABLE svn_version
	                ERROR_VARIABLE svn_version_error
	        )
	
			IF(svn_version_error)
				set(SVN_REVISION "\"svnversion failed\"")				
			ELSE(svn_version_error)
	        	# write a file with the SVN_REVISION define
	        	set(SVN_REVISION "\"${svn_version}\"")
	        ENDIF(svn_version_error)
	
	ELSE(Subversion_FOUND)
		MESSAGE(STATUS "WARNING: Looking for Subversion - not found (disable EMBED_SVN_REV with -DEMBED_SVN_REV=OFF)")
	        # write a file with the SVN_REVISION define
	        set(SVN_REVISION "\"unknown revision\"")
	ENDIF(Subversion_FOUND)
	
ELSE(EMBED_SVN_REV)
	#MESSAGE(STATUS "INFO: Not looking for Subversion - (Enable EMBED_SVN_REV with -DEMBED_SVN_REV=ON)")
        # write a file with the SVN_REVISION define
        set(SVN_REVISION  "\"unknown revision\"")
ENDIF(EMBED_SVN_REV)

IF(UNIX)
    execute_process(COMMAND date
        OUTPUT_VARIABLE compile_date OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # write a file with the COMPILE_DATE define
    set(COMPILE_DATE "\"${compile_date}\"")

    execute_process(COMMAND hostname
        OUTPUT_VARIABLE hostname OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # write a file with the BUILD_HOST define
    set(HOSTNAME "\"${CMAKE_SYSTEM_NAME} ${hostname}\"")
    
    execute_process(COMMAND touch "-t" "999901010000.00" "${CMAKE_BINARY_DIR}/CMakeCache.txt")   

ELSEIF(WIN32)
    #MESSAGE(STATUS "INFO: On Windows: time stamp generation unsupported!")
    # write a file with the COMPILE_DATE define
    set(COMPILE_DATE "__DATE__")
    #MESSAGE(STATUS "INFO: On Windows: hostname generation unsupported!")
    # write a file with the BUILD_HOST define
    set(HOSTNAME "\"${CMAKE_SYSTEM_NAME} WIN32\"")
ENDIF(UNIX)

file(WRITE ${CMAKE_BINARY_DIR}/compile_info.cpp 
"namespace ug\n"
"{\n"
"const char *SVN_REVISION=${SVN_REVISION};\n"
"const char *BUILD_HOST=${HOSTNAME};\n"
"const char *COMPILE_DATE=${COMPILE_DATE};\n"
"}\n")

if(buildEmbeddedPlugins)
	EXPORTSOURCES(/ ${CMAKE_BINARY_DIR}/compile_info.cpp)
else(buildEmbeddedPlugins)
	if(buildDynamicLib)
		add_library(CompileInfo SHARED ${CMAKE_BINARY_DIR}/compile_info.cpp)
	else(buildDynamicLib)
		add_library(CompileInfo_s STATIC ${CMAKE_BINARY_DIR}/compile_info.cpp)
	endif(buildDynamicLib)
endif(buildEmbeddedPlugins)
