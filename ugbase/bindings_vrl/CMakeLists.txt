cmake_minimum_required(VERSION 2.6)

####
# bindings_vrl Library
####

project(P_BINDINGS_VRL)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(${P_UGBASE_SOURCE_DIR})

include("../../ug_cmake_includes.txt")


set(srcBindingsVRL
                    bindings_vrl.cpp
                    type_converter.cpp
                    canvas.cpp
                    messaging.cpp
                    vrl_user_number.cpp
                    invocation.cpp
                    playground.cpp
                    user_data.cpp
)


set(Sources ${srcBindingsVRL})
EXPORTSOURCES(bindings_vrl  ${Sources})

if(NOT UG_SHARED)
add_library(bindings_vrl  ${Sources})
endif(NOT UG_SHARED)



####
# make global svn revision available as C++ macro
####

include(FindSubversion)

IF(Subversion_FOUND)
        MESSAGE(STATUS "Looking for Subversion - found")
        Subversion_WC_INFO(${CMAKE_SOURCE_DIR} "svn")
        EXECUTE_PROCESS(
                COMMAND "svnversion" "--no-newline" ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE svn_version
                ERROR_VARIABLE svn_version
        )
        #message( ${svn_version} )
        #add_definitions(-DSVN_REVISION=${svn_version})

        # write a file with the SVN_REVISION define
        file(WRITE svnrevision.h "#define SVN_REVISION \"${svn_version}\"\n")

ELSE(Subversion_FOUND)
        MESSAGE(STATUS "Looking for Subversion - not found")
        # write a file with the SVN_REVISION define
        file(WRITE svnrevision.h "#define SVN_REVISION \"unknown revision\"\n")
ENDIF(Subversion_FOUND)

####
# make compile date stamp available as C++ macro
####

IF(UNIX)
    MESSAGE(STATUS "On Unix: using date command to generate time stamp")
    execute_process(COMMAND date
        OUTPUT_VARIABLE compile_date OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    # write a file with the COMPILE_DATE define
    file(WRITE compiledate.h "#define COMPILE_DATE \"${compile_date}\"\n")
ELSEIF(WIN32)
    MESSAGE(STATUS "On Windows: time stamp generation unsupported!")
    # write a file with the COMPILE_DATE define
    file(WRITE compiledate.h "#define COMPILE_DATE \"unknown date\"\n")
ENDIF(UNIX)