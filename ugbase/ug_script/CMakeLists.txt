cmake_minimum_required(VERSION 2.6)

####
# uginterface Library
####

project(P_UGSCRIPT)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# This include directory is required by lib_grid.
include_directories(${P_UGBASE_SOURCE_DIR})

include_directories(${P_UGSCRIPT_SOURCE_DIR}/externals/lua)

include("../../ug_cmake_includes.txt")

set(srcLua	externals/lua/lapi.c
			externals/lua/lapi.h
			externals/lua/lauxlib.c
			externals/lua/lauxlib.h
			externals/lua/lbaselib.c
			externals/lua/lcode.c
			externals/lua/lcode.h
			externals/lua/ldblib.c
			externals/lua/ldebug.c
			externals/lua/ldebug.h
			externals/lua/ldo.c
			externals/lua/ldo.h
			externals/lua/ldump.c
			externals/lua/lfunc.c
			externals/lua/lfunc.h
			externals/lua/lgc.c
			externals/lua/lgc.h
			externals/lua/linit.c
			externals/lua/liolib.c
			externals/lua/llex.c
			externals/lua/llex.h
			externals/lua/llimits.h
			externals/lua/lmathlib.c
			externals/lua/lmem.c
			externals/lua/lmem.h
			externals/lua/loadlib.c
			externals/lua/lobject.c
			externals/lua/lobject.h
			externals/lua/lopcodes.c
			externals/lua/lopcodes.h
			externals/lua/loslib.c
			externals/lua/lparser.c
			externals/lua/lparser.h
			externals/lua/lstate.c
			externals/lua/lstate.h
			externals/lua/lstring.c
			externals/lua/lstring.h
			externals/lua/lstrlib.c
			externals/lua/ltable.c
			externals/lua/ltable.h
			externals/lua/ltablib.c
			externals/lua/ltm.c
			externals/lua/ltm.h
			externals/lua/lua.h
			externals/lua/luaconf.h
			externals/lua/lualib.h
			externals/lua/lundump.c
			externals/lua/lundump.h
			externals/lua/lvm.c
			externals/lua/lvm.h
			externals/lua/lzio.c
			externals/lua/lzio.h)


set(srcUGScript		ug_script.cpp
					info_commands.cpp
					bindings/bindings_lua.cpp)

# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(NOT NO_ALGEBRA)
	# add the source files
	set(srcUGScript ${srcUGScript}
					user_data/user_data.cpp
					user_data/boundary_number.cpp)
else(NOT NO_ALGEBRA)

endif(NOT NO_ALGEBRA)


set(Sources ${srcUGScript} ${srcLua})
EXPORTSOURCES(ug_script  ${Sources})

if(NOT UG_SHARED)
add_library(ugscript  ${Sources})
endif(NOT UG_SHARED)
