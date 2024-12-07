# Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
# Author: Sebastian Reiter
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

cmake_minimum_required(VERSION 2.8.12...3.27.1)

####
# uginterface Library
####

project(P_UGSCRIPT)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(${P_UGSCRIPT_SOURCE_DIR}/externals/lua)

include("../../../cmake/ug_includes.cmake")

if(UNIX)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLUA_USE_POSIX")
endif(UNIX)

if (NOT USE_LUAJIT) 
# use default version (shipped with ug4) 
set(srcLua	externals/lua/src/lapi.c
			externals/lua/src/lauxlib.c
			externals/lua/src/lbaselib.c
			externals/lua/src/lcode.c
			externals/lua/src/lctype.c
			externals/lua/src/lcorolib.c
			externals/lua/src/lutf8lib.c
			externals/lua/src/ldblib.c
			externals/lua/src/ldebug.c
			externals/lua/src/ldo.c
			externals/lua/src/ldump.c
			externals/lua/src/lfunc.c
			externals/lua/src/lgc.c
			externals/lua/src/linit.c
			externals/lua/src/liolib.c
			externals/lua/src/llex.c
			externals/lua/src/lmathlib.c
			externals/lua/src/lmem.c
			externals/lua/src/loadlib.c
			externals/lua/src/lobject.c
			externals/lua/src/lopcodes.c
			externals/lua/src/loslib.c
			externals/lua/src/lparser.c
			externals/lua/src/lstate.c
			externals/lua/src/lstring.c
			externals/lua/src/lstrlib.c
			externals/lua/src/ltable.c
			externals/lua/src/ltablib.c
			externals/lua/src/ltm.c
			externals/lua/src/lundump.c
			externals/lua/src/lvm.c
			externals/lua/src/lzio.c)

set(headersLua	externals/lua/src/lapi.h
				externals/lua/src/lauxlib.h
				externals/lua/src/lcode.h
				externals/lua/src/ldebug.h
				externals/lua/src/ldo.h
				externals/lua/src/lfunc.h
				externals/lua/src/lgc.h
				externals/lua/src/llex.h
				externals/lua/src/llimits.h
				externals/lua/src/lmem.h
				externals/lua/src/lobject.h
				externals/lua/src/lopcodes.h
				externals/lua/src/lparser.h
				externals/lua/src/lstate.h
				externals/lua/src/lstring.h
				externals/lua/src/ltable.h
				externals/lua/src/ltm.h
				externals/lua/src/lua.h
				externals/lua/src/luaconf.h
				externals/lua/src/lualib.h
				externals/lua/src/lundump.h
				externals/lua/src/lvm.h
				externals/lua/src/lzio.h)
else (NOT USE_LUAJIT)
# not need to compile when using luajit library 
set(srcLua )
endif(NOT USE_LUAJIT)

set(srcUGScript		lua_util.cpp
					info_commands.cpp
					bindings_lua.cpp
					lua_debug.cpp
					lua_stack.cpp
					lua_serialization.cpp
					lua_table_handle.cpp)

set(srcLUACompiler	compiler/parser.cpp
					compiler/lexer.cpp
					compiler/lua_parser_class.cpp
					compiler/lua_parser_class_create_c.cpp
					compiler/lua_parser_class_create_jitsg.cpp
					compiler/lua_parser_class_create_lua.cpp
					compiler/lua_parser_class_create_vm.cpp
					compiler/lua_parser_class_reduce.cpp
					compiler/converter.cpp
					compiler/parser.y
					compiler/lexer.l
					compiler/lua_compiler.cpp
					compiler/system_call.cpp)



set(Sources ${srcUGScript} ${srcLua})

# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(buildAlgebra)
	set(Sources ${Sources} lua_user_data.cpp)
endif(buildAlgebra)

if(USE_LUA2C)
	set(Sources ${Sources} ${srcLUACompiler})
endif(USE_LUA2C)


if(BUILD_ONE_LIB)
	message("one lib")
	EXPORTSOURCES(ugbase/bindings/lua  ${Sources})
	add_compile_definitions( USE_MKSTEMP=1)
	add_compile_definitions( LUA_USE_POSIX)
	add_compile_definitions(LUA_TMPNAMTEMPLATE)
else(BUILD_ONE_LIB)
	message("not one lib")
	add_library(bindings_lua  ${Sources})
	target_compile_definitions(bindings_lua PRIVATE USE_MKSTEMP=1)
	target_compile_definitions(bindings_lua PRIVATE LUA_TMPNAMTEMPLATE=1)
endif(BUILD_ONE_LIB)
