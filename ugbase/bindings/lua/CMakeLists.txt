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

cmake_minimum_required(VERSION 2.6)

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
set(srcLua	externals/lua/lapi.c
			externals/lua/lauxlib.c
			externals/lua/lbaselib.c
			externals/lua/lcode.c
			externals/lua/ldblib.c
			externals/lua/ldebug.c
			externals/lua/ldo.c
			externals/lua/ldump.c
			externals/lua/lfunc.c
			externals/lua/lgc.c
			externals/lua/linit.c
			externals/lua/liolib.c
			externals/lua/llex.c
			externals/lua/lmathlib.c
			externals/lua/lmem.c
			externals/lua/loadlib.c
			externals/lua/lobject.c
			externals/lua/lopcodes.c
			externals/lua/loslib.c
			externals/lua/lparser.c
			externals/lua/lstate.c
			externals/lua/lstring.c
			externals/lua/lstrlib.c
			externals/lua/ltable.c
			externals/lua/ltablib.c
			externals/lua/ltm.c
			externals/lua/lundump.c
			externals/lua/lvm.c
			externals/lua/lzio.c)

set(headersLua	externals/lua/lapi.h
				externals/lua/lauxlib.h
				externals/lua/lcode.h
				externals/lua/ldebug.h
				externals/lua/ldo.h
				externals/lua/lfunc.h
				externals/lua/lgc.h
				externals/lua/llex.h
				externals/lua/llimits.h
				externals/lua/lmem.h
				externals/lua/lobject.h
				externals/lua/lopcodes.h
				externals/lua/lparser.h
				externals/lua/lstate.h
				externals/lua/lstring.h
				externals/lua/ltable.h
				externals/lua/ltm.h
				externals/lua/lua.h
				externals/lua/luaconf.h
				externals/lua/lualib.h
				externals/lua/lundump.h
				externals/lua/lvm.h
				externals/lua/lzio.h)
else (NOT USE_LUAJIT)
# not need to compile when using luajit library 
set(srcLua )
endif(NOT USE_LUAJIT)

#find_package(BISON)
#find_package(FLEX)
#BISON_TARGET(MyParser parser.y parser.cc)
#FLEX_TARGET(MyScanner lexer.l  lexer.cc)
#ADD_FLEX_BISON_DEPENDENCY(MyScanner MyParser)
# Create target for the parser
#ADD_CUSTOM_TARGET(FooParser echo “Creating parser.c”)
# Create custom command for flex/lex (note the outputs)



set(srcUGScript		lua_util.cpp
					info_commands.cpp
					bindings_lua.cpp
					lua_debug.cpp
					lua_stack.cpp
					lua_serialization.cpp)

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
	EXPORTSOURCES(ugbase/bindings/lua  ${Sources})
else(BUILD_ONE_LIB)
	add_library(bindings_lua  ${Sources})
endif(BUILD_ONE_LIB)
