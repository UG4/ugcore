cmake_minimum_required(VERSION 2.6)

####
# uginterface Library
####

project(P_UGSCRIPT)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(${P_UGSCRIPT_SOURCE_DIR}/externals/lua)

include("../../../ug_cmake_includes.txt")

if(UNIX)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLUA_USE_POSIX")
endif(UNIX)

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
					lua_debug.cpp)

set(srcLUA2C		converter/parser.cpp
					converter/lexer.cpp
					converter/pclass.cpp
					converter/converter.cpp
					converter/parser.y
					converter/lexer.l
					converter/parser.cpp
					converter/lexer.cpp
					converter/lua2c.cpp)


# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(buildAlgebra)
	# add the source files
	set(srcUGScript ${srcUGScript}
					lua_user_data.cpp)
else(buildAlgebra)

endif(buildAlgebra)

if(USE_LUA2C)
	set(srcUGScript ${srcUGScript}
		${srcLUA2C})
endif()


set(Sources ${srcUGScript} ${srcLua})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/bindings/lua  ${Sources})
else(BUILD_ONE_LIB)
	add_library(bindings_lua  ${Sources})
endif(BUILD_ONE_LIB)
