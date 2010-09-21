// created by Andreas Vogel

#include "luabind/luabind.hpp"
#include "ug_script.h"
#include "lib_discretization/lib_discretization.h"


using namespace std;
using namespace ug;

#define CHECK_POINTER(grid, retVal) {if(!grid){LOG("invalid pointer\n"); return retVal;}}

namespace ldscript
{
/*
P1ConformDoFManager<MultiGrid>* new_P1ConformDoFManager(MultiGrid* sh)
{
	P1ConformDoFManager<MultiGrid>* manager = new P1ConformDoFManager<MultiGrid>(*sh);
	LOG("  P1ConformDoFManager for MultiGrid created\n");
	return manager;
}
*/
}//	end of namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	register the script-module at the module manager.
//	do this in this empty namespace
//	please note that this trick will not work if you link this file
//	into any form of library. Only if it is compiled and directly linked
//	to the executable, this code will actually be executed.
//	This is due the automatic removal of unreferenced code by the linker.
namespace ug{
namespace script
{
///	initialises lib-grids script methods
bool InitLibDiscScript(lua_State* L)
{
	using namespace ldscript;
	using namespace luabind;
/*
	module(L)[
//	bind classes
		class_<P1ConformDoFManager<MultiGrid> >("P1ConformDoFManagerMultiGrid"),

//	bind methods to lua
		def("new_P1ConformDoFManager", new_P1ConformDoFManager)

	];
*/
	return true;
}

void FinalizeLibDiscScript(lua_State* L)
{
}

}}//	end of namespace

