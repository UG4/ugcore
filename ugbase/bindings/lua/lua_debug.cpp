/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

//extern libraries
#include <cassert>
#include <cstring>
#include <string>
#include <stack>

// ug libraries
#include "ug.h"

#include "bridge/bridge.h"

#include "common/util/string_util.h"
#include "registry/class_helper.h"
#include "registry/registry.h"

#ifdef USE_LUAJIT
#include <lua.hpp>
#else
#include "externals/lua/src/lua.hpp"
#endif

#ifdef UG_DISC
#include "lua_user_data.h"
#endif

#ifdef UG_PROFILER
#include "common/profiler/runtime_profile_info.h"
#endif


#include "info_commands.h"
#include "lua_debug.h"
#include "lua_util.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace ug {
namespace script {

//////////////////////////////////////////////////////////////////////////////////////////////////////
// globals


static bool bDebugging = false;
static int debugMode  = DEBUG_CONTINUE;
static bool bProfiling = false;
static std::map<std::string, std::map<int, bool> > breakpoints;
static debug_return (*pDebugShell)() = nullptr;
static std::string lastsource;
static int lastline = -1;
static int currentDepth = -1;
#ifdef UG_PROFILER
static bool bStartProfiling=false;
static bool bEndProfiling=false;
static int profilingEndDepth=0;
static std::map<const char*, std::map<int, pRuntimeProfileInfo> >pis;
static std::stack<pRuntimeProfileInfo> pisStack;
#endif

bool IsLUADebug()
{
	return bProfiling;
}
int curHookMask = 0;
//////////////////////////////////////////////////////////////////////////////////////////////////////

extern stack<string> stkPathes;

//////////////////////////////////////////////////////////////////////////////////////////////////////
void LuaCallHook(lua_State *L, lua_Debug *ar);

//////////////////////////////////////////////////////////////////////////////////////////////////////

void FinalizeLUADebug()
{
	breakpoints.clear();
#ifdef UG_PROFILER
	pis.clear();
#endif
	lastsource.clear();
}

int SetDebugShell(debug_return (*s)())
{
	pDebugShell=s;
	return 0;
}

void CheckHook()
{
	if(bDebugging == false && bProfiling == false)
	{
		if(curHookMask)
		{
			lua_sethook (GetDefaultLuaState(), nullptr, 0, 0);
			curHookMask=0;
		}
	}
	else
	{
		constexpr int DEBUG_HOOK_MASK = LUA_MASKCALL | LUA_MASKRET | LUA_MASKLINE;
		constexpr int PROFILER_HOOK_MASK = LUA_MASKCALL | LUA_MASKRET;
		
		if(bDebugging && curHookMask != DEBUG_HOOK_MASK)
		{
			lua_sethook (GetDefaultLuaState(), LuaCallHook, DEBUG_HOOK_MASK, 0);
			curHookMask = DEBUG_HOOK_MASK;
		}
		else if(bProfiling && curHookMask != PROFILER_HOOK_MASK)
		{
			lua_sethook (GetDefaultLuaState(), LuaCallHook, PROFILER_HOOK_MASK, 0);
			curHookMask = PROFILER_HOOK_MASK;
		}
	}

}


void AddBreakpoint(const char*source, const int line)
{
	if(pDebugShell==nullptr)
	{
		UG_LOG("No Debug Shell set!\n");
		return;
	}
	const char *s=nullptr;

	const string relativeFilename=source;
	string absoluteFilename;

	if(GetAbsoluteUGScriptFilename(relativeFilename, absoluteFilename))
	{
		breakpoints[s][line]=true;
		bDebugging = true;
		CheckHook();
		UG_LOG("breakpoint at " << s << ":" << line << "\n")
	}
	else
	{
		UG_LOG("file " << s << " not found\n");
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void PrintBreakpoints()
{
	std::map<std::string, std::map<int, bool> >::iterator it1;
	std::map<int, bool>::iterator it2;
	for(it1 = breakpoints.begin(); it1 != breakpoints.end(); ++it1)
	{
		std::map<int, bool> &m = it1->second;
		for(it2 = m.begin(); it2 != m.end(); ++it2)
		{
			UG_LOG(it1->first << ":" << it2->first << ( it2->second?" enabled":" disabled") << "\n")
		}
	}
}

void DebugHold()
{
	lastsource="";
	debugMode=DEBUG_STEP;
	bDebugging=true;
	CheckHook();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void breakpoint()
{
	if(pDebugShell==nullptr)
	{
		UG_LOG("Breakpoint reached, no Debug Shell set.\n");
		return;
	}

	debug_return r=pDebugShell();
	if(r == DEBUG_EXIT)
		UGForceExit();
	else if(r == DEBUG_CONTINUE)
	{
		debugMode = DEBUG_CONTINUE;
		bDebugging = !breakpoints.empty();
		CheckHook();
	}
	else if(r == DEBUG_NEXT || r == DEBUG_STEP || r == DEBUG_FINISH)
	{
		debugMode=r;
		bDebugging=true;
		CheckHook();
	}
}

int getDepth()
{
	lua_State *L = GetDefaultLuaState();
	int depth=0;
	lua_Debug entry;
	for(int i = 0; lua_getstack(L, i, &entry); i++)
	{
		lua_getinfo(L, "Sln", &entry);
		if(entry.currentline >= 0)
			depth++;
	}
	return depth;
}

void breakpoint_in_script()
{
	lua_Debug entry;
	lua_State *L = GetDefaultLuaState();
	currentDepth =getDepth();
	for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
	{
		lua_getinfo(L, "Sln", &entry);
		if(entry.currentline >= 0)
		{
			lastsource = entry.source+1;
			lastline = entry.currentline;
			breakpoint();
			return;
		}
	}

	breakpoint();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void luaDebug(lua_State *L, const char *source, int line)
{
	if(source == nullptr || line < 0 || pDebugShell==nullptr) return;
	if(source[0]=='@') source++;

	bool bfound=false;
	if(debugMode == DEBUG_NEXT || debugMode == DEBUG_FINISH)
	{
		int d = getDepth();
		if( ((debugMode == DEBUG_NEXT && d <= currentDepth)
				|| (debugMode == DEBUG_FINISH && d < currentDepth))
				&& (lastsource==source && lastline == line) == false)
		{
			lastsource = source;
			lastline = line;
			currentDepth = d;
			bfound =true;
		}
	}
	else if(debugMode == DEBUG_STEP)
	{
		if((lastsource==source && lastline == line) == false)
		{
			lastsource = source;
			lastline = line;
			currentDepth = getDepth();
			bfound =true;
		}
	}


	if(!bfound && !breakpoints.empty())
	{
		lua_Debug entry;
		for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
		{
			lua_getinfo(L, "Sln", &entry);
			if(entry.currentline >= 0)
			{
				std::map<int, bool> &m = breakpoints[entry.source+1];
				auto it = m.find(entry.currentline);
				//UG_LOG(entry.source+1 << ":" << entry.currentline << "\n");
				if(it != m.end() && it->second == true &&
						(lastline != entry.currentline || lastsource!=entry.source+1))
				{
					lastsource = entry.source+1;
					lastline = entry.currentline;
					bfound=true;
					currentDepth = getDepth();
					break;
				}
			}
			if(bfound) {break;}
		}


	}

	if(!bfound) return;

	breakpoint();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void LuaCallHook(lua_State *L, lua_Debug *ar)
{
	if(bDebugging)
	{			
		if(ar->event == LUA_HOOKCALL || ar->event ==LUA_HOOKLINE)
		{
			
			lua_getinfo(L, "Sln", ar);
			auto source = "unknown";
			int line = 0;
			if(ar->currentline < 0)
			{
				lua_Debug entry;
				for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
				{
					lua_getinfo(L, "Sln", &entry);
					if(entry.currentline >= 0)
					{
						source = entry.source;
						line = entry.currentline;
						break;
					}
				}
			}
			else
			{
				source = ar->source;
				line = ar->currentline;
			}
			luaDebug(L, source, line);
		}
	}
#if UG_PROFILER	
	const static bool bDebugLuaProfiler = false;

	// Lua profiling is done the following way:
	// We only want to profile those lines, that actually do something and not
	// all lines, since they may be in inner loops and cause significant overhead.
	// Therefore, we only react on lua-calls and returns. In addition, we don't
	// want to profile lua-callbacks that are invoked from the C++ code, since
	// those are commonly used in the very inner loops - this would be to much
	// profiling overhead.
	//
	// NOTE: we currently use a map with const char* identifier. This will
	//		 only work if lua does not change the strings in between. It this
	//		 is observed we must change that to std::string
	if(bProfiling)
	{	
		if(ar->event == LUA_HOOKLINE || ar->event == LUA_HOOKCOUNT){
			return; // nothing to do: but may be enabled for debugging
		}

		// profileDepthInCcall is used to distinguish where a function call/return
		// is executed. We only want to profile those calls that are done within
		// lua, but not from C++ to Lua. If a call is to a C-function, the
		// depth will be initialized with 1 and no profiling is performed until
		// the depth is back to zero (i.e. when the C-call returns). It may however
		// happen, that several events happen, since lua code is executed from
		// C. In order to detect, if the hook-return is due to a lua-return or
		// from the initalizing C-call, we count the depth.
		static int profileDepthInCcall = 0;

		// check if within a c-call
		if(profileDepthInCcall != 0)
		{		
			// if another call, this must be inside lua but invoked from c
			if(ar->event == LUA_HOOKCALL) {
				profileDepthInCcall++; return;
			}
			// if a return, reduce depth. If depth == 0, we also have to
			// finish the current profile node, since this must be the return
			// from the original c-call
			else if(ar->event == LUA_HOOKRET || ar->event == LUA_HOOKTAILCALL){
				if(--profileDepthInCcall != 0) return;
			}
			else{
				UG_LOG("WARNING: LuaProfiling: Wrong hook event passed (" << ar->event << ")\n");
			}
		}
		
		// if a call is given, but end profiling requested, we can leave
		if(bEndProfiling == true && ar->event == LUA_HOOKCALL)
			return;

		// fill information about the event
		lua_getinfo(L, "Sln", ar);

		// this is the call/ret invoked, however not the called !
		if(bDebugLuaProfiler){
			std::string type = "call";
			if(ar->event == LUA_HOOKRET) type = "ret ";
			else if(ar->event == LUA_HOOKTAILCALL) type = "tailret ";
			UG_LOG(repeat(' ', pisStack.size()) << "## lua profile: source: "<<ar->source<<", line: "
			       <<ar->currentline<<" "<<type<<" " << ar->what);
		}


		// we are only interested in 'Lua' or 'C' call/return events
		if(ar->what[0] != 'L' && ar->what[0] != 'C' && ar->what[0] != 't')
		{
			if(bDebugLuaProfiler){ UG_LOG("\n"); }
			return;
		}

		// we get the debug info of the calling file and line
		lua_Debug entry;
		if(lua_getstack(L, 1, &entry) == 0)
			UG_THROW("LuaProfiling: Cannot get debug info from stack.")
		if(lua_getinfo(L, "Sln", &entry) == 0)
			UG_THROW("LuaProfiling: Cannot read debug info.")

		// get source an line of the event
		const char* source = entry.source;
		int line = entry.currentline;

		if(bDebugLuaProfiler){
			std::string type = "call";
			if(entry.event == LUA_HOOKRET) type = "ret ";
			else if(entry.event == LUA_HOOKTAILCALL) type = "tailret ";
			UG_LOG(",  corr: source: "<<source<<", line: "<<line<<" "<<type<<" " << entry.what << "\n");
		}


		// may still be unavailable, then we ignore this issue
		//if(line < 0) return;

		// a call
		//UG_LOG("ar->event  = " << ar->event  << endl)
		if(ar->event == LUA_HOOKCALL)
		{				
			 // get info from map
			 pRuntimeProfileInfo &pi = pis[source][line];

			 // if not yet initialized, create new node
			 if(pi == nullptr){
				 char buf[1024];
				 if(source[0]=='@') source++;
				 if(strncmp(source, "./../scripts/", 13)==0)
					 sprintf(buf, "!%s:%d ", source+13, line);
				 else
					 sprintf(buf, "@%s:%d ", source, line);

				 pi = new RuntimeProfileInfo(buf, true, "lua", false, source, true, line);
			 }

			 if(bDebugLuaProfiler){UG_LOG(repeat(' ', pisStack.size()) << "PUSH " << *pi << "\n");}
			 // start profiling
			 pi->beginNode();
			 pisStack.push(pi);

			 // if this is a call to C, we disable all successive call-event
			 // tracings until this call returns. This is done by setting
			 // profileDepthInCcall > 1. The only exception is the call to
			 // 'ug_load_script', that is a c-call but should be profiled
			 //	also internally since it just loads other lua-scripts
			 if(ar->name != nullptr && strcmp(ar->name, "ug_load_script") == 0)
				 return;
			 if(ar->what[0] == 'C')
				 	 profileDepthInCcall = 1;
		}

		// a return
		else if(ar->event == LUA_HOOKRET || ar->event == LUA_HOOKTAILCALL)
		{		
			 // if only starting the profiling, do not react on returns
			 if(bStartProfiling) { bStartProfiling = false; return; }

			 // if profiling is to be ended, so something (what?)
			 if(profilingEndDepth > 0) {
				 profilingEndDepth--;
			 }
			 else
			 {
				 // get profile node
				 pRuntimeProfileInfo pi;
				 if(line < 0 || ar->event == LUA_HOOKTAILCALL)
					 pi = pisStack.top();
				 else
				 {
					 pRuntimeProfileInfo p2 = pis[source][line];
					 pi = pisStack.top();
					 if(pi != p2) {UG_LOG("\nPIs not matching.\n"); UG_ASSERT(0, "lua debug bug."); }
				 }
				 pisStack.pop();

				 // release the node
				 // if this part causes trouble, the std::map should be used
				 // with std::string identifier, since I don't know if the
				 // lua const char* are persistent
				 pi->endNode();
				 if(bDebugLuaProfiler){UG_LOG(repeat(' ', pisStack.size()) << "POP  " << *pi << "\n"); }
				 // if we are ending profiling we have to check if we are in
				 // a c call
				 if(bEndProfiling && profileDepthInCcall == 0)
				 {
					 UG_LOG("Profiling ended.\n");
					 bProfiling=false;
					 bEndProfiling=false;
					 CheckHook();
				 }
			 }
		}
		else{
			UG_ASSERT(0, "wrong event already parsed: " << ar->event);
		}
	}
#endif

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void ProfileLUA(bool bProfile)
{
#ifdef UG_PROFILER
	if(bProfiling==false && bProfile==true)
	{	
		bStartProfiling=true;
		bEndProfiling=false;
		bProfiling=true;
		CheckHook();
		
	}
	else if(bProfiling == true && bProfile==false)
	{
		bEndProfiling=true;
	}
#else
	UG_LOG("No profiler available.\n");
#endif
}

void DebugList()
{
	lua_Debug entry;
	int depth=getDepth();
	lua_State *L = GetDefaultLuaState();
	for(int i = 0; lua_getstack(L, i, &entry); i++)
	{
		lua_getinfo(L, "Sln", &entry);
		if(entry.currentline <0) continue;
		if(depth==currentDepth)
		{
			lastsource = entry.source+1;
			lastline = entry.currentline;
			UG_LOG(entry.source+1 << ":" << entry.currentline << "\n");
			if(entry.currentline-3<0)
				{UG_LOG(GetFileLines(entry.source+1, 0, entry.currentline+5, true) << "\n");}
			else
				{UG_LOG(GetFileLines(entry.source+1, entry.currentline-3, entry.currentline+5, true) << "\n");}
			UG_LOG("\n");
		}
		depth--;
	}
}
void DebugBacktrace(int fromLevel)
{
	bridge::LuaStackTrace(fromLevel);
}



void UpdateDepth()
{
	lua_Debug entry;
	int depth=getDepth();
	lua_State *L = GetDefaultLuaState();
	for(int i = 0; lua_getstack(L, i, &entry); i++)
	{
		lua_getinfo(L, "Sln", &entry);
		if(entry.currentline <0) continue;
		if(depth==currentDepth)
		{
			lastsource = entry.source+1;
			lastline = entry.currentline;
			UG_LOG(entry.source+1 << ":" << entry.currentline);
			UG_LOG(" " << GetFileLine(entry.short_src, entry.currentline));
			UG_LOG("\n");
		}
		depth--;
	}
}
void DebugUp()
{
	if(currentDepth>0)
	{
		currentDepth--;
		UpdateDepth();
	}
	else
		{UG_LOG("already at base level.\n");}
}
void DebugDown()
{
	if(currentDepth < getDepth())
	{
		currentDepth++;
		UpdateDepth();
	}
	else
	{UG_LOG("already at max level.\n");}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
bool RegisterLuaDebug(bridge::Registry &reg)
{
	reg.add_function("breakpoint", &AddBreakpoint, "/ug4/lua");
	reg.add_function("breakpoint", &breakpoint_in_script, "/ug4/lua");
	reg.add_function("print_breakpoints", &PrintBreakpoints, "/ug4/lua");
	reg.add_function("DebugBacktrace", &DebugBacktrace, "/ug4/lua");

	reg.add_function("ProfileLUA", &ProfileLUA, "/ug4/lua");
	return true;
}



void SetLuaDebug(lua_State* L, const string &id)
{
	const string& name = id;
	string rest = name;
	string pre;

	while(true)
	{
		const int dotPos = rest.find('.');
		if(dotPos == -1) break;
		string sub = rest.substr(0, dotPos);
		rest = rest.substr(dotPos+1, rest.size());

		//ug::bridge::SetLuaNamespace(pre+sub+".group", pre+sub+".*");
		//string p = std::string("") + "debugID." + pre+sub + " = " + "debugID." + pre+sub + " or {}\n" + "function debugID." + pre+sub + ".set_group_level(level) GetLogAssistant():set_debug_level(\"" + pre+sub + ".*\", level) end\n";
		string p = std::string("") + "debugID." + pre+sub + " = " + "debugID." + pre+sub + " or {}\n" + "debugID." + pre+sub + ".id = \"" + pre+sub + "\"\n"; // todo avoid +
		//UG_LOGN(p);
		ParseAndExecuteBuffer(p.c_str(), "SetLuaDebug");
		pre += sub+".";

	}

	//ug::bridge::SetLuaNamespace(pre+rest+".id", pre+rest);
//	string p = std::string("") + "debugID." + name + " = " + "debugID." + name + " or {}\n" + "function debugID." + name + ".set_level(level) GetLogAssistant():set_debug_level(\"" + name + "\", level) end\n";
	string p = std::string("") + "debugID." + name + " = " + "debugID." + name + " or {}\n" + "debugID." + name + ".id = \"" + name + "\"\n";
	//UG_LOGN(p);
	ParseAndExecuteBuffer(p.c_str(), "SetLuaDebug");
}

void SetLuaDebugIDs(lua_State* L)
{
	ParseAndExecuteBuffer(
			"debugID = {}\n"
			"function SetDebugLevel(did, level)\n"
			"if(did == nil) then\n"
			"print(\"ERROR: Debug Node not existing. Perhaps you did not include the plugin?\")\n"
			"DebugBacktrace(0)"
			"else GetLogAssistant():set_debug_level((did.id) ..\"*\", level) end end",
			"SetLuaDebugIDs");
	const vector<string> &s = DebugIDManager::instance().get_registered_debug_IDs_arr();
	for(size_t i=0; i<s.size(); i++)
		SetLuaDebug(L, s[i]);
}
}
}