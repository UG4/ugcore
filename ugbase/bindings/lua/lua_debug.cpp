// created by Martin Rupp
// martin.rupp@gcsc.uni-frankfurt.de
// y12 m03 d17

//////////////////////////////////////////////////////////////////////////////////////////////////////
//extern libraries
#include <cassert>
#include <cstring>
#include <string>
#include <stack>

//////////////////////////////////////////////////////////////////////////////////////////////////////
// ug libraries
#include "ug.h"
#include "lua_util.h"
#include "common/util/file_util.h"
#include "bindings_lua.h"
#include "bridge/bridge.h"
#include "registry/class_helper.h"
#include "lua_user_data.h"
#include "registry/registry.h"
#include "info_commands.h"
#include "lua_debug.h"
#include "common/profiler/runtime_profile_info.h"
#include "common/util/string_util.h"

extern "C" {
#include "externals/lua/lstate.h"
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace ug
{

namespace script
{

//////////////////////////////////////////////////////////////////////////////////////////////////////
// globals


static bool bDebugging = false;
static int debugMode  = DEBUG_CONTINUE;
static bool bProfiling = false;
static std::map<std::string, std::map<int, bool> > breakpoints;
static debug_return (*pDebugShell)() = NULL;
static std::string lastsource;
static int lastline = -1;
static int currentDepth = -1;
#ifdef UG_PROFILER
static bool bStartProfiling=false;
static bool bEndProfiling=false;
static bool bProfileLUALines=true;
static int profilingDepth=0;
static int profilingEndDepth=0;
static std::map<const char*, std::map<int, pRuntimeProfileInfo> >pis;
#endif


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
			lua_sethook (GetDefaultLuaState(), NULL, 0, 0);
			curHookMask=0;
		}
	}
	else
	{
		const int DEBUG_HOOK_MASK = LUA_MASKCALL | LUA_MASKRET | LUA_MASKLINE;
		const int PROFILER_HOOK_MASK = LUA_MASKCALL | LUA_MASKRET;
		
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


void AddBreakpoint(const char*source, int line)
{
	if(pDebugShell==NULL)
	{
		UG_LOG("No Debug Shell set!\n");
		return;
	}
	const char *s=NULL;

	string relativeFilename=source;
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
		std::map<int, bool> &m = (*it1).second;
		for(it2 = m.begin(); it2 != m.end(); ++it2)
		{
			UG_LOG((*it1).first << ":" << (*it2).first << ((*it2).second?" enabled":" disabled") << "\n")
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void breakpoint()
{
	if(pDebugShell==NULL)
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
		bDebugging = breakpoints.size() > 0;
		CheckHook();
		return;
	}
	else if(r == DEBUG_NEXT || r == DEBUG_STEP || r == DEBUG_FINISH)
	{
		debugMode=r;
		bDebugging=true;
		CheckHook();
		return;
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
	if(source == NULL || line < 0 || pDebugShell==NULL) return;
	if(source[0]=='@') source++;

	bool bfound=false;
	if(debugMode == DEBUG_NEXT || debugMode == DEBUG_FINISH)
	{
		int d = getDepth();
		if( ((debugMode == DEBUG_NEXT && d <= currentDepth)
				|| (debugMode == DEBUG_FINISH && d < currentDepth))
				&& (lastsource.compare(source)==0 && lastline == line) == false)
		{
			lastsource = source;
			lastline = line;
			currentDepth = d;
			bfound =true;
		}
	}
	else if(debugMode == DEBUG_STEP)
	{
		if((lastsource.compare(source)==0 && lastline == line) == false)
		{
			lastsource = source;
			lastline = line;
			currentDepth = getDepth();
			bfound =true;
		}
	}


	if(!bfound && breakpoints.size() > 0)
	{
		lua_Debug entry;
		for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
		{
			lua_getinfo(L, "Sln", &entry);
			if(entry.currentline >= 0)
			{
				std::map<int, bool> &m = breakpoints[entry.source+1];
				std::map<int, bool>::iterator it = m.find(entry.currentline);
				//UG_LOG(entry.source+1 << ":" << entry.currentline << "\n");
				if(it != m.end() && (*it).second == true &&
						(lastline != entry.currentline || lastsource.compare(entry.source+1)!=0))
				{
					lastsource = entry.source+1;
					lastline = entry.currentline;
					bfound=true;
					currentDepth = getDepth();
					break;
				}
			}
			if(bfound) break;
		}


	}

	if(!bfound) return;

	breakpoint();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void LuaCallHook(lua_State *L, lua_Debug *ar)
{
#if 0
	UG_LOG("depth = " << profilingDepth << "\n");
	{
	UG_LOG("------------------------\n");
	{
		if(ar->event == LUA_HOOKCALL) UG_LOG("HOOKCALL\n");
		if(ar->event == LUA_HOOKLINE) UG_LOG("HOOKLINE\n");
		if(ar->event == LUA_HOOKRET) UG_LOG("HOOKRET\n");
	    lua_Debug entry;
	    for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
		{
	    	int status = lua_getinfo(L, "Sln", &entry);
	    	//if(entry.currentline <0) continue;
	    	if(entry.short_src && entry.currentline>0)
	    		UG_LOG(entry.short_src << ":" << entry.currentline);
	    	if(entry.what)
	    	{
	    		UG_LOG(" what=" << entry.what);
	    	}
	    	if(entry.name)
	    	{
	    		UG_LOG(" entry.name=" << entry.name);
	    	}
	    	UG_LOG("\n");
	    }
	}
	}
#endif
	if(bDebugging)
	{			
		if(ar->event == LUA_HOOKCALL || ar->event ==LUA_HOOKLINE)
		{
			
			lua_getinfo(L, "Sln", ar);
			const char *source = "unknown";
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
	if(bProfiling)
	{	
		if(profilingDepth != 0)
		{		
			if(ar->event == LUA_HOOKCALL) { profilingDepth++; return; }
			if(ar->event == LUA_HOOKRET)
			{
				if(--profilingDepth != 0)
					return;		
			}	
		}
		
		
		if(bEndProfiling==false && ar->event == LUA_HOOKCALL)
		{				
			lua_getinfo(L, "Sln", ar);
			const char *source = ar->source;
			int line = ar->currentline;

			if(ar->what[0] == 'L' || ar->what[0] == 'C')
			{
				lua_Debug entry;
				if(line < 0 && bProfileLUALines && lua_getstack(L, 1, &entry))
				{
					lua_getinfo(L, "Sl", &entry);
					source = entry.source;
					line = entry.currentline;
				}

				if(line >= 0)
				{
					// be sure that this is const char*
					//if(line>0) line--;
					pRuntimeProfileInfo &pi = pis[source][line];
					//UG_LOG("start profile node " << source << ":" << line << "\n");

					// if null, create new node
					if(pi == NULL)
					{
						char buf[1024] = "LUAunknown ";
						if(source[0]=='@') source++;
						if(strncmp(source, "./../scripts/", 13)==0)
							sprintf(buf, "!%s:%d ", source+13, line);
						else
							sprintf(buf, "@%s:%d ", source, line);
						//const char*p = GetFileLine(source, line).c_str();
						//strncat(buf, p+strspn(p, " \t"), 254);

						pis[source][line] = new RuntimeProfileInfo(buf, true, "lua", false, source, true, line);
						// UG_LOG(buf);
					 }

					 pRuntimeProfileInfo &pi2 = pis[source][line];
					 pi2->beginNode();
				}
			}
			if(ar->what[0] == 'C' && (ar->name == NULL || strcmp(ar->name, "ug_load_script") != 0))
				profilingDepth=1;
		}
		else if(ar->event == LUA_HOOKRET)
		{		
			if(bStartProfiling) { bStartProfiling = false; return; }
			lua_getinfo(L, "Sln", ar);
			if(ar->what[0] == 'L' || ar->what[0] == 'C')
			{
				/*int line = ar->currentline;
				lua_Debug entry;
				if(line < 0 && bProfileLUALines && lua_getstack(L, 1, &entry))
				{
					lua_getinfo(L, "Sln", &entry);
					line = entry.currentline;
				}*/
				if(profilingEndDepth>0)
					profilingEndDepth--;
				//UG_ASSERT(pis[ar->source][ar->linedefined]->isCurNode(), "profiler nodes not matching. forgot a PROFILE_END?");
				else
				{		
					const char *source = ar->source;
					int line = ar->currentline;
					pRuntimeProfileInfo &pi = pis[source][line];

					//UG_LOG("end profile node\n");
					pi->endNode();
					if(bEndProfiling && profilingDepth==0)
					{
						UG_LOG("Profiling ended.\n");
						bProfiling=false;
						bEndProfiling=false;
						CheckHook();
					}
				}		
			}
		}
	}
#endif

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void ProfileLUA(bool b)
{
#ifdef UG_PROFILER
	if(bProfiling==false && b==true)
	{	
		bStartProfiling=true;
		bEndProfiling=false;
		bProfiling=true;
		CheckHook();
		
	}
	else if(bProfiling == true && b==false)
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
	// todo
}
void DebugBacktrace()
{
	ug::bridge::LuaStackTrace(GetDefaultLuaState());
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
bool RegisterLuaDebug(ug::bridge::Registry &reg)
{
	reg.add_function("breakpoint", &AddBreakpoint, "/ug4/lua");
	reg.add_function("breakpoint", &breakpoint_in_script, "/ug4/lua");
	reg.add_function("print_breakpoints", &PrintBreakpoints, "/ug4/lua");
	reg.add_function("DebugBacktrace", &DebugBacktrace, "/ug4/lua");

	reg.add_function("ProfileLUA", &ProfileLUA, "/ug4/lua");
	return true;
}


}
}
