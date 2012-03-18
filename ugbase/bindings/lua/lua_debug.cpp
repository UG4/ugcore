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
#include "common/os_dependent/file_util.h"
#include "bindings_lua.h"
#include "bridge/bridge.h"
#include "registry/class_helper.h"
#include "info_commands.h"
#include "lua_user_data.h"
#include "registry/registry.h"

#include "lua_debug.h"

extern "C" {
#include "externals/lua/lstate.h"
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace ug
{

namespace bridge
{
string GetFileLines(const char *filename, size_t fromline, size_t toline, bool includeLineNumbers=false);
}
namespace script
{

//////////////////////////////////////////////////////////////////////////////////////////////////////
// globals


static bool bDebugging = false;
static int debugMode  = DEBUG_CONTINUE;
static bool bProfiling = false;
static std::map<std::string, std::map<int, bool> > breakpoints;
static debug_return (*pDebugShell)() = NULL;

//////////////////////////////////////////////////////////////////////////////////////////////////////

extern stack<string> stkPathes;

//////////////////////////////////////////////////////////////////////////////////////////////////////
int SetDebugShell(debug_return (*s)())
{
	pDebugShell=s;
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
void LuaCallHook(lua_State *L, lua_Debug *ar);

//////////////////////////////////////////////////////////////////////////////////////////////////////
inline void SetDebugHook()
{
	static bool hookset = false;
	if(!hookset)
	{
		lua_sethook (GetDefaultLuaState(), LuaCallHook, LUA_MASKCALL | LUA_MASKRET | LUA_MASKLINE, 0);
		hookset=true;
	}
}

// Lua Profiling (mrupp)
#ifdef PROFILE_BRIDGE
struct s_profileInformation
{
	s_profileInformation()
	{
		Shiny::ProfileZone pi = {NULL, Shiny::ProfileZone::STATE_HIDDEN, NULL, { { 0, 0 }, { 0, 0 }, { 0, 0 } }};
		profileInformation = pi;
		profilerCache =	&Shiny::ProfileNode::_dummy;
	}
	Shiny::ProfileZone profileInformation;
	Shiny::ProfileNodeCache profilerCache;
	char m_name[255];

	bool is_initialised()
	{
		return profileInformation.name != NULL;
	}

	void init(const char*name, char id)
	{
		if(id)
		{
			m_name[0]=id;
			strcpy(m_name+1, name);
		}
		else
			strcpy(m_name, name);
		profileInformation.name = m_name;
	}
};

typedef s_profileInformation* ps_profileInformation ;
#endif

void AddBreakpoint(const char*source, int line)
{
	if(pDebugShell==NULL)
	{
		UG_LOG("No Debug Shell set!\n");
		return;
	}
	string file;
	const char *s=NULL;
	if(!stkPathes.empty())
	{
		file =stkPathes.top();
		file.append("/").append(source);
		if(FileExists(file.c_str()))
			s = file.c_str();
	}
	if(FileExists(source)) s = source;
	if(s)
	{
		SetDebugHook();
		bDebugging=true;
		UG_LOG("breakpoint at " << s << ":" << line << "\n")
		breakpoints[s][line]=true;
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

static int level=0;
static int numLevels=0;

static std::string lastsource;
static 	int lastline = -1;
static int currentDepth = -1;

//////////////////////////////////////////////////////////////////////////////////////////////////////
void breakpoint()
{
	if(pDebugShell==NULL)
	{
		UG_LOG("Breakpoint reached, no Debug Shell set.\n");
		return;
	}
	debug_return r=pDebugShell();
	if(r == DEBUG_EXIT || r == DEBUG_CONTINUE)
	{
		debugMode = DEBUG_CONTINUE;
		return;
	}
	else if(r == DEBUG_NEXT || r == DEBUG_STEP)
	{
		bDebugging=true;
		debugMode=r;
		SetDebugHook();
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

	/*UG_LOG("------------------------\n");
	{
	lua_Debug entry;
	for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
		{
	    	int status = lua_getinfo(L, "Sln", &entry);
	    	if(entry.currentline <0) continue;
	    	if(!status || !entry.short_src || entry.currentline < 0) break;
			UG_LOG(depth << ": " << entry.short_src << ":" << entry.currentline << " " << entry.what);
			//UG_LOG(" " << GetFileLine(entry.short_src, entry.currentline));
			UG_LOG("\n");
	    }
	}/**/
	//fill up the debug structure with information from the lua stack
	lua_getinfo(L, "Sln", ar);
	//if(ar->what[0] == 'L' || ar->what[0] == 'C')
	{
		if(ar->event == LUA_HOOKCALL || (bDebugging && ar->event ==LUA_HOOKLINE))
		{
			lua_Debug entry;
			const char *source = "unknown";
			int line = 0;
			bool found=false;
			if(ar->currentline < 0)
			{
				for(int depth = 0; lua_getstack(L, depth, &entry); depth++)
				{
					lua_getinfo(L, "Sln", &entry);
					if(entry.currentline >= 0)
					{
						source = entry.source;
						line = entry.currentline;
						found = true;
						break;
					}
				}
			}
			else
			{
				source = ar->source;
				line = ar->currentline;
			}

			if(bDebugging)	luaDebug(L, source, line);
			if(bProfiling)
			if(ar->what[0] == 'L' || ar->what[0] == 'C')
			{
				static std::map<const char *, std::map<int, ps_profileInformation> >pis;

				if(ar->event == LUA_HOOKCALL)
				{
					ps_profileInformation &pi = pis[source][line];
					if(pi == NULL)
					{
						pi = new s_profileInformation;

						if(found)
						{
							char buf[255];
							sprintf(buf, "%s:%d", source, line);
							pi->init(buf, 0);
						}
						else
							pi->init("LUAunknown", '@');


					 }
					 Shiny::ProfileManager::instance._beginNode(&pi->profilerCache, &pi->profileInformation);
				}
			}
		}
		else if(ar->event == LUA_HOOKRET)
		{
			if(bProfiling)
			{
				//UG_ASSERT(Shiny::ProfileManager::instance._curNode == pis[ar->source][ar->linedefined]->profilerCache, "profiler nodes not matching. forgot a PROFILE_END?");
				Shiny::ProfileManager::instance._endCurNode();
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void ProfileLUA()
{
	SetDebugHook();
	bProfiling=true;
}



void DebugList()
{
	lua_Debug entry;
	int depth=getDepth();
	lua_State *L = GetDefaultLuaState();
	for(int i = 0; lua_getstack(L, i, &entry); i++)
	{
		int status = lua_getinfo(L, "Sln", &entry);
		if(entry.currentline <0) continue;
		if(depth==currentDepth)
		{
			lastsource = entry.source+1;
			lastline = entry.currentline;
			UG_LOG(entry.source+1 << ":" << entry.currentline << "\n");
			if(entry.currentline-3<0)
				{UG_LOG(ug::bridge::GetFileLines(entry.source+1, 0, entry.currentline+5, true) << "\n");}
			else
				{UG_LOG(ug::bridge::GetFileLines(entry.source+1, entry.currentline-3, entry.currentline+5, true) << "\n");}
			UG_LOG("\n");
		}
		depth--;
	}
	// todo
}
void DebugBacktrace()
{
	ug::bridge::lua_stacktrace(GetDefaultLuaState());
}



void UpdateDepth()
{
	lua_Debug entry;
	int depth=getDepth();
	lua_State *L = GetDefaultLuaState();
	for(int i = 0; lua_getstack(L, i, &entry); i++)
	{
		int status = lua_getinfo(L, "Sln", &entry);
		if(entry.currentline <0) continue;
		if(depth==currentDepth)
		{
			lastsource = entry.source+1;
			lastline = entry.currentline;
			UG_LOG(entry.source+1 << ":" << entry.currentline);
			UG_LOG(" " << ug::bridge::GetFileLine(entry.short_src, entry.currentline));
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

	reg.add_function("ProfileLUA", &ProfileLUA, "/ug4/lua");
	return true;
}


}
}
