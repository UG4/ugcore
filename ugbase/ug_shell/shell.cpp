/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#include <string>
#include <iostream>
#include <sstream>

#include "common/util/path_provider.h"
#include "common/profiler/profiler.h"
#include "bindings/lua/info_commands.h"
#include "bindings/lua/lua_util.h"
#include "shell.h"
#include "completion.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
//ø #include "pcl/pcl_util.h"
#endif

#if defined(UG_USE_READLINE)
	#include <stdio.h>
	#include <readline/readline.h>
	#include <readline/history.h>
#else
	#if defined(UG_USE_LINENOISE)
		#include "externals/linenoise/linenoise.h"
	#endif
#endif

////////////////////////////////////////////////////////////////////////////////
// interactive shells

#define UG_PROMPT "ug:> "
#define UG_DEBUG_PROMPT "debug:> "

using namespace std;

namespace ug{


using namespace script;
namespace bridge
{

#if defined(UG_USE_READLINE)
	string ug_readline(const char *prompt, bool &quit)
	{
		quit=false;
		string ret;
		char *str=readline(prompt);
		if(str==nullptr)
			ret="";
		else ret=str;
		free(str);
		return ret;
	}
	void ug_cacheline(string str)
	{
		add_history(str.c_str());
	}
#else
	#if defined(UG_USE_LINENOISE)
		string ug_readline(const char *prompt, bool &quit)
		{
			quit=false;
			string ret;
			char *str=linenoise(prompt);
			if(str==nullptr)
			{
				quit=true;
				ret="";
			}
			else ret=str;
			free(str);
			return ret;
		}
		void ug_cacheline(string str)
		{
			linenoiseHistoryAdd(str.c_str());
		}
	#else
		string ug_readline(const char* prompt, bool &quit)
		{
			quit=false;
			cout << prompt;
			string strBuffer;
			getline(cin, strBuffer);
			return strBuffer;
		}
		void ug_cacheline(string str) {}
	#endif
#endif

////////////////////////////////////////////////////////////////////////////////
// normal shell
void InitShell()
{
#if defined(UG_USE_LINENOISE)
	linenoiseSetCompletionFunction(CompletionFunction);
	SetDebugShell(DebugShell);
#endif
}

////////////////////////////////////////////////////////////////////////////////
// normal shell
int RunShell(const char *prompt)
{
	if(prompt == nullptr) prompt=UG_PROMPT;

	//	run the shell
	const char *completions[] ={"quit", "exit"};
	while(true)
	{
#if defined(UG_USE_LINENOISE)
		SetOtherCompletions(completions, sizeof(completions)/sizeof(completions[0])); // ø todo std::size ?
#endif

		bool quit;
		string buf;
		{
			PROFILE_BEGIN(ug_readline);
			buf = ug_readline(prompt, quit);
		}


#if defined(UG_USE_LINENOISE)
		SetOtherCompletions(nullptr, 0);
#endif
		if(quit) break;
		size_t len = buf.length();
		if(len)
		{
			if(buf=="exit" || buf=="quit")
				break;

			ug_cacheline(buf);

			if(buf[len-1]=='?')
			{
				buf.resize(len-1);
				UGTypeInfo(buf.c_str());
				continue;
			}
			if(len > 6 && strncmp(buf.c_str(), "print ", 6)==0 && buf[6] != '(')
			{
				ParseAndExecuteBuffer((std::string("PrintTemporaryObject=")+buf.substr(6)).c_str(), "debug shell");
				UGTypeInfo("PrintTemporaryObject");
				continue;
			}

			try
			{
				ParseAndExecuteBuffer(buf.c_str(), "interactive shell");
			}
			catch(LuaError& err)
			{
				PathProvider::clear_current_path_stack();
				if(err.show_msg()){
					if(!err.get_msg().empty()){
						UG_LOG("LUA-ERROR: \n");
						for(size_t i=0;i<err.num_msg();++i)
							UG_LOG(err.get_msg(i)<<endl);
					}
				}
			}
			catch(UGError &err)
			{
				PathProvider::clear_current_path_stack();
				UG_LOG("UGError:\n");
				for(size_t i=0; i<err.num_msg(); i++)
					UG_LOG(err.get_file(i) << ":" << err.get_line(i) << " : " << err.get_msg(i));
				UG_LOG("\n");
			}
		}
	}
//todo:	clear the history (add ug_freelinecache)
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
// debug shell
debug_return DebugShell()
{

	static string last;

	LuaStackTrace(0);
	//ug::bridge::LuaPrintCurrentLine(GetDefaultLuaState());

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
	{
		UG_LOG("Parallel Shell not available currently.");
		return DEBUG_CONTINUE;
	}
#endif

	//	run the shell
	const char *completions[]={"quit", "exit", "step", "next", "cont", "continue", "finish", "list", "backtrace", "bt",
			"up", "down"};
	while(true)
	{
#if defined(UG_USE_LINENOISE)
		SetOtherCompletions(completions, sizeof(completions)/sizeof(completions[0])); // ø todo std::size ?
#endif
		bool quit;
		string buf;
		{
			PROFILE_BEGIN(ug_readline);
			buf = ug_readline(UG_DEBUG_PROMPT, quit);
		}


#if defined(UG_USE_LINENOISE)
		SetOtherCompletions(nullptr, 0);
#endif
		if(quit) return DEBUG_EXIT;
		size_t len = buf.length();
		if(len)
		{
			ug_cacheline(buf);
			last = buf;
		}
		else
		{
			if(!last.empty())
			{
				buf = last;
				UG_LOG("debug:> " << last << "\n");
			}
		}

		if(buf=="exit" || buf=="quit")
			return DEBUG_EXIT;
		else if(buf=="continue" || buf=="cont")
			return DEBUG_CONTINUE;
		else if(buf=="next")
			return DEBUG_NEXT;
		else if(buf=="step")
			return DEBUG_STEP;
		else if(buf=="finish")
			return DEBUG_FINISH;
		else if(buf=="list")
		{
			DebugList();
			continue;
		}
		else if(buf=="backtrace" || buf=="bt")
		{
			DebugBacktrace(0);
			continue;
		}
		else if(buf=="up")
		{
			DebugUp();
			continue;
		}
		else if(buf=="down")
		{
			DebugDown();
			continue;
		}
		else if(len > 6 && strncmp(buf.c_str(), "print ", 6)==0 && buf[6] != '(')
		{
			ParseAndExecuteBuffer((std::string("PrintTemporaryObject=")+buf.substr(6)).c_str(), "debug shell");
			UGTypeInfo("PrintTemporaryObject");
			continue;
		}

		if(len)
		{
			if(buf[len-1]=='?')
			{
				buf.resize(len-1);
				UGTypeInfo(buf.c_str());
				continue;
			}

			try
			{
				ParseAndExecuteBuffer(buf.c_str(), "debug shell");
			}
			catch(LuaError& err)
			{
				PathProvider::clear_current_path_stack();
				if(err.show_msg()){
					if(!err.get_msg().empty()){
						UG_LOG("LUA-ERROR: \n");
						for(size_t i=0;i<err.num_msg();++i)
							UG_LOG(err.get_msg(i)<<endl);
					}
				}
			}
		}
	}
//todo:	clear the history (add ug_freelinecache)
	return DEBUG_EXIT;
}

}
}
