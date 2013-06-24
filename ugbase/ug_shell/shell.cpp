/*
 * debug_shell.cpp
 *
 *  Created on: 21.06.2013
 *      Author: mrupp
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
#include "pcl/pcl_util.h"
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
//#define UG_PROMPT "\e[1mug:>\e[0m"

#define UG_PROMPT "ug:> "

using namespace std;

namespace ug{

using namespace script;
namespace bridge
{

#if defined(UG_USE_READLINE)
	string ug_readline(const char *prompt)
	{
		string ret;
		char *str=readline(prompt);
		if(str==NULL)
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
		string ug_readline(const char *prompt)
		{
			string ret;
			char *str=linenoise(prompt);
			if(str==NULL)
				ret="";
			else ret=str;
			free(str);
			return ret;
		}
		void ug_cacheline(string str)
		{
			linenoiseHistoryAdd(str.c_str());
		}
	#else
		string ug_readline(const char* prompt)
		{
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
	if(prompt == NULL) prompt=UG_PROMPT;

	//	run the shell
	const char *completitions[] ={"quit", "exit"};
	while(1)
	{
#if defined(UG_USE_LINENOISE)
		SetOtherCompletions(completitions, sizeof(completitions)/sizeof(completitions[0]));
#endif

		PROFILE_BEGIN(ug_readline);
			string buf = ug_readline(prompt);
		PROFILE_END();

#if defined(UG_USE_LINENOISE)
		SetOtherCompletions(NULL, 0);
#endif
		size_t len = buf.length();
		if(len)
		{
			if(buf.compare("exit")==0 || buf.compare("quit")==0)
				break;

			if(len > 6 && strncmp(buf.c_str(), "print ", 6)==0 && buf[6] != '(')
			{
				bridge::UGTypeInfo(buf.c_str()+6);
				continue;
			}
			if(buf[len-1]=='?')
			{
				buf.resize(len-1);
				bridge::UGTypeInfo(buf.c_str());
				continue;
			}

			ug_cacheline(buf);

			try
			{
				script::ParseBuffer(buf.c_str(), "interactive shell");
			}
			catch(LuaError& err)
			{
				PathProvider::clear_current_path_stack();
				if(!err.get_msg().empty()){
					UG_LOG("LUA-ERROR: \n");
					for(size_t i=0;i<err.num_msg();++i)
						UG_LOG(err.get_msg(i)<<endl);
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
	static string last="";
	ug::bridge::LuaStackTrace(GetDefaultLuaState());
	//ug::bridge::LuaPrintCurrentLine(GetDefaultLuaState());

#ifdef UG_PARALLEL
	if(pcl::GetNumProcesses() > 1)
	{
		UG_LOG("Parallel Shell not available currently.");
		return DEBUG_CONTINUE;
	}
#endif

	//	run the shell
	const char *completitions[]={"quit", "exit", "step", "next", "cont", "continue", "finish", "list", "backtrace", "bt",
			"up", "down"};
	while(1)
	{
#if defined(UG_USE_LINENOISE)
		SetOtherCompletions(completitions, sizeof(completitions)/sizeof(completitions[0]));
#endif
		PROFILE_BEGIN(ug_readline);
			string buf = ug_readline("debug:> ");
		PROFILE_END();

#if defined(UG_USE_LINENOISE)
		SetOtherCompletions(NULL, 0);
#endif
		size_t len = buf.length();
		if(len)
		{
			ug_cacheline(buf);
			last = buf;
		}
		else
		{
			if(last.length() > 0)
			{
				buf = last;
				UG_LOG("debug:> " << last << "\n");
			}
		}

		if(buf.compare("exit")==0 || buf.compare("quit")==0)
			return DEBUG_EXIT;
		else if(buf.compare("continue")==0 || buf.compare("cont")==0)
			return DEBUG_CONTINUE;
		else if(buf.compare("next")==0)
			return DEBUG_NEXT;
		else if(buf.compare("step")==0)
			return DEBUG_STEP;
		else if(buf.compare("finish")==0)
			return DEBUG_FINISH;
		else if(buf.compare("list")==0)
		{
			DebugList();
			continue;
		}
		else if(buf.compare("backtrace")==0 || buf.compare("bt")==0)
		{
			DebugBacktrace();
			continue;
		}
		else if(buf.compare("up")==0)
		{
			DebugUp();
			continue;
		}
		else if(buf.compare("down")==0)
		{
			DebugDown();
			continue;
		}
		else if(len > 6 && strncmp(buf.c_str(), "print ", 6)==0 && buf[6] != '(')
		{
			bridge::UGTypeInfo(buf.c_str()+6);
			continue;
		}

		if(len)
		{
			if(buf[len-1]=='?')
			{
				buf.resize(len-1);
				bridge::UGTypeInfo(buf.c_str());
				continue;
			}

			try
			{
				script::ParseBuffer(buf.c_str(), "debug shell");
			}
			catch(LuaError& err)
			{
				PathProvider::clear_current_path_stack();
				if(!err.get_msg().empty()){
					UG_LOG("LUA-ERROR: \n");
					for(size_t i=0;i<err.num_msg();++i)
						UG_LOG(err.get_msg(i)<<endl);
				}
			}
		}
	}
//todo:	clear the history (add ug_freelinecache)
	return DEBUG_EXIT;
}

}
}
