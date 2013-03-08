// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d19

#include <iostream>
#include <sstream>

#include "ug.h"
#include "bindings/lua/lua_util.h"
#include "bindings/lua/lua_debug.h"
#include "bridge/bridge.h"
#include "bindings/lua/info_commands.h"
#include "common/util/parameter_parsing.h"
#include "common/util/file_util.h"

#include "common/util/os_info.h"
#include "common/util/path_provider.h"
#ifdef UG_PLUGINS
	#include "common/util/plugin_util.h"
	#ifdef UG_EMBEDDED_PLUGINS
		#include "embedded_plugins.h"
	#endif
#endif

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif


#include "compile_info/compile_info.h"


using namespace std;
using namespace ug;
using namespace script;


////////////////////////////////////////////////////////////////////////////////
// interactive shells
//#define UG_PROMPT "\e[1mug:>\e[0m"
#define UG_PROMPT "ug:> "

#if defined(UG_USE_READLINE)
	#include <stdio.h>
	#include <readline/readline.h>
	#include <readline/history.h>
	void ug_cacheline(string str)
	{
		add_history(str.c_str());
	}
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
#else
	#if defined(UG_USE_LINENOISE)
		#include "externals/linenoise/linenoise.h"

		int CompletionFunction(char *buf, int len, int buflen, int iPrintCompletionList);

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

		extern int iOtherCompletitionsLength;
		extern const char **pOtherCompletitions;
	#else
		#define ug_cacheline(str)	
		string ug_readline(const char* prompt)
		{
			cout << prompt;
			string strBuffer;
			getline(cin, strBuffer);
			return strBuffer;
		}

	#endif
#endif



////////////////////////////////////////////////////////////////////////////////
// normal shell
int runShell(const char *prompt=UG_PROMPT)
{
	//	run the shell
	const char *completitions[] ={"quit", "exit"};
	while(1)
	{
#if defined(UG_USE_LINENOISE)
		iOtherCompletitionsLength=sizeof(completitions)/sizeof(completitions[0]);
		pOtherCompletitions=completitions;
#endif

		PROFILE_BEGIN(ug_readline);
			string buf = ug_readline(prompt);
		PROFILE_END();

#if defined(UG_USE_LINENOISE)
		iOtherCompletitionsLength=0;
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
debug_return debugShell()
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
		iOtherCompletitionsLength=sizeof(completitions)/sizeof(completitions[0]);
		pOtherCompletitions=completitions;
#endif
		PROFILE_BEGIN(ug_readline);
			string buf = ug_readline("debug:> ");
		PROFILE_END();

#if defined(UG_USE_LINENOISE)
		iOtherCompletitionsLength=0;
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

#ifdef UG_DEBUG
void SharedLibrariesLoaded()
{
}
#endif


////////////////////////////////////////////////////////////////////////////////
// main
int main(int argc, char* argv[])
{			
	#ifdef UG_PARALLEL
		pcl::Init(&argc, &argv);
	#endif

	PROFILE_FUNC();
	PROFILE_BEGIN(ugshellInit);


//	check if an output-process has been specified
	int outputProc = 0;
//	if "-outproc" is not found, outputProc won't be changed.
	ParamToInt(outputProc, "-outproc", argc, argv);
	GetLogAssistant().set_output_process(outputProc);
	
	const char* logFileName = NULL;
	if(ParamToString(&logFileName, "-logtofile", argc, argv))
		GetLogAssistant().enable_file_output(true, logFileName);

	if(FindParam("-noterm", argc, argv))
		GetLogAssistant().enable_terminal_output(false);
		
	if(FindParam("-profile", argc, argv))
		UGOutputProfileStatsOnExit(true);

//	ATTENTION
//	Make sure to initialize ug before accessing the registry or performing any
//	UG_LOGs. This is important, since mpi has to be initialized in a parallel
//	environment.

//	initialize
//	NOTE: UGInit is not called!!! We'll do this ourselfs, to allow early output.
//	init ug. bindings are created here
/*
	if(UGInit(&argc, &argv, outputProc) != 0)
	{
		std::cout << "ERROR occured in 'UGInit'. Unable to start UG Shell. "
				"Aborting programm execution.\n";
		UGFinalize();
		exit(1);
	}
*/
	bool bAbort=false;
	int ret = 0;

	LOG("********************************************************************************\n");
	std::string aux_str(""); // for nicer output we need some padding with spaces ...

	aux_str.append("* ugshell - ug").append(UGGetVersionString()).append(", head revision '").append(UGSvnRevision()).append("',");
	LOG(AppendSpacesToString(aux_str,80-1).append("*\n"));

	aux_str = "";
	aux_str.append("*                    compiled '").append(UGCompileDate()).append("'");
	LOG(AppendSpacesToString(aux_str,80-1).append("*\n"));

	aux_str = "";
	aux_str.append("*                    on '").append(UGBuildHost()).append("'.");
	LOG(AppendSpacesToString(aux_str,80-1).append("*\n"));

	LOG("*                                                                              *\n");
	LOG("* arguments:                                                                   *\n");
	LOG("*   -outproc id:         Sets the output-proc to id. Default is 0.             *\n");
	LOG("*   -ex scriptname:      Executes the specified script.                        *\n");
	LOG("*   -noquit:             Runs the interactive shell after specified script.    *\n");
	LOG("*   -noterm:             Terminal logging will be disabled.                    *\n");
	LOG("*   -logtofile filename: Output will be written to the specified file.         *\n");
#ifdef UG_PROFILER
	LOG("*   -profile:            Shows profile-output when the application terminates. *\n");
#endif
	LOG("* Additional parameters are passed to the script through ugargc and ugargv.    *\n");
	LOG("*                                                                              *\n");

//	initialize the rest of ug

//	INIT PATH
	try{
		LOG("* Initializing: paths... ");
		if(InitPaths((argv)[0]))	{UG_LOG("done");}
	//	if only false is returned, the error is non-fatal. Still continue shell
		else{
			UG_LOG("fail");
			UG_ERR_LOG("Path Initialization failed. Expect file access problem.\n")
			UG_ERR_LOG("Check environment variable UG4_ROOT.\n")
		}
	}
	catch(UGError& err)
	{
	//	if ugerror is throw, an internal fatal error occured, we termiate shell
		UG_ERR_LOG("UGError occurred during Path Initialization:\n");
		for(size_t i=0; i<err.num_msg(); i++)
			UG_ERR_LOG(err.get_file(i) << ":" << err.get_line(i) << " : " << err.get_msg(i) << "\n");
		bAbort = true;
	}
	catch(...){
		UG_ERR_LOG("Unknown error received during Path Initialization.\n");
		bAbort = true;
	}

//	INIT STANDARD BRIDGE
	try{
		UG_LOG(", bridge... ");
		bridge::InitBridge();
		UG_LOG("done");
	}
	catch(UGError& err)
	{
		UG_LOG("fail");
		UG_ERR_LOG("UGError occurred during Standard Bridge Initialization:\n");
		for(size_t i=0; i<err.num_msg(); i++)
			UG_ERR_LOG(err.get_file(i) << ":" << err.get_line(i) << " : " << err.get_msg(i) << "\n");
		bAbort = true;
	}
	catch(...){
		UG_LOG("fail");
		UG_ERR_LOG("Unknown error received during Standard Bridge Initialization.\n");
		bAbort = true;
	}

//	INIT PLUGINS
	try{
		#ifdef UG_PLUGINS
			UG_LOG(", plugins... ");
			#ifdef UG_EMBEDDED_PLUGINS
				InitializeEmbeddedPlugins(&bridge::GetUGRegistry(), "ug4/");
				UG_LOG("done");
			#else
				if(LoadPlugins(PathProvider::get_path(PLUGIN_PATH).c_str(), "ug4/"))	{UG_LOG("done");}
				else																	{UG_LOG("fail");}
			#endif
			UG_LOG("                 *\n");
		#else
			UG_LOG("                                  *\n");
		#endif
	}
	catch(UGError &err)
	{
		UG_LOG("fail");
		UG_LOG("                 *\n");
	//	if registration of plugin fails, we do abort the shell
	//	this could be skipped if registering of plugin would be done more
	//	carefully. (try/catch in load plugins)
		PathProvider::clear_current_path_stack();
		UG_ERR_LOG("UGError occurred during Plugin Initialization:\n");
		for(size_t i=0; i<err.num_msg(); i++)
			UG_ERR_LOG(err.get_file(i) << ":" << err.get_line(i) << " : " << err.get_msg(i) << "\n");
		bAbort = true;
	}
	catch(...){
		UG_LOG("fail");
		UG_LOG("                 *\n");
		UG_ERR_LOG("Unknown error received during Plugin Initialization.\n");
		bAbort = true;
	}

	LOG("********************************************************************************\n");


#ifdef UG_DEBUG
	SharedLibrariesLoaded();
#endif

	if(!bAbort)
	{
		try{
		//	check that registry is consistent. Else abort.
			if(bridge::GetUGRegistry().check_consistency() == false){
				UG_ERR_LOG("#### Registry ERROR: Registering of Standard Bridges and Plugins failed.")
				bAbort = true;
			}

			if(!bAbort){
			//	register the lua only functonality at the registry
				RegisterDefaultLuaBridge(&bridge::GetUGRegistry());

			//	check that registry is consistent. Else abort.
				if(bridge::GetUGRegistry().check_consistency() == false) {
					UG_ERR_LOG("#### Registry ERROR: Registering of Standard LUA Bridges failed.")
					bAbort =true;
				}
			}
		}
		catch(UGError &err)
		{
			UG_ERR_LOG("UGError occurred on registry check:\n");
			for(size_t i=0; i<err.num_msg(); i++)
				UG_ERR_LOG(err.get_file(i) << ":" << err.get_line(i) << " : " << err.get_msg(i) << "\n");
			bAbort = true;
		}
		catch(...){
			UG_ERR_LOG("Unknown error occurred on registry check.\n");
			bAbort = true;
		}
	}

	if(!bAbort){
		bool runInteractiveShell = true;

	//	we want to forward argc and argv to the lua-environment.
	//	we'll create a table for that.
	//	Please note that ugargv will neither contain the name of the program, nor
	//	the script, if one was specified.
		lua_State* L = script::GetDefaultLuaState();

	//	if a script shall be executed we do so now.
		int scriptParamIndex = GetParamIndex("-ex", argc, argv);
		int firstParamIndex = 1;
		const char* scriptName = NULL;

		if(scriptParamIndex != -1){
			if(argc > scriptParamIndex + 1){
			//	offset to the parameterlist
				firstParamIndex = scriptParamIndex + 2;
			//	get the name of the script
				scriptName = argv[scriptParamIndex + 1];
			}
		}

		int iNoQuit = GetParamIndex("-noquit", argc, argv);

	//	create ugargc and ugargv in lua

		lua_newtable(L);
		int ugargc=0;
		for(int i = firstParamIndex; i < argc; ++i){
			if(i == iNoQuit) continue;
		//	push the index to the table
			lua_pushnumber(L, ++ugargc);
		//	push the value to the table
			lua_pushstring(L, argv[i]);
		//	create the entry
			lua_settable(L, -3);
		}
		//	set the tables name
		lua_setglobal(L, "ugargv");

		lua_pushnumber(L, ugargc);
		lua_setglobal(L, "ugargc");

		// replace LUAs print function with our own, to use UG_LOG
		lua_register(L, "print", UGLuaPrint );
		lua_register(L, "write", UGLuaWrite );
	
		PROFILE_END(); // ugshellInit
	
	#if defined(UG_USE_LINENOISE)
		linenoiseSetCompletionFunction(CompletionFunction);
		SetDebugShell(debugShell);
	#endif
	
	
	//	if a script has been specified, then execute it now
	//	if a script is executed, we won't execute the interactive shell.
		if(scriptName)
		{
			try{
				if(!LoadUGScript_Parallel(scriptName))
				{
					UG_LOG("Can not find specified script ('" << scriptName << "'). Aborting.\n");
					bAbort=true;
					ret=1;
				}
			}
			catch(LuaError& err) {
				PathProvider::clear_current_path_stack();
				if(!err.get_msg().empty()){
					UG_LOG("LUA-ERROR: \n");
					for(size_t i=0;i<err.num_msg();++i)
						UG_LOG(err.get_msg(i)<<endl);
				}
				UG_LOG("aborting script parsing...\n");
			}
			catch(UGError &err)
			{
				PathProvider::clear_current_path_stack();
				UG_LOG("UGError:\n");
				for(size_t i=0; i<err.num_msg(); i++)
					UG_LOG(err.get_file(i) << ":" << err.get_line(i) << " : " << err.get_msg(i));
				UG_LOG("\n");
				UG_LOG("aborting script parsing...\n");
			}

			if(FindParam("-noquit", argc, argv))
				runInteractiveShell = true;
			else
				runInteractiveShell = false;
		}

		//	interactive shell may only be executed in serial environment
		#ifdef UG_PARALLEL
			if(runInteractiveShell){
				if(pcl::GetNumProcesses() > 1){
					UG_LOG("Parallel environment detected. Disabling interactive shell.\n");
					runInteractiveShell = false;
				}

			}
		#endif

		if(runInteractiveShell)
			ret = runShell();

	}//bAbort
	else
	{
		UG_LOG("Errors occured on startup (see error log). Terminating shell.\n");
	}

	PROFILE_BEGIN(ugshellFinalize);

	// Until fully tested for dependency problems, I commented out the following lines.
	// For using valgrind, you'll have to uncomment them.
	//ReleaseDefaultLuaState();
	#ifdef UG_EMBEDDED_PLUGINS
		//FinalizeEmbeddPlugins();
	#else
		//UnloadPlugins();
	#endif

	PROFILE_END();
	UGFinalize();

	UG_LOG(endl);

	return ret;
}

