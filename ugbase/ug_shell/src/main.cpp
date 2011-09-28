// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d19

#include <iostream>
#include <sstream>

#include "ug.h"
#include "bindings/lua/ug_script.h"
#include "ug_bridge/ug_bridge.h"
#include "bindings/lua/info_commands.h"
#include "common/util/parameter_parsing.h"

using namespace std;
using namespace ug;
using namespace script;

#define UG_PROMPT "ug:> "

#if defined(UG_USE_READLINE)
	#include <stdio.h>
	#include <readline/readline.h>
	#include <readline/history.h>
	#define ug_readline()		(readline(UG_PROMPT))
	#define ug_cacheline(str)	(add_history(str))
	#define ug_freeline(str)	(free(str))
#else
	#if defined(UG_USE_LINENOISE)
		#include "externals/linenoise/linenoise.h"
		#define ug_readline()		(linenoise(UG_PROMPT))
		#define ug_cacheline(str)	(linenoiseHistoryAdd(str))
		#define ug_freeline(str)	(free(str))		
	#else
		#define ug_readline()		(ReadlinePseudo(UG_PROMPT))
		#define ug_cacheline(str)	
		#define ug_freeline(str)	(delete[](str))	
	#endif
#endif


int CompletionFunction(char *buf, int len, int buflen, int iPrintCompletionList);


///	a fallback if the readline lib is not available.
/**	\todo this method could be optimized.*/
char* ReadlinePseudo(const char* prompt)
{
	cout << prompt;
	string strBuffer;
	getline(cin, strBuffer);
	
	char* buffer = new char[strBuffer.size() + 1];
	memcpy(buffer, strBuffer.c_str(), strBuffer.size() + 1);

	return buffer;
}

int main(int argc, char* argv[])
{
	PROFILE_FUNC();
//	check if an output-process has been specified
	int outputProc = 0;
//	if "-outproc" is not found, outputProc won't be changed.
	ParamToInt(outputProc, "-outproc", argc, argv);
	
	char* logFileName = NULL;
	if(ParamToString(&logFileName, "-logtofile", argc, argv))
		GetLogAssistant().enable_file_output(true, logFileName);

	if(FindParam("-noterm", argc, argv))
		GetLogAssistant().enable_terminal_output(false);
	
//	ATTENTION
//	Make sure to initialize ug before accessing the registry or performing any
//	UG_LOGs. This is important, since mpi has to be initialized in a parallel
//	environment.

//	init ug. bindings are created here
	if(UGInit(&argc, &argv, outputProc) != 0)
	{
		std::cout << "ERROR occured in 'UGInit'. Unable to start UG Shell. "
				"Aborting programm execution.\n";
		UGFinalize();
		exit(1);
	}

//	register the load script method.
//	this currently can't be done by the bridge, since it would
//	introduce a dependency to scripting (this has to be avoided at all costs!!!)
//todo: move it to ug_script
	bridge::GetUGRegistry()
		.add_function("ug_load_script", &LoadUGScript, "/ug4/shell",
					"success", "", "Loads and parses a script and returns whether it succeeded.")
		.add_function("ug_file_exists", &FileExists, "/ug4/shell",
					 "exists", "", "Returns true if a path exists, false if not.")
		.add_function("exit", &UGForceExit, "/ug4/shell",
					 "", "", "Immediatly terminates the application.");

	bool runInteractiveShell = true;

	LOG("********************************************************************************\n");
	LOG("* ugshell - v4.0.1                                                             *\n");
	LOG("*                                                                              *\n");
	LOG("* arguments:                                                                   *\n");
	LOG("*   -outproc id:         Sets the output-proc to id. Default is 0.             *\n");
	LOG("*   -ex scriptname:      Executes the specified script.                        *\n");
	LOG("*   -noquit:             Runs the interactive shell after specified script.    *\n");
	LOG("*   -noterm:             Terminal logging will be disabled.                    *\n");
	LOG("*   -logtofile filename: Output will be written to the specified file.         *\n");
	LOG("* Additional parameters are passed to the script through ugargc and ugargv.    *\n");
	LOG("********************************************************************************\n");

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

	if(FindParam("-noquit", argc, argv))
		firstParamIndex++;

//	create ugargc and ugargv in lua
	int ugargc = argc - firstParamIndex;
	lua_pushnumber(L, ugargc);
	lua_setglobal(L, "ugargc");
	
	lua_newtable(L);
	for(int i = 0; i < ugargc; ++i){
	//	push the index to the table
		lua_pushnumber(L, i+1);
	//	push the value to the table
		lua_pushstring(L, argv[firstParamIndex + i]);
	//	create the entry
		lua_settable(L, -3);
	}
//	set the tables name
	lua_setglobal(L, "ugargv");
	
	// replace LUAs print function with our own, to use UG_LOG
	lua_register(L, "print", UGLuaPrint );
	lua_register(L, "write", UGLuaWrite );

//	if a script has been specified, then execute it now
//	if a script is executed, we won't execute the interactive shell.
	if(scriptName)
	{
		try{
			if(!LoadUGScript(scriptName))
			{
				UG_LOG("Can not find specified script ('" << scriptName << "'). Aborting.\n");
				UGFinalize();
				return 1;
			}
		}
		catch(LuaError& err) {
			UG_LOG("PARSE ERROR: \n");
			for(size_t i=0;i<err.num_msg();++i)
				UG_LOG(err.get_msg(i)<<endl);

			if(err.terminate()){
				UGFinalize();
				return 0; // exit with code 0
			}
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
	{
#if defined(UG_USE_LINENOISE)
		linenoiseSetCompletionFunction(CompletionFunction);
#endif

	//	run the shell
		while(1)
		{
			char* buffer = ug_readline();
			if(buffer){
				if(!(strcmp(buffer, "exit") && strcmp(buffer, "quit")))
					break;
					
				size_t len = strlen(buffer);
				if(len)
				{
					if(buffer[len-1]=='?')
					{
						buffer[len-1] = 0x00;
						bridge::UGTypeInfo(buffer);
						continue;
					}
					
					ug_cacheline(buffer);

					try
					{
						script::ParseBuffer(buffer, "interactive shell");
					}
					catch(LuaError& err)
					{
						UG_LOG("PARSE ERROR: \n");
						for(size_t i=0;i<err.num_msg();++i)
							UG_LOG(err.get_msg(i)<<endl);
						if(err.terminate())
						{
							UGFinalize();
							return 0; // exit with code 0
						}
					}
					
					
				}
				ug_freeline(buffer);
			}
		}
	//todo:	clear the history (add ug_freelinecache)
	}

	LOG(endl);

	UGFinalize();

	return 0;
}

