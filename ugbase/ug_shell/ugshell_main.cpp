/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include <iostream>
#include <sstream>

#include "ug.h"
#include "bindings/lua/lua_util.h"
#include "bridge/bridge.h"

#include "common/util/parameter_parsing.h"
#include "common/util/file_util.h"
#include "common/profiler/memtracker.h"
#include "common/util/os_info.h"
#include "common/util/path_provider.h"
#include "common/profiler/profile_node.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif
#include "common/catch_std.h"

void get_undeleted();
#include "compile_info/compile_info.h"
#include "bindings/lua/lua_debug.h"
#include "shell.h"
#include "common/util/stringify.h"

using namespace std;
using namespace ug;
using namespace script;


namespace ug
{

namespace bridge
{

void SetLuaNamespace(string name, string value);
}
}
/**
 * \defgroup ugbase_ugshell UGShell
 * \ingroup ugbase
 * \brief the shell for ug4
 * \{
 */


#ifdef UG_DEBUG
void SharedLibrariesLoaded()
{
}
#endif


//	a symbol preceding error messages
static const char* errSymb = " % ";
void quit_all_mpi_procs_in_parallel()
{
#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1){
		UG_LOG(errSymb<<"ABORTING all mpi processes."<<std::endl)
		pcl::Abort();
	}
#endif
}

////////////////////////////////////////////////////////////////////////////////
// main
int main(int argc, char* argv[])
{			
	PROFILE_FUNC();
	PROFILE_BEGIN(ugshellInit);

	#ifdef UG_PARALLEL
		pcl::Init(&argc, &argv);
	#endif


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
	aux_str.append("*                    on '").append(UGBuildHost()).append("'");
#ifdef UG_DEBUG
	aux_str.append(" (DEBUG).");
#else
	aux_str.append(".");
#endif
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
	//	if an UGError is thrown, an internal fatal error occured, we terminate shell
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
		UG_LOG(", plugins... ");
		if(UGInitPlugins() == true)
		{	UG_LOG("done"); }
		else																	
		{	UG_LOG("fail");	}
		
		UG_LOG("                 *\n");
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

		SetLuaDebugIDs(L);
		SetLuaUGArgs(L, argc, argv, firstParamIndex, iNoQuit);

		// replace LUAs print function with our own, to use UG_LOG
		RegisterStdLUAFunctions(L);
	
		ug::bridge::InitShell();

		PROFILE_END(); // ugshellInit

		EnableMemTracker(true);
	//	if a script has been specified, then execute it now
	//	if a script is executed, we won't execute the interactive shell.
		if(scriptName)
		{
			try{
				if(!LoadUGScript_Parallel(scriptName))
				{
					UG_LOG("Cannot find specified script ('" << scriptName << "'). Aborting.\n");
					bAbort=true;
				}
			}
			catch(SoftAbort& err){
				UG_LOG("Execution of script-buffer aborted with the following message:\n")
				UG_LOG(err.get_msg() << std::endl);
			}
			catch(LuaError& err) {
				PathProvider::clear_current_path_stack();
				if(err.show_msg()){
					if(!err.get_msg().empty()){
						UG_LOG("LUA-ERROR: \n");
						for(size_t i=0;i<err.num_msg();++i)
							UG_LOG(err.get_msg(i)<<endl);
					}
				}
				UG_LOG(errSymb<<"ABORTING script parsing.\n");
				quit_all_mpi_procs_in_parallel();
				bAbort=true;
			}
			catch(UGError &err)
			{
				PathProvider::clear_current_path_stack();
				UG_LOG("UGError:\n");
				for(size_t i=0; i<err.num_msg(); i++)
					UG_LOG(err.get_file(i) << ":" << err.get_line(i) << " : " << err.get_msg(i));
				UG_LOG("\n");
				UG_LOG(errSymb<<"ABORTING script parsing.\n");
				quit_all_mpi_procs_in_parallel();
				bAbort=true;
			}
			CATCH_STD_EXCEPTIONS();

//#ifdef UG_DEBUG
			try{
				script::ParseAndExecuteBuffer("if util ~= nil and util.CheckAndPrintHelp ~= nil then util.CheckAndPrintHelp(\"\\n-help: Available Command Line Arguments:\") end\n"
							"if util ~= nil and util.PrintIgnoredArguments ~= nil then print(\"\") util.PrintIgnoredArguments() end\n"
									, "ugshell_main");
			}
			catch(SoftAbort& err){
				UG_LOG("Execution of script-buffer aborted with the following message:\n");
				UG_LOG(err.get_msg() << std::endl);
			}

			if(FindParam("-noquit", argc, argv))
				runInteractiveShell = true;
			else
				runInteractiveShell = false;
//#endif
		}
		EnableMemTracker(false);
		//DisplayVacantMemory();

		//	interactive shell may only be executed in serial environment
		#ifdef UG_PARALLEL
			if(runInteractiveShell){
				if(pcl::NumProcs() > 1){
					UG_LOG("Parallel environment detected. Disabling interactive shell.\n");
					runInteractiveShell = false;
				}

			}
		#endif

		if(runInteractiveShell){
			try{
				ret = ug::bridge::RunShell();
			}
			catch(SoftAbort& err){
				UG_LOG("Execution of interactive shell aborted with the following message:\n");
				UG_LOG(err.get_msg() << std::endl);
			}
		}

	}//bAbort
	else
	{
		UG_LOG("Errors occured on startup (see error log). Terminating shell.\n");
	}

	PROFILE_BEGIN(ugshellFinalize);

	// Until fully tested for dependency problems, I commented out the following lines.
	// For using valgrind, you'll have to uncomment them.
	//ReleaseDefaultLuaState();
	#ifdef UG_PLUGINS
		//UnloadPlugins();
	#endif

	PROFILE_END();
	UGFinalize();

	ug::UGProfileNode::CheckForTooSmallNodes();
	UG_LOG(endl);

	// If shell aborted and no custom ret-value defined then return '1'; return 'ret' otherwise
	return bAbort&&ret==0?1:ret;
}

// end group ugbase_ugshell
/// \}
