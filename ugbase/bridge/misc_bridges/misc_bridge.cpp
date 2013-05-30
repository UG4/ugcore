/*
 * misc_bridge.cpp
 *
 *  Created on: 21.03.2011
 *      Author: Martin Rupp
 */

#include "registry/registry.h"
#include "bridge/bridge.h"
#include <string>
#include <sstream>
#include <ctime>
#include "common/util/string_util.h"
#include "compile_info/compile_info.h"
#include "common/util/crc32.h"
#include "../plugins/experimental/jitsg/stopwatch.h"
using namespace std;

void ug_backtrace();
namespace ug
{
	
void PrintLUA();
namespace bridge
{

/// \addtogroup misc_bridge
/// \{

uint32 GetLogAssistantTag(const char *s)
{
	/*if(strcmp(s, "MAIN") == 0) return LogAssistant::MAIN;
	if(strcmp(s, "APP") == 0) return LogAssistant::APP;
	if(strcmp(s, "LIB_GRID") == 0) return LogAssistant::LIB_GRID;
	if(strcmp(s, "LIB_GRID_REFINER") ==0) return LogAssistant::LIB_GRID_REFINER;
	if(strcmp(s, "LIB_DISC") == 0) return LogAssistant::LIB_DISC;
	if(strcmp(s, "LIB_DISC_ASSEMBLE") == 0) return LogAssistant::LIB_DISC_ASSEMBLE;
	if(strcmp(s, "LIB_DISC_D3F") == 0) return LogAssistant::LIB_DISC_D3F;
	if(strcmp(s, "LIB_DISC_MULTIGRID") == 0) return LogAssistant::LIB_DISC_MULTIGRID;
	if(strcmp(s, "LIB_DISC_NEWTON") == 0) return LogAssistant::LIB_DISC_NEWTON;
	if(strcmp(s, "LIB_DISC_LINKER") == 0) return LogAssistant::LIB_DISC_LINKER;
	if(strcmp(s, "LIB_DISC_TRANSFER") == 0) return LogAssistant::LIB_DISC_TRANSFER;
	if(strcmp(s, "LIB_DISC_DISCRETE_FUNCTION") == 0) return LogAssistant::LIB_DISC_DISCRETE_FUNCTION;
	if(strcmp(s, "LIB_DISC_OUTPUT") == 0) return LogAssistant::LIB_DISC_OUTPUT;
	if(strcmp(s, "LIB_DISC_OPERATOR_INVERSE") == 0) return LogAssistant::LIB_DISC_OPERATOR_INVERSE;
	if(strcmp(s, "LIB_ALG_LINEAR_OPERATOR") == 0) return LogAssistant::LIB_ALG_LINEAR_OPERATOR;
	if(strcmp(s, "LIB_ALG_LINEAR_SOLVER") == 0) return LogAssistant::LIB_ALG_LINEAR_SOLVER;
	if(strcmp(s, "LIB_ALG_VECTOR") == 0) return LogAssistant::LIB_ALG_VECTOR;
	if(strcmp(s, "LIB_ALG_MATRIX") == 0) return LogAssistant::LIB_ALG_MATRIX;
	if(strcmp(s, "LIB_ALG_AMG") == 0) return LogAssistant::LIB_ALG_AMG;
	if(strcmp(s, "LIB_PCL") == 0) return LogAssistant::LIB_PCL;*/
	if(GetDebugIDManager().debug_id_registered(s) == false)
	{
		UG_LOG("DebugID " << s << " not found!");
		return 0;
	}
	else return crc32(s);
}

////////////////////////////////////////////////////////////////////////////////
/*
  A bunch of methods used for printing build configuration
  (info for every primary and some dependent cmake setting(s)),
  order as by a plain 'cmake ../'

  (To generate a list of macro definitions:
  grep add_definitions cmake/ug_includes.cmake
  grep "\-D" cmake/ug_includes.cmake.)
*/

// DIM:


#ifdef UG_DIM_1
bool IsDefinedUG_DIM_1() { return true; }
#else
bool IsDefinedUG_DIM_1() { return false; }
#endif

#ifdef UG_DIM_2
bool IsDefinedUG_DIM_2() { return true; }
#else
bool IsDefinedUG_DIM_2() { return false; }
#endif

#ifdef UG_DIM_3
bool IsDefinedUG_DIM_3() { return true; }
#else
bool IsDefinedUG_DIM_3() { return false; }
#endif

// CPU:
#ifdef UG_CPU_1
bool IsDefinedUG_CPU_1() { return true; }
#else
bool IsDefinedUG_CPU_1() { return false; }
#endif

#ifdef UG_CPU_2
bool IsDefinedUG_CPU_2() { return true; }
#else
bool IsDefinedUG_CPU_2() { return false; }
#endif

#ifdef UG_CPU_3
bool IsDefinedUG_CPU_3() { return true; }
#else
bool IsDefinedUG_CPU_3() { return false; }
#endif

#ifdef UG_CPU_4
bool IsDefinedUG_CPU_4() { return true; }
#else
bool IsDefinedUG_CPU_4() { return false; }
#endif

#ifdef UG_CPU_VAR
bool IsDefinedUG_CPU_VAR() { return true; }
#else
bool IsDefinedUG_CPU_VAR() { return false; }
#endif

// STATIC:
#ifdef UG_STATIC
bool IsDefinedUG_STATIC() { return true; }
#else
bool IsDefinedUG_STATIC() { return false; }
#endif

// DEBUG, DEBUG_LOGS:
#ifdef UG_DEBUG
bool IsDefinedUG_DEBUG() { return true; }
#else
bool IsDefinedUG_DEBUG() { return false; }
#endif

#ifdef UG_ENABLE_DEBUG_LOGS
bool IsDefinedUG_ENABLE_DEBUG_LOGS() { return true; }
#else
bool IsDefinedUG_ENABLE_DEBUG_LOGS() { return false; }
#endif

// PARALLEL:
#ifdef UG_PARALLEL
bool IsDefinedUG_PARALLEL() { return true; }
#else
bool IsDefinedUG_PARALLEL() { return false; }
#endif

// PCL_DEBUG_BARRIER:
#ifdef PCL_DEBUG_BARRIER_ENABLED
bool IsDefinedPCL_DEBUG_BARRIER_ENABLED() { return true; }
#else
bool IsDefinedPCL_DEBUG_BARRIER_ENABLED() { return false; }
#endif

// PROFILE:
#ifdef UG_PROFILER
bool IsDefinedUG_PROFILER() { return true; }
#else
bool IsDefinedUG_PROFILER() { return false; }
#endif

// PROFILE_PCL:
#ifdef PROFILE_PCL
bool IsDefinedPROFILE_PCL() { return true; }
#else
bool IsDefinedPROFILE_PCL() { return false; }
#endif

// ALGEBRA - derived, no output by 'cmake' until now:
#ifdef UG_ALGEBRA
bool IsDefinedUG_ALGEBRA() { return true; }
#else
bool IsDefinedUG_ALGEBRA() { return false; }
#endif

// LAPACK_AVAILABLE - derived:
#ifdef LAPACK_AVAILABLE
bool IsDefinedLAPACK_AVAILABLE() { return true; }
#else
bool IsDefinedLAPACK_AVAILABLE() { return false; }
#endif

// BLAS_AVAILABLE - derived:
#ifdef BLAS_AVAILABLE
bool IsDefinedBLAS_AVAILABLE() { return true; }
#else
bool IsDefinedBLAS_AVAILABLE() { return false; }
#endif

// External libraries:
#ifdef UG_METIS
bool IsDefinedUG_METIS() { return true; }
#else
bool IsDefinedUG_METIS() { return false; }
#endif

#ifdef UG_PARMETIS
bool IsDefinedUG_PARMETIS() { return true; }
#else
bool IsDefinedUG_PARMETIS() { return false; }
#endif

#ifdef UG_TETGEN
bool IsDefinedUG_TETGEN() { return true; }
#else
bool IsDefinedUG_TETGEN() { return false; }
#endif

#ifdef UG_HYPRE
bool IsDefinedUG_HYPRE() { return true; }
#else
bool IsDefinedUG_HYPRE() { return false; }
#endif

#ifdef UG_HLIBPRO
bool IsDefinedUG_HLIBPRO() { return true; }
#else
bool IsDefinedUG_HLIBPRO() { return false; }
#endif

// VRL - derived (depends on target):
#ifdef UG_FOR_VRL
bool IsDefinedUG_FOR_VRL() { return true; }
#else
bool IsDefinedUG_FOR_VRL() { return false; }
#endif

// PLUGIN - derived (depends on target)
#ifdef UG_PLUGINS
bool IsDefinedUG_PLUGINS() { return true; }
#else
bool IsDefinedUG_PLUGINS() { return false; }
#endif

// BRIDGE - derived (depends on target)
#ifdef UG_BRIDGE
bool IsDefinedUG_BRIDGE() { return true; }
#else
bool IsDefinedUG_BRIDGE() { return false; }
#endif


/// prints CMake build parameters in a quite compact (pairwise) form
void PrintBuildConfiguration()
{
	std::string aux_str("");

	UG_LOG("--------------------------------------------------------------------------------\n");
	UG_LOG("Build configuration:\n\n");

	// TODO: Maybe there is a nicer order for displaying the parameters!?

	// 1. Primary (direct, independent) cmake parameters:
	UG_LOG("1. Primary 'cmake' parameters:\n");

	// first pair
	aux_str = "";
	aux_str.append("TARGET:            ").append(UG_TARGET);
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));
	aux_str = "";
	aux_str.append("STATIC:            ").append( (IsDefinedUG_STATIC() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append("\n"));

	// next pair
	aux_str = "";
	aux_str.append("PARALLEL:          ").append( (IsDefinedUG_PARALLEL() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));
	aux_str = "";
	aux_str.append("PCL_DEBUG_BARRIER: ").append( (IsDefinedPCL_DEBUG_BARRIER_ENABLED() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append("\n"));

	// next pair
	aux_str = "";
	aux_str.append("DEBUG:             ").append( (IsDefinedUG_DEBUG() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));
	aux_str = "";
	aux_str.append("DEBUG_LOGS:        ").append( (IsDefinedUG_ENABLE_DEBUG_LOGS() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append("\n"));

	// next pair
	aux_str = "";
	aux_str.append("PROFILER:          ").append( (IsDefinedUG_PROFILER() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));
	aux_str = "";
	aux_str.append("PROFILE_PCL:       ").append( (IsDefinedPROFILE_PCL() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append("\n"));

	// next pair
	// Please note that there are also independent cmake parameters 'LAPACK' and 'BLAS',
	// but we are only interested if LAPACK/BLAS was found or not (if requested).
	aux_str = "";
	aux_str.append("LAPACK available:  ").append( (IsDefinedLAPACK_AVAILABLE() ? "YES" : "NO ") );
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));
	aux_str = "";
	aux_str.append("BLAS   available:  ").append( (IsDefinedBLAS_AVAILABLE() ? "YES" : "NO ") );
	UG_LOG(AppendSpacesToString(aux_str,40).append("\n"));

	// DIM gets its own line
	UG_LOG("DIM:               ");
	if (IsDefinedUG_DIM_1() && IsDefinedUG_DIM_2() && IsDefinedUG_DIM_3()) {
		UG_LOG("ALL");
	} else {
		UG_LOG( (IsDefinedUG_DIM_1() ? "1 " : "") ); // or maybe a "-" for dim not defined!?
		UG_LOG( (IsDefinedUG_DIM_2() ? "2 " : "") );
		UG_LOG( (IsDefinedUG_DIM_3() ? "3 " : "") );
	}
	UG_LOG("\n");

	// next pair
	aux_str = "";
	aux_str.append("CPU:               ");
	if (IsDefinedUG_CPU_1() && IsDefinedUG_CPU_2() &&
		IsDefinedUG_CPU_3() && IsDefinedUG_CPU_4() &&
		IsDefinedUG_CPU_VAR()) {
		aux_str.append("ALL");
	} else {
		aux_str.append( (IsDefinedUG_CPU_1() ? "1 " : "") );
		aux_str.append( (IsDefinedUG_CPU_2() ? "2 " : "") );
		aux_str.append( (IsDefinedUG_CPU_3() ? "3 " : "") );
		aux_str.append( (IsDefinedUG_CPU_4() ? "4 " : "") );
		aux_str.append( (IsDefinedUG_CPU_VAR() ? "VAR" : "") );
	}
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));

	// We've decided so far not to display the following derived parameters!

	// 2. External stuff:
	UG_LOG("\n2. External libraries:\n");
	aux_str = "";
	aux_str.append("METIS:             ").append( (IsDefinedUG_METIS() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));

	aux_str = "";
	aux_str.append("PARMETIS:          ").append( (IsDefinedUG_PARMETIS() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append("\n"));

	// next pair
	aux_str = "";
	aux_str.append("TETGEN:            ").append( (IsDefinedUG_TETGEN() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));

	aux_str = "";
	aux_str.append("HYPRE:             ").append( (IsDefinedUG_HYPRE() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append("\n"));

	// next pair
	aux_str = "";
	aux_str.append("HLIBPRO:           ").append( (IsDefinedUG_HLIBPRO() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append(""));

	aux_str = "";
	//aux_str.append("?:               ").append( (IsDefinedUG_?() ? "ON " : "OFF") );
	UG_LOG(AppendSpacesToString(aux_str,40).append("\n"));
	UG_LOG("--------------------------------------------------------------------------------\n");
}

void PrintBuildConfigurationExtended()
{
	// We've decided so far not to display the following derived parameters:
	// 2. Derived parameters (no direct parameters to cmake):
	UG_LOG("Derived parameters:\n");
	UG_LOG("Build for VRL:     ");
	UG_LOG( (IsDefinedUG_FOR_VRL() ? "ON " : "OFF") );
	UG_LOG("\n");

	UG_LOG("PLUGIN:            ");
	UG_LOG( (IsDefinedUG_PLUGINS() ? "ON " : "OFF") );
	UG_LOG("\n");

	UG_LOG("BRIDGE:            ");
	UG_LOG( (IsDefinedUG_BRIDGE() ? "ON " : "OFF") );
	UG_LOG("\n");

	UG_LOG("ALGEBRA:           ");
	UG_LOG( (IsDefinedUG_ALGEBRA() ? "ON " : "OFF") );
	UG_LOG("\n\n");

// Derived parameters END


	UG_LOG("--------------------------------------------------------------------------------\n");
}

////////////////////////////////////////////////////////////////////////////////

void SetDebugLevel(const char* strTag, int level)
{
	GetLogAssistant().set_debug_level(strTag, level);
}

string GetSVNRevision()
{
	return string(UGSvnRevision());
}

string GetCompileDate()
{
	return string(UGCompileDate());
}

string GetBuildHostname()
{
	return string(UGBuildHost());
}


double GetClockS()
{
	return clock()/((double)CLOCKS_PER_SEC);
}



void RegisterBridge_Misc(Registry &reg, string parentGroup)
{

	{
		stringstream ss; ss << parentGroup << "/Util/Log";
		string grp = ss.str();
		reg.add_class_<LogAssistant>("LogAssistant", grp)
			.add_method("enable_file_output", &LogAssistant::enable_file_output,
					"", "bEnable#filename", "Please note that only the filename given at the first call is considered")
			.add_method("rename_log_file", &LogAssistant::rename_log_file,
					"", "filename", "Renames previously opened logfile to the name given")
			.add_method("enable_terminal_output", &LogAssistant::enable_terminal_output,
					"", "bEnable", "enables or disables terminal output (enabled by default)")
			.add_method("set_debug_levels", &LogAssistant::set_debug_levels, "", "lev", "sets the debug level of all tags to 'lev'")

			.add_method("set_debug_level", static_cast<bool (LogAssistant::*)(const char *, int)>(&LogAssistant::set_debug_level_noninline), "", "tag#lev", "sets the debug level of Tag 'tag' to level 'lev'")
			.add_method("get_debug_level", static_cast<int (LogAssistant::*)(const char *) const>(&LogAssistant::get_debug_level_noninline), "", "tag#lev", "gets the debug level of Tag 'tag'")

			.add_method("get_debug_IDs", &LogAssistant::get_registered_debug_IDs, "debug IDs")

			.add_method("set_output_process", &LogAssistant::set_output_process, "", "procRank", "Sets the process which shall output its data.")
			.add_method("get_output_process", &LogAssistant::get_output_process, "", "", "Returns the process which outputs its data.")
			.add_method("is_output_process", &LogAssistant::set_output_process, "", "procRank", "Returns whether the process outputs its data.")
			.add_method("flush", &LogAssistant::flush, "", "", "flush all buffers, especially the file buffer.")
			;
		reg.add_function("GetLogAssistant", &GetLogAssistant, grp, "the log assistant");
		//reg.add_function("GetLogAssistantTag", &GetLogAssistantTag, grp, "integer tag for use int set_debug_level");
		//reg.add_function("SetDebugLevel", &SetDebugLevel, grp);
		reg.add_function("GetClockS", &GetClockS, grp);
		reg.add_function("PrintLUA", &PrintLUA, grp);


		reg.add_class_<stopwatch>("stopwatch", grp)
				.add_constructor()
				.add_method("start", &stopwatch::start)
				.add_method("diff", &stopwatch::diff);
	}

	{
		stringstream ss; ss << parentGroup << "/Util/Internal";
		string grp = ss.str();
		

		reg.add_function("DefinedUG_DEBUG", &IsDefinedUG_DEBUG, grp, ""); // TODO: only for backward compatibility - remove if not / no longer used!
		reg.add_function("DefinedUG_ENABLE_DEBUG_LOGS", &IsDefinedUG_ENABLE_DEBUG_LOGS, grp, ""); // TODO: only for backward compatibility - remove if not / no longer used!
		
#define ADD_DEFINED_FUNC(VAR) reg.add_function("IsDefined"#VAR, IsDefined##VAR, grp)
		ADD_DEFINED_FUNC(UG_DEBUG);
		ADD_DEFINED_FUNC(UG_PARALLEL);
		ADD_DEFINED_FUNC(UG_DIM_1);
		ADD_DEFINED_FUNC(UG_DIM_2);
		ADD_DEFINED_FUNC(UG_DIM_3);
		ADD_DEFINED_FUNC(UG_CPU_1);
		ADD_DEFINED_FUNC(UG_CPU_2);
		ADD_DEFINED_FUNC(UG_CPU_3);		
		ADD_DEFINED_FUNC(UG_ENABLE_DEBUG_LOGS);
		ADD_DEFINED_FUNC(LAPACK_AVAILABLE);
		ADD_DEFINED_FUNC(BLAS_AVAILABLE);
		ADD_DEFINED_FUNC(UG_METIS);
		ADD_DEFINED_FUNC(UG_HYPRE);
		ADD_DEFINED_FUNC(UG_HLIBPRO);

		reg.add_function("PrintBuildConfiguration", &PrintBuildConfiguration, grp, "");
		reg.add_function("PrintBuildConfigurationExtended", &PrintBuildConfigurationExtended, grp, "");
		
		

		reg.add_function("GetSVNRevision", &GetSVNRevision, grp);
		reg.add_function("GetCompileDate", &GetCompileDate, grp);
		reg.add_function("GetBuildHostname", &GetBuildHostname, grp);

		reg.add_function("LevenshteinDistance", &LevenshteinDistance, grp);
		
		reg.add_function("ug_backtrace", &ug_backtrace, grp);
		
		
	}
	

}

// end group misc_bridge
/// \}

}	// end namespace bridge

} // end namespace ug
