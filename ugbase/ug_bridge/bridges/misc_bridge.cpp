/*
 * misc_bridge.cpp
 *
 *  Created on: 21.03.2011
 *      Author: Martin Rupp
 */


#include "ug_script/ug_script.h"
#include "../registry.h"
#include "../ug_bridge.h"

namespace ug
{
namespace bridge
{

// todo: support enums natively and remove this
template <>
struct PLStack<LogAssistant::Tags>
{
	static void push(ParameterStack& ps)
	{
		ps.push_integer();
	}
	static void write(ParameterStack& ps, LogAssistant::Tags data, int index)
	{
		ps.set_integer(index, data);
	}
	static LogAssistant::Tags read(const ParameterStack& ps, int index)
	{
		return (LogAssistant::Tags) ps.to_integer(index);
	}
};


LogAssistant::Tags GetLogAssistantTag(const char *s)
{
	if(strcmp(s, "MAIN") == 0) return LogAssistant::MAIN;
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
	if(strcmp(s, "LIB_PCL") == 0) return LogAssistant::LIB_PCL;
	UG_LOG("Tag " << s << " not found!");
	return LogAssistant::MAIN;
}

#ifdef UG_DEBUG
bool DefinedUG_DEBUG() { return true; }
#else
bool DefinedUG_DEBUG() { return false; }
#endif


#ifdef UG_ENABLE_DEBUG_LOGS
bool DefinedUG_ENABLE_DEBUG_LOGS() { return true; }
#else
bool DefinedUG_ENABLE_DEBUG_LOGS() { return false; }
#endif

void SetDebugLevel(const char* strTag, int level)
{
	GetLogAssistant().set_debug_level(GetLogAssistantTag(strTag), level);
}

void int_srand(int seed)
{
	srand((unsigned int)seed);
}

bool RegisterMiscFunctions(Registry &reg, const char* parentGroup)
{
	{
		std::stringstream group; group << parentGroup << "/log";

		reg.add_class_<LogAssistant>("LogAssistant", group.str().c_str())
			.add_method("enable_file_output", &LogAssistant::enable_file_output,
					"", "bEnable#filename", "Please note that only the filename given at the first call is considered")
			.add_method("enable_terminal_output", &LogAssistant::enable_terminal_output,
					"", "bEnable", "enables or disables terminal output (enabled by default)")
			.add_method("set_debug_levels", &LogAssistant::set_debug_levels,
					"", "lev", "sets the debug level of all tags to 'lev'")
			.add_method("set_debug_level", &LogAssistant::set_debug_level, "", "tag#lev", "sets the debug level of Tag 'tag' to level 'lev'")
			.add_method("get_debug_level", &LogAssistant::get_debug_level, "", "tag", "use GetLogAssistantTag to get a tag from a tag-string.")
			;
		reg.add_function("GetLogAssistant", &GetLogAssistant, group.str().c_str(), "the log assistant");
		reg.add_function("GetLogAssistantTag", &GetLogAssistantTag, group.str().c_str(), "integer tag for use int set_debug_level");
		reg.add_function("SetDebugLevel", &SetDebugLevel, group.str().c_str());
	}

	{

		std::stringstream group; group << parentGroup << "/internal";

		reg.add_function("DefinedUG_DEBUG", &DefinedUG_DEBUG, group.str().c_str(), "");
		reg.add_function("DefinedUG_ENABLE_DEBUG_LOGS", &DefinedUG_ENABLE_DEBUG_LOGS, group.str().c_str(), "");
		reg.add_function("srand", int_srand, "", "seed", "The pseudo-random number generator is initialized using the argument passed as seed.");
	}

	return true;
}


}

}
