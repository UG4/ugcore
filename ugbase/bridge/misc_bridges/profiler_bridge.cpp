/*
 * profiler_bridge.cpp
 *
 *  Created on: 11.02.2011
 *      Author: Martin Rupp
 */

#include "registry/registry.h"
#include "bridge/bridge.h"
#include "common/profiler/profiler.h"
#include "common/profiler/profile_node.h"
#include "ug.h" // Required for UGOutputProfileStatsOnExit.
#include <string>
#include <sstream>
#include "../util_overloaded.h"
using namespace std;

namespace ug
{

static void UpdateProfiler_BridgeImpl(number damping){
	PROFILER_UPDATE(damping);
}

static void SetShinyCallLoggingMaxFrequency(size_t maxFreq)
{
#ifdef SHINY_CALL_LOGGING
	ug::g_ShinyCallLoggingMaxFreq = maxFreq;
#else
	UG_LOG("SHINY CALL LOGGING NOT ENABLED! Enable with 'cmake -DSHINY_CALL_LOGGING=ON ..'")
#endif
}

  //void PrintLUA();
namespace bridge
{


/// \defgroup profiler_bridge Profiler Bridge
/// \ingroup misc_bridge
/// \{
//////////////////////////////////////////////////////////////////////////////////////////

void RegisterBridge_Profiler(Registry &reg, string parentGroup)
{
	stringstream ss; ss << parentGroup << "/Util/Profiler";
	string grp = ss.str();

	reg.add_class_<UGProfileNode>("UGProfileNode", grp)
		// call tree
		.add_method("call_tree",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, call_tree, ()),
				"string with call tree")
		.add_method("call_tree",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, call_tree, (double dSkipMarginal)),
				"string with call tree",
				"dSkipMarginal")

		// self time
		.add_method("child_self_time_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, child_self_time_sorted, ()),
				"string with sorted childs", "", "childs are sorted by self time")
		.add_method("child_self_time_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, child_self_time_sorted, (double dSkipMarginal)),
				"string with sorted childs", "dSkipMarginal", "childs are sorted by self time")

		// total time
		.add_method("total_time_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, total_time_sorted, ()),
				"string with sorted childs", "", "childs are sorted by total time")
		.add_method("total_time_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, total_time_sorted, (double dSkipMarginal)),
				"string with sorted childs", "dSkipMarginal", "childs are sorted by total time")


		// self memory
		.add_method("child_self_memory_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, child_self_memory_sorted, ()),
				"string with sorted childs", "", "childs are sorted by self memory")
		.add_method("child_self_memory_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, child_self_memory_sorted, (double dSkipMarginal)),
				"string with sorted childs", "dSkipMarginal", "childs are sorted by self memory")

		// total memory
		.add_method("total_memory_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, total_memory_sorted, ()),
				"string with sorted childs", "", "childs are sorted by total memory")
		.add_method("total_memory_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, total_memory_sorted, (double dSkipMarginal)),
				"string with sorted childs", "dSkipMarginal", "childs are sorted by total memory")

		// entry count
		.add_method("entry_count_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, entry_count_sorted, ()),
				"string with sorted childs", "", "childs are sorted by entry count")
		.add_method("entry_count_sorted",
				OVERLOADED_CONST_METHOD_PTR(string, UGProfileNode, entry_count_sorted, (double dSkipMarginal)),
				"string with sorted childs", "dSkipMarginal", "childs are sorted by entry count")

		// misc
		.add_method("get_avg_entry_count", &UGProfileNode::get_avg_entry_count,
				"number of entries in this profiler node", "")
		.add_method("get_avg_self_time_ms", &UGProfileNode::get_avg_self_time_ms,
				"time in milliseconds spend in this node excluding subnodes", "")
		.add_method("get_avg_total_time_ms", &UGProfileNode::get_avg_total_time_ms,
				"time in milliseconds spend in this node including subnodes", "")
		.add_method("is_valid", &UGProfileNode::valid, "true if node has been found", "")

	  		.add_method("groups", &UGProfileNode::groups, "", "")

		;
		/*.add_method("__tostring", &UGProfileNode::tostring, "tostring")
		.add_method("__unm", &UGProfileNode::unm, "unm")
		.add_method("__add", &UGProfileNode::add, "add");*/

	//	reg.add_function("PrintLUA", &PrintLUA, grp);


	typedef const UGProfileNode* TT;
	reg.add_function("GetProfileNode", OVERLOADED_FUNCTION_PTR(const UGProfileNode *, GetProfileNode, (const char*name)),
			grp, "a profile node", "name", "if root = null, return");
	reg.add_function("GetProfileNode", OVERLOADED_FUNCTION_PTR(const UGProfileNode *, GetProfileNode, (const char*name, const UGProfileNode*)),
			grp, "a profile node", "name", "if root = null, return");
	reg.add_function("GetProfilerAvailable", &GetProfilerAvailable, grp,
	                 "true if profiler available");
	reg.add_function("SetOutputProfileStats", &UGOutputProfileStatsOnExit, grp,  
	                 "", "bOutput", "if set to true and profiler available, profile stats are printed at the end of the program. true is default");
	reg.add_function("WriteProfileData",
					OVERLOADED_FUNCTION_PTR(void, WriteProfileDataXML, (const char*)),
					 grp,
	                 "", "filename|save-dialog|endings=[\"pdxml\"]", "writes a XML-file with profile data viewable with the ShinyProfileViewer. Pick a filename ending with .pdxml");
	reg.add_function("WriteProfileData",
					OVERLOADED_FUNCTION_PTR(void, WriteProfileDataXML, (const char*, int)),
					 grp,
	                 "", "filename|save-dialog|endings=[\"pdxml\"]", "writes a XML-file with profile data viewable with the ShinyProfileViewer. Pick a filename ending with .pdxml");
	reg.add_function("WriteCallLog",
						OVERLOADED_FUNCTION_PTR(void, WriteCallLog, (const char*)),
						 grp,
		                 "", "filename|save-dialog|endings=[\"txt\"]", "writes txt file with call log");
	reg.add_function("WriteCallLog",
					OVERLOADED_FUNCTION_PTR(void, WriteCallLog, (const char*, int)),
					 grp,
	                 "", "filename|save-dialog|endings=[\"txt\"]", "writes txt file with call log");

	reg.add_function("UpdateProfiler", &UpdateProfiler_BridgeImpl, grp);

	reg.add_function("SetShinyCallLoggingMaxFrequency", &SetShinyCallLoggingMaxFrequency, grp, "", "maxFreq");

}


// end group profiler_bridge
/// \}

} // namespace bridge

} // namespace ug
