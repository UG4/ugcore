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

using namespace std;

namespace ug
{
  //void PrintLUA();
namespace bridge
{

//////////////////////////////////////////////////////////////////////////////////////////

void RegisterBridge_Profiler(Registry &reg, string parentGroup)
{
	stringstream ss; ss << parentGroup << "/Util/Profiler";
	string grp = ss.str();

	reg.add_class_<UGProfileNode>("UGProfileNode", grp)
		.add_method("call_tree", static_cast<string(UGProfileNode::*)() const>(&UGProfileNode::call_tree), "string with call tree")
		.add_method("call_tree",
				static_cast<string(UGProfileNode::*)(double dSkipMarginal) const>(&UGProfileNode::call_tree), "string with call tree",
				"dSkipMarginal")

		.add_method("child_self_time_sorted",
				static_cast<string(UGProfileNode::*)() const>(&UGProfileNode::child_self_time_sorted),
				"string with sorted childs", "", "childs are sorted by self time")
		.add_method("child_self_time_sorted",
				static_cast<string(UGProfileNode::*)(double dSkipMarginal) const>(&UGProfileNode::child_self_time_sorted),
				"string with sorted childs", "dSkipMarginal", "childs are sorted by self time")

		.add_method("total_time_sorted",
				static_cast<string(UGProfileNode::*)() const>(&UGProfileNode::total_time_sorted),
				"string with sorted childs", "", "childs are sorted by total time")
		.add_method("total_time_sorted",
				static_cast<string(UGProfileNode::*)(double dSkipMarginal) const>(&UGProfileNode::total_time_sorted),
				"string with sorted childs", "dSkipMarginal", "childs are sorted by total time")

		.add_method("entry_count_sorted",
				static_cast<string(UGProfileNode::*)() const>(&UGProfileNode::entry_count_sorted),
				"string with sorted childs", "", "childs are sorted by entry count")
		.add_method("entry_count_sorted",
				static_cast<string(UGProfileNode::*)(double dSkipMarginal) const>(&UGProfileNode::entry_count_sorted),
				"string with sorted childs", "dSkipMarginal", "childs are sorted by entry count")

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

	reg.add_function("GetProfileNode", &GetProfileNode, grp);
	reg.add_function("WriteProfileData", &WriteProfileData, grp);
	reg.add_function("GetProfilerAvailable", &GetProfilerAvailable, grp, "true if profiler available");
	reg.add_function("SetOutputProfileStats", &UGOutputProfileStatsOnExit, grp, "", "bOutput",
			"if set to true and profiler available, profile stats are printed at the end of the program. true is default");
}


} // namespace bridge

} // namespace ug
