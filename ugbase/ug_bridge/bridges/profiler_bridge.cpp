/*
 * profiler_bridge.cpp
 *
 *  Created on: 11.02.2011
 *      Author: Martin Rupp
 */


#include "ug.h"
#include "ug_script/ug_script.h"

namespace ug
{
namespace bridge
{
#if SHINY_PROFILER

// note: for some really strange reason, shiny multiplies every time by 0.9 when you call PROFILER_UPDATE
// and since update(0.9) is called at least once at the end of UGFinalize, we need to compensate for that
// (WE do call update with damping = 1.0 of course)
#define SHINY_DAMPING_FACTOR 0.9

class UGProfilerNode : public Shiny::ProfileNode
{
public:
	double get_avg_entry_count() const
	{
		return data.entryCount.avg * SHINY_DAMPING_FACTOR;
	}

	double get_avg_self_time_ms() const
	{
		return data.selfTicks.avg / 1000.0 * SHINY_DAMPING_FACTOR;
	}

	double get_avg_total_time_ms() const
	{
		return data.totalTicksAvg() / 1000.0 * SHINY_DAMPING_FACTOR;
	}

	bool is_valid() const
	{
		return this != NULL;
	}
};


const UGProfilerNode *GetProfileNode(const char *name)
{
	Shiny::ProfileManager::instance.update(1.0); // WE call with damping = 1.0

	const Shiny::ProfileNode *node = &Shiny::ProfileManager::instance.rootNode;
	do
	{
		if(strcmp(node->zone->name, name) == 0)
			return reinterpret_cast<const UGProfilerNode*> (node);
		node = node->findNextInTree();
	} while (node);

	UG_LOG("Profiler Node \"" << name << "\" not found\n");
	return NULL;
}

bool RegisterProfileFunctions(Registry &reg, const char* parentGroup)
{
	std::stringstream group; group << parentGroup << "/Profiler";

	reg.add_class_<UGProfilerNode>("UGProfilerNode", group.str().c_str())
		.add_method("get_avg_entry_count", &UGProfilerNode::get_avg_entry_count, "", "")
		.add_method("get_avg_self_time_ms", &UGProfilerNode::get_avg_self_time_ms, "", "")
		.add_method("get_avg_total_time_ms", &UGProfilerNode::get_avg_total_time_ms, "", "")
		.add_method("is_valid", &UGProfilerNode::is_valid, "", "");
	reg.add_function("GetProfileNode", &GetProfileNode);

	return true;
}

#else
bool AddProfileFunctionsToRegistry(Registry &reg)
{
	return true;
}
#endif

}

}
