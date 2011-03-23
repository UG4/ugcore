/*
 * profiler_bridge.cpp
 *
 *  Created on: 11.02.2011
 *      Author: Martin Rupp
 */


#include "ug_script/ug_script.h"
#include "../registry.h"
#include "../ug_bridge.h"

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
	/// \return number of entries in this profiler node
	double get_avg_entry_count() const
	{
		if(!is_valid()) return 0.0;
		return data.entryCount.avg * SHINY_DAMPING_FACTOR;
	}

	/// \return time in milliseconds spend in this node excluding subnodes
	double get_avg_self_time_ms() const
	{
		if(!is_valid()) return 0.0;
		return data.selfTicks.avg / 1000.0 * SHINY_DAMPING_FACTOR;
	}

	/// \return time in milliseconds spend in this node including subnodes
	double get_avg_total_time_ms() const
	{
		if(!is_valid()) return 0.0;
		return data.totalTicksAvg() / 1000.0 * SHINY_DAMPING_FACTOR;
	}

	/// \return true if node has been found
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

bool GetProfilerAvailable()
{
	return true;
}

#else

// dummy profiler node
class UGProfilerNode
{
public:
	/// \return number of entries in this profiler node
	double get_avg_entry_count() const
	{
		return 0;
	}

	/// \return time in milliseconds spend in this node excluding subnodes
	double get_avg_self_time_ms() const
	{
		return 0.0;
	}

	/// \return time in milliseconds spend in this node including subnodes
	double get_avg_total_time_ms() const
	{
		return 0.0;
	}

	/// \return true if node has been found
	bool is_valid() const
	{
		return false;
	}

	/*const char * tostring()
	{
		return "hello world";
	}

	void unm()
	{
		UG_LOG("unm!\n");
	}

	void add(const UGProfilerNode *other)
	{
		UG_LOG("oha oha!\n");
	}*/
};


const UGProfilerNode *GetProfileNode(const char *name)
{
	return NULL;
}

bool GetProfilerAvailable()
{
	return false;
}

#endif


bool RegisterProfileFunctions(Registry &reg, const char* parentGroup)
{
	std::stringstream group; group << parentGroup << "/Profiler";

	reg.add_class_<UGProfilerNode>("UGProfilerNode", group.str().c_str())
		.add_method("get_avg_entry_count", &UGProfilerNode::get_avg_entry_count,
				"number of entries in this profiler node", "")
		.add_method("get_avg_self_time_ms", &UGProfilerNode::get_avg_self_time_ms,
				"time in milliseconds spend in this node excluding subnodes", "")
		.add_method("get_avg_total_time_ms", &UGProfilerNode::get_avg_total_time_ms,
				"time in milliseconds spend in this node including subnodes", "")
		.add_method("is_valid", &UGProfilerNode::is_valid, "true if node has been found", "");
		/*.add_method("__tostring", &UGProfilerNode::tostring, "tostring")
		.add_method("__unm", &UGProfilerNode::unm, "unm")
		.add_method("__add", &UGProfilerNode::add, "add");*/
	reg.add_function("GetProfileNode", &GetProfileNode, group.str().c_str());
	reg.add_function("GetProfilerAvailable", &GetProfilerAvailable, group.str().c_str(), "true if profiler available");

	return true;
}


}

}
