/*
 * profiler_bridge.cpp
 *
 *  Created on: 11.02.2011
 *      Author: Martin Rupp
 */


#include "ug_script/ug_script.h"
#include "../registry.h"
#include "../ug_bridge.h"
#include "common/profiler/profiler.h"
#include <iomanip>
#include <math.h> // für floor

namespace ug
{
namespace bridge
{
#if SHINY_PROFILER

std::string cut(std::string &s, size_t L)
{
	return s.substr(0, L);
}

std::string cut(const char *p, size_t L)
{
	std::string s(p);
	return s.substr(0, L);
}

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
		return data.entryCount.avg; // * SHINY_DAMPING_FACTOR;
	}

	/// \return time in milliseconds spend in this node excluding subnodes
	double get_avg_self_time_ms() const
	{
		return get_avg_self_time() / 1000.0;
	}

	/// \return time in milliseconds spend in this node including subnodes
	double get_avg_total_time_ms() const
	{
		return get_avg_total_time() / 1000.0;
	}


	/// \return time in seconds spend in this node excluding subnodes
	double get_avg_self_time() const
	{
		if(!is_valid()) return 0.0;
		return data.selfTicks.avg * SHINY_DAMPING_FACTOR;
	}

	/// \return time in seconds spend in this node including subnodes
	double get_avg_total_time() const
	{
		if(!is_valid()) return 0.0;
		return data.totalTicksAvg() * SHINY_DAMPING_FACTOR;
	}




	std::string print_node(double full, size_t offset=0) const
	{
		if(!is_valid()) return "";
		double totalTicksAvg = get_avg_total_time();
		const Shiny::TimeUnit *selfUnit = Shiny::GetTimeUnit(get_avg_self_time());
		const Shiny::TimeUnit *totalUnit = Shiny::GetTimeUnit(totalTicksAvg);
		std::stringstream s;
		if(offset)	s << std::setw(offset) << " ";

		s <<	std::left << std::setw(Shiny::OUTPUT_WIDTH_NAME-offset) << cut(zone->name, Shiny::OUTPUT_WIDTH_NAME-offset) <<
				std::right << std::setw(Shiny::OUTPUT_WIDTH_HIT) << floor(get_avg_entry_count()) << " " <<
				std::setprecision(Shiny::OUTPUT_WIDTH_TIME-1) <<
				std::setw(Shiny::OUTPUT_WIDTH_TIME) << get_avg_self_time() * selfUnit->invTickFreq << " " << selfUnit->suffix << " " <<
				std::setw(Shiny::OUTPUT_WIDTH_PERC) << floor(get_avg_self_time() / full * 100) << "% " <<
				std::setw(Shiny::OUTPUT_WIDTH_TIME) << totalTicksAvg * totalUnit->invTickFreq << " " << totalUnit->suffix << " " <<
				std::setw(Shiny::OUTPUT_WIDTH_PERC) << floor(totalTicksAvg / full * 100) << "% ";
		return s.str();
	}

	const UGProfilerNode *get_first_child() const
	{
		return reinterpret_cast<const UGProfilerNode*>(firstChild);
	}

	const UGProfilerNode *get_last_child() const
	{
		return reinterpret_cast<const UGProfilerNode*>(lastChild);
	}

	const UGProfilerNode *get_next_sibling() const
	{
		return reinterpret_cast<const UGProfilerNode*>(nextSibling);
	}

	std::string call_tree() const
	{
		if(!is_valid()) return "Profile Node not valid!";

		std::stringstream s;
		UGProfilerNode::log_header(s, "call tree");

		rec_print(get_avg_total_time(), s);

		return s.str();
	}

	std::string child_self_time_sorted() const
	{
		return child_sorted("self time sorted", UGProfilerNode::self_time_sort);
	}

	std::string total_time_sorted() const
	{
		return child_sorted("total time sorted", UGProfilerNode::total_time_sort);
	}

	std::string entry_count_sorted() const
	{
		return child_sorted("entry count sorted", UGProfilerNode::entry_count_sort);
	}

	// \return true if node has been found
	bool is_valid() const
	{
		return this != NULL;
	}

private:
	void rec_print(double full, std::stringstream &s, size_t offset=0) const
	{
		if(!is_valid()) return;
		s << print_node(full, offset) << "\n";
		for(const UGProfilerNode *p=get_first_child(); p != NULL; p=p->get_next_sibling())
		{
			p->rec_print(full, s, offset+1);
			if(p==get_last_child())
				break;
		}
	}

	void add_nodes(std::vector<const UGProfilerNode*> &nodes) const
	{
		nodes.push_back(this);
		for(const UGProfilerNode *p=get_first_child(); p != NULL; p=p->get_next_sibling())
		{
			p->add_nodes(nodes);
			if(p==get_last_child())
				break;
		}
	}

	std::string child_sorted(const char *name, bool sortFunction(const UGProfilerNode *a, const UGProfilerNode *b)) const
	{
		if(!is_valid()) return "";
		std::stringstream s;
		std::vector<const UGProfilerNode*> nodes;
		add_nodes(nodes);
		sort(nodes.begin(), nodes.end(), sortFunction);

		UGProfilerNode::log_header(s, name);
		for(size_t i=0; i<nodes.size(); i++)
			s << nodes[i]->print_node(get_avg_total_time()) << "\n";
		return s.str();
	}

	static void log_header(std::stringstream &s, const char *name)
	{

		s << 	std::left << std::setw(Shiny::OUTPUT_WIDTH_NAME) << name << " " <<
				std::right << std::setw(Shiny::OUTPUT_WIDTH_HIT) << "hits" << " " <<
				std::setw(Shiny::OUTPUT_WIDTH_TIME+4+Shiny::OUTPUT_WIDTH_PERC+1) << "self time" << " " <<
				std::setw(Shiny::OUTPUT_WIDTH_TIME+4+Shiny::OUTPUT_WIDTH_PERC+1) << "total time"  << " \n";
	}

	static bool self_time_sort(const UGProfilerNode *a, const UGProfilerNode *b)
	{
		return a->get_avg_self_time() < b->get_avg_self_time();
	}

	static bool total_time_sort(const UGProfilerNode *a, const UGProfilerNode *b)
	{
		return a->get_avg_total_time() < b->get_avg_total_time();
	}

	static bool entry_count_sort(const UGProfilerNode *a, const UGProfilerNode *b)
	{
		return a->get_avg_entry_count() < b->get_avg_entry_count();
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

	std::string call_tree()
	{
		return "Profiler not availabel!";
	}

	std::string child_self_time_sorted() const
	{
		return "Profiler not availabel!";
	}

	std::string total_time_sorted() const
	{
		return "Profiler not availabel!";
	}

	std::string entry_count_sorted() const
	{
		return "Profiler not availabel!";
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
		.add_method("call_tree", &UGProfilerNode::call_tree, "string with call tree")
		.add_method("child_self_time_sorted", &UGProfilerNode::child_self_time_sorted, "string with sorted childs", "", "childs are sorted by self time")
		.add_method("total_time_sorted", &UGProfilerNode::total_time_sorted, "string with sorted childs", "", "childs are sorted by total time")
		.add_method("entry_count_sorted", &UGProfilerNode::entry_count_sorted, "string with sorted childs", "", "childs are sorted by entry count")
		.add_method("get_avg_entry_count", &UGProfilerNode::get_avg_entry_count,
				"number of entries in this profiler node", "")
		.add_method("get_avg_self_time_ms", &UGProfilerNode::get_avg_self_time_ms,
				"time in milliseconds spend in this node excluding subnodes", "")
		.add_method("get_avg_total_time_ms", &UGProfilerNode::get_avg_total_time_ms,
				"time in milliseconds spend in this node including subnodes", "")
		.add_method("is_valid", &UGProfilerNode::is_valid, "true if node has been found", "")

		;
		/*.add_method("__tostring", &UGProfilerNode::tostring, "tostring")
		.add_method("__unm", &UGProfilerNode::unm, "unm")
		.add_method("__add", &UGProfilerNode::add, "add");*/
	reg.add_function("GetProfileNode", &GetProfileNode, group.str().c_str());
	reg.add_function("GetProfilerAvailable", &GetProfilerAvailable, group.str().c_str(), "true if profiler available");

	return true;
}


}

}
