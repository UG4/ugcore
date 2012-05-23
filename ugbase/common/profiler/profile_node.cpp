/*
 * profile_node.cpp
 *
 *  Created on: May 22, 2012
 *      Author: Martin Rupp
 */


#include "profiler.h"
#include <iomanip>
#include <cmath> // for floor
#include <string>
#include "common/log.h"
#include "profile_node.h"

using namespace std;


namespace ug
{

#if SHINY_PROFILER

static const int PROFILER_BRIDGE_OUTPUT_WIDTH_NAME = 50; // Shiny::OUTPUT_WIDTH_NAME
static const int PROFILER_BRIDGE_OUTPUT_WIDTH_HIT  = 13; // Shiny::OUTPUT_WIDTH_HIT
static const int PROFILER_BRIDGE_OUTPUT_WIDTH_TIME =  7; // Shiny::OUTPUT_WIDTH_TIME
static const int PROFILER_BRIDGE_OUTPUT_WIDTH_PERC =  4; // Shiny::OUTPUT_WIDTH_PERC


/*static string cut(string &s, size_t L)
{
	return s.substr(0, L);
}*/

static string cut(const char *p, size_t L)
{
	string s(p);
	return s.substr(0, L);
}


// note: for some really strange reason, shiny multiplies every time by 0.9 when you call PROFILER_UPDATE
// and since update(0.9) is called at least once at the end of UGFinalize, we need to compensate for that
// (WE do call update with damping = 1.0 of course)
#define SHINY_DAMPING_FACTOR 0.9

double UGProfileNode::get_avg_entry_count() const
{
	if(!valid()) return 0.0;
	return data.entryCount.avg; // * SHINY_DAMPING_FACTOR;
}

double UGProfileNode::get_avg_self_time_ms() const
{
	return get_avg_self_time() / 1000.0;
}

double UGProfileNode::get_avg_total_time_ms() const
{
	return get_avg_total_time() / 1000.0;
}

double UGProfileNode::get_avg_self_time() const
{
	if(!valid()) return 0.0;
	return data.selfTicks.avg * SHINY_DAMPING_FACTOR;
}

double UGProfileNode::get_avg_total_time() const
{
	if(!valid()) return 0.0;
	return data.totalTicksAvg() * SHINY_DAMPING_FACTOR;
}

string UGProfileNode::call_tree(double dSkipMarginal) const
{
	if(!valid()) return "Profile Node not valid!";

	stringstream s;
	UGProfileNode::log_header(s, "call tree");

	rec_print(get_avg_total_time(), s, 0, dSkipMarginal);

	return s.str();
}

string UGProfileNode::call_tree() const
{
	return call_tree(0.0);
}

string UGProfileNode::child_self_time_sorted(double dSkipMarginal) const
{
	return child_sorted("self time sorted", UGProfileNode::self_time_sort, dSkipMarginal);
}
string UGProfileNode::child_self_time_sorted() const
{
	return child_self_time_sorted(0.0);
}

string UGProfileNode::total_time_sorted(double dSkipMarginal) const
{
	return child_sorted("total time sorted", UGProfileNode::total_time_sort, dSkipMarginal);
}
string UGProfileNode::total_time_sorted() const
{
	return total_time_sorted(0.0);
}

string UGProfileNode::entry_count_sorted(double dSkipMarginal) const
{
	return child_sorted("entry count sorted", UGProfileNode::entry_count_sort, dSkipMarginal);
}
string UGProfileNode::entry_count_sorted() const
{
	return entry_count_sorted(0.0);
}

bool UGProfileNode::valid() const
{
	return this != NULL;
}

// private functions


string UGProfileNode::print_node(double full, size_t offset) const
{
	if(!valid()) return "";
	double totalTicksAvg = get_avg_total_time();
	const Shiny::TimeUnit *selfUnit = Shiny::GetTimeUnit(get_avg_self_time());
	const Shiny::TimeUnit *totalUnit = Shiny::GetTimeUnit(totalTicksAvg);
	stringstream s;
	if(offset)	s << setw(offset) << " ";

	s <<	left << std::setw(PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset) << cut(zone->name, PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset) <<
			right << std::setw(PROFILER_BRIDGE_OUTPUT_WIDTH_HIT) << floor(get_avg_entry_count()) << " " <<
			setprecision(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME-1) <<
			setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME) << get_avg_self_time() * selfUnit->invTickFreq << " " <<
			left << setw(2) << selfUnit->suffix << " " <<
			right << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_PERC) << floor(get_avg_self_time() / full * 100) << "%  " <<
			setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME) << totalTicksAvg * totalUnit->invTickFreq << " " <<
			left << setw(2) << totalUnit->suffix << " " <<
			right << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_PERC) << floor(totalTicksAvg / full * 100) << "%  ";
	return s.str();
}

const UGProfileNode *UGProfileNode::get_first_child() const
{
	return reinterpret_cast<const UGProfileNode*>(firstChild);
}

const UGProfileNode *UGProfileNode::get_last_child() const
{
	return reinterpret_cast<const UGProfileNode*>(lastChild);
}

const UGProfileNode *UGProfileNode::get_next_sibling() const
{
	return reinterpret_cast<const UGProfileNode*>(nextSibling);
}

void UGProfileNode::rec_print(double full, stringstream &s, size_t offset, double dSkipMarginal) const
{
	if(!valid()) return;
	if(dSkipMarginal==0.0 || full*dSkipMarginal < get_avg_total_time())
	{
		s << print_node(full, offset) << "\n";
		for(const UGProfileNode *p=get_first_child(); p != NULL; p=p->get_next_sibling())
		{
			p->rec_print(full, s, offset+1, dSkipMarginal);
			if(p==get_last_child())
				break;
		}
	}
}

void UGProfileNode::add_nodes(vector<const UGProfileNode*> &nodes) const
{
	nodes.push_back(this);
	for(const UGProfileNode *p=get_first_child(); p != NULL; p=p->get_next_sibling())
	{
		p->add_nodes(nodes);
		if(p==get_last_child())
			break;
	}
}

string UGProfileNode::child_sorted(const char *name, bool sortFunction(const UGProfileNode *a, const UGProfileNode *b),
		double dSkipMarginal) const
{
	if(!valid()) return "";
	stringstream s;
	vector<const UGProfileNode*> nodes;
	add_nodes(nodes);
	sort(nodes.begin(), nodes.end(), sortFunction);

	UGProfileNode::log_header(s, name);
	double full = get_avg_total_time();
	for(size_t i=0; i<nodes.size(); i++)
	{
		if(dSkipMarginal==0.0 || full*dSkipMarginal < nodes[i]->get_avg_total_time())
			s << nodes[i]->print_node(full) << "\n";
	}
	return s.str();
}

void UGProfileNode::log_header(stringstream &s, const char *name)
{

	s << 	left << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_NAME) << name << " " <<
			right << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_HIT) << "hits" << " " <<
			setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME+4+PROFILER_BRIDGE_OUTPUT_WIDTH_PERC+1 -4) << "self time" << " " <<
			setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME+4+PROFILER_BRIDGE_OUTPUT_WIDTH_PERC+1) << "total time"  << " \n";
}

bool UGProfileNode::self_time_sort(const UGProfileNode *a, const UGProfileNode *b)
{
	return a->get_avg_self_time() < b->get_avg_self_time();
}

bool UGProfileNode::total_time_sort(const UGProfileNode *a, const UGProfileNode *b)
{
	return a->get_avg_total_time() < b->get_avg_total_time();
}

bool UGProfileNode::entry_count_sort(const UGProfileNode *a, const UGProfileNode *b)
{
	return a->get_avg_entry_count() < b->get_avg_entry_count();
}


const UGProfileNode *GetProfileNode(const char *name)
{
	Shiny::ProfileManager::instance.update(1.0); // WE call with damping = 1.0

	const Shiny::ProfileNode *node = &Shiny::ProfileManager::instance.rootNode;
	do
	{
		if(strcmp(node->zone->name, name) == 0)
			return reinterpret_cast<const UGProfileNode*> (node);
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


double UGProfileNode::get_avg_entry_count() const
{
	return 0;
}

/// \return time in milliseconds spend in this node excluding subnodes
double UGProfileNode::get_avg_self_time_ms() const
{
	return 0.0;
}

/// \return time in milliseconds spend in this node including subnodes
double UGProfileNode::get_avg_total_time_ms() const
{
	return 0.0;
}

string UGProfileNode::call_tree(double dSkipMarginal) const
{
	return "Profiler not available!";
}

string UGProfileNode::call_tree() const
{
	return call_tree(0.0);
}

string UGProfileNode::child_self_time_sorted(double dSkipMarginal) const
{
	return "Profiler not available!";
}

string UGProfileNode::child_self_time_sorted() const
{
	return child_self_time_sorted(0.0);
}

string UGProfileNode::total_time_sorted(double dSkipMarginal) const
{
	return "Profiler not available!";
}

string UGProfileNode::total_time_sorted() const
{
	return total_time_sorted(0.0);
}

string UGProfileNode::entry_count_sorted(double dSkipMarginal) const
{
	return "Profiler not available!";
}

string UGProfileNode::entry_count_sorted() const
{
	return entry_count_sorted(0.0);
}

/// \return true if node has been found
bool UGProfileNode::valid() const
{
	return false;
}



const UGProfileNode *GetProfileNode(const char *name)
{
	return NULL;
}

bool GetProfilerAvailable()
{
	return false;
}

#endif // SHINY

} // namespace ug


//////////////////////////////////////////////////////////////////////////////////////////
