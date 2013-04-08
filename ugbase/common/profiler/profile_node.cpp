/*
 * profile_node.cpp
 *
 *  Created on: May 22, 2012
 *      Author: Martin Rupp
 */


#include "profiler.h"
#include <iomanip>
#include <cmath> // for floor
#include <algorithm>
#include <string>
#include <string.h>
#include "common/log.h"
#include "common/util/string_util.h"
#include "profile_node.h"
#include "common/util/string_util.h"
#include "common/util/path_provider.h"
#include <map>
#include <fstream>
#include "compile_info/compile_info.h"
#include "pcl/pcl_base.h"


#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

using namespace std;


namespace ug
{

#if SHINY_PROFILER

static const int PROFILER_BRIDGE_OUTPUT_WIDTH_NAME = 70; // Shiny::OUTPUT_WIDTH_NAME
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

string XMLStringEscape(string s)
{
	ReplaceAll(s, "&", "&amp;");
	ReplaceAll(s, "\"", "&quot;");
	ReplaceAll(s, "\'", "&apos;");
	ReplaceAll(s, "<", "&lt;");
	ReplaceAll(s, ">", "&gt;");
	return s;
}

void UGProfileNode::write_node(ostream &s) const
{
	if(!valid()) return;
	
	/*
	 * <node>
	 * name
	 * group
	 * file
	 * line
	 * hits
	 * self
	 * total	 * 
	 * </node>	 
	 */
	
	s << "<node>\n";
	if(zone->name)
	{
		if(strcmp(zone->name, "<root>") == 0)
			s << "<name>root</name>\n";
		else
			s << "<name>" << XMLStringEscape(zone->name) << "</name>\n";
	}
	if(zone->groups)
		s << "<groups>" << XMLStringEscape(zone->groups) << "</groups>\n";
	
	if(zone->file != NULL)
	{
		s << "<file>";
		const string &ug4root = PathProvider::get_path(ROOT_PATH);
		if(StartsWith(zone->file, ug4root))
			s << "$" << (zone->file+ug4root.length());
		else
			s << zone->file;
		s << "</file>\n<line>" << zone->line << "</line>\n";
	}
	 
	s << "<hits>" << floor(get_avg_entry_count()) << "</hits>\n"
	  << "<self>" << get_avg_self_time() << "</self>\n"
	  << "<total>" << get_avg_total_time() << "</total>\n";
			
	for(const UGProfileNode *p=get_first_child(); p != NULL; p=p->get_next_sibling())
	{
		p->write_node(s);
		if(p==get_last_child())
			break;
	}

	s << "</node>\n";
}


string UGProfileNode::print_node(double full, size_t offset) const
{
	if(!valid()) return "";
	double totalTicksAvg = get_avg_total_time();
	const Shiny::TimeUnit *selfUnit = Shiny::GetTimeUnit(get_avg_self_time());
	const Shiny::TimeUnit *totalUnit = Shiny::GetTimeUnit(totalTicksAvg);
	stringstream s;
	if(offset)	s << setw(offset) << " ";

	if(zone->name[0] == '@')
	{
		if(strlen(zone->name) > (PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset) )
			s << "@... " << zone->name+strlen(zone->name)-(PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset-5);
		else
			s << left << std::setw(PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset) << cut(zone->name, PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset);

		const char *name = zone->name+1;
		const char *p = strchr(name, ':'); // search for line number
		if(!(p == NULL || p[0] == 0x00 || p[1] == 0x00))
		{
			int line = strtol(p+1, NULL, 10);
			char file[255];
			strncpy(file, name, p-name);
			file[p-name]=0x00;
			string str = GetFileLine(file, line);
			for(size_t i=0; i<str.size(); i++) if(str[i] == '\t') str[i] = ' ';
			s << "\n";
			if(offset)	s << setw(offset) << " ";
			s << left << std::setw(PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset) << cut(str.c_str(), PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset);
		}
	}
	else
		s << left << std::setw(PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset) << cut(zone->name, PROFILER_BRIDGE_OUTPUT_WIDTH_NAME-offset);

	s <<	right << std::setw(PROFILER_BRIDGE_OUTPUT_WIDTH_HIT) << floor(get_avg_entry_count()) << " " <<
			setprecision(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME-1) <<
			setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME) << get_avg_self_time() * selfUnit->invTickFreq << " " <<
			left << setw(2) << selfUnit->suffix << " " <<
			right << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_PERC) << floor(get_avg_self_time() / full * 100) << "%  " <<
			setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME) << totalTicksAvg * totalUnit->invTickFreq << " " <<
			left << setw(2) << totalUnit->suffix << " " <<
			right << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_PERC) << floor(totalTicksAvg / full * 100) << "%  "
			<< zone->groups;
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

string UGProfileNode::groups() const
{
	vector<const UGProfileNode*> nodes;
	add_nodes(nodes);

	map<string, double> mapGroups;
	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes[i]->zone->groups == NULL) continue;
		vector<string> g;
		TokenizeString(nodes[i]->zone->groups, g, ' ');
		for(size_t j=0; j<g.size(); j++)
			mapGroups[g[j]] += nodes[i]->get_avg_self_time();
	}

	vector<string> gs;
#ifdef UG_PARALLEL
	if(pcl::GetProcRank() == 0)
#endif
	for(map<string, double>::iterator it = mapGroups.begin(); it != mapGroups.end();++it)
		gs.push_back(it->first);

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator pc;
	pc.broadcast(gs);
	vector<double> t(gs.size(), 0.0), tMax, tMin;
	for(size_t i=0; i<gs.size(); i++)
		t[i] = mapGroups[gs[i]];
	pc.allreduce(t, tMax, PCL_RO_MAX);
	pc.allreduce(t, tMin, PCL_RO_MIN);
#endif
	stringstream s;
	for(size_t i=0; i<gs.size(); i++)
	{
		string name = gs[i];
		double time = mapGroups[name];
		const Shiny::TimeUnit *unit = Shiny::GetTimeUnit(time);
		s << left << std::setw(20) << name
		  << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME) << time * unit->invTickFreq << " " <<
			left << setw(2) << unit->suffix;
#ifdef UG_PARALLEL
		double maxTime = tMax[i];
		double minTime = tMin[i];
		double diffTime = maxTime - minTime;
		const Shiny::TimeUnit *maxUnit  = Shiny::GetTimeUnit(maxTime);
		const Shiny::TimeUnit *minUnit  = Shiny::GetTimeUnit(minTime);
		const Shiny::TimeUnit *diffUnit = Shiny::GetTimeUnit(diffTime);
		s << left << " max: " << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME)
				<< maxTime * maxUnit->invTickFreq << " " <<	left << setw(2) << maxUnit->suffix;
		s << left << " min: " << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME)
						<< minTime * maxUnit->invTickFreq << " " <<	left << setw(2) << minUnit->suffix;
		s << left << " diff: " << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME)
								<< diffTime * diffUnit->invTickFreq << " " <<	left << setw(2) << diffUnit->suffix
								<< " (" << diffTime/maxTime*100 << " %)";
#endif
		s << "\n";
	}
	return s.str();
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

/*
void WriteCompressedProfileData(const char *filename)
{
	// compressed:
	// filedata: file0, file1, file2, file3
	// zonedata: zone1 {name, fileid, line}, zone2, zone3	
    // nodes: node0 {self, total, children}, node1, node2, node3.
 
    // combine file and zone data, so that profile data is
 * 
 * // <files>...</files>
 * <zones><zone id="732"><name>bla</name><group>algebra</group><file>8</file><line>344</line> </zone> ... </zones>
 * <core id=0>
 *   <node zoneId=732><self>8734.2</self><total>8333.2</total></node>
    
}
*/

void WriteProfileData(const char *filename)
{
	
	Shiny::ProfileManager::instance.update(1.0); // WE call with damping = 1.0
	const Shiny::ProfileNode *node = &Shiny::ProfileManager::instance.rootNode;
	const UGProfileNode *pnRoot = reinterpret_cast<const UGProfileNode*> (node);
#ifdef UG_PARALLEL
	bool bProfileAll = false; // for the moment
	typedef pcl::SingleLevelLayout<pcl::OrderedInterface<size_t, std::vector> >
		IndexLayout;

	pcl::InterfaceCommunicator<IndexLayout> ic;
	if(pcl::GetProcRank()==0)
	{
#endif
		fstream f(filename, ios::out);
		f << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
		f << "<!-- ug4 created profile data -->\n";
		f << "<ProfileData>\n";
		f << "<AdditionalInfo>\n";
		f << "<SVNRevision>" << UGSvnRevision() << "</SVNRevision>\n";
		f << "<BuildHostname>" << UGBuildHost() << "</BuildHostname>\n";
		f << "<CompileDate>" << UGCompileDate() << "</CompileDate>\n";
		f << "</AdditionalInfo>\n";
		f << "<basepath>" << PathProvider::get_path(ROOT_PATH) << "</basepath>\n";

		f << "<core id=\"0\">\n";
		
		pnRoot->write_node(f);
		f << "</core>\n";	
		
#ifdef UG_PARALLEL
		if(bProfileAll)
		{
			std::vector<ug::BinaryBuffer> buffers(pcl::GetNumProcesses()-1);
			for(int i=1; i<pcl::GetNumProcesses(); i++)
				ic.receive_raw(i, buffers[i-1]);
			ic.communicate();

			for(int i=1; i<pcl::GetNumProcesses(); i++)
			{
				f << "\n<core id=\"" << i << "\">";
				string s;
				Deserialize(buffers[i-1], s);
				f << s;			
				f << "\n</core>";
			}		
		}
#endif
		f << "</ProfileData>\n";
		
#ifdef UG_PARALLEL		
	}	
	else		
	{
		if(bProfileAll)
		{
			stringstream ss;
			pnRoot->write_node(ss);
			BinaryBuffer buf;
			Serialize(buf, ss.str());
			ic.send_raw(0, buf.buffer(), buf.write_pos(), false);
			ic.communicate();
		}
	}	
#endif	
	
}

const UGProfileNode *GetProfileNode(const char *name)
{
	Shiny::ProfileManager::instance.update(1.0); // WE call with damping = 1.0
	const Shiny::ProfileNode *node = &Shiny::ProfileManager::instance.rootNode;
	if(name == NULL)
		return reinterpret_cast<const UGProfileNode*> (node);
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

void PrintLUA()
{
	const UGProfileNode *rootNode = GetProfileNode(NULL);
	vector<const UGProfileNode*> nodes;
	rootNode->add_nodes(nodes);
	double full = rootNode->get_avg_total_time();

	map<string, vector<double> > files;
	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes[i]->zone->groups == NULL ||
				strcmp(nodes[i]->zone->groups, "lua") != 0)
			continue;
		const char *name = nodes[i]->zone->name;
		cout << name << "\n";
		if(name[0]==0x00 || name[1]==0x00) continue;
		name++; // skip @
		const char *p = strchr(name, ':'); // search for line number
		if(p == NULL || p[0] == 0x00 || p[1] == 0x00) continue;
		int line = strtol(p+1, NULL, 10);
		if(line > 10000 || line < 0) continue;
		char file[255];
		strncpy(file, name, p-name);
		file[p-name]=0x00;
		vector<double> &v = files[file];
		if(v.size() < (size_t)line+1)
		{
			size_t s=v.size();
			v.resize(line+1);
			for(; s<(size_t)line+1; s++)
				v[s]=0.0;
		}
		v[line] = nodes[i]->get_avg_total_time();
	}


	for(map<string, vector<double> >::iterator it = files.begin(); it != files.end();
		++it)
	{
		string name = it->first;
		vector<double> &v = it->second;
		cout << "\n" << name << ":\n\n";

		char buf[512];
		fstream file(name.c_str(), ios::in);
		if(file.is_open() == false) continue;
		size_t lineNr=0;
		while(!file.eof())
		{
			file.getline(buf, 512);
			double time = lineNr >= v.size() ? 0.0 : v[lineNr];
			const Shiny::TimeUnit *unit = Shiny::GetTimeUnit(time);
			cout << std::resetiosflags( ::std::ios::scientific ) <<
					left << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_TIME)
					<< time * unit->invTickFreq << " " <<
				left << setw(2) << unit->suffix << " " <<
				right << setw(PROFILER_BRIDGE_OUTPUT_WIDTH_PERC) << floor(time / full * 100) << "%  "
				<< setw(3) << lineNr << "| " << buf << "\n";
			lineNr ++;
		}
	}
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

string UGProfileNode::groups() const
{
	return "Profiler not available!";
}

const UGProfileNode *GetProfileNode(const char *name)
{
	return NULL;
}

void WriteProfileData(const char *filename)
{
	return;
}

bool GetProfilerAvailable()
{
	return false;
}

void PrintLUA()
{
	UG_LOG("LUA Profiler not available.");
}

#endif // SHINY

} // namespace ug


//////////////////////////////////////////////////////////////////////////////////////////
