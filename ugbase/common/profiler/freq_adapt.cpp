/*
 * freq_adapt.cpp
 *
 *  Created on: 10.12.2014
 *      Author: andreasvogel
 */

#include <cpufreq.h>
#include "freq_adapt.h"
#include <iostream>
#include <fstream>
#include <cstdlib>	//used for atoi. Better: boost::lexical_cast
#include "common/common.h"
#include "common/util/string_util.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// FreqAdaptNodeManager
////////////////////////////////////////////////////////////////////////////////

void FreqAdaptNodeManager::add(AutoFreqAdaptNode* node)
{
	inst().m_nodes.push(node);
}

void FreqAdaptNodeManager::release_latest()
{
	if(!inst().m_nodes.empty()){
		AutoFreqAdaptNode* node = inst().m_nodes.top();
		inst().m_nodes.pop();
		node->release();
	}
}

FreqAdaptNodeManager::FreqAdaptNodeManager() {}

FreqAdaptNodeManager::~FreqAdaptNodeManager()
{
	//	release and deactivate all nodes
		while(!inst().m_nodes.empty())
			inst().release_latest();
}

FreqAdaptNodeManager& FreqAdaptNodeManager::inst()
{
	static FreqAdaptNodeManager pnm;
	return pnm;
}


////////////////////////////////////////////////////////////////////////////////
// AutoFreqAdaptNode
////////////////////////////////////////////////////////////////////////////////

AutoFreqAdaptNode::AutoFreqAdaptNode(unsigned long freq) : m_bActive(true) {
	// remember for auto-release
	FreqAdaptNodeManager::add(this);

	if(freq){
		cout << "CPU_FREQ: Changing frequency to " << freq << "\n";

		// remember curr freq
		if ( (m_prevFreq = cpufreq_get_freq_kernel(0)) == 0)
			UG_THROW("Error while getting frequency");

		// set new freq \todo: adjust to different architectures
		if (cpufreq_modify_policy_governor(0, "userspace") != 0)
			UG_THROW("Error while setting governor");
		if (cpufreq_set_frequency(0, freq) != 0)
			UG_LOG("error while setting frequency");
	}
	else{
		// no freq-change issued -> no change back
		m_prevFreq = 0;
	}
}

AutoFreqAdaptNode::~AutoFreqAdaptNode()
{
	if(m_bActive){
		FreqAdaptNodeManager::release_latest();
	}
}

void AutoFreqAdaptNode::release(){

	// check if freq was changed on entry
	if(m_prevFreq){

		cout << "CPU_FREQ: Resetting frequency to " << m_prevFreq << "\n";

		// set previous freq \todo: adjust to different architectures
		if (cpufreq_modify_policy_governor(0, "userspace") != 0)
			UG_THROW("Error while setting governor");
		if (cpufreq_set_frequency(0, m_prevFreq) != 0)
			UG_LOG("error while setting frequency");
	}
}


////////////////////////////////////////////////////////////////////////////////
// FreqAdaptValues
////////////////////////////////////////////////////////////////////////////////

FreqAdaptValues& FreqAdaptValues::inst()
{
	static FreqAdaptValues myInst;
	return myInst;
};

unsigned long FreqAdaptValues::find_freq(const char* file, const int line){

	cout << "CPU_FREQ: Node for: "<<file<<", line: " << line << endl;

	for(size_t i = 0; i < m_pos.size(); ++i){
		if(m_pos[i].line == line && m_pos[i].file == file){
			return m_pos[i].freq;
		}
	}

	return 0;
}

void FreqAdaptValues::set_freqs(std::string csvFile){
	ifstream in(csvFile.c_str());
	if(!in){
		UG_THROW("FreqAdaptValues::set_freqs: Couldn't open: '" << csvFile << "'\n");
	}

	m_pos.clear();
	cout << "FreqAdaptValues: parsing file: "<<csvFile << "\n";
	string line;
	while(!in.eof()){
		getline(in, line);

		size_t i2 = line.rfind(",");
		if(i2 == string::npos)
			continue;

		size_t i1 = line.rfind(",", i2 - 1);
		if(i1 == string::npos)
			continue;

		string file = ug::TrimString(line.substr(0, i1));
		string s2 = line.substr(i1 + 1, i2 - 1 - i1);
		string s3 = line.substr(i2 + 1, line.size() - i1 + 1);

		int line = atoi(s2.c_str());
		unsigned long freq = atoi(s3.c_str());
	//	instead of atoi, better use boost::lexical_cast
		// int lineNumber = boost::lexical_cast<int>(s2);
		// int frequency = boost::lexical_cast<int>(s3);

		cout << file << "===" << line << "===" << freq << "===" << endl;

		m_pos.push_back( FreqAdaptPoint(file,line,freq) );
	}
}


std::vector<FreqAdaptValues::FreqAdaptPoint> FreqAdaptValues::m_pos =
		std::vector<FreqAdaptValues::FreqAdaptPoint>();
