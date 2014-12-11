/*
 * freq_adapt.cpp
 *
 *  Created on: 10.12.2014
 *      Author: andreasvogel
 */

#include <cpufreq.h>
#include "freq_adapt.cpp"

////////////////////////////////////////////////////////////////////////////////
// AutoFreqAdaptNode
////////////////////////////////////////////////////////////////////////////////

void AutoFreqAdaptNode::add(AutoFreqAdaptNode* node)
{
	inst().m_nodes.push(node);
}

void AutoFreqAdaptNode::release_latest()
{
	if(!inst().m_nodes.empty()){
		AutoFreqAdaptNode* node = inst().m_nodes.top();
		inst().m_nodes.pop();
		node->release();
	}
}

AutoFreqAdaptNode::FreqAdaptNodeManager() {}

AutoFreqAdaptNode::~FreqAdaptNodeManager()
{
	//	release and deactivate all nodes
		while(!inst().m_nodes.empty())
			inst().release_latest();
}

static FreqAdaptNodeManager& AutoFreqAdaptNode::inst()
{
	static FreqAdaptNodeManager pnm;
	return pnm;
}


////////////////////////////////////////////////////////////////////////////////
// AutoFreqAdaptNode
////////////////////////////////////////////////////////////////////////////////

AutoFreqAdaptNode::AutoFreqAdaptNode(unsigned long freq) : m_bActive(true) {
	FreqAdaptNodeManager::add(this);

	if(freq){
		// remember curr freq
		if ( (m_prevFreq = cpufreq_get_freq_kernel(0)) == 0)
			UG_THROW("Error while getting frequency");

		// set new freq \todo: adjust to different architectures
		for (int i=0; i<40; i++) {
			if (cpufreq_modify_policy_governor(i, "userspace") != 0)
				UG_THROW("Error while setting governor");

			if (cpufreq_set_frequency(i, freq) != 0)
				UG_LOG("error while setting frequency");
		}
	}
	else{
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
	if(m_prevFreq){
		// set previous freq \todo: adjust to different architectures
		for (int i=0; i<40; i++) {
			if (cpufreq_modify_policy_governor(i, "userspace") != 0)
				UG_THROW("Error while setting governor");

			if (cpufreq_set_frequency(i, m_prevFreq) != 0)
				UG_LOG("error while setting frequency");
		}
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
	for(int i = 0; i < m_pos; ++i){
		if(m_pos[i].line == line && m_pos[i].file == file){
			return m_pos[i].freq;
		}
	}

	return 0;
}

static void FreqAdaptValues::read_marks(){
	// todo: Add code here
}


std::vector<FreqAdaptValues::FreqAdaptPoint> FreqAdaptValues::m_pos =
		std::vector<FreqAdaptValues::FreqAdaptPoint>();
