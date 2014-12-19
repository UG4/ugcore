/*
 * freq_adapt.cpp
 *
 *  Created on: 10.12.2014
 *      Author: andreasvogel
 */

#include <cpufreq.h>  // frequency adaption library
#include <iostream>   // usual i/o
#include <fstream>    // file i/o
#include <cstdlib>	  //used for atoi. Better: boost::lexical_cast
#include "freq_adapt.h"
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

		// set new frequency
		FreqAdaptValues::adjust_freq(freq);
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

		// set previous freq
		FreqAdaptValues::adjust_freq(freq);
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


void* FreqAdaptValues::freqAdaptWorker(void* This) {

	unsigned long freq;

	while (1) {

		// wait on frequency transition request
		pthread_mutex_lock(&freqAdapt_mutex);
		while (newFreq == 0) {
			pthread_cond_wait(&freqAdapt_condVar, &freqAdapt_mutex);
		}
		pthread_mutex_unlock(&freqAdapt_mutex);

		/* sleep some time to handle only the very last of a series of
                   requests as seen in case of multiple nested routines */
		usleep(10); // TODO: find an appropriate value

		pthread_mutex_lock(&freqAdapt_mutex);
		freq = newFreq;
		newFreq = 0;
		pthread_mutex_unlock(&freqAdapt_mutex);

		// do the actual frequency transition
		if (cpufreq_set_frequency(0, freq) != 0) {
			UG_THROW("FreqAdaptValues::freqAdaptWorker: Error while setting frequency");
		}

	}

}

void FreqAdaptValues::adjust_freq(unsigned long freq){
    pthread_mutex_lock(&freqAdapt_mutex);
    newFreq = freq;
    pthread_cond_signal(&freqAdapt_condVar);
    pthread_mutex_unlock(&freqAdapt_mutex);
}

void FreqAdaptValues::set_freqs(std::string csvFile){

	///////////////////////////
	// read file, line, freqs
	///////////////////////////
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


	///////////////////////////
	// start adaption thread
	///////////////////////////

	// TODO: check return values?

	newFreq = 0;

	pthread_mutex_init(&freqAdapt_mutex, NULL);
	pthread_cond_init(&freqAdapt_condVar, NULL);
	pthread_attr_init(&freqAdaptWorkerThreadAttr);
	pthread_attr_setdetachstate(&freqAdaptWorkerThreadAttr, PTHREAD_CREATE_DETACHED);
	// pin the new thread to the corresponding virtual core
	// TODO: can we achieve a better performance by pinning it to another core?
	CPU_ZERO(&processor_mask);
	CPU_SET(20,&processor_mask);
	pthread_attr_setaffinity_np(&freqAdaptWorkerThreadAttr, sizeof(cpu_set_t), &processor_mask);

	pthread_create(&freqAdaptWorkerThread, &freqAdaptWorkerThreadAttr, freqAdaptWorker, this);

}


std::vector<FreqAdaptValues::FreqAdaptPoint> FreqAdaptValues::m_pos =
		std::vector<FreqAdaptValues::FreqAdaptPoint>();

unsigned long FreqAdaptValues::newFreq = 0;
pthread_mutex_t FreqAdaptValues::freqAdapt_mutex = pthread_mutex_t();
pthread_cond_t FreqAdaptValues::freqAdapt_condVar = pthread_cond_t();

pthread_attr_t FreqAdaptValues::freqAdaptWorkerThreadAttr = pthread_attr_t();
cpu_set_t FreqAdaptValues::processor_mask = cpu_set_t();
pthread_t FreqAdaptValues::freqAdaptWorkerThread = pthread_t();

