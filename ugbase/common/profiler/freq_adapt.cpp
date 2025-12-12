/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "freq_adapt.h"

#include <cpufreq.h>  // frequency adaption library
#include <iostream>   // usual i/o
#include <fstream>    // file i/o
#include <cstdlib>	  //used for atoi. Better: boost::lexical_cast
#include <unistd.h>   // usleep

#include "common/common.h"
#include "common/util/string_util.h"

//#include <sys/time.h> // for testing purposes only, can be removed

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

	//timeval tv; // for testing purposes only, can be removed

	if(freq){
		// remember curr freq
		// we now use the most recently requested value instead of the actual current frequency
		m_prevFreq = FreqAdaptValues::newFreq;

		// for testing purposes only, can be removed
		//if ( gettimeofday(&tv, nullptr) != 0 ) {
		//	printf("error while getting time\n");
		//}
		//cout << "CPU_FREQ: time = " << tv.tv_sec*1000000+tv.tv_usec << " , Changing frequency to " << freq << " ...\n";

		// set new frequency
		FreqAdaptValues::adjust_freq(freq);

		// for testing purposes only, can be removed
		//if ( gettimeofday(&tv, nullptr) != 0 ) {
		//	printf("error while getting time\n");
		//}
		//cout << "CPU_FREQ: time = " << tv.tv_sec*1000000+tv.tv_usec << " , Changing frequency to " << freq << " ... done\n";
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

	//timeval tv; // for testing purposes only, can be removed

	if(m_bActive){
		// check if freq was changed on entry
		if(m_prevFreq){

			// for testing purposes only, can be removed
			//if ( gettimeofday(&tv, nullptr) != 0 ) {
			//	printf("error while getting time\n");
			//}
			//cout << "CPU_FREQ: time = " << tv.tv_sec*1000000+tv.tv_usec << " , Resetting frequency to " << m_prevFreq << " ...\n";

			// set previous freq
			FreqAdaptValues::adjust_freq(m_prevFreq);

			// for testing purposes only, can be removed
			//if ( gettimeofday(&tv, nullptr) != 0 ) {
			//	printf("error while getting time\n");
			//}
			//cout << "CPU_FREQ: time = " << tv.tv_sec*1000000+tv.tv_usec << " , Resetting frequency to " << m_prevFreq << " ... done\n";
		}

		m_bActive = false;
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

	// for testing purposes only, can be removed
	//cout << "CPU_FREQ: Node for: "<<file<<", line: " << line << endl;

	for(size_t i = 0; i < m_pos.size(); ++i){
		if(m_pos[i].line == line && m_pos[i].file == file){
			return m_pos[i].freq;
		}
	}

	return 0;
}


void* FreqAdaptValues::freqAdaptWorker(void* This) {

	unsigned long freq;

	//timeval tv; // for testing purposes only, can be removed

	while (true) {

		// wait on frequency transition request
		if ( pthread_mutex_lock(&freqAdapt_mutex) != 0 )
			UG_THROW("Error while locking freqAdapt_mutex");
		while (newFreq == 0) {
			if ( pthread_cond_wait(&freqAdapt_condVar, &freqAdapt_mutex) != 0 )
				UG_THROW("Error while waiting on freqAdapt_condVar");
		}
		if ( pthread_mutex_unlock(&freqAdapt_mutex) != 0 )
			UG_THROW("Error while unlocking freqAdapt_mutex");

		/* sleep some time to handle only the very last of a series of
                   requests as seen in case of multiple nested routines */
		usleep(10); // TODO: find an appropriate value

		if ( pthread_mutex_lock(&freqAdapt_mutex) != 0 )
			UG_THROW("Error while locking freqAdapt_mutex");
		freq = newFreq;
		newFreq = 0;
		if ( pthread_mutex_unlock(&freqAdapt_mutex) != 0 )
			UG_THROW("Error while unlocking freqAdapt_mutex");

		// for testing purposes only, can be removed
		//if ( gettimeofday(&tv, nullptr) != 0 ) {
		//	printf("error while getting time\n");
		//}
		//cout << "CPU_FREQ: time = " << tv.tv_sec*1000000+tv.tv_usec << " , calling cpufreq_set_frequency() with " << freq << " ...\n";

		// do the actual frequency transition
		if (cpufreq_set_frequency(0, freq) != 0) {
			UG_THROW("FreqAdaptValues::freqAdaptWorker: Error while setting frequency");
		}

		// for testing purposes only, can be removed
		//if ( gettimeofday(&tv, nullptr) != 0 ) {
		//	printf("error while getting time\n");
		//}
		//cout << "CPU_FREQ: time = " << tv.tv_sec*1000000+tv.tv_usec << " , calling cpufreq_set_frequency() with " << freq << " ... done\n";
	}

}

void FreqAdaptValues::adjust_freq(unsigned long freq){
	if ( pthread_mutex_lock(&freqAdapt_mutex) != 0 )
	        UG_THROW("Error while locking freqAdapt_mutex");
	newFreq = freq;
	if ( pthread_cond_signal(&freqAdapt_condVar) != 0 )
	        UG_THROW("Error while signaling freqAdapt_condVar");
	if ( pthread_mutex_unlock(&freqAdapt_mutex) != 0 )
	        UG_THROW("Error while unlocking freqAdapt_mutex");
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

	// for testing purposes only, can be removed
	//cout << "FreqAdaptValues: parsing file: "<<csvFile << "\n";

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
		// instead of atoi, better use boost::lexical_cast
		// int lineNumber = boost::lexical_cast<int>(s2);
		// int frequency = boost::lexical_cast<int>(s3);

		// for testing purposes only, can be removed
		//cout << file << "===" << line << "===" << freq << "===" << endl;

		m_pos.push_back( FreqAdaptPoint(file,line,freq) );
	}


	///////////////////////////
	// start adaption thread
	///////////////////////////

	newFreq = 0;

	if ( pthread_mutex_init(&freqAdapt_mutex, nullptr) != 0 )
	        UG_THROW("Error while initializing freqAdapt_mutex");
	if ( pthread_cond_init(&freqAdapt_condVar, nullptr) != 0 )
	        UG_THROW("Error while initializing freqAdapt_condVar");
	if ( pthread_attr_init(&freqAdaptWorkerThreadAttr) != 0 )
	        UG_THROW("Error while initializing freqAdaptWorkerThreadAttr");
	if ( pthread_attr_setdetachstate(&freqAdaptWorkerThreadAttr, PTHREAD_CREATE_DETACHED) != 0 )
	        UG_THROW("Error while setting detach state of freqAdaptWorkerThreadAttr");
	// pin the new thread to the corresponding virtual core
	// TODO: can we achieve a better performance by pinning it to another core?
	CPU_ZERO(&processor_mask);
	CPU_SET(20,&processor_mask);
	if ( pthread_attr_setaffinity_np(&freqAdaptWorkerThreadAttr, sizeof(cpu_set_t), &processor_mask) != 0 )
	        UG_THROW("Error while setting affinity of freqAdaptWorkerThreadAttr");

	if ( pthread_create(&freqAdaptWorkerThread, &freqAdaptWorkerThreadAttr, freqAdaptWorker, nullptr) != 0 )
	        UG_THROW("Error while creating thread freqAdaptWorkerThread");



	///////////////////////////
	// initialize newFreq
	///////////////////////////

	if ( (newFreq = cpufreq_get_freq_kernel(0)) == 0)
	        UG_THROW("Error while getting frequency");
}


std::vector<FreqAdaptValues::FreqAdaptPoint> FreqAdaptValues::m_pos =
		std::vector<FreqAdaptValues::FreqAdaptPoint>();

unsigned long FreqAdaptValues::newFreq = 0;
pthread_mutex_t FreqAdaptValues::freqAdapt_mutex = pthread_mutex_t();
pthread_cond_t FreqAdaptValues::freqAdapt_condVar = pthread_cond_t();

pthread_attr_t FreqAdaptValues::freqAdaptWorkerThreadAttr = pthread_attr_t();
cpu_set_t FreqAdaptValues::processor_mask = cpu_set_t();
pthread_t FreqAdaptValues::freqAdaptWorkerThread = pthread_t();
