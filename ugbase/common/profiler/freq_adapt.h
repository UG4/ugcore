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

#ifndef __H__UG__COMMON__PROFILER__FREQ_ADAPT__
#define __H__UG__COMMON__PROFILER__FREQ_ADAPT__

#include <stack>
#include <vector>
#include <string>
#include <pthread.h>  // used to start separate thread controlling freqency

class AutoFreqAdaptNode;

////////////////////////////////////////////////////////////////////////////////
// FreqAdaptNodeManager
////////////////////////////////////////////////////////////////////////////////

class FreqAdaptNodeManager
{
	public:
		static void add(AutoFreqAdaptNode* node);
		static void release_latest();

	private:
		FreqAdaptNodeManager();
		~FreqAdaptNodeManager();

	public:
		static FreqAdaptNodeManager& inst();

//	private:
		std::stack<AutoFreqAdaptNode*>	m_nodes;
};

////////////////////////////////////////////////////////////////////////////////
// AutoFreqAdaptNode
////////////////////////////////////////////////////////////////////////////////

class AutoFreqAdaptNode{
	friend class FreqAdaptNodeManager;

	public:
		AutoFreqAdaptNode(unsigned long freq);
		~AutoFreqAdaptNode();

	private:
		void release();
		inline bool is_active()		{return m_bActive;}

	protected:
		bool m_bActive;
		unsigned long m_prevFreq;
};


////////////////////////////////////////////////////////////////////////////////
// FreqAdaptValues
////////////////////////////////////////////////////////////////////////////////

class FreqAdaptValues {
	friend class AutoFreqAdaptNode;

	private:
	// disallow constructor, destructor; copy and assignment (intentionally left unimplemented)
		FreqAdaptValues() {};
		FreqAdaptValues(const FreqAdaptValues&);
		FreqAdaptValues& operator=(const FreqAdaptValues&);
		~FreqAdaptValues() {};

	// 	Singleton provider
		static FreqAdaptValues& inst();

	private:
	//	return freq if adjusted; if not contained, return 0
		unsigned long find_freq(const char* file, const int line);

	private:
		struct FreqAdaptPoint{
			std::string file;
			int line;
			unsigned long freq;
			FreqAdaptPoint(std::string _file, int _line, unsigned long _freq) :
				file(_file), line(_line), freq(_freq) {}
		};

	//	returns the continuous information
		static std::vector<FreqAdaptPoint> m_pos;

	public:
		// returns requested frequency (or 0 if freq not adjusted) at file and line
		static unsigned long freq(const char* file, const int line){
			return inst().find_freq(file, line);
		}

		// reads the database of (file, line, required frequency)
		static void set_freqs(std::string csvFile);

		// adjust the frequency of the cpu
		static void adjust_freq(unsigned long freq);

	protected:
		static void* freqAdaptWorker(void* This);

		static unsigned long newFreq;
		static pthread_mutex_t freqAdapt_mutex;
		static pthread_cond_t freqAdapt_condVar;

		static pthread_attr_t freqAdaptWorkerThreadAttr;
		static cpu_set_t processor_mask;
		static pthread_t freqAdaptWorkerThread;

};


#endif /* __H__UG__COMMON__PROFILER__FREQ_ADAPT__ */
