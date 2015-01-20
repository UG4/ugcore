/*
 * freq_adapt.h
 *
 *  Created on: 11.12.2014
 *      Author: andreasvogel
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
