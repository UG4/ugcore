/*
 * freq_adapt.h
 *
 *  Created on: 11.12.2014
 *      Author: andreasvogel
 */

#ifndef FREQ_ADAPT_H_
#define FREQ_ADAPT_H_

#include <stack>

class AutoFreqAdaptNode;

////////////////////////////////////////////////////////////////////////////////
// AutoFreqAdaptNode
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
	private:
	// disallow constructor, destructor; copy and assignment (intentionally left unimplemented)
		FreqAdaptValues() {};
		FreqAdaptValues(const FreqAdaptValues&);
		FreqAdaptValues& operator=(const FreqAdaptValues&);
		~FreqAdaptValues() {};

	// 	Singleton provider
		inline static FreqAdaptValues& inst();

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
		// returns requested frequency (or 0 if freq not adjusted)
		unsigned long freq(const char* file, const int line){
			return inst().find_freq();
		}

		static void read_marks( /* some database for reading */);

};

std::vector<FreqAdaptValues::FreqAdaptPoint> FreqAdaptValues::m_pos =
		std::vector<FreqAdaptValues::FreqAdaptPoint>();


#endif /* FREQ_ADAPT_H_ */
