/*
 * progress.h
 *
 *  Created on: 05.10.2012
 *      Author: mrupp
 */

#ifndef PROGRESS_H_
#define PROGRESS_H_

#include <iostream>
#include "common/util/string_util.h"
#include "common/log.h"
#include <ctime>
#include <string>
#include <sstream>

namespace ug
{

class Progress
{
	int m_minSecondsUntilProgress;
	double startS;
public:
	Progress(int minSecondsUntilProgress=1);

	~Progress()
	{
		stop();
		totalDepth--;
	}
	void set_length(int l)
	{
		m_length = l;
	}
	inline void start(double total, std::string msg="")
	{
		m_msg = msg;
		m_now = 0;
		m_total = total;
		bStarted=false;
		startS = clock_s();
	}
	double clock_s()
	{
		return clock() / CLOCKS_PER_SEC;
	}
	inline void set(double now)
	{
		if(now < 0 && now > m_total) return;
		if(bStarted && myDepth == lastUpdateDepth)
		{
			int i2=(int)(m_length*now/m_total);
			for(; posNow<i2; posNow++) { UG_LOG("-"); }
			m_now = now;
		}
		else if(clock_s() - startS > m_minSecondsUntilProgress)
		{
			if(m_msg.length() > 0)
				{UG_LOG("\n" << m_msg);}
			UG_LOG("\n." << repeat('_', m_length) << ".\n");
			UG_LOG("[");
			bStarted = true;
			posNow = 0;
			lastUpdateDepth = myDepth;
			set(now);

		}
	}
	inline void stop()
	{
		if(bStarted && myDepth == lastUpdateDepth)
		{
			int i=(int)(m_length*m_now/m_total);
			for(; i<m_length; i++) { UG_LOG("-"); }
			UG_LOG("]\n");
			lastUpdateDepth = myDepth;
			bStarted = false;
		}
	}
private:
	int posNow;
	bool bStarted;
	double m_total;
	double m_now;
	int m_length;
	std::string m_msg;

	static int totalDepth;
	static int lastUpdateDepth;
	int myDepth;
};

}

#define PROGRESS_START(progVarName, dSize, msg) \
	ug::Progress progVarName; { std::stringstream ss; ss << msg; progVarName.start(dSize, ss.str()); }

#define PROGRESS_UPDATE(progVarName, d) progVarName.set(d);
#define PROGRESS_FINISH(progVarName) progVarName.stop();


#endif /* PROCESS_H_ */
