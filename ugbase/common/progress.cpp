/*
 * progress.cpp
 *
 *  Created on: 17.05.2013
 *      Author: mrupp
 */

#include "progress.h"
namespace ug{

int Progress::totalDepth = 0;
int Progress::lastUpdateDepth = -1;

Progress::Progress(int minSecondsUntilProgress)
{
	m_minSecondsUntilProgress = minSecondsUntilProgress;
	m_length=100;
	bStarted=false;
	myDepth = totalDepth++;
}
/*

void Progress::setD(double now)
{
	if(now < 0 && now > m_total) return;
	if(bStarted)
	{
		if(myDepth == lastUpdateDepth)
		{
			int i2=(int)(m_length*now/m_total);
			if(i2 != posNow)
			{
				for(; posNow<i2; posNow++) { UG_LOG("-"); }
				GetLogAssistant().flush();
				calc_next_value();
			}
			m_now = now;
		}
	}
	else if(get_clock_s() - startS > m_minSecondsUntilProgress)
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
	else
	{
		posNow++;
		calc_next_value();
	}
}

*/

static inline void print_mark(int posNow, int length)
{
	if(posNow%(length/4)==0 )
		{UG_LOG("+");}
	else {UG_LOG("-");}
}

void Progress::setD(double now)
{
	if(now < 0 && now > m_total) return;
	if(bStarted)
	{
		if(myDepth == lastUpdateDepth)
		{
			int i2=(int)(m_length*now/m_total);
			if(i2 != posNow)
			{
				for(; posNow<i2; posNow++)
					print_mark(posNow, m_length);
				GetLogAssistant().flush();
				calc_next_value();
			}
			m_now = now;
		}
	}
	else if(get_clock_s() - startS > m_minSecondsUntilProgress)
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
	else
	{
		posNow++;
		calc_next_value();
	}
}
void Progress::stop()
{
	if(bStarted && myDepth == lastUpdateDepth)
	{
		int i=(int)(m_length*m_now/m_total);
		for(; i<m_length; i++)
			print_mark(i, m_length);
		UG_LOG("]");
		UG_LOG(" took " << reset_floats << get_clock_s()-startS << " s.\n");
		lastUpdateDepth = myDepth;
		bStarted = false;
	}
}
void Progress::start(double total, std::string msg)
{
	m_msg = msg;
	m_now = 0;
	m_total = total;
	bStarted=false;
	startS = get_clock_s();
	posNow = 0;
	calc_next_value();
}

}
