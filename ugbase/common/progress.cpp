/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#include "progress.h"

#include <utility>

#include "log.h"
#include "util/ostream_util.h"

namespace ug {

int Progress::totalDepth = 0;
int Progress::lastUpdateDepth = -1;

size_t g_minSecondUntilProgress = 3;

Progress::Progress(int minSecondsUntilProgress)
{
	if(minSecondsUntilProgress == -1)
		minSecondsUntilProgress = g_minSecondUntilProgress;
	m_minSecondsUntilProgress = minSecondsUntilProgress;
	m_length=100;
	bStarted=false;
	posNow=-1;
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

static inline void print_mark(int &posNow, int length)
{
	if(posNow%(length/4)==0 )
	{
		switch(posNow / (length/4))
		{
		case 0: UG_LOG("-"); break;
		case 1: UG_LOG("25%");posNow+=2; break;
		case 2: UG_LOG("50%");posNow+=2; break;
		case 3: UG_LOG("75%");posNow+=2; break;
		}
	}
	else {UG_LOG("-");}
}

void Progress::setD(double now)
{
	if(now < 0 && now > m_total) return;
	if(bStarted)
	{
		if(myDepth == lastUpdateDepth)
		{
			int i2=static_cast<int>(m_length * now / m_total);
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
		if(!m_msg.empty())
			{UG_LOG("\n" << m_msg);}
		UG_LOG("\n");
//		UG_LOG("." << repeat('_', m_length) << ".\n");
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
		posNow = -1;
	}
	if(!bStarted && posNow != -1 && get_clock_s() - startS > m_minSecondsUntilProgress*2)
	{
		UG_LOG(m_msg << ". took " << reset_floats << get_clock_s()-startS << " s.\n");
		bStarted = false;
		posNow = -1;
	}
}
void Progress::start(double total, std::string msg)
{
	m_msg = std::move(msg);
	m_now = 0;
	m_total = total;
	bStarted=false;
	startS = get_clock_s();
	posNow = 0;
	calc_next_value();
}

}
