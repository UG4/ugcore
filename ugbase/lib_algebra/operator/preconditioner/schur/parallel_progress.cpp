/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifdef UG_PARALLEL

#include "parallel_progress.h"
#include "pcl/pcl.h"
namespace ug{


ParallelProgress::ParallelProgress()
{
	m_length=100;
}

static inline void print_mark(int &posNow, int length)
{
	ug::LogAssistant& la = ug::GetLogAssistant();
	int op = la.get_output_process();
	la.set_output_process(-1);\
	UG_LOG("-");
	la.flush();
	la.set_output_process(op);
}

void ParallelProgress::setD(double now)
{
	if(now < 0 && now > m_total) return;


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
void ParallelProgress::stop()
{
	setD(m_total);
	UG_LOG("] took " << reset_floats << get_clock_s()-startS << " s.\n");
}
void ParallelProgress::start(double total, std::string msg, size_t numProc)
{
	if(pcl::ProcRank() == 0)
	{
		UG_LOG(msg << "\n[");
		for(int i=0; i<m_length; i++)
			UG_LOG("-");
		UG_LOG("]\n[");

		m_length = m_length / numProc + m_length%numProc;
	}
	else
		m_length /= numProc;

	m_msg = msg;
	m_now = 0;
	m_total = total;
	startS = get_clock_s();
	posNow = 0;
	calc_next_value();
}

}
#endif
