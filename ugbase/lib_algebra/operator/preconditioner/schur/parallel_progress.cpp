/*
 * progress.cpp
 *
 *  Created on: 17.05.2013
 *      Author: mrupp
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
