/*
 * progress.h
 *
 *  Created on: 05.10.2012
 *      Author: mrupp
 */

#ifndef PARALLEL_PROGRESS_H_
#define PARALLEL_PROGRESS_H_

#ifdef UG_PARALLEL

#include <iostream>
#include "common/util/string_util.h"
#include "common/log.h"
#include <ctime>
#include <string>
#include <sstream>

#include "common/stopwatch.h"
#include "common/util/ostream_util.h"

namespace ug
{

class ParallelProgress
{
	int m_minSecondsUntilProgress;
	double startS;
public:
	ParallelProgress();

	~ParallelProgress()
	{
		stop();
	}
	void set_length(int l)
	{
		m_length = l;
	}

	void calc_next_value()
	{
		dNextValueToUpdate = (posNow+1)*m_total/m_length;
		iNextValueToUpdate = (int) dNextValueToUpdate;
	}
	void start(double total, std::string msg, size_t numProc);

	inline void set(size_t now)
	{
		if(now < iNextValueToUpdate) return;
		setD(now);
	}
	inline void set(int now)
	{
		if(now < 0 || (size_t)now < iNextValueToUpdate) return;
		setD(now);
	}

	inline void set(double now)
	{
		if(now < dNextValueToUpdate) return;
		setD(now);
	}
	void setD(double now);

	void stop();
private:
	double dNextValueToUpdate;
	size_t iNextValueToUpdate;
	int posNow;
	double m_total;
	double m_now;
	int m_length;
	std::string m_msg;
};

}

#define PARALLEL_PROGRESS_START(progVarName, dSize, msg, numProcs) \
	ug::ParallelProgress progVarName; { std::stringstream ss; ss << msg; progVarName.start(dSize, ss.str(), numProcs); }

#define PARALLEL_PROGRESS_UPDATE(progVarName, d) progVarName.set(d);
#define PARALLEL_PROGRESS_FINISH(progVarName) progVarName.stop();


#endif

#endif /* PARALLEL_PROGRESS_H_ */
