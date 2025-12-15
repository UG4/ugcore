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

#ifndef PROGRESS_H_
#define PROGRESS_H_

#include <iostream>
//#include "common/util/string_util.h"
//#include "common/log.h"
//#include <ctime>
#include <string>
//#include <sstream>

#include "stopwatch.h"
//#include "util/ostream_util.h"

namespace ug {

class Progress
{
public:
	explicit Progress(int minSecondsUntilProgress=-1);

	~Progress()
	{
		stop();
		totalDepth--;
	}
	void set_length(int l)
	{
		m_length = l;
	}

	void calc_next_value()
	{
		dNextValueToUpdate = (posNow+1)*m_total/m_length;
		iNextValueToUpdate = static_cast<int>(dNextValueToUpdate);
	}
	void start(double total, std::string msg="");

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

	int m_minSecondsUntilProgress;
	double startS;
	double dNextValueToUpdate;
	size_t iNextValueToUpdate;
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

#define PROGRESS_START_WITH(progVarName, dSize, msg) \
	{ std::stringstream ss; ss << msg; progVarName.start(dSize, ss.str()); }

#define PROGRESS_UPDATE(progVarName, d) progVarName.set(d);
#define PROGRESS_FINISH(progVarName) progVarName.stop();


#endif