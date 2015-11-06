/*
 * Copyright (c) 2011-2014:  G-CSC, Goethe University Frankfurt
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

/**
 * \file parallel_coloring.cpp
 *
 * \author Martin Rupp
 *
 * \date 8.2.2011
 *
 * implementation file for coloring of processes
 *
 * Goethe-Center for Scientific Computing 2011.
 *
 */


#include "common/common.h"
#include "pcl/pcl.h"
#include "parallel_coloring.h"
#include "lib_algebra/common/stl_debug.h"


namespace ug
{
DebugID DBG_COLORING("ParallelColoring");

int ParallelColoring::color(pcl::InterfaceCommunicator<IndexLayout> &com)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	// 1. send all neighbors our number of neighbors
	size_t myDegree = pids.size();
	typedef std::set<int>::iterator setiterator;
	std::map<int, size_t> othersDegree;
	UG_DLOG(DBG_COLORING, 1, "my degree is " << myDegree);
	UG_DLOG(DBG_COLORING, 1, "sending degree data to ");
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid = *iter;
		com.send_raw(pid, &myDegree, sizeof(size_t), true);
		com.receive_raw(pid, &othersDegree[pid], sizeof(size_t));
		UG_DLOG(DBG_COLORING, 1, pid << " ");
	}
	UG_DLOG(DBG_COLORING, 1, "\n");

	com.communicate();

	int myPID = pcl::ProcRank();

/*#ifdef UG_ENABLE_DEBUG_LOGS
	IF_DEBUG(LIB_ALG_MATRIX, 1)
	{
		UG_DLOG(DBG_COLORING, 1, "process " <<  myPID << " here. I have degree " << myDegree << ".\n");
		for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			int pid = *iter;
			UG_DLOG(DBG_COLORING, 1, " neighbor pid " << pid << " has degree " << othersDegree[pid] << ".\n");
		}
	}
#endif*/

	stdvector<int> colors(myDegree, -1);

	UG_DLOG(DBG_COLORING, 1, "got degrees:\n");
	// get processes with higher/lower weight
	stdvector<int> processesWithHigherWeight, processesWithLowerWeight;
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid2 = *iter; size_t degree2 = othersDegree[pid2];
		UG_DLOG(DBG_COLORING, 1, " process " << pid2 << " has degree " << degree2 << "\n");
		if(myDegree < degree2 || (myDegree == degree2 && myPID < pid2))
			processesWithHigherWeight.push_back(pid2);
		else
			processesWithLowerWeight.push_back(pid2);
	}

	// 2. for all neighbors which have a higher weight than this process, issue a color receive
	UG_DLOG(DBG_COLORING, 1, "issue color receive from... ");
	colors.resize(processesWithHigherWeight.size(), -1);
	for(size_t i=0; i<processesWithHigherWeight.size(); i++)
	{
		int pid2 = processesWithHigherWeight[i];
		com.receive_raw(pid2, &colors[i], sizeof(int));
		UG_DLOG(DBG_COLORING, 1, pid2 << " ");
	}
	UG_DLOG(DBG_COLORING, 1, "\n");

	int myColor = 0;

	// 3. communicate & get a free color (only if actually got color from neighbors)
	if(processesWithHigherWeight.size() > 0)
	{
		//UG_DLOG(DBG_COLORING, 1, "communicate...");
		com.communicate();
		//UG_DLOG(DBG_COLORING, 1, "done.\n");

/*#ifdef UG_ENABLE_DEBUG_LOGS
		IF_DEBUG(LIB_ALG_MATRIX, 1)
		{
			UG_DLOG(DBG_COLORING, 1, "received colors from " << processesWithHigherWeight.size() << " processes!\n");
			for(size_t i=0; i<processesWithHigherWeight.size(); i++)
			{
				int pid2 = processesWithHigherWeight[i];
				if(colors[i] == -1)
					UG_DLOG(DBG_COLORING, 1, "pid " << pid2 << " did not send color?");
				UG_DLOG(DBG_COLORING, 1, " neighbor pid " << pid2 << " has color " << colors[i] << ".\n");
			}
		}
#endif*/

		// determine the smallest unused color for me
		stdvector<int> sortedcolors = colors;
		sort(sortedcolors.begin(), sortedcolors.end());
		myColor = -1; // dummy guess
		size_t i;
		for(i=0; i < sortedcolors.size(); i++)
		{
			if(myColor != sortedcolors[i])
			{
				myColor++; // next guess
				if(myColor != sortedcolors[i])
					break;
			}
		}
		if(i == sortedcolors.size())	// no gap found
			myColor++;			// introduce new color

		if(p_colors || p_processesWithLowerColor || p_processesWithHigherColor)
		{
			if(p_processesWithLowerColor) p_processesWithLowerColor->clear();
			if(p_processesWithHigherColor) p_processesWithHigherColor->clear();
			if(p_colors) p_colors->clear();
			UG_DLOG(DBG_COLORING, 1, "got colors from processes with higher weight: \n");
			for(size_t i=0; i<processesWithHigherWeight.size(); i++)
			{
				int pid2 = processesWithHigherWeight[i];
				UG_DLOG(DBG_COLORING, 1, " process " << pid2 << " has color " << colors[i] << "\n");
				if(p_processesWithLowerColor && colors[i] < myColor)
					p_processesWithLowerColor->push_back(pid2);
				else if (p_processesWithHigherColor && colors[i] > myColor)
					p_processesWithHigherColor->push_back(pid2);
				if(p_colors)
					(*p_colors)[pid2] = colors[i];
			}
		}
	}

	UG_DLOG(DBG_COLORING, 1, "GOT COLOR! My color is " << myColor << "\n");

	// 4. send color to all neighbors

	//UG_DLOG(DBG_COLORING, 1, "send color to all neighbors\n");
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid2 = *iter;
		com.send_raw(pid2, &myColor, sizeof(int), true);
	}
	UG_DLOG(DBG_COLORING, 1, "communicate...");
	com.communicate();
	UG_DLOG(DBG_COLORING, 1, "done.\n");

	// 5. get color from others
	// issue a color request from processors with lower rank
	if(processesWithLowerWeight.size() > 0)
	{
		//UG_LOG("issue a color request from processors with lower rank\n");
		UG_DLOG(DBG_COLORING, 1, "got colors from processes with lower weight: \n");
		colors.resize(processesWithLowerWeight.size(), -1);
		for(size_t i=0; i<processesWithLowerWeight.size(); i++)
		{
			int pid2 = processesWithLowerWeight[i];
			colors[i] = -1;
			com.receive_raw(pid2, &colors[i], sizeof(int));
		}

		//UG_DLOG(DBG_COLORING, 1, "communicate...");
		com.communicate();
		//UG_DLOG(DBG_COLORING, 1, "done.\n");

		if(p_colors || p_processesWithLowerColor || p_processesWithHigherColor)
		{
			for(size_t i=0; i<processesWithLowerWeight.size(); i++)
			{
				int pid2 = processesWithLowerWeight[i];
				UG_DLOG(DBG_COLORING, 1, " process " << pid2 << " has color " << colors[i] << "\n");
				if(p_processesWithLowerColor && colors[i] < myColor)
					p_processesWithLowerColor->push_back(pid2);
				else if (p_processesWithHigherColor && colors[i] > myColor)
					p_processesWithHigherColor->push_back(pid2);
				if(p_colors)
					(*p_colors)[pid2] = colors[i];
			}
		}
	}

	//UG_DLOG(DBG_COLORING, 1, "coloring done.\n");
	return myColor;
}


int ColorProcessorGraph(pcl::InterfaceCommunicator<IndexLayout> &com, std::set<int> &pids,
		std::vector<int> &processesWithLowerColor,
		std::vector<int> &processesWithHigherColor)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	ParallelColoring coloring;
	coloring.set_connections(pids);
	coloring.save_processes_with_lower_color_in(&processesWithLowerColor);
	coloring.save_processes_with_higher_color_in(&processesWithHigherColor);
	return coloring.color(com);
}


}
