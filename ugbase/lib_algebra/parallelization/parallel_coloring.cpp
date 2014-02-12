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

int ParallelColoring::color(pcl::InterfaceCommunicator<IndexLayout> &com)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	// 1. send all neighbors our number of neighbors
	size_t myDegree = pids.size();
	typedef std::set<int>::iterator setiterator;
	std::map<int, size_t> othersDegree;
	UG_DLOG(LIB_ALG_AMG, 1, "my degree is " << myDegree);
	UG_DLOG(LIB_ALG_AMG, 1, "sending degree data to ");
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid = *iter;
		com.send_raw(pid, &myDegree, sizeof(size_t), true);
		com.receive_raw(pid, &othersDegree[pid], sizeof(size_t));
		UG_DLOG(LIB_ALG_AMG, 1, pid << " ");
	}
	UG_DLOG(LIB_ALG_AMG, 1, "\n");

	com.communicate();

	int myPID = pcl::ProcRank();

/*#ifdef UG_ENABLE_DEBUG_LOGS
	IF_DEBUG(LIB_ALG_MATRIX, 1)
	{
		UG_DLOG(LIB_ALG_AMG, 1, "process " <<  myPID << " here. I have degree " << myDegree << ".\n");
		for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			int pid = *iter;
			UG_DLOG(LIB_ALG_AMG, 1, " neighbor pid " << pid << " has degree " << othersDegree[pid] << ".\n");
		}
	}
#endif*/

	stdvector<int> colors(myDegree, -1);

	UG_DLOG(LIB_ALG_AMG, 1, "got degrees:\n");
	// get processes with higher/lower weight
	stdvector<int> processesWithHigherWeight, processesWithLowerWeight;
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid2 = *iter; size_t degree2 = othersDegree[pid2];
		UG_DLOG(LIB_ALG_AMG, 1, " process " << pid2 << " has degree " << degree2 << "\n");
		if(myDegree < degree2 || (myDegree == degree2 && myPID < pid2))
			processesWithHigherWeight.push_back(pid2);
		else
			processesWithLowerWeight.push_back(pid2);
	}

	// 2. for all neighbors which have a higher weight than this process, issue a color receive
	UG_DLOG(LIB_ALG_AMG, 1, "issue color receive from... ");
	colors.resize(processesWithHigherWeight.size(), -1);
	for(size_t i=0; i<processesWithHigherWeight.size(); i++)
	{
		int pid2 = processesWithHigherWeight[i];
		com.receive_raw(pid2, &colors[i], sizeof(int));
		UG_DLOG(LIB_ALG_AMG, 1, pid2 << " ");
	}
	UG_DLOG(LIB_ALG_AMG, 1, "\n");

	int myColor = 0;

	// 3. communicate & get a free color (only if actually got color from neighbors)
	if(processesWithHigherWeight.size() > 0)
	{
		//UG_DLOG(LIB_ALG_AMG, 1, "communicate...");
		com.communicate();
		//UG_DLOG(LIB_ALG_AMG, 1, "done.\n");

/*#ifdef UG_ENABLE_DEBUG_LOGS
		IF_DEBUG(LIB_ALG_MATRIX, 1)
		{
			UG_DLOG(LIB_ALG_AMG, 1, "received colors from " << processesWithHigherWeight.size() << " processes!\n");
			for(size_t i=0; i<processesWithHigherWeight.size(); i++)
			{
				int pid2 = processesWithHigherWeight[i];
				if(colors[i] == -1)
					UG_DLOG(LIB_ALG_AMG, 1, "pid " << pid2 << " did not send color?");
				UG_DLOG(LIB_ALG_AMG, 1, " neighbor pid " << pid2 << " has color " << colors[i] << ".\n");
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
			UG_DLOG(LIB_ALG_AMG, 1, "got colors from processes with higher weight: \n");
			for(size_t i=0; i<processesWithHigherWeight.size(); i++)
			{
				int pid2 = processesWithHigherWeight[i];
				UG_DLOG(LIB_ALG_AMG, 1, " process " << pid2 << " has color " << colors[i] << "\n");
				if(p_processesWithLowerColor && colors[i] < myColor)
					p_processesWithLowerColor->push_back(pid2);
				else if (p_processesWithHigherColor && colors[i] > myColor)
					p_processesWithHigherColor->push_back(pid2);
				if(p_colors)
					(*p_colors)[pid2] = colors[i];
			}
		}
	}

	UG_DLOG(LIB_ALG_AMG, 1, "GOT COLOR! My color is " << myColor << "\n");

	// 4. send color to all neighbors

	//UG_DLOG(LIB_ALG_AMG, 1, "send color to all neighbors\n");
	for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
	{
		int pid2 = *iter;
		com.send_raw(pid2, &myColor, sizeof(int), true);
	}
	UG_DLOG(LIB_ALG_AMG, 1, "communicate...");
	com.communicate();
	UG_DLOG(LIB_ALG_AMG, 1, "done.\n");

	// 5. get color from others
	// issue a color request from processors with lower rank
	if(processesWithLowerWeight.size() > 0)
	{
		//UG_LOG("issue a color request from processors with lower rank\n");
		UG_DLOG(LIB_ALG_AMG, 1, "got colors from processes with lower weight: \n");
		colors.resize(processesWithLowerWeight.size(), -1);
		for(size_t i=0; i<processesWithLowerWeight.size(); i++)
		{
			int pid2 = processesWithLowerWeight[i];
			colors[i] = -1;
			com.receive_raw(pid2, &colors[i], sizeof(int));
		}

		//UG_DLOG(LIB_ALG_AMG, 1, "communicate...");
		com.communicate();
		//UG_DLOG(LIB_ALG_AMG, 1, "done.\n");

		if(p_colors || p_processesWithLowerColor || p_processesWithHigherColor)
		{
			for(size_t i=0; i<processesWithLowerWeight.size(); i++)
			{
				int pid2 = processesWithLowerWeight[i];
				UG_DLOG(LIB_ALG_AMG, 1, " process " << pid2 << " has color " << colors[i] << "\n");
				if(p_processesWithLowerColor && colors[i] < myColor)
					p_processesWithLowerColor->push_back(pid2);
				else if (p_processesWithHigherColor && colors[i] > myColor)
					p_processesWithHigherColor->push_back(pid2);
				if(p_colors)
					(*p_colors)[pid2] = colors[i];
			}
		}
	}

	//UG_DLOG(LIB_ALG_AMG, 1, "coloring done.\n");
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
