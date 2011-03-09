/**
 * \file famg_coloring.cpp
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


#include "ug.h"
#include "common/common.h"
#include "pcl/pcl.h"
#include "lib_algebra/common/stl_debug.h"
#include "lib_algebra/parallelization/parallel_index_layout.h"
#include <set>

namespace ug
{

class PCLProcessorColoring
{

	void set_connections(std::set<int> &pids)
	{

	}

	void set_connections(pcl::ParallelCommunicator<IndexLayout> &master,
			pcl::ParallelCommunicator<IndexLayout> &slave, bool bAddSlaveSlave)
	{
	}

	int color(pcl::ParallelCommunicator<IndexLayout> &com)
	{
		// 1. send all neighbors our number of neighbors
		size_t myDegree = pids.size();
		typedef std::set<int>::iterator setiterator;
		std::map<int, size_t> othersDegree;

		UG_DLOG(LIB_ALG_AMG, 2, "sending degree data to ");
		for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			int pid = *iter;
			com.send_raw(pid, &myDegree, sizeof(size_t), true);
			com.receive_raw(pid, &othersDegree[pid], sizeof(size_t));
			UG_DLOG(LIB_ALG_AMG, 2, pid << " ");
		}
		UG_DLOG(LIB_ALG_AMG, 2, "\n");

		com.communicate();

		int myPID = pcl::GetProcRank();

	#ifdef UG_ENABLE_DEBUG_LOGS
		IF_DEBUG(LIB_ALG_MATRIX, 2)
		{
			UG_DLOG(LIB_ALG_AMG, 2, "process " <<  myPID << " here. I have degree " << myDegree << ".\n");
			for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
			{
				int pid = *iter;
				UG_DLOG(LIB_ALG_AMG, 2, " neighbor pid " << pid << " has degree " << othersDegree[pid] << ".\n");
			}
		}
	#endif

		stdvector<int> colors(myDegree, -1);

		// get processes with higher/lower weight
		stdvector<int> processesWithHigherWeight, processesWithLowerWeight;
		for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			int pid2 = *iter; size_t degree2 = othersDegree[pid2];
			if(myDegree < degree2 || (myDegree == degree2 && myPID < pid2))
				processesWithHigherWeight.push_back(pid2);
			else
				processesWithLowerWeight.push_back(pid2);
		}

		// 2. for all neighbors which have a higher weight than this process, issue a color receive
		UG_DLOG(LIB_ALG_AMG, 2, "issue receive from... ");
		colors.resize(processesWithHigherWeight.size(), -1);
		for(size_t i=0; i<processesWithHigherWeight.size(); i++)
		{
			int pid2 = processesWithHigherWeight[i];
			com.receive_raw(pid2, &colors[i], sizeof(int));
			UG_DLOG(LIB_ALG_AMG, 2, pid2 << " ");
		}
		UG_DLOG(LIB_ALG_AMG, 2, "\n");

		int myColor = 0;

		// 3. communicate & get a free color (only if actually got color from neighbors)
		if(processesWithHigherWeight.size() > 0)
		{
			UG_DLOG(LIB_ALG_AMG, 2, "communicate...");
			com.communicate();
			UG_DLOG(LIB_ALG_AMG, 2, "done.\n");

	#ifdef UG_ENABLE_DEBUG_LOGS
			IF_DEBUG(LIB_ALG_MATRIX, 2)
			{
				UG_DLOG(LIB_ALG_AMG, 2, "received colors from " << processesWithHigherWeight.size() << " processes!\n");
				for(size_t i=0; i<processesWithHigherWeight.size(); i++)
				{
					int pid2 = processesWithHigherWeight[i];
					if(colors[i] == -1)
						UG_DLOG(LIB_ALG_AMG, 2, "pid " << pid2 << " did not send color?");
					UG_DLOG(LIB_ALG_AMG, 2, " neighbor pid " << pid2 << " has color " << colors[i] << ".\n");
				}
			}
	#endif

			// determine the smallest unused color for me
			sort(colors.begin(), colors.end());
			myColor = -1; // dummy guess
			size_t i;
			for(i=0; i < colors.size(); i++)
			{
				if(myColor != colors[i])
				{
					myColor++; // next guess
					if(myColor != colors[i])
						break;
				}
			}
			if(i == colors.size())	// no gap found
				myColor++;			// introduce new color
		}

		UG_DLOG(LIB_ALG_AMG, 2, "GOT COLOR! My color is " << myColor << "\n");

		// 4. send color to all neighbors

		UG_LOG("send color to all neighbors\n");
		for(setiterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			int pid2 = *iter;
			com.send_raw(pid2, &myColor, sizeof(int), true);
		}
		UG_DLOG(LIB_ALG_AMG, 2, "communicate...");
		com.communicate();
		UG_DLOG(LIB_ALG_AMG, 2, "done.\n");

		// 5. get color from others
		if(p_colors || p_processesWithLowerColor || p_processesWithHigherColor)
		{
			if(p_processesWithLowerColor) p_processesWithLowerColor->clear();
			if(p_processesWithHigherColor) p_processesWithHigherColor->clear();
			if(p_colors) p_colors->clear();
			for(size_t i=0; i<processesWithHigherWeight.size(); i++)
			{
				int pid2 = processesWithHigherWeight[i];
				if(p_processesWithLowerColor && colors[i] < myColor)
					p_processesWithLowerColor->push_back(pid2);
				else if (p_processesWithHigherColor && colors[i] > myColor)
					p_processesWithHigherColor->push_back(pid2);
				if(p_colors)
					p_colors->at(pid2) = colors[i];
			}
		}

		// issue a color request from processors with lower rank
		if(processesWithLowerWeight.size() > 0)
		{
			UG_LOG("issue a color request from processors with lower rank\n");
			colors.resize(processesWithLowerWeight.size(), -1);
			for(size_t i=0; i<processesWithLowerWeight.size(); i++)
			{
				int pid2 = processesWithLowerWeight[i];
				colors[i] = -1;
				com.receive_raw(pid2, &colors[i], sizeof(int));
			}

			UG_DLOG(LIB_ALG_AMG, 2, "communicate...");
			com.communicate();
			UG_DLOG(LIB_ALG_AMG, 2, "done.\n");

			if(p_colors || p_processesWithLowerColor || p_processesWithHigherColor)
			{
				for(size_t i=0; i<processesWithLowerWeight.size(); i++)
				{
					int pid2 = processesWithLowerWeight[i];
					if(p_processesWithLowerColor && colors[i] < myColor)
						p_processesWithLowerColor->push_back(pid2);
					else if (p_processesWithHigherColor && colors[i] > myColor)
						p_processesWithHigherColor->push_back(pid2);
					if(p_colors)
						p_colors->at(pid2) = colors[i];
				}
			}
		}

		UG_DLOG(LIB_ALG_AMG, 2, "coloring done.\n");
		return myColor;
	}




	std::map<int, int> *p_colors;
	std::vector<int> *p_processesWithLowerColor;
	std::vector<int> *p_processesWithHigherColor;
	std::set<int> pids;
};


}
