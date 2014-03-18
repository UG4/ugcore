/**
 * \file parallel_coloring.h
 *
 * \author Martin Rupp
 *
 * \date 9.3.2011
 *
 * header file for coloring of processes
 *
 * Goethe-Center for Scientific Computing 2011.
 *
 */


#include <map>
#include <vector>
#include <set>
#include "parallelization_util.h"

namespace ug
{
extern DebugID DBG_COLORING;

class ParallelColoring
{
public:
	ParallelColoring() : p_colors(NULL), p_processesWithLowerColor(NULL), p_processesWithHigherColor(NULL)
	{

	}
	void set_connections(std::set<int> &pids_)
	{
		pids = pids_;
	}
	/*void set_connections(pcl::InterfaceCommunicator<IndexLayout> &master,
			pcl::InterfaceCommunicator<IndexLayout> &slave, bool bAddSlaveSlave)
	{
	}*/

	int color(pcl::InterfaceCommunicator<IndexLayout> &com);

	void save_processes_with_lower_color_in(std::vector<int> *p)
	{
		p_processesWithLowerColor = p;
	}

	void save_processes_with_higher_color_in(std::vector<int> *p)
	{
		p_processesWithHigherColor = p;
	}

	void save_colors_in(std::map<int, int> *p)
	{
		p_colors = p;
	}

private:
	std::map<int, int> *p_colors;
	std::vector<int> *p_processesWithLowerColor;
	std::vector<int> *p_processesWithHigherColor;
	std::set<int> pids;
};

int ColorProcessorGraph(pcl::InterfaceCommunicator<IndexLayout> &com, std::set<int> &pids,
		std::vector<int> &processesWithLowerColor,
		std::vector<int> &processesWithHigherColor);

}
