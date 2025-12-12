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
#ifndef IG_UGBASE_LIB_ALGEBRA_PRALLELIZATION_PARALLEL_COLORING_H
#define IG_UGBASE_LIB_ALGEBRA_PRALLELIZATION_PARALLEL_COLORING_H

#include <map>
#include <vector>
#include <set>

#include "parallelization_util.h"

namespace ug {

extern DebugID DBG_COLORING;

class ParallelColoring
{
public:
	ParallelColoring() : p_colors(nullptr), p_processesWithLowerColor(nullptr), p_processesWithHigherColor(nullptr)
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

#endif