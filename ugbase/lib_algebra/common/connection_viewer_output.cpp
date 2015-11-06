/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include <string>
#include <fstream>
#include <vector>
#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

using namespace std;

namespace ug
{
namespace ConnectionViewer
{
bool g_parallelConnectionViewer=false;

#ifdef UG_PARALLEL
string GetParallelName(string name, const pcl::ProcessCommunicator &pc, bool bWriteHeader)
{
	if(pcl::NumProcs() == 1)
		return name;

	char buf[20];
	int rank = pcl::ProcRank();

	size_t iExtPos = name.find_last_of(".");
	string ext = name.substr(iExtPos+1);
	name.resize(iExtPos);
	if(bWriteHeader && pc.size() > 1 && rank == pc.get_proc_id(0))
	{
		size_t iSlashPos = name.find_last_of("/");
		if(iSlashPos == string::npos) iSlashPos = 0; else iSlashPos++;
		string name2 = name.substr(iSlashPos);
		fstream file((name+".p"+ext).c_str(), ios::out);
		file << pc.size() << "\n";
		for(size_t i=0; i<pc.size(); i++)
		{
			sprintf(buf, "_p%04d.%s", pc.get_proc_id(i), ext.c_str());
			file << name2 << buf << "\n";
		}
	}

	sprintf(buf, "_p%04d.%s", rank, ext.c_str());
	name.append(buf);
	return name;
}
#endif


bool AddMarkers(string filename, string markfilename)
{
	fstream fmat(filename.c_str(), ios::out | ios::app);
	fmat << "c " << markfilename << "\n";
    return true;
}

bool WriteMarkers(string markfilename, vector<bool> markers, float r, float g, float b, float alpha, int size)
{
	fstream fmark(markfilename.c_str(), ios::out);
	fmark << r << " " << g << " " << b  << " " << alpha << " " << size << "\n";
	for(size_t i=0; i<markers.size(); i++)
		if(markers[i])
			fmark << i << "\n";
    return true;
}

/*
bool GetNextGrid(vector<MathVector<3> > &coarsegrid, const vector<MathVector<3> > &grid, const vector<int> &newIndex, int nCoarse)
{
	coarsegrid.resize(nCoarse);
	assert(newIndex.size() == grid.size());
	for(size_t i=0; i<grid.size(); i++)
		if(newIndex[i] != -1)
		{
			assert(newIndex[i] < nCoarse && newIndex[i] >= 0);
			coarsegrid[newIndex[i]] = grid[i];
		}
    return true;
}*/

} // CV
} // ug
