/*
 * connection_viewer_output.cpp
 *
 *  Created on: 03.08.2011
 *      Author: mrupp
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

#ifdef UG_PARALLEL
string GetParallelName(string name, const pcl::ProcessCommunicator &pc, bool bWriteHeader)
{
	if(pcl::GetNumProcesses() == 1)
		return name;

	char buf[20];
	int rank = pcl::GetProcRank();

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
