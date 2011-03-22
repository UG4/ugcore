// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 17.03.2011 (m,d,y)
 
#include "pcl_util.h"

namespace pcl{

void CommunicateInvolvedProcesses(std::vector<int>& vReceiveFromRanksOut,
								  std::vector<int>& vSendToRanks,
								  const ProcessCommunicator& procComm)
{
	using namespace std;

	vReceiveFromRanksOut.clear();

//	if there is only one process in the communicator, there's
//	nothing to do.
	if(procComm.size() < 2)
		return;

	const int localProcRank = GetProcRank();

//	we'll use an array in which we'll store the number of
//	processes, to which each process wants to talk.
	vector<int> vNumAssProcs(procComm.size());

//	perform allgather with proc-numbers
	int procCount = (int)vSendToRanks.size();

	procComm.allgather(&procCount, 1, PCL_DT_INT,
					   &vNumAssProcs.front(),
					   1,
					   PCL_DT_INT);

//	sum the number of processes so that the absolute list size can
//	be determined.
	size_t listSize = 0;
	vector<int> vDisplacements(vNumAssProcs.size());
	for(size_t i = 0; i < vNumAssProcs.size(); ++i){
		vDisplacements[i] = (int)listSize;
		listSize += vNumAssProcs[i];
	}

//	perform allgather with proc-lists
//	this list will later on contain an adjacency-list for each proc.
//	adjacency in this case means, that processes want to communicate.
	vector<int> vGlobalProcList(listSize);

	procComm.allgatherv(&vSendToRanks.front(),
						procCount,
						PCL_DT_INT,
						&vGlobalProcList.front(),
						&vNumAssProcs.front(),
						&vDisplacements.front(),
						PCL_DT_INT);

//	we can now check for each process whether it wants to
//	communicate with this process. Simply iterate over
//	the adjacency list to do so.
	for(int i = 0; i < (int)vNumAssProcs.size(); ++i)
	{
		if(i != localProcRank){
		//	check whether the i-th proc wants to communicate with the local proc
			for(int j = 0; j < vNumAssProcs[i]; ++j)
			{
				if(vGlobalProcList[vDisplacements[i] + j] == localProcRank)
				{
					vReceiveFromRanksOut.push_back(procComm.get_proc_id(i));
				//	the j-th proc is handled completly. resume with the next.
					break;
				}
			}
		}
	}

//	vProcRanksInOut should now contain all process-ranks with which
//	the local proc should communicate.
}

}// end of namespace
