// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 27.01.2011 (m,d,y)
 
#include "parallelization_util.h"

namespace ug{

void BuildOneToManyLayout(IndexLayout& masterLayoutOut,
						  IndexLayout& slaveLayoutOut,
						  int rootProcID,
						  IndexLayout& masterLayout,
						  IndexLayout& slaveLayout,
						  pcl::ProcessCommunicator procComm)
{
	int numLocalMasterNodes = 0;
	int localRank = procComm.get_local_proc_id();
	std::vector<IndexLayout::Element> masterNodes;

//	create a new slave interface in slaveLayoutOut
	if(!masterLayout.empty()){
//		IndexLayout::Interface& interface = slaveLayoutOut.interface(rootProcID);
	//	now collect all master nodes
		CollectUniqueElements(masterNodes, masterLayout);

	//	write the number of elements to the numMasterNodes vector
		numLocalMasterNodes = (int)masterNodes.size();
	}

//todo:	gather numMasterNodes on process 0
//		create a master layout to all those old master nodes
	std::vector<int> numOldMasterNodes(procComm.size(), 0);
	procComm.gather(&numLocalMasterNodes, 1, PCL_DT_INT,
					&numOldMasterNodes.front(),
					numOldMasterNodes.size(),
					PCL_DT_INT, rootProcID);
}

}// end of namespace
