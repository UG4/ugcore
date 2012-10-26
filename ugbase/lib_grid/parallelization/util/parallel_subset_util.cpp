// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 16.02.2012 (m,d,y)
 
#include <vector>
#include "parallel_subset_util.h"

using namespace std;

namespace ug{

void UpdateGlobalMaxDimensionOfSubset(ISubsetHandler& sh,
									  const std::string propertyName,
									  pcl::ProcessCommunicator com)
{
//	since we only want to perform one communication, we gather all dims in
//	an array and perform the reduce operation on that.
//	first however we have to communicate the number of subsets, since it may
//	differ between processes (e.g. if a process didn't receive any grid yet).
	int numSubsets = sh.num_subsets();
	numSubsets = com.allreduce(numSubsets, PCL_RO_MAX);

	vector<int> dims(numSubsets, -1);

	for(int i = 0; i < sh.num_subsets(); ++i){
		int dim = -1;
		if(sh.contains_volumes(i))
			dim = 3;
		else if(sh.contains_faces(i))
			dim = 2;
		else if(sh.contains_edges(i))
			dim = 1;
		else if(sh.contains_vertices(i))
			dim = 0;

		dims[i] = dim;
	}

//	globally communicate dim
	vector<int> globalDims;
	globalDims.resize(dims.size());

	com.allreduce(&dims.front(), &globalDims.front(), dims.size(), PCL_RO_MAX);

//	now set the properties
	for(int i = 0; i < sh.num_subsets(); ++i)
		sh.subset_info(i).set_property(propertyName, globalDims[i]);
}

}// end of namespace
