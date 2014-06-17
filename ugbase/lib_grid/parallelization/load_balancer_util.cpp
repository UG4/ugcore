// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Oct 17, 2013 (d,m,y)

#include <cmath>
#include "load_balancer_util.h"

using namespace std;

namespace ug{

SPProcessHierarchy
CreateProcessHierarchy(size_t* numElemsOnLvl, size_t numLvls,
					   size_t minNumElemsPerProcPerLvl, size_t maxNumRedistProcs,
					   size_t maxNumProcs, int minDistLvl,
					   int maxLvlsWithoutRedist)
{
	using std::min;

	SPProcessHierarchy procH = ProcessHierarchy::create();

	if(minDistLvl < 0)
		minDistLvl = 0;

	if(minNumElemsPerProcPerLvl < 1)
		minNumElemsPerProcPerLvl = 1;

	if(maxNumRedistProcs > maxNumProcs)
		maxNumRedistProcs = maxNumProcs;

	if(((int)numLvls <= minDistLvl) || (maxNumProcs == 0) || (maxNumRedistProcs == 0))
	{
		procH->add_hierarchy_level(0, 1);
		return procH;
	}

//	find the level with the most elements first
	size_t largestLvl = (size_t)minDistLvl;

	for(size_t i = minDistLvl + 1; i < numLvls; ++i){
		if(numElemsOnLvl[i] > numElemsOnLvl[largestLvl])
			largestLvl = i;
	}

//	get the number of processes onto which the largest level shall be distributed
	size_t maxNumElemsTotal = numElemsOnLvl[largestLvl];
	if(maxNumElemsTotal == 0){
		procH->add_hierarchy_level(0, 1);
		return procH;
	}

	size_t numProcsInvolved = min<size_t>(maxNumElemsTotal / minNumElemsPerProcPerLvl,
										  maxNumProcs);

//	The following line only ensures that v-interfaces stay on fixed levels
//todo: remove the following line as soon as v-interfaces may change levels.
	//numProcsInvolved = (numProcsInvolved / maxNumRedistProcs) * maxNumRedistProcs;

//	create the process hierarchy
	size_t curNumProcs = 1;
	int lastDistLvl = minDistLvl - 2;
	for(size_t lvl = (size_t)minDistLvl; lvl < numLvls; ++lvl){
		int numRedistProcs = maxNumRedistProcs;
		if(curNumProcs * numRedistProcs > numProcsInvolved){
			numRedistProcs = numProcsInvolved / curNumProcs;
			if(numRedistProcs <= 1)
				break;
		}

		if((curNumProcs * numRedistProcs * minNumElemsPerProcPerLvl <= numElemsOnLvl[lvl])
			&& ((int)lvl - lastDistLvl > 1))
		{
			lastDistLvl = (int)lvl;
			procH->add_hierarchy_level(lvl, numRedistProcs);
			curNumProcs *= numRedistProcs;
		}
	}

	if(lastDistLvl < 0){
		procH->add_hierarchy_level(0, 1);
		return procH;
	}

//	check whether we have to perform redistributions to ensure quality
	const int qualityRedistMinNumLvls = maxLvlsWithoutRedist + 1;
	int numQualityRedists = (numLvls - lastDistLvl) / qualityRedistMinNumLvls;

	if(numQualityRedists > 0){
//		int redistStepSize = (numLvls - lastDistLvl) / (numQualityRedists + 1);
//		if(redistStepSize < 2)
//			redistStepSize = 2;
	//	the new version avoids unnecessary redistributions
		for(int i = 0; i < numQualityRedists; ++i){
			procH->add_hierarchy_level(lastDistLvl + (i+1) * qualityRedistMinNumLvls, 1);
		}
	}

	return procH;
}

}// end of namespace
