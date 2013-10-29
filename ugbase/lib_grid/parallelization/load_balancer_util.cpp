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
					   size_t maxNumProcs)
{
	using std::min;

//	we currently can't distribute on lvl 1 (well - we can - but LU won't work then)
	const int minDistLvl = 1;

	SPProcessHierarchy procH = ProcessHierarchy::create();

	if(minNumElemsPerProcPerLvl < 1)
		minNumElemsPerProcPerLvl = 1;

	if(maxNumRedistProcs > maxNumProcs)
		maxNumRedistProcs = maxNumProcs;

	if((numLvls == 0) || (maxNumProcs == 0) || (maxNumRedistProcs == 0)){
		procH->add_hierarchy_level(0, 1);
		return procH;
	}

//	find the level with the most elements first
	size_t largestLvl = 0;

	for(size_t i = 1; i < numLvls; ++i){
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
	const int qualityRedistMinNumLvls = 5;
	int numQualityRedists = (numLvls - lastDistLvl) / qualityRedistMinNumLvls;

	if(numQualityRedists > 0){
		int redistStepSize = (numLvls - lastDistLvl) / (numQualityRedists + 1);
		if(redistStepSize < 2)
			redistStepSize = 2;

		for(int i = 1; i <= numQualityRedists; ++i){
			int lvl = lastDistLvl + i * redistStepSize;
			procH->add_hierarchy_level(lvl, 1);
		}
	}

	return procH;
}

}// end of namespace
