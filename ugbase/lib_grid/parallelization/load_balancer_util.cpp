// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Oct 17, 2013 (d,m,y)

#include <cmath>
#include "load_balancer_util.h"

namespace ug{

SPProcessHierarchy
CreateProcessHierarchy(size_t* numElemsOnLvl, size_t numLvls,
					   size_t minNumElemsPerProcPerLvl, size_t maxNumRedistProcs,
					   size_t maxNumProcs)
{
	using std::min;

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

//	find out how many distribution-levels are required
	size_t numDistLvls = 0;
	size_t curNumProcs = 1;
	while((curNumProcs < numProcsInvolved) && (numDistLvls < numLvls)){
		curNumProcs *= maxNumRedistProcs;
		++numDistLvls;
	}

	if(numDistLvls == 0){
		procH->add_hierarchy_level(0, 1);
		return procH;
	}

//	find out to how many processes we'll distribute on each distribution level
	size_t numRedistProcs = std::pow(numProcsInvolved, 1./(double)numDistLvls) + SMALL;

//	finally, create the process hierarchy
	curNumProcs = 1;
	int lastDistLvl = -2;
	for(size_t lvl = 0; lvl < numLvls; ++lvl){
		if(procH->num_hierarchy_levels() > numDistLvls)
			break;

		if((curNumProcs * numRedistProcs * minNumElemsPerProcPerLvl <= numElemsOnLvl[lvl])
			&& ((int)lvl - lastDistLvl > 1))
		{
			lastDistLvl = (int)lvl;
			procH->add_hierarchy_level(lvl, numRedistProcs);
			curNumProcs *= numRedistProcs;
		}
	}

	return procH;
}

}// end of namespace
