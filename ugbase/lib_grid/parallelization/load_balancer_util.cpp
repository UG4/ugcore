/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

	if((static_cast<int>(numLvls) <= minDistLvl) || (maxNumProcs == 0) || (maxNumRedistProcs == 0))
	{
		procH->add_hierarchy_level(0, 1);
		return procH;
	}

//	find the level with the most elements first
	size_t largestLvl = static_cast<size_t>(minDistLvl);

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
	numProcsInvolved = (numProcsInvolved / maxNumRedistProcs) * maxNumRedistProcs;

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
			&& (static_cast<int>(lvl) - lastDistLvl > 1))
		{
			lastDistLvl = static_cast<int>(lvl);
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
