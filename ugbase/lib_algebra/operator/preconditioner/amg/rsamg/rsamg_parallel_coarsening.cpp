/*
 * rsamg_parallel_coarsening.cpp
 *
 *  Created on: 12.10.2011
 *      Author: mrupp
 */

#ifdef UG_PARALLEL

#include <iostream>

#include "common/assert.h"
#include "common/log.h"

using namespace std;

#include "rsamg_coarsening.h"
#include "rsamg_parallel_coarsening.h"

#include "lib_algebra/common/stl_debug.h"
#include "../amg_profiling.h"


#include "lib_algebra/parallelization/parallel_coloring.h"
#include "../send_interface.h"

namespace ug{

void FullSubdomainBlocking(const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();
	size_t N = nodes.size();
	for(size_t i=0; i<N; i++)
	{
		if(!nodes.is_inner(i))
			nodes.set_coarse(i);
	}

	// adjust ratings
	for(size_t i=0; i<N; i++)
	{
		if(!nodes.is_inner(i))
		{
			RemoveUnassignedNeighbors(graphST, PQ, nodes, i);
			MarkUnassignedNeighborsFine(graphST, PQ, nodes, i);
		}
	}
	Coarsen(graphST, PQ, nodes);
}


void ColoringCoarsen(pcl::ParallelCommunicator<IndexLayout> &communicator,
		IndexLayout &OLCoarseningSendLayout, IndexLayout &OLCoarseningReceiveLayout,
		const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes)
{
	AMG_PROFILE_FUNC();

	// Color
	std::set<int> pidsOL;
	std::vector<int> processesWithLowerColor, processesWithHigherColor;
	for(IndexLayout::iterator iter = OLCoarseningSendLayout.begin(); iter != OLCoarseningSendLayout.end(); ++iter)
	{
		if(OLCoarseningSendLayout.interface(iter).size())
			pidsOL.insert(OLCoarseningSendLayout.proc_id(iter));
	}
	for(IndexLayout::iterator iter = OLCoarseningReceiveLayout.begin(); iter != OLCoarseningReceiveLayout.end(); ++iter)
	{
		if(OLCoarseningReceiveLayout.interface(iter).size())
			pidsOL.insert(OLCoarseningReceiveLayout.proc_id(iter));
	}

	// int m_myColor =
	ColorProcessorGraph(communicator, pidsOL, processesWithLowerColor, processesWithHigherColor);

	// receive from processes with lower color

	StdArrayCommunicationScheme<AMGNodes> nodesCommunication(nodes);
	ReceiveOnInterfaces(communicator, processesWithHigherColor, OLCoarseningReceiveLayout, nodesCommunication);

	// update ratings
	for(size_t i=0; i<nodes.size(); i++)
	{
		if(nodes[i] == AMG_COARSE_RATING)
		{
			RemoveUnassignedNeighbors(graphST, PQ, nodes, i);
			MarkUnassignedNeighborsFine(graphST, PQ, nodes, i);
		}
		else if(nodes[i] == AMG_FINE_RATING)
			UpdateNeighborsOfFineNode(graphST, PQ, nodes, i);
	}

	// coarsen
	Coarsen(graphST, PQ, nodes);

	// send coarsening data
	SendOnInterfaces(communicator, processesWithLowerColor, OLCoarseningSendLayout, nodesCommunication);
}


}


#endif /* UG_PARALLEL */
