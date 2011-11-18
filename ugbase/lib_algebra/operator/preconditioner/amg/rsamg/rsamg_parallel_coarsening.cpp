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
#include "rsamg_coarsening.h"

namespace ug{

void CreateAllToAllFromMasterSlave(pcl::ParallelCommunicator<IndexLayout> &communicator,
		IndexLayout &OLCoarseningSendLayout, IndexLayout &OLCoarseningReceiveLayout,
		IndexLayout &OL1MasterLayout, IndexLayout &OL1SlaveLayout);



class FullSubdomainBlocking : public IParallelCoarsening
{
public:
	void coarsen(ParallelNodes &PN, stdvector<IndexLayout> vMasterLayouts, stdvector<IndexLayout> vSlaveLayouts,
		const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric)
	{
		AMG_PROFILE_FUNC();
		size_t N = nodes.size();
		for(size_t i=0; i<N; i++)
		{
			if(nodes[i].is_assigned()) continue;

			if(nodes.is_master(i))
			{
				nodes.set_coarse(i);
				PQ.remove(i);
			}
			else if(nodes.is_slave(i))
				nodes.set_coarse(i);
			else if(!nodes.is_inner(i))
				nodes[i].set_parallel_dont_care();
		}
		/* do NOT set neighboring nodes fine.
		 * otherwise we will get a ton of FF-Connections
		 * this COULD be OK if used with aggressive coarsening.
		 * */
		Coarsen(graphS, graphST, PQ, nodes, bUnsymmetric, false);
	}

	int overlap_depth_master()
	{
		return 0;
	}
	int overlap_depth_slave()
	{
		return 0;
	}

	const char *tostring()
	{
		return "FullSubdomainBlockingCoarsening";
	}
};


IParallelCoarsening *GetFullSubdomainBlockingCoarsening()
{
	return new FullSubdomainBlocking;
}




class RSAMGCoarseningCommunicationScheme :
	public CommunicationScheme<RSAMGCoarseningCommunicationScheme>
{
public:
	typedef char value_type;
	RSAMGCoarseningCommunicationScheme(AMGNodes &nodes,
			std::map<size_t, char> &recv) : m_nodes(nodes), m_recv(recv)
	{
	}

	char send(int pid, int index) const
	{
		UG_LOG("send node " << index << " = " << m_nodes[index] << "\n");
		return (char) m_nodes[index].rating;
	}

	void receive(int pid, int index, char b)
	{
		m_recv[index] = b;
	}

	int get_element_size()
	{
		return sizeof(char);
	}

private:
	AMGNodes &m_nodes;
	std::map<size_t, char> &m_recv;
};


class SimpleParallelCoarsening : public IParallelCoarsening
{
	void coarsen(ParallelNodes &PN, stdvector<IndexLayout> vMasterLayouts, stdvector<IndexLayout> vSlaveLayouts,
			const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric)
	{
		AMG_PROFILE_FUNC();

		// Coarsen
		Coarsen(graphS, graphST, PQ, nodes, bUnsymmetric, false);

		std::map<size_t, char> recv;
		StdArrayCommunicationScheme<AMGNodes> scheme(nodes);
		IndexLayout masterOL1, slaveOL1;

		AddLayout(masterOL1, vMasterLayouts[0]);
		AddLayout(masterOL1, vMasterLayouts[1]);

		AddLayout(slaveOL1, vSlaveLayouts[0]);
		AddLayout(slaveOL1, vSlaveLayouts[1]);

		CommunicateOnInterfaces(PN.get_communicator(), slaveOL1, masterOL1, scheme);
	}

	int overlap_depth_master()
	{
		return 1;
	}

	int overlap_depth_slave()
	{
		return 0;
	}

	const char *tostring()
	{
		return "SimpleParallelCoarsening";
	}
};

IParallelCoarsening *GetSimpleParallelCoarsening()
{
	return new SimpleParallelCoarsening;
}


class RS3Coarsening: public IParallelCoarsening
{
	void coarsen(ParallelNodes &PN, stdvector<IndexLayout> vMasterLayouts, stdvector<IndexLayout> vSlaveLayouts,
			const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric)
	{
		AMG_PROFILE_FUNC();

		// Coarsen
		Coarsen(graphS, graphST, PQ, nodes, bUnsymmetric, false);

		std::map<size_t, char> recv;
		StdArrayCommunicationScheme<AMGNodes> scheme(nodes);
		IndexLayout masterOL1, slaveOL1;

		AddLayout(masterOL1, vMasterLayouts[0]);
		AddLayout(masterOL1, vMasterLayouts[1]);

		AddLayout(slaveOL1, vSlaveLayouts[0]);
		AddLayout(slaveOL1, vSlaveLayouts[1]);

		CommunicateOnInterfaces(PN.get_communicator(), slaveOL1, masterOL1, scheme);

		// prevent FF
		int nrOfFFCoarseNodes=0;
		size_t N = graphS.size();
		vector<bool> marks(N, false);

		for(size_t i=0; i< N; i++)
		{
			if(nodes.is_master_or_inner(i)==false)
				continue;
			PreventFFConnection(graphS, graphST, nodes, i, marks, nrOfFFCoarseNodes);
		}

		CommunicateOnInterfaces(PN.get_communicator(), vSlaveLayouts[0], vMasterLayouts[0], scheme);
	}

	int overlap_depth_master()
	{
		return 1;
	}

	int overlap_depth_slave()
	{
		return 0;
	}

	const char *tostring()
	{
		return "SimpleParallelCoarsening";
	}
};


IParallelCoarsening *GetRS3Coarsening()
{
	UG_ASSERT(0, "not implemented yet");
	return NULL;
}


IParallelCoarsening *GetCLJPCoarsening()
{
	UG_ASSERT(0, "not implemented yet");
	return NULL;
}

IParallelCoarsening *GetFalgoutCoarsening()
{
	UG_ASSERT(0, "not implemented yet");
	return NULL;
}


IParallelCoarsening *GetMinimumSubdomainBlockingCoarsening()
{
	UG_ASSERT(0, "not implemented yet");
	return NULL;
}


IParallelCoarsening *GetCoarseGridClassificationCoarsening()
{
	UG_ASSERT(0, "not implemented yet");
	return NULL;
}


class ColorCoarsening : public IParallelCoarsening
{
	void coarsen(ParallelNodes &PN, stdvector<IndexLayout> vMasterLayouts, stdvector<IndexLayout> vSlaveLayouts,
			const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric)
	{
		AMG_PROFILE_FUNC();

		IndexLayout OL1MasterLayout, OL1SlaveLayout;

		AddLayout(OL1MasterLayout, vMasterLayouts[0]);
		AddLayout(OL1SlaveLayout, vSlaveLayouts[0]);
		AddLayout(OL1MasterLayout, vMasterLayouts[1]);
		AddLayout(OL1SlaveLayout, vSlaveLayouts[1]);


		// create slave-slave-connections
		IndexLayout OLCoarseningSendLayout, OLCoarseningReceiveLayout;
		CreateAllToAllFromMasterSlave(PN.get_communicator(),
				OLCoarseningSendLayout, OLCoarseningReceiveLayout, OL1MasterLayout, OL1SlaveLayout);


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
		ColorProcessorGraph(PN.get_communicator(), pidsOL, processesWithLowerColor, processesWithHigherColor);

		// receive from processes with lower color

		std::map<size_t, char> recv;
		RSAMGCoarseningCommunicationScheme nodesCommunication(nodes, recv);
		ReceiveOnInterfaces(PN.get_communicator(), processesWithHigherColor, OLCoarseningReceiveLayout, nodesCommunication);

		/*UG_LOG("nodes.get_unassigned() = " << nodes.get_unassigned() << "\n");
		UG_LOG("PQ.height() = " << PQ.height() << "\n");

		PN.print();
		nodes.print();*/
		PN.print();
		for(std::map<size_t, char>::iterator it = recv.begin(); it != recv.end(); ++it)
		{
			char b = it->second;
			int i = it->first;

			if(b == AMG_COARSE_RATING)
			{
				nodes.set_coarse(i);
				PQ.remove(i);
			}
			else if(b == AMG_FINE_RATING)
			{
				nodes.set_fine(i);
				PQ.remove(i);
			}
			UG_LOG("recv node " << i << " (" << PN.local_to_global(i) << ") = " << nodes[i] << "\n");
		}
		for(std::map<size_t, char>::iterator it = recv.begin(); it != recv.end(); ++it)
		{
			char b = it->second;
			int i = it->first;
			if(b == AMG_COARSE_RATING)
			{
				RemoveUnassignedNeighbors(graphST, PQ, nodes, i);
				MarkUnassignedNeighborsFine(graphST, PQ, nodes, i, false);

				if(bUnsymmetric)
					ChangeRatingOfUnassignedNeighbors(graphS, PQ, nodes, i, -1);
			}
			else if(b == AMG_FINE_RATING)
			{
				ChangeRatingOfUnassignedNeighbors(graphST, PQ, nodes, i, +1);
			}
		}

		/*UG_LOG("nodes.get_unassigned() = " << nodes.get_unassigned() << "\n");
		UG_LOG("PQ.height() = " << PQ.height() << "\n");*/

		/*for(size_t i=0; i<nodes.size(); i++)
		{
			if(!nodes.is_slave(i) && !nodes.is_master(i) && !nodes.is_inner(i))
			{
				pq.remove(i);
				nodes[i].set_parallel_dont_care();
			}
		}*/

		// coarsen
		Coarsen(graphS, graphST, PQ, nodes, bUnsymmetric, false);

		// send coarsening data

		SendOnInterfaces(PN.get_communicator(), processesWithLowerColor, OLCoarseningSendLayout, nodesCommunication);

		UG_LOG("phase2\n");
		recv.clear();
		SendOnInterfaces(PN.get_communicator(), processesWithHigherColor, vMasterLayouts[1], nodesCommunication);
		ReceiveOnInterfaces(PN.get_communicator(), processesWithLowerColor, vSlaveLayouts[1], nodesCommunication);

		for(std::map<size_t, char>::iterator it = recv.begin(); it != recv.end(); ++it)
		{
			char b = it->second;
			int i = it->first;

			if(b == AMG_COARSE_RATING)
			{
				nodes.set_coarse(i);
				UG_LOG("set coarse " << i << "\n");
			}
			else if(b == AMG_FINE_RATING)
			{
				nodes.set_fine(i);
				UG_LOG("set coarse " << i << "\n");
			}
			UG_LOG("recv node " << i << " (" << PN.local_to_global(i) << ") = " << nodes[i] << "\n");
		}
	}

	int overlap_depth_master()
	{
		return 1;
	}

	int overlap_depth_slave()
	{
		return 0;
	}

	const char *tostring()
	{
		return "ColorCoarsening";
	}
};

IParallelCoarsening *GetColorCoarsening()
{
	return new ColorCoarsening;
}

}


#endif /* UG_PARALLEL */
