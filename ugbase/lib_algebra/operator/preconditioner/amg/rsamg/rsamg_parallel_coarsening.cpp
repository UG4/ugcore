/*
 * rsamg_parallel_coarsening.cpp
 *
 *  Created on: 12.10.2011
 *      Author: mrupp
 */

#ifndef UG_PARALLEL
#error "This only works with a UG_PARALLEL define."
#endif

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


void PreventFFConnections(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes);
/**
* \class FullSubdomainBlocking
* \brief Full Subdomain Blocking Coarsening Scheme.
*
* This parallel coarsening scheme sets all border nodes coarse.
* With border nodes we mean nodes which are master or slave (but not "inner")
* \note This coarsening does not need overlap.
*/
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
		PreventFFConnections(graphS, graphST, nodes);
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
	public CommunicationScheme<RSAMGCoarseningCommunicationScheme, char>
{
public:
	typedef char value_type;
	RSAMGCoarseningCommunicationScheme(AMGNodes &nodes) : m_nodes(nodes)
	{
	}

	char send(int pid, int index) const
	{
		return (char) m_nodes[index].rating;
	}

	void receive(int pid, int index, char b)
	{
		if(b == AMG_COARSE_RATING || b == AMG_FINE_RATING
				|| b == AMG_DIRICHLET_RATING)
		{
			if((char) m_nodes[index].rating != b)
				m_recv[index] = b;
		}
	}

	int get_element_size()
	{
		return sizeof(char);
	}

public:
	void clear()
	{
		m_recv.clear();
	}
	void set_coarse_and_fine()
	{
		for(std::map<size_t, char>::iterator it = m_recv.begin(); it != m_recv.end(); ++it)
		{
			char b = it->second;
			int i = it->first;

			if(b == AMG_COARSE_RATING)
				m_nodes.set_coarse(i);
			else if(b == AMG_FINE_RATING)
				m_nodes.set_fine(i);
			else if(b == AMG_DIRICHLET_RATING)
				m_nodes.set_dirichlet(i);
			//UG_LOG("m_recv node " << i << " (" << PN.local_to_global(i) << ") = " << m_nodes[i] << "\n");
		}
	}
	void set_coarse_and_fine(nodeinfo_pq_type &PQ)
	{
		for(std::map<size_t, char>::iterator it = m_recv.begin(); it != m_recv.end(); ++it)
		{
			char b = it->second;
			int i = it->first;

			if(b == AMG_COARSE_RATING)
			{
				m_nodes.set_coarse(i);
				PQ.remove(i);
			}
			else if(b == AMG_FINE_RATING)
			{
				m_nodes.set_fine(i);
				PQ.remove(i);
			}
			//UG_LOG("m_recv node " << i << " (" << PN.local_to_global(i) << ") = " << m_nodes[i] << "\n");
		}
	}

	void change_ratings(const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, bool bUnsymmetric)
	{
		for(std::map<size_t, char>::iterator it = m_recv.begin(); it != m_recv.end(); ++it)
		{
			char b = it->second;
			int i = it->first;
			if(b == AMG_COARSE_RATING)
			{
				RemoveUnassignedNeighbors(graphST, PQ, m_nodes, i);
				MarkUnassignedNeighborsFine(graphST, PQ, m_nodes, i, false);

				if(bUnsymmetric)
					ChangeRatingOfUnassignedNeighbors(graphS, PQ, m_nodes, i, -1);
			}
			else if(b == AMG_FINE_RATING)
			{
				ChangeRatingOfUnassignedNeighbors(graphST, PQ, m_nodes, i, +1);
			}
		}
	}
private:
	AMGNodes &m_nodes;
	std::map<size_t, char> m_recv;
};

/**
* \class SimpleParallelCoarsening
* \brief Parallel Coarsening without border handling
* \note This coarsening does need masterOverlap 1.
*/
class SimpleParallelCoarsening : public IParallelCoarsening
{
	void coarsen(ParallelNodes &PN, stdvector<IndexLayout> vMasterLayouts, stdvector<IndexLayout> vSlaveLayouts,
			const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric)
	{
		AMG_PROFILE_FUNC();

		IndexLayout masterOL1, slaveOL1;
		AddLayout(masterOL1, vMasterLayouts[0]);
		AddLayout(masterOL1, vMasterLayouts[1]);
		AddLayout(slaveOL1, vSlaveLayouts[0]);
		AddLayout(slaveOL1, vSlaveLayouts[1]);

		Coarsen(graphS, graphST, PQ, nodes, bUnsymmetric, false);

		std::map<size_t, char> recv;
		StdArrayCommunicationScheme<AMGNodes> scheme(nodes);
		CommunicateOnInterfaces(PN.get_communicator(), masterOL1, slaveOL1, scheme);
		PreventFFConnections(graphS, graphST, nodes);
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


bool PreventFFConnectionPar(const cgraph &graphS, const cgraph &graphST, AMGNodes &nodes, size_t i, vector<bool> &marks,
		ParallelNodes &PN)
{
	if(nodes[i].is_coarse() || graphS.num_connections(i)==0) return false;

	// mark coarse nodes interpolating this fine node
	for(cgraph::const_row_iterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
	{
		if(nodes[(*it)].is_coarse())
			marks[(*it)] = true;
	}

	bool bCoarse=false;
	// prevent strong F-F connections without common Interpolation node
	for(cgraph::const_row_iterator it = graphST.begin_row(i); it != graphST.end_row(i); ++it)
	{
		if(nodes[*it].is_coarse()) //  || graphS.is_isolated(*it))
			continue;
		const AlgebraID &id = PN.local_to_global(*it);
		if(id.master_proc() > pcl::GetProcRank())
			continue;

		cgraph::const_row_iterator it2 = graphST.begin_row(*it);
		cgraph::const_row_iterator it2end = graphST.end_row(*it);
		for(; it2 != it2end; ++it2)
		{
			if(nodes[*it2].is_coarse() && marks[*it2])
				break;
		}
		if(it2 == it2end)
		{
			nodes.set_coarse(i);
			bCoarse=true;
			break;
		}
	}

	// remove marks
	for(cgraph::const_row_iterator conn = graphST.begin_row(i); conn != graphST.end_row(i); ++conn)
	{
		size_t index = (*conn);
		if(nodes[index].is_coarse())
			marks[index] = false;
	}
	return bCoarse;
}

/**
* \class RS3Coarsening
* \brief RS3 Parallel Coarsening
* First we do a normal coarsening on all elements (inner+master).
* After this step we know about all coarse/fine on OL0 and OL1
* Then we do a post-processing step: We eleminate all strong Fine-Fine-Connections
* then we exchange data on OL0 and OL1.
* \note This coarsening does need masterOverlap 1.
*/
class RS3Coarsening : public IParallelCoarsening
{
	void coarsen(ParallelNodes &PN, stdvector<IndexLayout> vMasterLayouts, stdvector<IndexLayout> vSlaveLayouts,
			const cgraph &graphS, const cgraph &graphST, nodeinfo_pq_type &PQ, AMGNodes &nodes, bool bUnsymmetric)
	{
		AMG_PROFILE_FUNC();

		// first do simple coarsening
		IndexLayout masterOL1, slaveOL1;
		AddLayout(masterOL1, vMasterLayouts[0]);
		AddLayout(masterOL1, vMasterLayouts[1]);
		AddLayout(slaveOL1, vSlaveLayouts[0]);
		AddLayout(slaveOL1, vSlaveLayouts[1]);

		size_t N = graphS.size();

		/*for(size_t i=0; i< N; i++)
		{
			if(nodes.is_master_or_inner(i)==false) continue;
			if(graphS.is_isolated(i))
				nodes.set_coarse(i);
		}*/
		nodes.print();

		Coarsen(graphS, graphST, PQ, nodes, bUnsymmetric, false);

		nodes.print();

		StdArrayCommunicationScheme<AMGNodes> scheme(nodes);
		CommunicateOnInterfaces(PN.get_communicator(), masterOL1, slaveOL1, scheme);

		nodes.print();

		// prevent strong FF-Connections

		int nrOfFFCoarseNodes=0;

		vector<bool> marks(N, false);

		for(size_t i=0; i< N; i++)
		{
			if(nodes.is_master_or_inner(i)==false) continue;
			bool b = PreventFFConnectionPar(graphS, graphST, nodes, i, marks, PN);
			if(b) nrOfFFCoarseNodes++;
		}

		CommunicateOnInterfaces(PN.get_communicator(), masterOL1, slaveOL1, scheme);
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
	return new RS3Coarsening;
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

/**
* \class ColorCoarsening
* \brief Parallel Coarsening by coloring
* First we create implicitely a graph G with (i,j) in G iff i and j are connected in OL0 or OL1.
* Then we color this graph
* Then all processors with color 0 may coarsen
* They send their data to processes with color 1
* then they do the coarsening and so on.
* \note This coarsening does need masterOverlap 1.
*/
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
		RSAMGCoarseningCommunicationScheme nodesCommunication(nodes);
		ReceiveOnInterfaces(PN.get_communicator(), processesWithHigherColor, OLCoarseningReceiveLayout, nodesCommunication);

		nodesCommunication.set_coarse_and_fine(PQ);
		nodesCommunication.change_ratings(graphS, graphST, PQ, bUnsymmetric);

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

		//UG_LOG("phase2\n");
		nodesCommunication.clear();
		SendOnInterfaces(PN.get_communicator(), processesWithHigherColor, vMasterLayouts[1], nodesCommunication);
		ReceiveOnInterfaces(PN.get_communicator(), processesWithLowerColor, vSlaveLayouts[1], nodesCommunication);
		nodesCommunication.set_coarse_and_fine();

		PreventFFConnections(graphS, graphST, nodes);
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

