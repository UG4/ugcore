/*
 * famg_communicate_prolongation.h
 *
 *  Created on: 27.04.2011
 *      Author: mrupp
 */

#ifndef FAMG_COMMUNICATE_PROLONGATION_H_
#define FAMG_COMMUNICATE_PROLONGATION_H_
#include "pcl/pcl.h"

namespace ug
{

/**
 *
 *
 * send prolongation from masters to slaves
 *
 */

// FAMGLevelCalculator::communicate_prolongation
//---------------------------------------------------------------------------
/** send prolongation from Master nodes to Slave nodes.
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	communicate_prolongation()
{
	pcl::ParallelCommunicator<IndexLayout> &communicator = (const_cast<matrix_type&>(A)).get_communicator();

	UG_DLOG(LIB_ALG_AMG, 4, "\n********** communicate_prolongation **************\n\n")
	AMG_PROFILE_FUNC();

	UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, m_famg.iDebugLevelCommunicateProlongation);

	//UG_DLOG(LIB_ALG_AMG, 1, "\n\n\nOL2Layout:\n")
	//PrintLayout(A_OL2.get_communicator(), OL2MasterLayout, OL2SlaveLayout);

	// 1. get global IDs (we would only need IDs on OL1)

	ParallelNodes PN(A_OL2.get_communicator(), A_OL2.get_master_layout(), A_OL2.get_slave_layout(), A_OL2.num_rows());

	// 2. send from master interface of A to slaves rows of P (with global IDs)

	for(IndexLayout::iterator iter = A.get_master_layout().begin(); iter != A.get_master_layout().end(); ++iter)
	{
		IndexLayout::Interface &interface = A.get_master_layout().interface(iter);
		int pid = A.get_master_layout().proc_id(iter);
		UG_DLOG(LIB_ALG_AMG, 4, "sending to pid " << pid << "\n");
		BinaryBuffer stream;

		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t localIndex = interface.get_element(iter);
			if(PoldIndices.num_connections(localIndex)==0)
				continue;
			SerializeRow(stream, PoldIndices, localIndex, PN);
		}
		communicator.send_raw(pid, stream.buffer(), stream.write_pos(), false);
	}

	// 3. set receive buffer
	typedef std::map<int, BinaryBuffer> BufferMap;
	BufferMap receivepack;

	for(IndexLayout::iterator iter = A.get_slave_layout().begin(); iter != A.get_slave_layout().end(); ++iter)
	{
		int pid = A.get_slave_layout().proc_id(iter);
		communicator.receive_raw(pid, receivepack[pid]);
	}

	// 4. communicate
	communicator.communicate();

	// 5. process receive buffers: add rows of P.
	for(IndexLayout::iterator iter = A.get_slave_layout().begin(); iter != A.get_slave_layout().end(); ++iter)
	{
		int pid = A.get_slave_layout().proc_id(iter);
		UG_DLOG(LIB_ALG_AMG, 4, "receiving from pid " << pid << "\n");
		BinaryBuffer &stream = receivepack[pid];

		stdvector<typename matrix_type::connection> cons;
		while(!stream.eof())
		{
		
			size_t localRowIndex = DeserializeRow(stream, cons, PN);
			UG_ASSERT(rating.i_must_assign(localRowIndex) != true, "i must assign is true for " << localRowIndex << "?");
			//UG_ASSERT(rating[localRowIndex].is_fine(), "node " << localRowIndex << " is not fine!");

			if(cons.size() == 1)
			{
				if(rating[localRowIndex].is_coarse() == false)
				{
					UG_DLOG(LIB_ALG_AMG, 4, "post-setted " << localRowIndex << " coarse.\n");
					rating.external_set_coarse(localRowIndex);
				}
			}
			else
			{
				if(rating[localRowIndex].is_fine() == false)
				{
					UG_ASSERT(rating.is_slave(localRowIndex), localRowIndex);
					rating.set_fine(localRowIndex);
					UG_DLOG(LIB_ALG_AMG, 4, "post-setted " << localRowIndex << " fine.\n");
				}

				for(size_t i=0; i<cons.size(); i++)
				{
					size_t localColIndex = cons[i].iIndex;
					//UG_ASSERT(rating[cons[i].iIndex].is_coarse(), "node " << cons[i].iIndex << " is not coarse");
					if( //localToGlobal[localColIndex].first == pcl::GetProcRank() &&
							rating[localColIndex].is_coarse() == false)
					{
						//UG_ASSERT(rating.is_slave(localColIndex)
							//	|| rating.is_slave(localColIndex, 1), localColIndex << ": " << rating.OL_type(localColIndex));
						rating.external_set_coarse(localColIndex);
						UG_DLOG(LIB_ALG_AMG, 4, "post-setted " << localColIndex << " coarse.\n");
					}
				}
			}
			PoldIndices.set_matrix_row(localRowIndex, &cons[0], cons.size());
		}
	}

	// create interfaces
	// 1. add all coarse masters to master interface
	for(IndexLayout::iterator iter = A.get_master_layout().begin(); iter != A.get_master_layout().end(); ++iter)
	{
		IndexLayout::Interface &interface = A.get_master_layout().interface(iter);
		int pid = A.get_master_layout().proc_id(iter);
		UG_DLOG(LIB_ALG_AMG, 4, "creating master interface to processor " << pid << ": ");

		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{

			size_t localIndex = interface.get_element(iter);
			if(rating[localIndex].is_coarse())
			{
				IndexLayout::Interface &nextLevelMasterInterface = nextLevelMasterLayout.interface(pid);
				nextLevelMasterInterface.push_back(localIndex);
				UG_DLOG(LIB_ALG_AMG, 4, localIndex << " ");
			}
		}
		UG_DLOG(LIB_ALG_AMG, 4, "\n");
	}

	// 2. add all coarse slaves to slave interface
	for(IndexLayout::iterator iter = A.get_slave_layout().begin(); iter != A.get_slave_layout().end(); ++iter)
	{
		IndexLayout::Interface &interface = A.get_slave_layout().interface(iter);
		int pid = A.get_slave_layout().proc_id(iter);
		UG_DLOG(LIB_ALG_AMG, 4, "creating Slave interface to processor " << pid << ": ");

		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t localIndex = interface.get_element(iter);
			if(rating[localIndex].is_coarse())
			{
				IndexLayout::Interface &nextLevelSlaveInterface = nextLevelSlaveLayout.interface(pid);
				nextLevelSlaveInterface.push_back(localIndex);
				UG_DLOG(LIB_ALG_AMG, 4, localIndex << " ");
			}
		}
		UG_DLOG(LIB_ALG_AMG, 4, "\n");
	}

	// 3. if we interpolate our node from a node which is master on another processor
	// we need to add this connection to the slave interface and inform this processor
	// this can happen on an number of ways:
	// 1. master0 node interpolating from slave0 node
	// 2. slave0 (original slave) node interpolating from slave1 node

	BufferMap sendpack;
	stdvector<bool> bInLayout(overlapSize[1], false);
	for(size_t i=0; i<overlapSize[0]; i++)
	{
		for(typename prolongation_matrix_type::row_iterator conn = PoldIndices.begin_row(i); conn != PoldIndices.end_row(i); ++conn)
		{
			size_t localIndex = conn.index();
			// master0 and
			if(localIndex < overlapSize[0] || bInLayout[localIndex]) continue;
			//UG_ASSERT(!rating.is_master(i) ||
				//				rating.is_slave(localIndex, 1), "master nodes may only be interpolated by slave0 nodes");

			const AlgebraID &globalID = PN.local_to_global(localIndex);
			UG_ASSERT(globalID.first != pcl::GetProcRank(), globalID);

			bInLayout[localIndex] = true;
			size_t s = globalID.second;
			Serialize(sendpack[globalID.first], s);
			nextLevelSlaveLayout.interface(globalID.first).push_back(localIndex);
			// ug_assert localIndex < overlapSize[1]
			UG_ASSERT(localIndex < overlapSize[1], localIndex << " " << overlapSize[1]);
			UG_DLOG(LIB_ALG_AMG, 2, "informing processor " << globalID.first << " that we have a connection to " << globalID << " (local " << localIndex << ")\n");
		}
	}

	// eigentlich OL1. ne doch OL2, da 
	// ein slave-knoten von einem knoten interpolieren kann,
	// der quasi im OL2 liegt (slave = OL1)
	for(IndexLayout::iterator iter = OL2SlaveLayout.begin(); iter != OL2SlaveLayout.end(); ++iter)
	{
		int pid = OL2SlaveLayout.proc_id(iter);
		BinaryBuffer &stream = sendpack[pid];
		communicator.send_raw(pid, stream.buffer(), stream.write_pos(), false);
	}

	BufferMap receivepack2;
	for(IndexLayout::iterator iter = OL2MasterLayout.begin(); iter != OL2MasterLayout.end(); ++iter)
	{
		int pid = OL2MasterLayout.proc_id(iter);
		communicator.receive_raw(pid, receivepack2[pid]);
	}


	// communicate
	communicator.communicate();

	// process receive buffers
	for(IndexLayout::iterator iter = OL2MasterLayout.begin(); iter != OL2MasterLayout.end(); ++iter)
	{
		int pid = OL2MasterLayout.proc_id(iter);
		BinaryBuffer &stream = receivepack2[pid];
		if(stream.eof()) continue;
		IndexLayout::Interface &nextLevelMasterInterface = nextLevelMasterLayout.interface(pid);
		while(!stream.eof())
		{
			size_t localIndex;
			Deserialize(stream, localIndex);
			UG_DLOG(LIB_ALG_AMG, 4, "got information from processor " << pid << " that he has a slave of my " << localIndex << "\n");
			nextLevelMasterInterface.push_back(localIndex);
		}
	}

	IF_DEBUG(LIB_ALG_AMG, 4)
	{
		UG_LOG("\n\nnextLevelMasterLayout:\n");
		PrintLayout(nextLevelMasterLayout);
		UG_LOG("\n\nnextLevelSlaveLayout:\n");
		PrintLayout(nextLevelSlaveLayout);
		UG_LOG("\n\n nextLevel Layouts:\n")
		UG_ASSERT(PrintLayout(communicator, nextLevelMasterLayout, nextLevelSlaveLayout), "layout broken");
	}
}

}
#endif /* FAMG_COMMUNICATE_PROLONGATION_H_ */
