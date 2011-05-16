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
template<typename matrix_type, typename TLocalToGlobal>
void SerializeRow(BinaryBuffer &stream, const matrix_type &mat, size_t localRowIndex, const TLocalToGlobal &localToGlobal)
{
	const AlgebraID &globalRowIndex = localToGlobal[localRowIndex];

	// serialize global row index
	Serialize(stream, globalRowIndex);

	size_t num_connections = mat.num_connections(localRowIndex);

	// serialize number of connections
	Serialize(stream, num_connections);
	UG_DLOG(LIB_ALG_AMG, 4, "Sending row " << localRowIndex << " (" << globalRowIndex << "), " << num_connections << " cons: ");

	for(typename matrix_type::const_row_iterator conn = mat.begin_row(localRowIndex);
						conn != mat.end_row(localRowIndex); ++conn)
	{
		size_t localColIndex = conn.index();
		const AlgebraID &globalColIndex = localToGlobal[localColIndex];
		UG_DLOG(LIB_ALG_AMG, 4, localColIndex << " (" << globalColIndex << ") -> " << conn.value() << " ");

		// serialize connection
		Serialize(stream, globalColIndex);
		Serialize(stream, conn.value());
	}
	UG_DLOG(LIB_ALG_AMG, 4, "\n");
}


template<typename TConnectionType, typename TGlobalToLocal>
size_t DeserializeRow(BinaryBuffer &stream, stdvector<TConnectionType> &cons, const TGlobalToLocal &globalToLocal)
{
	AlgebraID globalRowIndex;

	// serialize global row index
	Deserialize(stream, globalRowIndex);
	size_t localRowIndex = globalToLocal[globalRowIndex];

	UG_DLOG(LIB_ALG_AMG, 4, "Got row " << localRowIndex << " (" << globalRowIndex << "), ");
	size_t num_connections;

	// serialize number of connections
	Deserialize(stream, num_connections);
	
	UG_DLOG(LIB_ALG_AMG, 4, num_connections << " connections: ")

	cons.resize(num_connections);
	for(size_t i =0; i<num_connections; i++)
	{
		AlgebraID globalColIndex;
		Deserialize(stream, globalColIndex);
		cons[i].iIndex = globalToLocal[globalColIndex];
		Deserialize(stream, cons[i].dValue);
		UG_DLOG(LIB_ALG_AMG, 4, cons[i].iIndex << " (" << globalColIndex << ") -> " << cons[i].dValue << " ");
	}
	UG_DLOG(LIB_ALG_AMG, 4, "\n");
	return localRowIndex;
}


// FAMGLevelCalculator::communicate_prolongation
//---------------------------------------------------------------------------
/** send prolongation from Master nodes to Slave nodes.
 */
template<typename matrix_type, typename prolongation_matrix_type, typename vector_type>
void
FAMGLevelCalculator<matrix_type, prolongation_matrix_type, vector_type>::
	communicate_prolongation()
{
	UG_DLOG(LIB_ALG_AMG, 4, "\n********** communicate_prolongation **************\n\n")
	AMG_PROFILE_FUNC();

	//UG_DLOG(LIB_ALG_AMG, 1, "\n\n\nOL2Layout:\n")
	//PrintLayout(A_OL2.get_communicator(), OL2MasterLayout, OL2SlaveLayout);

	// 1. get global IDs (we would only need IDs on OL1)
	std::vector<AlgebraID> localToGlobal;
	GenerateGlobalAlgebraIDs(localToGlobal, A_OL2.num_rows(), OL2MasterLayout, OL2SlaveLayout);
	IF_DEBUG(LIB_ALG_AMG, 4)
	{
		for(size_t i=0; i<localToGlobal.size(); i++)
		{	UG_DLOG(LIB_ALG_AMG, 4, "local " << i << " = " << localToGlobal[i] << "\n");	}
	}
	NewNodesNummerator globalToLocal(localToGlobal);

	// 2. send from master interface of A to slaves rows of P (with global IDs)
	pcl::ParallelCommunicator<IndexLayout> &communicator = (const_cast<matrix_type&>(A)).get_communicator();

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
			SerializeRow(stream, PoldIndices, localIndex, localToGlobal);
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
		
			size_t localRowIndex = DeserializeRow(stream, cons, globalToLocal);
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
					UG_ASSERT(rating.is_slave(localRowIndex, 0), localRowIndex);
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
						UG_ASSERT(rating.is_slave(localColIndex, 0)	|| rating.is_slave(localColIndex, 1), localColIndex << ": " << rating.OL_type(localColIndex));
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

	// 2. add all slave masters to slave interface
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

	// 3. if we interpolate our node from a node which is on another processor (in OL1),
	// we need to add this connection to the slave interface and inform this processor
	BufferMap sendpack;
	stdvector<bool> bInLayout(overlapSize[1], false);
	for(size_t i=0; i<overlapSize[0]; i++)
	{
		for(typename prolongation_matrix_type::row_iterator conn = PoldIndices.begin_row(i); conn != PoldIndices.end_row(i); ++conn)
		{
			size_t localIndex = conn.index();
			UG_ASSERT(localIndex < overlapSize[1], "interpolating from " << localIndex << " which is on OL2 (" << localToGlobal[localIndex] << ")?");
			if(localIndex < overlapSize[0] || bInLayout[localIndex]) continue;

			AlgebraID &globalID = localToGlobal[localIndex];
			UG_ASSERT(globalID.first != pcl::GetProcRank(), globalID);

			bInLayout[localIndex] = true;
			size_t s = globalID.second;
			Serialize(sendpack[globalID.first], s);
			nextLevelSlaveLayout.interface(globalID.first).push_back(localIndex);
			//UG_LOG("informing processor " << globalID.first << " that we have a connection to " << globalID << " (local " << localIndex << ")");
		}
	}


	for(IndexLayout::iterator iter = OL1SlaveLayout.begin(); iter != OL1SlaveLayout.end(); ++iter)
	{
		int pid = OL1SlaveLayout.proc_id(iter);
		BinaryBuffer &stream = sendpack[pid];
		communicator.send_raw(pid, stream.buffer(), stream.write_pos(), false);
	}

	BufferMap receivepack2;
	for(IndexLayout::iterator iter = OL1MasterLayout.begin(); iter != OL1MasterLayout.end(); ++iter)
	{
		int pid = OL1MasterLayout.proc_id(iter);
		communicator.receive_raw(pid, receivepack2[pid]);
	}


	// communicate
	communicator.communicate();

	// process receive buffers
	for(IndexLayout::iterator iter = OL1MasterLayout.begin(); iter != OL1MasterLayout.end(); ++iter)
	{
		int pid = OL1MasterLayout.proc_id(iter);
		BinaryBuffer &stream = receivepack2[pid];
		if(stream.eof()) continue;
		IndexLayout::Interface &nextLevelMasterInterface = nextLevelMasterLayout.interface(pid);
		while(!stream.eof())
		{
			size_t localIndex;
			Deserialize(stream, localIndex);
			//UG_LOG("got information from processor " << pid << " that he has a slave of my " << localIndex << "\n");
			nextLevelMasterInterface.push_back(localIndex);
		}
	}

	IF_DEBUG(LIB_ALG_AMG, 4)
	{
		UG_LOG("\n\nnextLevelMasterLayout:\n");
		PrintLayout(nextLevelMasterLayout);
		UG_LOG("\n\n nextLevel Layouts:\n")
		UG_ASSERT(PrintLayout(communicator, nextLevelMasterLayout, nextLevelSlaveLayout), "layout broken");
	}
}

}
#endif /* FAMG_COMMUNICATE_PROLONGATION_H_ */
