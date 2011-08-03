/*
 * collect_matrix.h
 *
 *  Created on: 02.05.2011
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__COLLECT_MATRIX_H_
#define __H__LIB_ALGEBRA__COLLECT_MATRIX_H_

#include "lib_algebra/parallelization/new_nodes_nummerator.h"

#include "parallelization.h"

namespace ug{

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


template<typename matrix_type>
void SendMatrix(const matrix_type &A, IndexLayout &verticalSlaveLayout,	int destproc, std::vector<AlgebraID> localToGlobal)
{
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** SendMatrix ************\n\n");

	pcl::ParallelCommunicator<IndexLayout> &communicator = (const_cast<matrix_type&>(A)).get_communicator();
	BinaryBuffer stream;

	Serialize(stream, A.num_rows());
	for(size_t i=0; i<A.num_rows(); i++)
		SerializeRow(stream, A, i, localToGlobal);

	SerializeLayout(stream, A.get_master_layout(), localToGlobal);
	SerializeLayout(stream, A.get_slave_layout(), localToGlobal);

	IndexLayout::Interface &verticalInterface = verticalSlaveLayout.interface(0);
	for(size_t i=0; i<A.num_rows(); i++)
		verticalInterface.push_back(i);

	UG_DLOG(LIB_ALG_AMG, 3, "Srcproc " << pcl::GetProcRank() << " is sending " << stream.write_pos() << " bytes of data to destproc " << destproc << "\n");
	communicator.send_raw(destproc, stream.buffer(), stream.write_pos(), false);
	communicator.communicate();
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

	UG_DLOG(LIB_ALG_AMG, 4, num_connections << " connections: ");

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

// ReceiveMatrix
//---------------------------------------------------------------------------
/**
 *	Receives a distributed matrix from several processors
 * \param A				(in) input matrix
 * \param M				(out) collected matrix
 * \param masterLayout	(out) created master layout to processors in srcprocs
 * \param
 * \param srcprocs		list of source processors
 *
 */
template<typename matrix_type>
void ReceiveMatrix(const matrix_type &A, matrix_type &M, IndexLayout &verticalMasterLayout,	const std::vector<int> &srcprocs,
		NewNodesNummerator &globalToLocal)
{
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** ReceiveMatrix ************\n\n");
	pcl::ParallelCommunicator<IndexLayout> &communicator = (const_cast<matrix_type&>(A)).get_communicator();

	M = A;
	IndexLayout *pNull=NULL;
	M.set_master_layout(*pNull);
	M.set_slave_layout(*pNull);
	typedef std::map<int, BinaryBuffer> BufferMap;
	BufferMap streams;

	UG_DLOG(LIB_ALG_AMG, 3, "DestProc " << pcl::GetProcRank() << " is waiting on data from ");
	for(size_t i=0; i<srcprocs.size(); i++)
	{
		UG_DLOG(LIB_ALG_AMG, 3, srcprocs[i] << " ");
		communicator.receive_raw(srcprocs[i], streams[srcprocs[i]]);
	}
	UG_LOG("\n");
	communicator.communicate();

	AlgebraID globalRowIndex, globalColIndex;;
	size_t num_connections, numRows;

	for(size_t i=0; i<srcprocs.size(); i++)
	{
		int pid = srcprocs[i];
		BinaryBuffer &stream = streams[pid];
		stream.set_read_pos(0);

		UG_DLOG(LIB_ALG_AMG, 4, "received " << stream.write_pos() << " bytes of data from process " << pid << "\n");
		IndexLayout::Interface &verticalInterface = verticalMasterLayout.interface(pid);
		typename matrix_type::connection con;

		Deserialize(stream, numRows);
		for(size_t i=0; i<numRows; i++)
		{
			// serialize global row index, number of connections
			Deserialize(stream, globalRowIndex);
			Deserialize(stream, num_connections);

			size_t localRowIndex = globalToLocal.get_index_or_create_new(globalRowIndex);
			verticalInterface.push_back(localRowIndex);
			UG_DLOG(LIB_ALG_AMG, 4, "Got row " << localRowIndex << " (" << globalRowIndex << "), ");
			UG_DLOG(LIB_ALG_AMG, 4, num_connections << " connections: ");

			for(size_t pid =0; pid<num_connections; pid++)
			{
				Deserialize(stream, globalColIndex);
				Deserialize(stream, con.dValue);

				con.iIndex = globalToLocal.get_index_or_create_new(globalColIndex);
				UG_DLOG(LIB_ALG_AMG, 4, con.iIndex << " (" << globalColIndex << ") -> " << con.dValue << " ");
			}
			UG_DLOG(LIB_ALG_AMG, 4, "\n");
		}
	}

	M.resize(globalToLocal.get_new_indices_size(),  globalToLocal.get_new_indices_size());

	for(size_t i=0; i<srcprocs.size(); i++)
	{
		int pid = srcprocs[i];
		BinaryBuffer &stream = streams[pid];
		stream.set_read_pos(0);
		stdvector<typename matrix_type::connection> cons;

		Deserialize(stream, numRows);
		for(size_t i=0; i<numRows; i++)
		{
			size_t localRowIndex = DeserializeRow(stream, cons, globalToLocal);
			if(cons.size())
				M.add_matrix_row(localRowIndex, &cons[0], cons.size());
		}
	}

	//UG_DLOG(LIB_ALG_AMG, 4, "\n** the matrix M: \n\n");
	//M.print();
	//UG_DLOG(LIB_ALG_AMG, 4, "\n");

	//UG_LOG("COLLECTED LAYOUT:\n");
	//PrintLayout(communicator, masterLayout, slaveLayout);
}

template<typename matrix_type>
void collect_matrix(matrix_type &A, matrix_type &M, IndexLayout &masterLayout, IndexLayout &slaveLayout)
{
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** SendMatrix ************\n\n");
	std::vector<int> srcprocs;
	pcl::ProcessCommunicator &pc = A.get_process_communicator();

	std::vector<AlgebraID> localToGlobal;
	GenerateGlobalAlgebraIDs(localToGlobal, A.num_rows(), A.get_master_layout(), A.get_slave_layout());

	if(pcl::GetProcRank() == pc.get_proc_id(0))
	{
		NewNodesNummerator globalToLocal(localToGlobal);
		srcprocs.resize(pc.size()-1);
		for(size_t i=1; i<pc.size(); i++) srcprocs[i-1] = pc.get_proc_id(i);
		ReceiveMatrix(A, M, masterLayout, srcprocs, globalToLocal);
	}
	else
		SendMatrix(A, slaveLayout, pc.get_proc_id(0), localToGlobal);
}

} // namespace ug

#endif /* __H__LIB_ALGEBRA__COLLECT_MATRIX_H_ */
