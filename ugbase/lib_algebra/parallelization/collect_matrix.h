/*
 * collect_matrix.h
 *
 *  Created on: 02.05.2011
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__COLLECT_MATRIX_H_
#define __H__LIB_ALGEBRA__COLLECT_MATRIX_H_

#include "parallel_nodes.h"
#include "serialize_interfaces.h"


namespace ug{

template<typename matrix_type>
void SerializeRow(BinaryBuffer &stream, const matrix_type &mat, size_t localRowIndex, ParallelNodes &PN)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	const AlgebraID &globalRowIndex = PN.local_to_global(localRowIndex);

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
		const AlgebraID &globalColIndex = PN.local_to_global(localColIndex);
		UG_DLOG(LIB_ALG_AMG, 4, localColIndex << " (" << globalColIndex << ") -> " << conn.value() << " ");

		// serialize connection
		Serialize(stream, globalColIndex);
		Serialize(stream, conn.value());
	}
	UG_DLOG(LIB_ALG_AMG, 4, "\n");
}


template<typename matrix_type>
void SendMatrix(const matrix_type &A, IndexLayout &verticalSlaveLayout,	int destproc, ParallelNodes &PN)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** SendMatrix ************\n\n");

	pcl::InterfaceCommunicator<IndexLayout> &communicator = A.layouts()->comm();
	BinaryBuffer stream;

	Serialize(stream, A.num_rows());
	for(size_t i=0; i<A.num_rows(); i++)
		SerializeRow(stream, A, i, PN);

	SerializeLayout(stream, A.layouts()->master(), PN);
	SerializeLayout(stream, A.layouts()->slave(), PN);

	IndexLayout::Interface &verticalInterface = verticalSlaveLayout.interface(destproc);
	for(size_t i=0; i<A.num_rows(); i++)
		verticalInterface.push_back(i);

	UG_DLOG(LIB_ALG_AMG, 3, "Srcproc " << pcl::GetProcRank() << " is sending " << stream.write_pos() << " bytes of data to destproc " << destproc << "\n");
	communicator.send_raw(destproc, stream.buffer(), stream.write_pos(), false);
	communicator.communicate();
}

template<typename TConnectionType>
size_t DeserializeRow(BinaryBuffer &stream, stdvector<TConnectionType> &cons, ParallelNodes &PN)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	AlgebraID globalRowIndex;

	// serialize global row index
	Deserialize(stream, globalRowIndex);
	size_t localRowIndex = PN.global_to_local(globalRowIndex);

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
		cons[i].iIndex = PN.global_to_local(globalColIndex);
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
		ParallelNodes &PN)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** ReceiveMatrix ************\n\n");
	pcl::InterfaceCommunicator<IndexLayout> &communicator = A.layouts()->comm();

	M = A;
	M.set_layouts(SmartPtr<HorizontalAlgebraLayouts>(new HorizontalAlgebraLayouts));
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

			size_t localRowIndex = PN.get_local_index_or_create_new(globalRowIndex, 0);
			verticalInterface.push_back(localRowIndex);
			UG_DLOG(LIB_ALG_AMG, 4, "Got row " << localRowIndex << " (" << globalRowIndex << "), ");
			UG_DLOG(LIB_ALG_AMG, 4, num_connections << " connections: ");

			for(size_t pid =0; pid<num_connections; pid++)
			{
				Deserialize(stream, globalColIndex);
				Deserialize(stream, con.dValue);

				con.iIndex = PN.get_local_index_or_create_new(globalColIndex, 0);
				UG_DLOG(LIB_ALG_AMG, 4, con.iIndex << " (" << globalColIndex << ") -> " << con.dValue << " ");
			}
			UG_DLOG(LIB_ALG_AMG, 4, "\n");
		}
	}

	M.resize(PN.local_size(), PN.local_size());

	for(size_t i=0; i<srcprocs.size(); i++)
	{
		int pid = srcprocs[i];
		BinaryBuffer &stream = streams[pid];
		stream.set_read_pos(0);
		stdvector<typename matrix_type::connection> cons;

		Deserialize(stream, numRows);
		for(size_t i=0; i<numRows; i++)
		{
			size_t localRowIndex = DeserializeRow(stream, cons, PN);
			if(cons.size())
				M.add_matrix_row(localRowIndex, &cons[0], cons.size());
		}
	}

	//UG_DLOG(LIB_ALG_AMG, 4, "\n** the matrix M: \n\n");
	//M.print();
	//UG_DLOG(LIB_ALG_AMG, 4, "\n");

	//UG_LOG("COLLECTED LAYOUT:\n");
	//PrintLayout(processCommunicator, communicator, masterLayout, slaveLayout);
}

template<typename matrix_type>
void collect_matrix(matrix_type &A, matrix_type &M, IndexLayout &masterLayout, IndexLayout &slaveLayout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** SendMatrix ************\n\n");
	std::vector<int> srcprocs;
	const pcl::ProcessCommunicator &pc = A.layouts()->proc_comm();

	ParallelNodes PN(A.layouts(), A.num_rows());

	if(pcl::GetProcRank() == pc.get_proc_id(0))
	{
		srcprocs.resize(pc.size()-1);
		for(size_t i=1; i<pc.size(); i++) srcprocs[i-1] = pc.get_proc_id(i);
		ReceiveMatrix(A, M, masterLayout, srcprocs, PN);
	}
	else
		SendMatrix(A, slaveLayout, pc.get_proc_id(0), PN);
}

} // namespace ug

#endif /* __H__LIB_ALGEBRA__COLLECT_MATRIX_H_ */
