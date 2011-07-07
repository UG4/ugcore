/*
 * collect_matrix.h
 *
 *  Created on: 02.05.2011
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__COLLECT_MATRIX_H_
#define __H__LIB_ALGEBRA__COLLECT_MATRIX_H_

namespace ug{

template<typename matrix_type, typename TLocalToGlobal>
void SerializeRow(BinaryBuffer &stream, const matrix_type &mat, size_t localRowIndex, const TLocalToGlobal &localToGlobal);
template<typename TConnectionType, typename TGlobalToLocal>
size_t DeserializeRow(BinaryBuffer &stream, stdvector<TConnectionType> &cons, const TGlobalToLocal &globalToLocal);

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


template<typename TConnectionType, typename TGlobalToLocal>
size_t DeserializeRow(BinaryStream &stream, stdvector<TConnectionType> &cons, const TGlobalToLocal &globalToLocal)
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
template<typename matrix_type>
void collect_matrix(const matrix_type &A, matrix_type &M, IndexLayout &masterLayout, IndexLayout &slaveLayout,
		int srcproc, const std::vector<int> &destprocs)
{
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** collect_matrix ************\n\n");

	pcl::ParallelCommunicator<IndexLayout> &communicator = (const_cast<matrix_type&>(A)).get_communicator();
	std::vector<AlgebraID> localToGlobal;
	GenerateGlobalAlgebraIDs(localToGlobal, A.num_rows(), A.get_master_layout(), A.get_slave_layout());
	NewNodesNummerator globalToLocal(localToGlobal);

	/*UG_DLOG(LIB_ALG_AMG, 4, "\n** the matrix A: \n\n");
	A.print();
	UG_LOG("\n");*/
	if(pcl::GetProcRank() != srcproc)
	{

		BinaryBuffer stream;
		for(size_t i=0; i<A.num_rows(); i++)
			SerializeRow(stream, A, i, localToGlobal);

		IndexLayout::Interface &interface = slaveLayout.interface(0);
		for(size_t i=0; i<A.num_rows(); i++)
			interface.push_back(i);

		UG_DLOG(LIB_ALG_AMG, 3, "Destproc " << pcl::GetProcRank() << " is sending " << stream.write_pos() << " bytes of data to sourceproc " <<  srcproc << "\n");
		communicator.send_raw(srcproc, stream.buffer(), stream.write_pos(), false);
		communicator.communicate();
	}
	else
	{
		M = A;
		typedef std::map<int, BinaryBuffer> BufferMap;
		BufferMap streams;

		UG_DLOG(LIB_ALG_AMG, 3, "SourceProc " << srcproc << " is waiting on data from ");
		for(size_t i=0; i<destprocs.size(); i++)
		{
			UG_DLOG(LIB_ALG_AMG, 3, destprocs[i] << " ");
			communicator.receive_raw(destprocs[i], streams[destprocs[i]]);
		}
		UG_LOG("\n");
		communicator.communicate();

		for(size_t i=0; i<destprocs.size(); i++)
		{
			int pid = destprocs[i];
			// copy stream (bug in BinaryStream, otherwise I would use stream.seekg(ios::begin)
			BinaryBuffer stream; stream.write((const char*)streams[pid].buffer(), streams[pid].write_pos());

			UG_DLOG(LIB_ALG_AMG, 4, "received " << stream.write_pos() << " bytes of data from process " << pid << "\n");
			IndexLayout::Interface &interface = masterLayout.interface(pid);
			stdvector<typename matrix_type::connection> cons;
			while(!stream.eof())
			{
				AlgebraID globalRowIndex;

				// serialize global row index
				Deserialize(stream, globalRowIndex);
				size_t localRowIndex = globalToLocal.get_index_or_create_new(globalRowIndex);

				interface.push_back(localRowIndex);

				UG_DLOG(LIB_ALG_AMG, 4, "Got row " << localRowIndex << " (" << globalRowIndex << "), ");
				size_t num_connections;

				// serialize number of connections
				Deserialize(stream, num_connections);

				UG_DLOG(LIB_ALG_AMG, 4, num_connections << " connections: ")

				cons.resize(num_connections);
				for(size_t pid =0; pid<num_connections; pid++)
				{
					AlgebraID globalColIndex;
					Deserialize(stream, globalColIndex);
					cons[pid].iIndex = globalToLocal.get_index_or_create_new(globalColIndex);
					Deserialize(stream, cons[pid].dValue);
					UG_DLOG(LIB_ALG_AMG, 4, cons[pid].iIndex << " (" << globalColIndex << ") -> " << cons[pid].dValue << " ");
				}
				UG_DLOG(LIB_ALG_AMG, 4, "\n");
			}
		}

		M.resize(globalToLocal.get_new_indices_size(),  globalToLocal.get_new_indices_size());

		for(size_t i=0; i<destprocs.size(); i++)
		{
			int pid = destprocs[i];
			BinaryBuffer &stream = streams[pid];
			stdvector<typename matrix_type::connection> cons;
			while(!stream.eof())
			{
				size_t localRowIndex = DeserializeRow(stream, cons, globalToLocal);
				if(cons.size())
					M.add_matrix_row(localRowIndex, &cons[0], cons.size());
			}
		}
		//UG_DLOG(LIB_ALG_AMG, 4, "\n** the matrix M: \n\n");
		//M.print();
		//UG_DLOG(LIB_ALG_AMG, 4, "\n");

	}

	//UG_LOG("COLLECTED LAYOUT:\n");
	//PrintLayout(communicator, masterLayout, slaveLayout);
}

template<typename matrix_type>
void collect_matrix(const matrix_type &A, matrix_type &M, IndexLayout &masterLayout, IndexLayout &slaveLayout)
{
	std::vector<int> destprocs;
	destprocs.resize(pcl::GetNumProcesses()-1);
	for(int i=1; i<pcl::GetNumProcesses(); i++) destprocs[i-1] = i;
	collect_matrix(A, M, masterLayout, slaveLayout, 0, destprocs);
}

} // namespace ug

#endif /* __H__LIB_ALGEBRA__COLLECT_MATRIX_H_ */
