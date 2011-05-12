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

template<typename matrix_type>
void collect_matrix(const matrix_type &A, matrix_type &M, IndexLayout &masterLayout, IndexLayout &slaveLayout)
{
	UG_DLOG(LIB_ALG_AMG, 4, "\n*********** collect_matrix ************\n\n");
	pcl::ParallelCommunicator<IndexLayout> &communicator = (const_cast<matrix_type&>(A)).get_communicator();
	std::vector<AlgebraID> localToGlobal;
	GenerateGlobalAlgebraIDs(localToGlobal, A.num_rows(), A.get_master_layout(), A.get_slave_layout());
	NewNodesNummerator globalToLocal(localToGlobal);

	//UG_DLOG(LIB_ALG_AMG, 4, "\n** the matrix A: \n\n");
	//A.print();
	//UG_LOG("\n");

	if(pcl::GetProcRank() != 0)
	{
		BinaryBuffer stream;
		for(size_t i=0; i<A.num_rows(); i++)
			SerializeRow(stream, A, i, localToGlobal);

		IndexLayout::Interface &interface = slaveLayout.interface(0);
		for(size_t i=0; i<A.num_rows(); i++)
			interface.push_back(i);

		UG_DLOG(LIB_ALG_AMG, 4, "sending " << stream.write_pos() << " size to 0\n");
		communicator.send_raw(0, stream.buffer(), stream.write_pos(), false);
		communicator.communicate();


	}
	else
	{
		M = A;
		typedef std::map<int, BinaryBuffer> BufferMap;
		BufferMap streams;

		for(int i=1; i<pcl::GetNumProcesses(); i++)
			communicator.receive_raw(i, streams[i]);
		communicator.communicate();

		for(int i=1; i<pcl::GetNumProcesses(); i++)
		{
			// copy stream (bug in BinaryStream, otherwise I would use stream.seekg(ios::begin)
			BinaryBuffer stream; stream.write((const char*)streams[i].buffer(), streams[i].write_pos());

			UG_DLOG(LIB_ALG_AMG, 4, "received " << stream.write_pos() << " bytes of data from process " << i << "\n");
			IndexLayout::Interface &interface = masterLayout.interface(i);
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
				for(size_t i =0; i<num_connections; i++)
				{
					AlgebraID globalColIndex;
					Deserialize(stream, globalColIndex);
					cons[i].iIndex = globalToLocal.get_index_or_create_new(globalColIndex);
					Deserialize(stream, cons[i].dValue);
					UG_DLOG(LIB_ALG_AMG, 4, cons[i].iIndex << " (" << globalColIndex << ") -> " << cons[i].dValue << " ");
				}
				UG_DLOG(LIB_ALG_AMG, 4, "\n");
			}
		}

		M.resize(globalToLocal.get_new_indices_size(),  globalToLocal.get_new_indices_size());

		for(int i=1; i<pcl::GetNumProcesses(); i++)
		{
			// copy stream (bug in BinaryStream, otherwise I would use stream.seekg(ios::begin)
			BinaryBuffer &stream = streams[i];
			stdvector<typename matrix_type::connection> cons;
			while(!stream.eof())
			{
				size_t localRowIndex = DeserializeRow(stream, cons, globalToLocal);
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

} // namespace ug

#endif /* __H__LIB_ALGEBRA__COLLECT_MATRIX_H_ */
