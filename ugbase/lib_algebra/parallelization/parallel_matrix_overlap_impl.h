/*
 * parallelization.h
 *
 *  Created on: 17.01.2010
 *      Author: Martin Rupp / Sebastian Reiter
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_OVERLAP_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_OVERLAP_IMPL__

#include "parallelization_util.h"

namespace ug
{

class NewNodeNummerator
{
private:
	std::map<AlgebraID, size_t> m_indicesMap;
	typedef std::map<AlgebraID,size_t>::iterator iterator;
	iterator m_it;
	size_t newIndex;

public:
	NewNodeNummerator(size_t new_indices_start) { newIndex = new_indices_start; }

	size_t get_index_or_create_new(AlgebraID &global_index)
	{
		UG_DLOG(LIB_ALG_MATRIX, 0, "get_index_or_create_new for global index " << global_index.first << " | " << global_index.second << ": ");
		if(global_index.first == pcl::GetProcRank())
		{
			UG_DLOG(LIB_ALG_MATRIX, 0, "is on this processor\n");
			return global_index.second;
		}
		else
		{
			std::pair<iterator, bool> ret = m_indicesMap.insert(pair<AlgebraID, size_t> (global_index, newIndex));
			if(ret.second)
			{
				UG_DLOG(LIB_ALG_MATRIX, 0, "created new index " << newIndex << "\n");
				newIndex++;
			}
			else
				UG_DLOG(LIB_ALG_MATRIX, 0, "already has index " << newIndex << "\n");
			return ret.first->second;
		}
	}

	size_t get_index_if_available(AlgebraID &global_index, bool &has_index)
	{
		UG_DLOG(LIB_ALG_MATRIX, 0, "get_index_if_available for global index " << global_index.first << " | " << global_index.second << ": ");
		if(global_index.first == pcl::GetProcRank())
		{
			UG_DLOG(LIB_ALG_MATRIX, 0, "is on this processor\n");
			has_index = true;
			return global_index.second;
		}
		iterator it = m_indicesMap.find(global_index);
		if(it == m_indicesMap.end())
		{
			UG_DLOG(LIB_ALG_MATRIX, 0, "index not found\n");
			has_index = false;
			return -1;
		}
		else
		{
			UG_DLOG(LIB_ALG_MATRIX, 0, "index is " << it->second << "\n");
			has_index = true;
			return it->second;
		}
	}

	size_t get_new_indices_size()
	{
		return newIndex;
	}

};

template<typename matrix_type>
void GenerateOverlap(const ParallelMatrix<matrix_type> &mat2, ParallelMatrix<matrix_type> &newMat)
{

	const int i=5;

	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (mat2);

	pcl::ParallelCommunicator<IndexLayout> &com = mat.get_communicator();
	UG_DLOG(LIB_ALG_MATRIX, 0, "GENERATE OVERLAP START\n");

	UG_DLOG(LIB_ALG_MATRIX, 0, "matrix is " << mat.num_rows() << " x " << mat.num_cols() << "\n");
	UG_ASSERT(mat.num_rows() == mat.num_cols(), "atm only for square matrices");

	IF_DEBUG(LIB_ALG_MATRIX, 3)
		mat.print();

	// create new matrix
	newMat.create_as_copy_of(mat);

	newMat.set_communicator(mat.get_communicator());
	newMat.set_process_communicator(mat.get_process_communicator());
	newMat.copy_storage_type(mat);

	IndexLayout &masterLayout = mat.get_master_layout();
	IndexLayout &slaveLayout = mat.get_slave_layout();

	// generate global algebra indices
	UG_DLOG(LIB_ALG_MATRIX, 0, "generate " << mat.num_rows() << " global_ids\n");
	std::vector<AlgebraID> global_ids;
	GenerateGlobalAlgebraIDs(global_ids, mat.num_rows(), masterLayout, slaveLayout);

	IF_DEBUG(LIB_ALG_MATRIX, 3)
	{
		for(size_t i = 0; i<global_ids.size(); i++)
			UG_DLOG(LIB_ALG_MATRIX, 0, "local id " << i << " is global id (" << global_ids[i].first << " | " << global_ids[i].second << ")\n");
	}


	StreamPack sendpack, receivepack;

	typedef IndexLayout::Interface Interface;

	// mark slaves
	vector<bool> slave(mat.num_rows(), false);

	for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		Interface &interface = slaveLayout.interface(iter);
		for(Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			slave[interface.get_element(iter2)] = true;
	}

	IndexLayout masterOLLayout;
	IndexLayout slaveOLLayout;

	IF_DEBUG(LIB_ALG_MATRIX, 3)
	{
		UG_DLOG(LIB_ALG_MATRIX, 0, "slave/master list:\n");
		for(size_t i = 0; i<slave.size(); i++)
			UG_DLOG(LIB_ALG_MATRIX, 0, i << " is " << (slave[i] ? "slave\n" : "master\n"));
	}

	// collect data
	std::vector<size_t> N1;
	for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		Interface &interface = slaveLayout.interface(iter);
		size_t pid = slaveLayout.proc_id(iter);

		Interface &interface2 = masterOLLayout.interface(pid);

		UG_DLOG(LIB_ALG_MATRIX, 0, "collecting data for slave interface with processor " << pid << "\n");

		BinaryStream &stream = *sendpack.get_stream(pid);

		vector<bool> node_done(mat.num_rows(), false);

		for(Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			// mark 1-neighbors, which are not slaves.
			// put in stream << globID << numCons << connection1 ... connection[numCons] << ..

			size_t slave_local_index = interface.get_element(iter2);

			GetNeighborhood1(mat, slave_local_index, N1);
			UG_DLOG(LIB_ALG_MATRIX, 0, "node " << slave_local_index << ".\n");

			for(size_t i=0; i<N1.size(); i++)
			{
				size_t local_index = N1[i];
				IF_DEBUG(LIB_ALG_MATRIX, 3)
				{
					UG_DLOG(LIB_ALG_MATRIX, 0, local_index << " (global id " << global_ids[local_index].first << " | " << global_ids[local_index].second << ") ");
					if(slave[local_index]) UG_DLOG(LIB_ALG_MATRIX, 0, "slave.\n");
					if(node_done[local_index]) UG_DLOG(LIB_ALG_MATRIX, 0, "already done.\n");
				}
				if(slave[local_index] || node_done[local_index]) continue;

				node_done[local_index] = true;
				interface2.push_back(local_index);

				AlgebraID &global_index = global_ids[local_index];

				Serialize(stream, global_index);
				Serialize(stream, mat.num_connections(local_index));
				UG_DLOG(LIB_ALG_MATRIX, 0, mat.num_connections(local_index) << " connections: \n");
				for(typename matrix_type::rowIterator conn = mat.beginRow(local_index); !conn.isEnd(); ++conn)
				{
					UG_DLOG(LIB_ALG_MATRIX, 0, "  " << conn.index() << " (global id " << global_ids[conn.index()].first << " | " << global_ids[conn.index()].second << ") -> " << conn.value() << "\n");
					Serialize(stream, global_ids[conn.index()]);
					Serialize(stream, conn.value());
				}
			}
		}

		com.send_raw(pid, stream.buffer(), stream.size(), false);
	}

	// receive data
	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		size_t pid = masterLayout.proc_id(iter);
		UG_DLOG(LIB_ALG_MATRIX, 0, "receiving data for master interface with processor " << pid << "\n");
		BinaryStream &stream = *receivepack.get_stream(pid);
		com.receive_raw(pid, stream);
	}

	UG_DLOG(LIB_ALG_MATRIX, 0, "communicate.\n");
	com.communicate();

	// process data

	NewNodeNummerator nodeNummerator(mat.num_rows());
	size_t num_connections;

	// get all new nodes and their indices
	UG_DLOG(LIB_ALG_MATRIX, 0, "get all new nodes and their indices\n");
	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		size_t pid = masterLayout.proc_id(iter);
		BinaryStream &stream2 = *receivepack.get_stream(pid);
		BinaryStream stream;

		Interface &interface2 = slaveOLLayout.interface(pid);

		stream.write((const char*)stream2.buffer(), stream2.size());

		UG_DLOG(LIB_ALG_MATRIX, 0, "processing data for master interface with processor " << pid << ". size is " << stream.size() << "\n");

		while(stream.can_read_more())
		{
			AlgebraID global_index;
			Deserialize(stream, global_index);
			UG_DLOG(LIB_ALG_MATRIX, 0, " global id " << global_index.first << " | " << global_index.second << ". ");

			// create local index for this global index
			size_t local_index = nodeNummerator.get_index_or_create_new(global_index);
			interface2.push_back(local_index);
			UG_DLOG(LIB_ALG_MATRIX, 0, " local index is " << local_index << ".\n");

			// skip connection information (matrix row) in stream
			Deserialize(stream, num_connections);
			stream.read_jump((sizeof(AlgebraID) + sizeof(typename matrix_type::value_type)) * num_connections);
		}


	}

	// resize matrix to be now new_indices_size x new_indices_size.
	size_t new_indices_size = nodeNummerator.get_new_indices_size();
	UG_DLOG(LIB_ALG_MATRIX, 0, "resize matrix to be now " << new_indices_size << " x " << new_indices_size << ".\n");
	newMat.resize(new_indices_size, new_indices_size);

	// iterate again over all streams to get the matrix lines
	UG_DLOG(LIB_ALG_MATRIX, 0, "iterate again over all streams to get the matrix lines\n");
	vector<typename matrix_type::connection> cons;
	bool has_index;
	AlgebraID global_row_index, global_col_index;
	size_t j;
	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		size_t pid = masterLayout.proc_id(iter);
		BinaryStream &stream = *receivepack.get_stream(pid);

		UG_DLOG(LIB_ALG_MATRIX, 0, "processing data for master interface with processor " << pid << ". size is " << stream.size() << "\n");

		while(stream.can_read_more())
		{
			Deserialize(stream, global_row_index);

			// get local index for this global index
			size_t local_row_index = nodeNummerator.get_index_if_available(global_row_index, has_index);
			UG_ASSERT(has_index, "global id " << global_row_index.first << " | " << global_row_index.second << " has no associated local id");

			UG_DLOG(LIB_ALG_MATRIX, 0, "Processing global id " << global_row_index.first << " | " << global_row_index.second << ", local id " << local_row_index << ". ");

			// get nr of connections
			Deserialize(stream, num_connections);
			if(cons.size() < num_connections) cons.resize(num_connections);
			UG_DLOG(LIB_ALG_MATRIX, 0, num_connections << " connections:\n");

			j=0;
			for(size_t i=0; i<num_connections; i++)
			{
				// get global column index, and associated local index
				Deserialize(stream, global_col_index);
				cons[j].iIndex = nodeNummerator.get_index_if_available(global_col_index, has_index);
				// if we have an associated local index, deserialize value and save connection
				if(has_index)
				{
					UG_DLOG(LIB_ALG_MATRIX, 0, cons[j].iIndex << " (global index " << global_col_index.first << " | " << global_col_index.second << ")");
					Deserialize(stream, cons[j++].dValue);
					UG_DLOG(LIB_ALG_MATRIX, 0, " -> " << cons[j-1].dValue << "\n");

				}
				// otherwise: skip this connection
				else
				{
					UG_DLOG(LIB_ALG_MATRIX, 0, " connection to global index " << global_col_index.first << " | " << global_col_index.second << ", but is not on this processor\n");
					stream.read_jump(sizeof(typename matrix_type::value_type));
				}
			}

			// set matrix row
			UG_DLOG(LIB_ALG_MATRIX, 0, "set matrix row: ");
			for(size_t i=0; i<j; i++)
				UG_DLOG(LIB_ALG_MATRIX, 0, "(" << cons[i].iIndex << "-> " << cons[i].dValue << ") ");
			UG_DLOG(LIB_ALG_MATRIX, 0, "\n");
			newMat.set_matrix_row(local_row_index, &cons[0], j);
		}
	}

	UG_DLOG(LIB_ALG_MATRIX, 0, "new matrix\n\n");
	IF_DEBUG(LIB_ALG_MATRIX, 3)
	{
		newMat.print();

		UG_DLOG(LIB_ALG_MATRIX, 0, "master:\n");
		for(IndexLayout::iterator iter = masterOLLayout.begin(); iter != masterOLLayout.end(); ++iter)
		{
			size_t pid = masterOLLayout.proc_id(iter);
			UG_DLOG(LIB_ALG_MATRIX, 0, "to processor " << pid << ": ");
			Interface &interface = masterOLLayout.interface(iter);
			for(Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				size_t i = interface.get_element(iter2);
				UG_DLOG(LIB_ALG_MATRIX, 0, i << " ");
			}
			UG_DLOG(LIB_ALG_MATRIX, 0, "\n");
		}

		UG_DLOG(LIB_ALG_MATRIX, 0, "slave:\n");
		for(IndexLayout::iterator iter = slaveOLLayout.begin(); iter != slaveOLLayout.end(); ++iter)
		{
			size_t pid = slaveOLLayout.proc_id(iter);
			UG_DLOG(LIB_ALG_MATRIX, 0, "to processor " << pid << ": ");
			Interface &interface = slaveOLLayout.interface(iter);
			for(Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				size_t i = interface.get_element(iter2);
				UG_DLOG(LIB_ALG_MATRIX, 0, i << " ");
			}
			UG_DLOG(LIB_ALG_MATRIX, 0, "\n");
		}
	}
	newMat.set_slave_layout(slaveOLLayout);
	newMat.set_slave_layout(masterOLLayout);
}


}
#endif
