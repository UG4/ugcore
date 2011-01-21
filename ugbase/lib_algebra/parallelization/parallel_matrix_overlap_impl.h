/*
 * parallelization.h
 *
 *  Created on: 17.01.2010
 *      Author: Martin Rupp / Sebastian Reiter
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_OVERLAP_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_OVERLAP_IMPL__

#include "new_node_nummerator.h"
#include "parallelization_util.h"

namespace ug
{

// ACHTUNG!!! Overlap 2 kann sich Ÿber 2 Prozessoren erstrecken!
//       p1     p2
//      o   o   o  | (o)  o   o
//      o   o   o  | (o)  o   o
//      o   o   A  | (A)  B   o
//     --------------------------
//     (o) (o) (A)       (B) (o)
// p3   o   o   o         C   o
//      o   o   o         C   o
//      o   o   o         C   o
//
// slave-knoten sind in (klammern).
// soll nun der Ÿberlap p1-p2 bestimmt werden, so kann p2 vom A slave nicht mehr die 2-Umgebung bestimmen, denn dazu gehšrt schon C:
//      o   o   o  (o) (o)
//      o   o   o  (o) (o)
//      o   o   A  (B) (o)
//     (o) (o) (o) (C)
//     (o) (o) (o)

// (5-punkt stern)
// es muss also iterativ vorgegangen werden

template<typename matrix_type>
void GenerateOverlap(const ParallelMatrix<matrix_type> &mat2, ParallelMatrix<matrix_type> &newMat, IndexLayout &masterOLLayout, IndexLayout &slaveOLLayout, size_t overlap_depth=1)
{
	UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (mat2);

	pcl::ParallelCommunicator<IndexLayout> &com = mat.get_communicator();
	UG_DLOG(LIB_ALG_MATRIX, 4, "GENERATE OVERLAP START\n");

	UG_DLOG(LIB_ALG_MATRIX, 4, "matrix is " << mat.num_rows() << " x " << mat.num_cols() << "\n");
	UG_ASSERT(mat.num_rows() == mat.num_cols(), "atm only for square matrices");

	IF_DEBUG(LIB_ALG_MATRIX, 4)
		mat.print();

	// create new matrix
	newMat.create_as_copy_of(mat);

	newMat.set_communicator(mat.get_communicator());
	newMat.set_process_communicator(mat.get_process_communicator());
	newMat.copy_storage_type(mat);

	IndexLayout &masterLayout = mat.get_master_layout();
	IndexLayout &slaveLayout = mat.get_slave_layout();

	// generate global algebra indices
	UG_DLOG(LIB_ALG_MATRIX, 4, "generate " << mat.num_rows() << " global_ids\n");
	std::vector<AlgebraID> global_ids;
	GenerateGlobalAlgebraIDs(global_ids, mat.num_rows(), masterLayout, slaveLayout);

	IF_DEBUG(LIB_ALG_MATRIX, 4)
	{
		for(size_t i = 0; i<global_ids.size(); i++)
			UG_DLOG(LIB_ALG_MATRIX, 4, "local id " << i << " is global id (" << global_ids[i].first << " | " << global_ids[i].second << ")\n");
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

	AddLayout(masterOLLayout, masterLayout);
	AddLayout(slaveOLLayout, slaveLayout);
	PrintLayout(slaveOLLayout);
	PrintLayout(masterOLLayout);


	IF_DEBUG(LIB_ALG_MATRIX, 4)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "slave/master list:\n");
		for(size_t i = 0; i<slave.size(); i++)
			UG_DLOG(LIB_ALG_MATRIX, 4, i << " is " << (slave[i] ? "slave\n" : "master\n"));
	}


	// collect data
	std::vector<size_t> N1;
	std::vector<bool> bVisited(mat.num_rows(), false);
	for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		Interface &interface = slaveLayout.interface(iter);
		size_t pid = slaveLayout.proc_id(iter);

		Interface &interface2 = masterOLLayout.interface(pid);

		UG_DLOG(LIB_ALG_MATRIX, 4, "collecting data for slave interface with processor " << pid << "\n");

		BinaryStream &stream = *sendpack.get_stream(pid);

		vector<bool> node_done(mat.num_rows(), false);

		for(Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			// mark 1-neighbors, which are not slaves.
			// put in stream << globID << numCons << connection1 ... connection[numCons] << ..

			size_t slave_local_index = interface.get_element(iter2);

			GetNeighborhood(mat, slave_local_index, overlap_depth, N1);
			UG_DLOG(LIB_ALG_MATRIX, 4, "node " << slave_local_index << ".\n");
			UG_DLOG(LIB_ALG_MATRIX, 4, "neighbors: ")
			for(size_t i=0; i<N1.size(); i++)
				UG_DLOG(LIB_ALG_MATRIX, 4, N1[i] << " ");
			UG_DLOG(LIB_ALG_MATRIX, 4, "\n");

			for(size_t i=0; i<N1.size(); i++)
			{
				size_t local_index = N1[i];
				UG_DLOG(LIB_ALG_MATRIX, 4, "* " << overlap_depth << " neighbor " << i << " is " << local_index << ":\n");
				IF_DEBUG(LIB_ALG_MATRIX, 4)
				{
					UG_DLOG(LIB_ALG_MATRIX, 4, local_index << " (global id " << global_ids[local_index].first << " | " << global_ids[local_index].second << ") ");
					if(slave[local_index]) UG_DLOG(LIB_ALG_MATRIX, 4, "slave.\n");
					if(node_done[local_index]) UG_DLOG(LIB_ALG_MATRIX, 4, "already done.\n");
				}
				if(slave[local_index] || node_done[local_index]) continue;

				// hier muss man nochmal ŸberprŸfen, welche verbindungen wir schicken mŸssen.
				// hier stand mal if(slave[local_index]) continue;, aber dann hŠtte zB. der master nicht seine volle zeile erhalten
				// allerdings mŸssen wir natŸrlich verhindern, dass wir verbindungen schicken, die schon da sind, zumindest bei der rŸckrichtung mŸssen wir aufpassen

				node_done[local_index] = true;
				interface2.push_back(local_index);

				AlgebraID &global_index = global_ids[local_index];

				Serialize(stream, global_index);
				Serialize(stream, mat.num_connections(local_index));
				UG_DLOG(LIB_ALG_MATRIX, 4, mat.num_connections(local_index) << " connections: \n");
				for(typename matrix_type::rowIterator conn = mat.beginRow(local_index); !conn.isEnd(); ++conn)
				{
					UG_DLOG(LIB_ALG_MATRIX, 4, "  " << conn.index() << " (global id " << global_ids[conn.index()].first << " | " << global_ids[conn.index()].second << ") -> " << conn.value() << "\n");
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
		UG_DLOG(LIB_ALG_MATRIX, 4, "receiving data for master interface with processor " << pid << "\n");
		BinaryStream &stream = *receivepack.get_stream(pid);
		com.receive_raw(pid, stream);
	}

	UG_DLOG(LIB_ALG_MATRIX, 4, "communicate.\n");
	com.communicate();

	// process data

	NewNodeNummerator nodeNummerator(mat.num_rows());
	size_t num_connections;
	AlgebraID global_row_index, global_col_index;
	typename matrix_type::value_type value;

	// get all new nodes and their indices
	UG_DLOG(LIB_ALG_MATRIX, 4, "get all new nodes and their indices\n");
	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		size_t pid = masterLayout.proc_id(iter);
		BinaryStream &stream2 = *receivepack.get_stream(pid);
		BinaryStream stream;

		Interface &interface2 = slaveOLLayout.interface(pid);

		stream.write((const char*)stream2.buffer(), stream2.size());

		UG_DLOG(LIB_ALG_MATRIX, 4, "processing data for master interface with processor " << pid << ". size is " << stream.size() << "\n");

		while(stream.can_read_more())
		{
			Deserialize(stream, global_row_index);
			UG_DLOG(LIB_ALG_MATRIX, 4, " global id " << global_row_index.first << " | " << global_row_index.second << ". ");

			// create local index for this global index
			size_t local_index = nodeNummerator.get_index_or_create_new(global_row_index);
			interface2.push_back(local_index);
			UG_DLOG(LIB_ALG_MATRIX, 4, " local index is " << local_index << ".\n");

			// skip connection information (matrix row) in stream
			Deserialize(stream, num_connections);

			// this is only possible with static values:
			//stream.read_jump((sizeof(AlgebraID) + sizeof(typename matrix_type::value_type)) * num_connections);

			for(size_t i=0; i<num_connections; i++)
			{
				// get global column index, and associated local index
				Deserialize(stream, global_col_index);
				Deserialize(stream, value);
			}
		}


	}

	// resize matrix to be now new_indices_size x new_indices_size.
	size_t new_indices_size = nodeNummerator.get_new_indices_size();
	UG_DLOG(LIB_ALG_MATRIX, 4, "resize matrix to be now " << new_indices_size << " x " << new_indices_size << ".\n");
	newMat.resize(new_indices_size, new_indices_size);

	// iterate again over all streams to get the matrix lines
	UG_DLOG(LIB_ALG_MATRIX, 4, "iterate again over all streams to get the matrix lines\n");
	vector<typename matrix_type::connection> cons;
	bool has_index;
	size_t j;
	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		size_t pid = masterLayout.proc_id(iter);
		BinaryStream &stream = *receivepack.get_stream(pid);

		UG_DLOG(LIB_ALG_MATRIX, 4, "processing data for master interface with processor " << pid << ". size is " << stream.size() << "\n");

		while(stream.can_read_more())
		{
			Deserialize(stream, global_row_index);

			// get local index for this global index
			size_t local_row_index = nodeNummerator.get_index_if_available(global_row_index, has_index);
			UG_ASSERT(has_index, "global id " << global_row_index.first << " | " << global_row_index.second << " has no associated local id");

			UG_DLOG(LIB_ALG_MATRIX, 4, "Processing global id " << global_row_index.first << " | " << global_row_index.second << ", local id " << local_row_index << ". ");

			// get nr of connections
			Deserialize(stream, num_connections);
			if(cons.size() < num_connections) cons.resize(num_connections);
			UG_DLOG(LIB_ALG_MATRIX, 4, num_connections << " connections:\n");

			j=0;
			for(size_t i=0; i<num_connections; i++)
			{
				// get global column index, and associated local index
				Deserialize(stream, global_col_index);
				Deserialize(stream, cons[j].dValue);
				cons[j].iIndex = nodeNummerator.get_index_if_available(global_col_index, has_index);

				UG_DLOG(LIB_ALG_MATRIX, 4, cons[j].iIndex << " (global index " << global_col_index.first << " | " << global_col_index.second << ")");
				UG_DLOG(LIB_ALG_MATRIX, 4, " -> " << cons[j].dValue << "\n");
				// if we have an associated local index, save connection
				if(has_index)
					j++;
				else
					UG_DLOG(LIB_ALG_MATRIX, 4, " connection to global index " << global_col_index.first << " | " << global_col_index.second << ", but is not on this processor\n");
			}

			// set matrix row
			UG_DLOG(LIB_ALG_MATRIX, 4, "set matrix row: ");
			for(size_t i=0; i<j; i++)
				UG_DLOG(LIB_ALG_MATRIX, 4, "(" << cons[i].iIndex << "-> " << cons[i].dValue << ") ");
			UG_DLOG(LIB_ALG_MATRIX, 4, "\n");
			newMat.add_matrix_row(local_row_index, &cons[0], j);
		}
	}

	UG_DLOG(LIB_ALG_MATRIX, 4, "new matrix\n\n");
	IF_DEBUG(LIB_ALG_MATRIX, 4)
	{
		newMat.print();

		UG_DLOG(LIB_ALG_MATRIX, 4, "master Layout:\n");
		PrintLayout(masterOLLayout);

		UG_DLOG(LIB_ALG_MATRIX, 4, "slave:\n");
		PrintLayout(slaveOLLayout);
	}
	newMat.set_slave_layout(slaveOLLayout);
	newMat.set_master_layout(masterOLLayout);
}


}
#endif
