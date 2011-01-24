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
//      o   o   o  (o) (o)
//      o   o   o  (o) (o)
//      o   o   A  (B) (o)
//     (o) (o) (o) (C)
//     (o) (o) (o)

// (5-punkt stern)
// frage: kann sich ein Ÿberlapp Ÿber mehrere Prozessoren erstrecken?
// o A | (A) B | (B) C


// GenerateOverlap_CollectData
//--------------------------------------------------
/**
 * \brief Subfunction for GenerateOverlap. Collects matrix rows and sends them via a ParallelCommunicator
 * \param com 				used ParallelCommunicator
 * \param mat				Matrix of which to send matrix rows
 * \param layout			we collect data for all interfaces of this layout
 * \param overlapLayout		where overlap nodes are saved
 * \param sendpack			used StreamPack for putting data stream to
 * \param overlap_depth		depth of overlap
 * \param global_ids		global algebra ids
 *
 * For all Interfaces I in layout:
 *   create interface I2 in overlapLayout with same destination pid as interface I
 *   For all nodes i, which have a nodal distance of at most overlap_depth from a node in the interface I
 *     put the complete matrix row via sendpack.get_stream(pid)
 *     is this node on this processor, its local index is added to interface I2 (these become overlap nodes)
 *    send data via com.
 *
 * format for sending matrix rows:
 * size_t global_row_id;
 * size_t num_connections;
 * {
 *  size_t global_col_id;
 *  matrix_type::value_type value;
 * } con[num_connections];
 *
 * always use Serialize/Deserialize for all this values, since matrix_type::value_type may not be of fixed size.
 *
 * \sa GenerateOverlap
 */
template<typename matrix_type>
void GenerateOverlap_CollectData(pcl::ParallelCommunicator<IndexLayout> &com, matrix_type &mat, IndexLayout &layout,
		IndexLayout &overlapLayout, StreamPack &sendpack, size_t overlap_depth, std::vector<AlgebraID> &global_ids)
{
	typedef IndexLayout::Interface Interface;

	std::vector<size_t> N1;
	std::vector<bool> bVisited(mat.num_rows(), false);

	// For all Interfaces in layout:
	for(IndexLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
	{
		Interface &interface = layout.interface(iter);
		size_t pid = layout.proc_id(iter);

		// create interface in overlapLayout with same destination pid
		Interface &overlapInterface = overlapLayout.interface(pid);

		UG_DLOG(LIB_ALG_MATRIX, 4, "collecting data for slave interface with processor " << pid << "\n");

		// get a stream to put data to process pid to
		BinaryStream &stream = *sendpack.get_stream(pid);

		vector<bool> node_done(mat.num_rows(), false);

		for(Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
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
					//if(slave[local_index]) UG_DLOG(LIB_ALG_MATRIX, 4, "slave.\n");
					if(node_done[local_index]) UG_DLOG(LIB_ALG_MATRIX, 4, "already done.\n");
				}
				if(node_done[local_index]) continue;
				node_done[local_index] = true;

				// if the node is on this processor, its local index is added to overlapInterface
				AlgebraID &global_index = global_ids[local_index];
				if(global_index.first == pcl::GetProcRank())
					overlapInterface.push_back(local_index);


				// serialize matrix row
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

		// use communicator to send data
		com.send_raw(pid, stream.buffer(), stream.size(), false);
	}
}

// GenerateOverlap_GetNewIndices
//--------------------------------------------------
/**
 * \brief Subfunction for GenerateOverlap. Extracts new indices from received data
 * \param receivepack		StreamPack where to extract new indices of
 * \param overlapLayout		node which are master on another processor are added here
 * \param nodeNummerator	handles creation of new indices on this processor
 *
 * when this function is finished, all global row indices send through receivepack / GenerateOverlap_CollectData
 * to this processor have unique indices available through nodeNummerator, and all corresponding interfaces for
 * the overlap nodes are set.
 * \sa GenerateOverlap
 * \sa GenerateOverlap_CollectData
 */
template<typename matrix_type>
void GenerateOverlap_GetNewIndices(StreamPack &receivepack, IndexLayout &overlapLayout, NewNodeNummerator &nodeNummerator)
{
	typedef IndexLayout::Interface Interface;
	UG_DLOG(LIB_ALG_MATRIX, 4, "get all new nodes and their indices\n");

	size_t num_connections;
	AlgebraID global_row_index, global_col_index;
	typename matrix_type::value_type value;

	for(StreamPack::iterator iter = receivepack.begin(); iter != receivepack.end(); ++iter)
	{
		size_t pid = iter->first;
		// copy stream (bug in BinaryStream, otherwise I would use stream.seekg(ios::begin)
		BinaryStream &stream2 = *iter->second;
		BinaryStream stream; stream.write((const char*)stream2.buffer(), stream2.size());

		Interface &interface2 = overlapLayout.interface(pid);

		while(stream.can_read_more())
		{
			Deserialize(stream, global_row_index);

			// create local index for this global index
			if(global_row_index.first != pcl::GetProcRank())
			{
				size_t local_index = nodeNummerator.get_index_or_create_new(global_row_index);
				interface2.push_back(local_index);
				UG_DLOG(LIB_ALG_MATRIX, 4, " global id " << global_row_index.first << " | " << global_row_index.second << ". local index is " << local_index << ".\n");
			}
			else
			{
				UG_DLOG(LIB_ALG_MATRIX, 4, " global id " << global_row_index.first << " | " << global_row_index.second << ". local index is " << global_row_index.second << " (on this processor).\n");
			}

			// skip connection information (matrix row) in stream
			Deserialize(stream, num_connections);

			// this is only possible with static values:
			//stream.read_jump((sizeof(AlgebraID) + sizeof(typename matrix_type::value_type)) * num_connections);

			for(size_t i=0; i<num_connections; i++)
			{
				Deserialize(stream, global_col_index);
				Deserialize(stream, value);
			}
		}

	}

}

// GenerateOverlap_AddMatrixRows
//--------------------------------------------------
/**
 * \brief Subfunction for GenerateOverlap. Extracts matrix rows from received data, and adds them a matrix
 * \param receivepack		StreamPack where to extract data from
 * \param newMat			matrix to add new rows to
 * \param nodeNummerator	maps global indices to local indices
 *
 *
 * \sa GenerateOverlap_CollectData
 * \sa GenerateOverlap
 * \sa GenerateOverlap_GetNewIndices
 */
template<typename matrix_type>
void GenerateOverlap_AddMatrixRows(StreamPack &receivepack, matrix_type &newMat, NewNodeNummerator &nodeNummerator)
{
	UG_DLOG(LIB_ALG_MATRIX, 4, "iterate again over all streams to get the matrix lines\n");

	size_t num_connections;
	AlgebraID global_row_index, global_col_index;
	vector<typename matrix_type::connection> cons;
	bool has_index;
	size_t j;
	for(StreamPack::iterator iter = receivepack.begin(); iter != receivepack.end(); ++iter)
	{
		size_t pid = iter->first;
		BinaryStream &stream = *iter->second;

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
}


// GenerateOverlap
//--------------------------------------------------
/**
 * \brief Generates a new matrix with overlap from another matrix
 * \param _mat				matrix to create overlap from
 * \param newMat			matrix to store overlap matrix in
 * \param masterOLLayout
 * \param masterOLLayout
 * \param overlap_depth
 *
 * \sa GenerateOverlap_CollectData
 * \sa GenerateOverlap_GetNewIndices
 * \sa GenerateOverlap_AddMatrixRows
 */
template<typename matrix_type>
void GenerateOverlap(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat, IndexLayout &masterOLLayout, IndexLayout &slaveOLLayout, size_t overlap_depth=1)
{
	typedef IndexLayout::Interface Interface;

	// pcl does not use const much
	UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	pcl::ParallelCommunicator<IndexLayout> &com = mat.get_communicator();
	UG_DLOG(LIB_ALG_MATRIX, 4, "GENERATE OVERLAP START\n");

	UG_DLOG(LIB_ALG_MATRIX, 4, "matrix is " << mat.num_rows() << " x " << mat.num_cols() << "\n");
	UG_ASSERT(mat.num_rows() == mat.num_cols(), "atm only for square matrices");

	IF_DEBUG(LIB_ALG_MATRIX, 4)
		mat.print();


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
	//StreamPack sendpack_master, receivepack_master;

	newMat.create_as_copy_of(mat);

	// collect data
	//-----------------

	GenerateOverlap_CollectData(com, mat, slaveLayout, masterOLLayout, sendpack, overlap_depth, global_ids);
	//CollectData(com, mat, masterLayout, masterOLLayout, sendpack_master, overlap_depth-1, global_ids);

	// receive data
	//-----------------
	for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		size_t pid = masterLayout.proc_id(iter);
		com.receive_raw(pid, *receivepack.get_stream(pid));
	}

	com.communicate();

	// process data
	//-----------------
	NewNodeNummerator nodeNummerator(mat.num_rows());

	// get all new nodes and their indices
	GenerateOverlap_GetNewIndices<matrix_type>(receivepack, slaveOLLayout, nodeNummerator);

	// create as copy & resize matrix to be now new_indices_size x new_indices_size.
	size_t new_indices_size = nodeNummerator.get_new_indices_size();
	UG_DLOG(LIB_ALG_MATRIX, 4, "resize matrix to be now " << new_indices_size << " x " << new_indices_size << ".\n");

	newMat.resize(new_indices_size, new_indices_size);

	// iterate again over all streams to get the matrix lines

	GenerateOverlap_AddMatrixRows(newMat, receivepack, nodeNummerator);

	// done!

	AddLayout(masterOLLayout, masterLayout);
	AddLayout(slaveOLLayout, slaveLayout);

	newMat.set_slave_layout(slaveOLLayout);
	newMat.set_master_layout(masterOLLayout);
	newMat.set_communicator(mat.get_communicator());
	newMat.set_process_communicator(mat.get_process_communicator());
	newMat.copy_storage_type(mat);


	UG_DLOG(LIB_ALG_MATRIX, 4, "new matrix\n\n");
	IF_DEBUG(LIB_ALG_MATRIX, 4)
	{
		newMat.print();

		UG_DLOG(LIB_ALG_MATRIX, 4, "master Layout:\n");
		PrintLayout(masterOLLayout);

		UG_DLOG(LIB_ALG_MATRIX, 4, "slave:\n");
		PrintLayout(slaveOLLayout);
	}


}


}
#endif
