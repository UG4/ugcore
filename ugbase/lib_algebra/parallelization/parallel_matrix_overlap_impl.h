/*
 * parallelization.h
 *
 *  Created on: 17.01.2010
 *      Author: Martin Rupp / Sebastian Reiter
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_OVERLAP_IMPL__
#define __H__LIB_ALGEBRA__PARALLELIZATION__PARALLEL_MATRIX_OVERLAP_IMPL__

#include <vector>
#include <set>
#include <map>

#include "new_nodes_nummerator.h"
#include "parallelization_util.h"
//#include "test_layout.h"
#include "common/util/sort_util.h"


namespace ug
{


template<typename TLayout>
void AddLayout(TLayout &destLayout, TLayout &sourceLayout)
{
	for(typename TLayout::iterator iter = sourceLayout.begin(); iter != sourceLayout.end(); ++iter)
	{
		typename TLayout::Interface &source_interface = sourceLayout.interface(iter);
		typename TLayout::Interface &dest_interface = destLayout.interface(sourceLayout.proc_id(iter));
		for(typename TLayout::Interface::iterator iter2 = source_interface.begin(); iter2 != source_interface.end(); ++iter2)
			dest_interface.push_back(source_interface.get_element(iter2));
	}
}

template<typename TLayout>
void PrintLayout(TLayout &layout)
{
	for(typename TLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
	{
		size_t pid = layout.proc_id(iter);
		UG_LOG("to processor " << pid << ": ");
		typename TLayout::Interface &interface = layout.interface(iter);
		for(typename TLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			UG_LOG(interface.get_element(iter2) << "  ");
		UG_LOG("\n");
	}
}



// neuer algorithmus:
// overlap 0
// 1. slave-knoten verschicken ihre Zeile an den Master.
// 2. Master nimmt diese entgegen, zeile wird fŸr vorhandene Knoten addiert
// 3. fertig.

// overlap 1
// 1. slave-knoten verschicken ihre Zeile an den Master.
//    /!\ werden verknŸpfungen zu anderen prozessoren verschickt, werden die prozessoren informiert
//    /!\ unter umstŠnden wird so ein prozessor "von 2 seiten" informiert. dann muss es eine
// 2. verschicke die matrixzeilen und benachrichtungen
// 3. nehme matrixzeilen und benachrichtungen entgegen
// 4. verarbeite benachrichtungen: erzeuge u.U. neue Master Knoten
// 5. verarbeite matrixzeilen und erzeugt u.U. neue Knoten (Slave). Diese werden als "neu" markiert
//
// \param mat matrix to take values from
// \param layout : for all interfaces I in layout, send matrix rows of nodes in I with com
// \param create_new_nodes : if true, do not add nodes to masterOLLayout
// \param masterOLLayout : for nodes which will be new on the other processors, add them to this layout
// \param pids ... unclear if this is needed



template<typename matrix_type>
class GenerateOverlapClass
{
private:
	typedef IndexLayout::Interface Interface;

	pcl::ParallelCommunicator<IndexLayout> &m_com; ///< communicator

	std::vector<AlgebraID> m_globalIDs; ///< global IDs

	/// map for marking nodes
	/**
	 * for each processor we need to have a list which of our master nodes exist on their processor
	 * this is important because we sometimes will need to add them to interfaces
	 *
	 * the map serves two functions:
	 * - knowing which processor has copies of our own master nodes, and so constructing correct interfaces
	 * - knowing which notifications are sent
	 *
	 * - if i is a master node (that is m_globalIDs[i].first == pcl::GetProcRank())
	 *   then mark[pid][i] = true means: process pid already knows that i is a slave node on his processor,
	 *   and this processor knows that i the associated master
	 *   that means: i is in a master interface on this processor and in a slave interface on pid.
	 * - if is is not a master node
	 *   then mark[pid][i] = true means: we already sent a notification to the owner of i (processor m_globalIDs[i].first)
	 *   that processor pid has a copy of i.
	 *

	 */
	std::map<int, std::vector<bool> > m_mapMark;


	matrix_type &m_mat;		///< the original matrix (should be const)
	matrix_type &m_newMat;	///< the new to create matrix

	IndexLayout &m_totalMasterLayout;	///< layout combining all master layouts from overlap 0 to overlap_depth-1
	IndexLayout &m_totalSlaveLayout;	///< layout combining all slave layouts from overlap 0 to overlap_depth-1

	std::vector<IndexLayout> &m_vMasterLayouts; 	///< m_vMasterLayout[i] is the master layout from overlap i without others
	std::vector<IndexLayout> &m_vSlaveLayouts;	///< m_vSlaveLayout[i] is the slave layout from overlap i without others
	size_t m_overlapDepth;						///< overlap depth to be achieved




	/// collect_matrix_row_data
	/**
	 * \param mat
	 * \param layout
	 * \param bAddNodesToLayout if true, create new nodes in interfaces/send notifications
	 * \param pids if bAddNodesToLayout=true, set of processor ids we need to communicate notifications to
	 * \param mNewLayout if bAddNodesToLayout=true, the function adds nodes which are not already on processor pid to the mNewLayout[pid]
	 *
	 * for all interfaces I in layout
	 * 		pid = associated processor of I
	 * 		for all indices i in I
	 * 			prepare to send matrix row mat(i, .) to processor pid
	 *
	 * 			if index j of this row is not on this processor
	 * 				prepare to send a notification to owner of i that i has a copy on pid if not already done so
	 * 			else
	 * 				add i to newInterface to pid if not already in some interface to pid.
	 * 		send matrix rows
	 *
	 * for all pids
	 * 		send notifications
	 *
	 *\sa m_mapMark
	 */
	void collect_matrix_row_data(const matrix_type &mat,
			IndexLayout &layout, bool bAddNodesToLayout, std::set<int> &pids,
			std::map<int, std::vector<size_t> > &mNewLayout)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlapClass::collect_matrix_row_data\n");

		StreamPack notifications_pack;

		// zeilen verschicken
		for(IndexLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
		{
			Interface &interface = layout.interface(iter);
			int pid = layout.proc_id(iter);

			// get a stream to put data to process pid to
			BinaryStream stream;
			std::vector<bool> &mark = m_mapMark[pid];

			UG_DLOG(LIB_ALG_MATRIX, 4, "PID " << pid << "\n");

			if(mark.size() == 0) mark.resize(m_mat.num_rows(), false);


			for(Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				size_t local_index = interface.get_element(iter2);

				AlgebraID &globalRowIndex = m_globalIDs[local_index];

				std::vector<size_t> *p_vNewInterface = NULL;

				// serialize global row index
				Serialize(stream, globalRowIndex);

				// get number of connections != 0.0
				size_t num_connections = 0;
				for(typename matrix_type::cRowIterator conn = mat.beginRow(local_index); !conn.isEnd(); ++conn)
				{
					if(conn.value() == 0.0) continue;
					num_connections++;
				}

				// serialize number of connections
				Serialize(stream, num_connections);

				UG_DLOG(LIB_ALG_MATRIX, 4, "GID " << globalRowIndex.first << " | " << globalRowIndex.second << "\n");

				// serialize connections
				for(typename matrix_type::cRowIterator conn = mat.beginRow(local_index); !conn.isEnd(); ++conn)
				{
					if(conn.value() == 0.0) continue;
					AlgebraID &globalColIndex = m_globalIDs[conn.index()];
					size_t localColIndex = conn.index();

					UG_DLOG(LIB_ALG_MATRIX, 4, "  " << conn.index() << " (global id " <<
							globalColIndex.first << " | " << globalColIndex.second << ") -> " << conn.value() << "\n");

					// serialize connection
					Serialize(stream, globalColIndex);
					Serialize(stream, conn.value());

					// check if nodes can be added
					if(bAddNodesToLayout == false)
					{
						UG_DLOG(LIB_ALG_MATRIX, 4, "not adding node to masterOLInterface because pmNewLayout = NULL\n");
					}
					else
					{
						if(globalColIndex.first == pcl::GetProcRank())
						{
							// add i to mNewLayout[pid] if not already in some interface to pid.

							UG_ASSERT(globalColIndex.second == localColIndex, "");

							// only add nodes once to interface
							// for master nodes: mark[pid][i] = true means: we already send some matrix row with this index, so there is already a interface.
							if(mark[localColIndex] == false)
							{
								mark[localColIndex] = true;
								UG_DLOG(LIB_ALG_MATRIX, 4, " adding node to masterOLInterface .interface(" << pid << ").push_back(" << localColIndex << ")\n");
								if(p_vNewInterface == NULL)
									p_vNewInterface = &mNewLayout[pid];
								p_vNewInterface->push_back(localColIndex);
							}
							else
								UG_DLOG(LIB_ALG_MATRIX, 4, "  not adding node to masterOLInterface because it is already\n");
							//node_marked[localColIndex] = true;
						}
						else
						{
							// if localColIndex is not on this processor
							//	prepare to send a notification to owner of (globalColIndex.first) that globalColIndex has a copy on pid if not already done so

							UG_DLOG(LIB_ALG_MATRIX, 4, "  not adding node to masterOLInterface because it is on another processor\n");

							if(globalColIndex.first != pid) // obviously the owner knows that his node is on his processor
							{
								UG_DLOG(LIB_ALG_MATRIX, 4, " processor " << globalColIndex.first << " gets informed that his local node "
																		<< globalColIndex.second << " has a copy on processor " << pid << "\n");

								// for slave nodes: mark[pid][i]=true means we already send notification to the owner that local index i has a copy on processor pid.
								std::vector<bool> &otherMark = m_mapMark[pid];
								otherMark.resize(m_globalIDs.size(), false);

								if(otherMark[localColIndex] == false)
								{
									otherMark[localColIndex] = true;
									BinaryStream &notifications_stream = *notifications_pack.get_stream(globalColIndex.first);

									// serialize notification: global_index.second, pid.
									Serialize(notifications_stream, globalColIndex.second);
									Serialize(notifications_stream, pid);
									//UG_ASSERT(find(pids.begin(), pids.end(), globalColIndex.first) != pids.end(), "sending notification to a processor "
									// << globalColIndex.first << ", not in pid list???");
								}
								else
									UG_DLOG(LIB_ALG_MATRIX, 4, " not sending notification because already send one\n");
							}
							else
								UG_DLOG(LIB_ALG_MATRIX, 4, "  not sending notification because owner knows that his node is on his\n");

						}
					}




				}
			}

			// use communicator to send data
			UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": sending matrixrow data to proc " << pid << " (size " << stream.size() << ")\n");
			m_com.send_raw(pid, stream.buffer(), stream.size(), false);
		}

		if(bAddNodesToLayout)
		{
			for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
			{
				int pid = *iter;
				BinaryStream &stream = *notifications_pack.get_stream(pid);
				m_com.send_raw(pid, stream.buffer(), stream.size(), false);
				UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": sending notifications data to proc " << pid << " (size " << stream.size() << ") \n");
			}
		}
	}

	/// receive_notifications
	/**
	 * processing notifications about our master nodes having slaves on other processors
	 * \param notifications_pack notifications stream in form (size_t) local_index, (int) pid.
	 * \param masterOLLayout layout where to add masters if appropriate
	 *
	 * for a notification that the index local_index has now a copy on processor pid, we add local_index to masterOLLayout[pid], if
	 * local_index is not already in some master interface to pid.
	 *
	 * \sa m_mapMark
	 */
	void receive_notifications(StreamPack &notifications_pack,	std::map<int, std::vector<size_t> > &masterOLLayout)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nProcessing notifications\n");

		std::vector<int> pids;
		for(StreamPack::iterator iter = notifications_pack.begin(); iter != notifications_pack.end(); ++iter)
			pids.push_back(iter->first);
		sort(pids.begin(), pids.end());

		for(size_t i=0; i<pids.size(); i++)
		{
			int pid = pids[i];

			UG_DLOG(LIB_ALG_MATRIX, 4, "PID " << pid << "\n");

			BinaryStream &stream = *notifications_pack.get_stream(pid);
			while(stream.can_read_more())
			{
				size_t local_index;
				int pid2;
				Deserialize(stream, local_index);
				Deserialize(stream, pid2);
				std::vector<bool> &mark = m_mapMark[pid2];

				UG_DLOG(LIB_ALG_MATRIX, 4, "Got a notification from processor " << pid << " that local index " << local_index << " has a copy on processor " << pid2 << "\n");

				UG_ASSERT(m_globalIDs[local_index].first == pcl::GetProcRank() && m_globalIDs[local_index].second == local_index, "got notification about a node not master on this proc?");
				UG_ASSERT(local_index < m_mat.num_rows(), "");

				mark.resize(m_globalIDs.size(), false);
				if(mark[local_index] == false)
				{
					UG_ASSERT(mark.size() == m_globalIDs.size(), "");
					mark[local_index]=true;

					// this is a notification that local_index is now slave on processor pid2.
					masterOLLayout[pid2].push_back(local_index);
				}
				else
					UG_DLOG(LIB_ALG_MATRIX, 4, " but is already.\n")
			}
		}
	}

	/// create_mark_map
	/**
	 * \param masterLayout the layout to mark nodes of
	 * for every interface I in masterLayout
	 * 		pid = associated processor to interface I
	 * 		mark all nodes in I so that we know they are already in some master interface to pid.
	 *
	 * \sa m_mapMark
	 */
	void create_mark_map(IndexLayout &masterLayout)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlap_CreateMarks\n");

		for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
		{
			IndexLayout::Interface &interface = masterLayout.interface(iter);
			int pid = masterLayout.proc_id(iter);

			UG_DLOG(LIB_ALG_MATRIX, 4, "PID " << pid << "\n");

			std::vector<bool> &mark = m_mapMark[pid];
			mark.resize(m_mat.num_rows(), false);
			UG_DLOG(LIB_ALG_MATRIX, 4, "created marks for interface to processor " << pid << ": ");

			for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				size_t local_index = interface.get_element(iter2);
				UG_DLOG(LIB_ALG_MATRIX, 4, local_index << " ");
				mark[local_index] = true;
			}
		}
	}


	/// get_new_indices
	/**
	 * \param matrixrow_pack data received from other processors
	 * \param overlapLayout layout where to add not yet existing nodes
	 * \param nodeNummerator tells us if globalID has a local ID or creates one
	 */
	void get_new_indices(StreamPack &matrixrow_pack, std::map<int, std::vector<size_t> > &overlapLayout,
			NewNodesNummerator &nodeNummerator)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlap_GetNewIndices\n");
		typedef IndexLayout::Interface Interface;

		size_t num_connections;
		AlgebraID globalRowIndex, globalColIndex;
		typename matrix_type::value_type value;

		std::vector<int> pids;
		for(StreamPack::iterator iter = matrixrow_pack.begin(); iter != matrixrow_pack.end(); ++iter)
			pids.push_back(iter->first);

		sort(pids.begin(), pids.end());
		bool has_index; bool bCreated;

		for(size_t i=0; i<pids.size(); i++)
		{
			int pid = pids[i];
			UG_DLOG(LIB_ALG_MATRIX, 4, "PID " << pid << "\n");

			// copy stream (bug in BinaryStream, otherwise I would use stream.seekg(ios::begin)
			BinaryStream &stream2 = *matrixrow_pack.get_stream(pids[i]);
			BinaryStream stream; stream.write((const char*)stream2.buffer(), stream2.size());

			while(stream.can_read_more())
			{
				Deserialize(stream, globalRowIndex);

				size_t localRowIndex = nodeNummerator.get_index_if_available(globalRowIndex, has_index);
				UG_ASSERT(has_index, "global id " << globalRowIndex.first << " | " << globalRowIndex.second << " has no associated local id");
				UG_DLOG(LIB_ALG_MATRIX, 4, "GID " << globalRowIndex.first << " | " << globalRowIndex.second << " has local ID " << localRowIndex << "\n");

				Deserialize(stream, num_connections);

				for(size_t i=0; i<num_connections; i++)
				{
					Deserialize(stream, globalColIndex);
					size_t localColIndex = nodeNummerator.get_index_or_create_new(globalColIndex, bCreated);
					UG_DLOG(LIB_ALG_MATRIX, 4, " connection GID " << globalColIndex.first << " | " << globalColIndex.second << " has local ID " << localColIndex <<
							(bCreated ? " and was created\n" : "\n"));
					if(bCreated)
						overlapLayout[globalColIndex.first].push_back(localColIndex);
					Deserialize(stream, value);
				}
			}

		}

		UG_DLOG(LIB_ALG_MATRIX, 4, "\n");
	}

	/// process_matrix_rows
	/**
	 * \param matrixrow_pack 	StreamPack where to extract data from
	 * \param nodeNummerator	maps global indices to local indices
	 * \param 					bSet if true, use overwrite (set) for received matrix rows. else add matrix rows.
	 *
	 */
	void process_matrix_rows(StreamPack &matrixrow_pack, NewNodesNummerator &nodeNummerator, bool bSet)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\niterate again over all streams to get the matrix lines\n");

		size_t num_connections;
		AlgebraID globalRowIndex, globalColIndex;
		std::vector<typename matrix_type::connection> cons;
		bool has_index;

		for(StreamPack::iterator iter = matrixrow_pack.begin(); iter != matrixrow_pack.end(); ++iter)
		{
			int pid = iter->first;
			BinaryStream &stream = *iter->second;

			UG_DLOG(LIB_ALG_MATRIX, 4, "processing data for master interface with processor " << pid << ". size is " << stream.size() << "\n");

			while(stream.can_read_more())
			{
				Deserialize(stream, globalRowIndex);

				// get local index for this global index
				size_t localRowIndex = nodeNummerator.get_index_if_available(globalRowIndex, has_index);
				UG_ASSERT(has_index, "global id " << globalRowIndex.first << " | " << globalRowIndex.second << " has no associated local id");

				UG_DLOG(LIB_ALG_MATRIX, 4, "Processing global id " << globalRowIndex.first << " | " << globalRowIndex.second << ", local id " << localRowIndex << ". ");

				// get nr of connections
				Deserialize(stream, num_connections);
				if(cons.size() < num_connections) cons.resize(num_connections);
				UG_DLOG(LIB_ALG_MATRIX, 4, num_connections << " connections:\n");

				size_t j=0;
				for(size_t i=0; i<num_connections; i++)
				{
					// get global column index, and associated local index
					Deserialize(stream, globalColIndex);
					Deserialize(stream, cons[j].dValue);
					cons[j].iIndex = nodeNummerator.get_index_if_available(globalColIndex, has_index);

					UG_DLOG(LIB_ALG_MATRIX, 4, cons[j].iIndex << " (global index " << globalColIndex.first << " | " << globalColIndex.second << ")");
					UG_DLOG(LIB_ALG_MATRIX, 4, " -> " << cons[j].dValue << "\n");

					//UG_ASSERT(has_index, "global id " << globalColIndex.first << " | " << globalColIndex.second << " has no associated local id");
					if(!has_index)
						UG_DLOG(LIB_ALG_MATRIX, 4, "has no index\n")
					else j++;
				}

				// set matrix row
				UG_DLOG(LIB_ALG_MATRIX, 4, "set matrix row: ");
				for(size_t i=0; i<j; i++)
					UG_DLOG(LIB_ALG_MATRIX, 4, "(" << cons[i].iIndex << "-> " << cons[i].dValue << ") ");
				UG_DLOG(LIB_ALG_MATRIX, 4, "\n");

				if(bSet)
					m_newMat.set_matrix_row(localRowIndex, &cons[0], j);
				else
				{
					m_newMat.add_matrix_row(localRowIndex, &cons[0], j);
				}

			}
		}
	}


	// GenerateOverlap
	//--------------------------------------------------
	/**
	 * \brief Generates a new matrix with overlap from another matrix
	 * \param _mat				matrix to create overlap from
	 * \param newMat			matrix to store overlap matrix in
	 * \param masterOLLayout	Layout
	 * \param masterOLLayout
	 * \param overlap_depth
	 *
	 * \sa GenerateOverlap_CollectData
	 * \sa GenerateOverlap_GetNewIndices
	 * \sa GenerateOverlap_AddMatrixRows

	 */

	/// insert_map_into_layout_sorted
	/**
	 * \param m the map mapping a processor id to a std::vector<size_t>
	 * \param layout the layout where to add to
	 *
	 * master and slave lists can be permuted due to the construction.
	 * we sort these lists by taking the globalIDs as sort criteria.
	 * then master and associated slave interface are correct
	 */
	void insert_map_into_layout_sorted(std::map<int, std::vector<size_t> > &m, IndexLayout &layout)
	{
		for(std::map<int, std::vector<size_t> >::iterator iter = m.begin(); iter != m.end(); ++iter)
		{
			int pid = iter->first;
			std::vector<size_t> &v = iter->second;

			IndicesSort(v.begin(), v.end(), m_globalIDs);
			Interface &interface = layout.interface(pid);
			for(size_t i=0; i<v.size(); i++)
				interface.push_back(v[i]);
		}
	}

	/// communicate
	/**
	 * \brief creates one overlap level one-sided
	 * \param sendingNodesLayout Layout of which nodes send their matrix rows
	 * \param receivingNodesLayout Layout of which nodes receive matrix rows
	 * \param bCreateNewNodes if true, create new indices/interfaces for globalIndices which are not yet on this or the other processor
	 * \param newSlavesLayout new slaves are added to this layout
	 * \param newMastersLayout new masters are added to this layout
	 * \param pids used to send/receive notifications
	 * \param nodeNummerator used to number (new) indices
	 * \param bSet if true, use overwrite (set) for received matrix rows. else add matrix rows.
	 * \param level current overlap level
	 */

	void communicate(IndexLayout &sendingNodesLayout, IndexLayout &receivingNodesLayout,
			bool bCreateNewNodes,
			IndexLayout &newSlavesLayout, IndexLayout &newMastersLayout,
			std::set<int> &pids, NewNodesNummerator &nodeNummerator, bool bSet, size_t level)
	{

		std::map<int, std::vector<size_t> > newSlaves;
		std::map<int, std::vector<size_t> > newMasters;

		collect_matrix_row_data(m_newMat, sendingNodesLayout, bCreateNewNodes, pids, newMasters);

		// receive data
		//-----------------


		StreamPack notifications_pack, matrixrow_pack;

		for(IndexLayout::iterator iter = receivingNodesLayout.begin(); iter != receivingNodesLayout.end();
				++iter)
		{
			int pid = receivingNodesLayout.proc_id(iter);
			UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": preparing buffer for matrixrow data from proc " << pid << "\n");
			m_com.receive_raw(pid, *matrixrow_pack.get_stream(pid));
		}

		if(bCreateNewNodes)
		{
			for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
			{
				size_t pid = *iter;
				UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": preparing buffer for notifications data from proc " << pid << "\n");
				m_com.receive_raw(pid, *notifications_pack.get_stream(pid));
			}
		}


		try
		{
			m_com.communicate();
		}
		catch(...)
		{
			UG_ASSERT(0, "unknown error here");
		}

		IF_DEBUG(LIB_ALG_MATRIX, 4)
		{
			for(IndexLayout::iterator iter = receivingNodesLayout.begin(); iter != receivingNodesLayout.end();
					++iter)
			{
				int pid = receivingNodesLayout.proc_id(iter);
				UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": received " << matrixrow_pack.get_stream(pid)->size() << " bytes of matrixrow data from proc " << pid << "\n");
			}

			if(bCreateNewNodes)
			{
				for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
				{
					size_t pid = *iter;
					UG_ASSERT(notifications_pack.has_stream(pid), "pid = " << pid);
					UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": received " << notifications_pack.get_stream(pid)->size() << " bytes of notifications data from proc " << pid << "\n");
				}
			}
		}
		// process data
		//-----------------

		if(bCreateNewNodes)
		{
			receive_notifications(notifications_pack, newMasters);
			//GenerateOverlap_CreateMarks(pids, slaveLayout, masterLayout, master_map, newMat.num_rows());

			// get all new nodes and their indices
			get_new_indices(matrixrow_pack, newSlaves, nodeNummerator);


			insert_map_into_layout_sorted(newMasters, newMastersLayout);
			insert_map_into_layout_sorted(newSlaves, newSlavesLayout);

			AddLayout(m_totalMasterLayout, newMastersLayout);
			AddLayout(m_totalSlaveLayout, newSlavesLayout);

			AddLayout(m_vMasterLayouts[level], newMastersLayout);
			AddLayout(m_vSlaveLayouts[level], newSlavesLayout);


			for(IndexLayout::iterator iter = newMastersLayout.begin(); iter != newMastersLayout.end(); ++iter)
				pids.insert(newMastersLayout.proc_id(iter));
			for(IndexLayout::iterator iter = newSlavesLayout.begin(); iter != newSlavesLayout.end(); ++iter)
				pids.insert(newSlavesLayout.proc_id(iter));
		}

		// create as copy & resize matrix to be now new_indices_size x new_indices_size.
		size_t new_indices_size = nodeNummerator.get_new_indices_size();
		UG_DLOG(LIB_ALG_MATRIX, 4, "resize matrix to be now " << new_indices_size << " x " << new_indices_size << ".\n");

		m_newMat.resize(new_indices_size, new_indices_size);

		// iterate again over all streams to get the matrix lines
		process_matrix_rows(matrixrow_pack, nodeNummerator, bSet);

		// send master rows to slaves
		IF_DEBUG(LIB_ALG_MATRIX, 4)
			if(bCreateNewNodes)
			{
				UG_DLOG(LIB_ALG_MATRIX, 4, "master OL Layout:\n");
				PrintLayout(newMastersLayout);
				UG_DLOG(LIB_ALG_MATRIX, 4, "slave OL Layout:\n");
				PrintLayout(newSlavesLayout);
				UG_DLOG(LIB_ALG_MATRIX, 4, "OL Layout:\n");
				PrintLayout(m_com, newMastersLayout, newSlavesLayout);
			}
	}

public:
	GenerateOverlapClass(matrix_type &mat, matrix_type &newMat, IndexLayout &totalMasterLayout, IndexLayout &totalSlaveLayout,
			std::vector<IndexLayout> vMasterLayouts, std::vector<IndexLayout> vSlaveLayouts, size_t overlapDepth=1) :
		m_com(mat.get_communicator()), m_mat(mat), m_newMat(newMat), m_totalMasterLayout(totalMasterLayout), m_totalSlaveLayout(totalSlaveLayout),
		m_vMasterLayouts(vMasterLayouts), m_vSlaveLayouts(vSlaveLayouts), m_overlapDepth(overlapDepth)
	{

	}

	/// calculate
	/**
	 * calculates overlap
	 */
	bool calculate()
	{
		IF_DEBUG(LIB_ALG_MATRIX, 4)
		{
			UG_DLOG(LIB_ALG_MATRIX, 4, "GENERATE OVERLAP START\n");

			UG_DLOG(LIB_ALG_MATRIX, 4, "matrix is " << m_mat.num_rows() << " x " << m_mat.num_cols() << "\n");
			m_mat.print();
		}
		UG_ASSERT(m_mat.num_rows() == m_mat.num_cols(), "atm only for square matrices");

		IndexLayout &masterLayout = m_mat.get_master_layout();
		IndexLayout &slaveLayout = m_mat.get_slave_layout();

		TestLayout(m_com, masterLayout, slaveLayout);

		// generate global algebra indices
		UG_DLOG(LIB_ALG_MATRIX, 4, "generate " << m_mat.num_rows() << " m_globalIDs\n");
		GenerateGlobalAlgebraIDs(m_globalIDs, m_mat.num_rows(), masterLayout, slaveLayout);

		IF_DEBUG(LIB_ALG_MATRIX, 4)
		{
			for(size_t i=0; i<m_globalIDs.size(); i++)
				UG_DLOG(LIB_ALG_MATRIX, 4, "  " << i << ": global id " << m_globalIDs[i].first << " | " << m_globalIDs[i].second << "\n")
		}

		m_newMat.create_as_copy_of(m_mat);

		NewNodesNummerator nodeNummerator(m_globalIDs);


		// collect data
		//-----------------

		std::vector<IndexLayout> masterOLLayouts, slaveOLLayouts;
		masterOLLayouts.clear();
		masterOLLayouts.resize(m_overlapDepth+1);
		slaveOLLayouts.clear();
		slaveOLLayouts.resize(m_overlapDepth+1);

		std::vector<IndexLayout> backward_masterOLLayouts, backward_slaveOLLayouts;
		backward_masterOLLayouts.clear();
		backward_masterOLLayouts.resize(m_overlapDepth+1);
		backward_slaveOLLayouts.clear();
		backward_slaveOLLayouts.resize(m_overlapDepth+1);

		// TODO: try to remove these pid numbers or reduce them by introducing receivePIDs, sendPIDs
		// these are necessary because notifications can occur to a processor not in the current layout
		std::set<int> pids;
		for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
			pids.insert(iter->first);
		for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
			pids.insert(iter->first);


		create_mark_map(masterLayout);
		m_vMasterLayouts.resize(m_overlapDepth+1);
		m_vSlaveLayouts.resize(m_overlapDepth+1);
		AddLayout(m_vMasterLayouts[0], masterLayout);
		AddLayout(m_vSlaveLayouts[0], slaveLayout);

		for(size_t current_overlap=0; current_overlap <= m_overlapDepth; current_overlap++)
		{
			bool bCreateNewNodes = (current_overlap == m_overlapDepth ? false : true);

			IF_DEBUG(LIB_ALG_MATRIX, 4)
			{
				UG_DLOG(LIB_ALG_MATRIX, 4, "\n---------------------\ncurrentOL: " << current_overlap << "\n");
				if(bCreateNewNodes)
					UG_DLOG(LIB_ALG_MATRIX, 4, "(creating New Nodes)\n");
				UG_DLOG(LIB_ALG_MATRIX, 4, "---------------------\n\n");
			}

			IndexLayout *send_layout;
			if(current_overlap == 0)
				send_layout = &slaveLayout;
			else
				send_layout = &masterOLLayouts[current_overlap-1];

			IndexLayout *receive_layout;
			if(current_overlap == 0)
				receive_layout = &masterLayout;
			else
				receive_layout = &slaveOLLayouts[current_overlap-1];


			communicate(*send_layout, *receive_layout, bCreateNewNodes,
					slaveOLLayouts[current_overlap], masterOLLayouts[current_overlap], pids,
					nodeNummerator, false, current_overlap+1);

			// backwards

			IndexLayout *backward_send_layout;
			if(current_overlap == 0)
				backward_send_layout = &masterLayout;
			else
				backward_send_layout = &backward_masterOLLayouts[current_overlap-1];

			IndexLayout *backward_receive_layout;
			if(current_overlap == 0)
				backward_receive_layout = &slaveLayout;
			else
				backward_receive_layout = &backward_slaveOLLayouts[current_overlap-1];

			communicate(*backward_send_layout, *backward_receive_layout, bCreateNewNodes,
					backward_slaveOLLayouts[current_overlap], backward_masterOLLayouts[current_overlap], pids,
								nodeNummerator, true, current_overlap+1);

			// done!
		}


		m_newMat.set_layouts(m_totalMasterLayout, m_totalSlaveLayout);
		m_newMat.set_communicator(m_mat.get_communicator());
		m_newMat.set_process_communicator(m_mat.get_process_communicator());
		m_newMat.copy_storage_type(m_mat);



		IF_DEBUG(LIB_ALG_MATRIX, 4)
		{
			UG_DLOG(LIB_ALG_MATRIX, 4, "new matrix\n\n");
			m_newMat.print();

			UG_DLOG(LIB_ALG_MATRIX, 4, "master OL Layout:\n");
			PrintLayout(m_totalMasterLayout);
			UG_DLOG(LIB_ALG_MATRIX, 4, "slave OL Layout:\n");
			PrintLayout(m_totalSlaveLayout);

			UG_DLOG(LIB_ALG_MATRIX, 4, "OL Layout:\n");
			PrintLayout(m_com, m_totalMasterLayout, m_totalSlaveLayout);
		}

		return true;

	}
};

// TODO: one "bug" remains: dirichlet nodes, which have only connection to themselfs = 1.0, get afterwards 2.0 (because rows are not additive there)
template<typename matrix_type>
bool GenerateOverlap(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat,
		IndexLayout &totalMasterLayout, IndexLayout &totalSlaveLayout, std::vector<IndexLayout> vMasterLayouts, std::vector<IndexLayout> vSlaveLayouts,
		size_t overlapDepth=1)
{
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat, totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts, overlapDepth);
	return c.calculate();
}

}
#endif
