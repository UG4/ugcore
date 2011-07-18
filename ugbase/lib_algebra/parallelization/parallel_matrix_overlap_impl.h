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
#include "common/util/binary_buffer.h"

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
	typedef std::map<int, BinaryBuffer>	BufferMap;

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
	//size_t m_overlapDepth;						///< overlap depth to be achieved

public:
	std::vector<size_t> m_overlapSize;
	size_t m_overlapDepthMaster;
	size_t m_overlapDepthSlave;
	bool m_masterDirichletLast;
	bool m_slaveDirichletLast;
private:

	/// collect_matrix_row_data
	/**
	 * problem here is that we have to keep track of the interfaces. when we send a connection to a node
	 * which is not on that processor, we need to add an interface from that node to its master:
	 * if we are master of the node, it is simple. if we are not, we need to send a notification to the
	 * masters' processor
	 *
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
		PROFILE_FUNC();
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlapClass::collect_matrix_row_data\n");

		BufferMap notificationsPack;

		// zeilen verschicken
		for(IndexLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
		{
			Interface &interface = layout.interface(iter);
			int pid = layout.proc_id(iter);

			// get a stream to put data to process pid to
			BinaryBuffer stream;
			std::vector<bool> &mark = m_mapMark[pid];

			UG_DLOG(LIB_ALG_MATRIX, 4, "PID " << pid << "\n");

			if(mark.size() == 0) mark.resize(m_mat.num_rows(), false);


			for(Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				size_t localIndex = interface.get_element(iter2);

				AlgebraID &globalRowIndex = m_globalIDs[localIndex];

				std::vector<size_t> *p_vNewInterface = NULL;

				// serialize global row index
				Serialize(stream, globalRowIndex);

				size_t numConnections = mat.num_connections(localIndex);

				// serialize number of connections
				Serialize(stream, numConnections);

				UG_DLOG(LIB_ALG_MATRIX, 4, "GID " << globalRowIndex << "\n");

				// serialize connections
				for(typename matrix_type::const_row_iterator conn = mat.begin_row(localIndex); conn != mat.end_row(localIndex); ++conn)
				{
					AlgebraID &globalColIndex = m_globalIDs[conn.index()];
					size_t localColIndex = conn.index();

					UG_DLOG(LIB_ALG_MATRIX, 4, "  " << conn.index() << " (global id " <<
							globalColIndex << ") -> " << conn.value() << "\n");

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
									BinaryBuffer &notificationsStream = notificationsPack[globalColIndex.first];

									// serialize notification: globalColIndex.second, pid.
									Serialize(notificationsStream, globalColIndex.second);
									Serialize(notificationsStream, pid);
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
			UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": sending matrixrow data to proc " << pid << " (size " << stream.write_pos() << ")\n");
			m_com.send_raw(pid, stream.buffer(), stream.write_pos(), false);
		}

		if(bAddNodesToLayout)
		{
			for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
			{
				int pid = *iter;
				BinaryBuffer &stream = notificationsPack[pid];
				m_com.send_raw(pid, stream.buffer(), stream.write_pos(), false);
				UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": sending notifications data to proc " << pid << " (size " << stream.write_pos() << ") \n");
			}
		}
	}

	/// receive_notifications
	/**
	 * processing notifications about our master nodes having slaves on other processors
	 * \param notificationsPack notifications stream in form (size_t) localIndex, (int) pid.
	 * \param masterOLLayout layout where to add masters if appropriate
	 *
	 * for a notification that the index localIndex has now a copy on processor pid, we add localIndex to masterOLLayout[pid], if
	 * localIndex is not already in some master interface to pid.
	 *
	 * \sa m_mapMark
	 */
	void receive_notifications(BufferMap &notificationsPack,	std::map<int, std::vector<size_t> > &masterOLLayout)
	{
		PROFILE_FUNC();
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nProcessing notifications\n");

		std::vector<int> pids;
		for(BufferMap::iterator iter = notificationsPack.begin(); iter != notificationsPack.end(); ++iter)
			pids.push_back(iter->first);
		sort(pids.begin(), pids.end());

		for(size_t i=0; i<pids.size(); i++)
		{
			int pid = pids[i];

			UG_DLOG(LIB_ALG_MATRIX, 4, "PID " << pid << "\n");

			BinaryBuffer &stream = notificationsPack[pid];
			while(!stream.eof())
			{
				size_t localIndex;
				int pid2;
				Deserialize(stream, localIndex);
				Deserialize(stream, pid2);
				std::vector<bool> &mark = m_mapMark[pid2];

				UG_DLOG(LIB_ALG_MATRIX, 4, "Got a notification from processor " << pid << " that local index " << localIndex << " has a copy on processor " << pid2 << "\n");

				UG_ASSERT(m_globalIDs[localIndex].first == pcl::GetProcRank() && m_globalIDs[localIndex].second == localIndex, "got notification about a node not master on this proc?");
				UG_ASSERT(localIndex < m_mat.num_rows(), "");

				mark.resize(m_globalIDs.size(), false);
				if(mark[localIndex] == false)
				{
					UG_ASSERT(mark.size() == m_globalIDs.size(), "");
					mark[localIndex]=true;

					// this is a notification that localIndex is now slave on processor pid2.
					masterOLLayout[pid2].push_back(localIndex);
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
		PROFILE_FUNC();
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlapClass::CreateMarks\n");

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
				size_t localIndex = interface.get_element(iter2);
				UG_DLOG(LIB_ALG_MATRIX, 4, localIndex << " ");
				mark[localIndex] = true;
			}
		}
	}


	/// get_new_indices
	/**
	 * \param matrixrowPack data received from other processors
	 * \param overlapLayout layout where to add not yet existing nodes
	 * \param nodeNummerator tells us if globalID has a local ID or creates one
	 */
	void get_new_indices(BufferMap &matrixrowPack, std::map<int, std::vector<size_t> > &overlapLayout,
			NewNodesNummerator &nodeNummerator)
	{
		PROFILE_FUNC();
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlapClass::GetNewIndices\n");
		typedef IndexLayout::Interface Interface;

		size_t numConnections;
		AlgebraID globalRowIndex, globalColIndex;
		typename matrix_type::value_type value;

		std::vector<int> pids;
		for(BufferMap::iterator iter = matrixrowPack.begin(); iter != matrixrowPack.end(); ++iter)
			pids.push_back(iter->first);

		sort(pids.begin(), pids.end());
		bool bCreated;

		for(size_t i=0; i<pids.size(); i++)
		{
			UG_DLOG(LIB_ALG_MATRIX, 4, "PID " << pids[i] << "\n");

			// copy stream (bug in BinaryStream, otherwise I would use stream.seekg(ios::begin)
			BinaryBuffer &stream2 = matrixrowPack[pids[i]];
			BinaryBuffer stream; stream.write((const char*)stream2.buffer(), stream2.write_pos());

			while(!stream.eof())
			{
				Deserialize(stream, globalRowIndex);

				Deserialize(stream, numConnections);

#ifdef UG_ENABLE_DEBUG_LOGS
				IF_DEBUG(LIB_ALG_MATRIX, 4)
				{
					size_t localRowIndex = nodeNummerator[globalRowIndex];
					UG_DLOG(LIB_ALG_MATRIX, 4, "GID " << globalRowIndex << " has local ID " << localRowIndex << ", and " << numConnections << " connections\n");
				}
#endif



				for(size_t i=0; i<numConnections; i++)
				{
					Deserialize(stream, globalColIndex);
					size_t localColIndex = nodeNummerator.get_index_or_create_new(globalColIndex, bCreated);
					UG_DLOG(LIB_ALG_MATRIX, 4, " connection GID " << globalColIndex << " has local ID " << localColIndex <<
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
	 * \param matrixrowPack 	StreamPack where to extract data from
	 * \param nodeNummerator	maps global indices to local indices
	 * \param 					bSet if true, use overwrite (set) for received matrix rows. else add matrix rows.
	 *
	 */
	void process_matrix_rows(BufferMap &matrixrowPack, NewNodesNummerator &nodeNummerator, bool bSet)
	{
		PROFILE_FUNC();
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\niterate again over all streams to get the matrix lines\n");

		size_t numConnections;
		AlgebraID globalRowIndex, globalColIndex;
		std::vector<typename matrix_type::connection> cons;
		bool has_index;

		for(BufferMap::iterator iter = matrixrowPack.begin(); iter != matrixrowPack.end(); ++iter)
		{
			BinaryBuffer &stream = iter->second;

			UG_DLOG(LIB_ALG_MATRIX, 4, "processing data for master interface with processor " << iter->first << ". size is " << stream.write_pos() << "\n");

			while(!stream.eof())
			{
				Deserialize(stream, globalRowIndex);

				// get local index for this global index
				size_t localRowIndex = nodeNummerator[globalRowIndex];

				UG_DLOG(LIB_ALG_MATRIX, 4, "Processing global id " << globalRowIndex << ", local id " << localRowIndex << ". ");

				// get nr of connections
				Deserialize(stream, numConnections);
				if(cons.size() < numConnections) cons.resize(numConnections);
				UG_DLOG(LIB_ALG_MATRIX, 4, numConnections << " connections:\n");

				size_t j=0;
				for(size_t i=0; i<numConnections; i++)
				{
					// get global column index, and associated local index
					Deserialize(stream, globalColIndex);
					Deserialize(stream, cons[j].dValue);
					cons[j].iIndex = nodeNummerator.get_index_if_available(globalColIndex, has_index);

					if(has_index) 	{ UG_DLOG(LIB_ALG_MATRIX, 4, cons[j].iIndex);   }
					else			{ UG_DLOG(LIB_ALG_MATRIX, 4, "n/a"); }
					UG_DLOG(LIB_ALG_MATRIX, 4,  " (global index " << globalColIndex << ")" << " -> " << cons[j].dValue << "\n");

					//UG_ASSERT(has_index, "global id " << globalColIndex << " has no associated local id");
					if(has_index) j++;
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
		PROFILE_FUNC();
		for(std::map<int, std::vector<size_t> >::iterator iter = m.begin(); iter != m.end(); ++iter)
		{
			int pid = iter->first;
			std::vector<size_t> &v = iter->second;

			sort(v.begin(), v.end(), CompareIndicesBy(m_globalIDs) );
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
		PROFILE_FUNC();

		std::map<int, std::vector<size_t> > newSlaves;
		std::map<int, std::vector<size_t> > newMasters;

		collect_matrix_row_data(m_newMat, sendingNodesLayout, bCreateNewNodes, pids, newMasters);

		// receive data
		//-----------------


		BufferMap notificationsPack, matrixrowPack;

		for(IndexLayout::iterator iter = receivingNodesLayout.begin(); iter != receivingNodesLayout.end();
				++iter)
		{
			int pid = receivingNodesLayout.proc_id(iter);
			UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": preparing buffer for matrixrow data from proc " << pid << "\n");
			m_com.receive_raw(pid, matrixrowPack[pid]);
		}

		if(bCreateNewNodes)
		{
			for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
			{
				size_t pid = *iter;
				UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": preparing buffer for notifications data from proc " << pid << "\n");
				m_com.receive_raw(pid, notificationsPack[pid]);
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
				UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": received " << matrixrowPack[receivingNodesLayout.proc_id(iter)].write_pos() << " bytes of matrixrow data from proc " <<
						receivingNodesLayout.proc_id(iter) << "\n");
			}

			if(bCreateNewNodes)
			{
				for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
				{
					//size_t pid = *iter;
					UG_ASSERT(notificationsPack.find(*iter) != notificationsPack.end(), "pid = " << *iter);
					UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": received " << notificationsPack[*iter].write_pos() << " bytes of notifications data from proc " << *iter << "\n");
				}
			}
		}
		// process data
		//-----------------

		if(bCreateNewNodes)
		{
			receive_notifications(notificationsPack, newMasters);
			//GenerateOverlap_CreateMarks(pids, slaveLayout, masterLayout, master_map, newMat.num_rows());

			// get all new nodes and their indices
			get_new_indices(matrixrowPack, newSlaves, nodeNummerator);


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
		process_matrix_rows(matrixrowPack, nodeNummerator, bSet);

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
		 */
	GenerateOverlapClass(matrix_type &mat, matrix_type &newMat, IndexLayout &totalMasterLayout, IndexLayout &totalSlaveLayout,
			std::vector<IndexLayout> &vMasterLayouts, std::vector<IndexLayout> &vSlaveLayouts) :
		m_com(mat.get_communicator()), m_mat(mat), m_newMat(newMat), m_totalMasterLayout(totalMasterLayout), m_totalSlaveLayout(totalSlaveLayout),
		m_vMasterLayouts(vMasterLayouts), m_vSlaveLayouts(vSlaveLayouts)
	{

	}

	/// calculate
	/**
	 * calculates overlap
	 */
	bool calculate()
	{
		PROFILE_FUNC();
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
				UG_DLOG(LIB_ALG_MATRIX, 4, "  " << i << ": global id " << m_globalIDs[i] << "\n")
		}

		m_newMat.set_as_copy_of(m_mat);

		NewNodesNummerator nodeNummerator(m_globalIDs);


		// collect data
		//-----------------

		size_t maxOverlap = std::max(m_overlapDepthMaster, m_overlapDepthSlave);

		std::vector<IndexLayout> masterOLLayouts, slaveOLLayouts;
		masterOLLayouts.clear();
		masterOLLayouts.resize(maxOverlap+1);
		slaveOLLayouts.clear();
		slaveOLLayouts.resize(maxOverlap+1);

		std::vector<IndexLayout> backward_masterOLLayouts, backward_slaveOLLayouts;
		backward_masterOLLayouts.clear();
		backward_masterOLLayouts.resize(maxOverlap+1);
		backward_slaveOLLayouts.clear();
		backward_slaveOLLayouts.resize(maxOverlap+1);

		// TODO: try to remove these pid numbers or reduce them by introducing receivePIDs, sendPIDs
		// these are necessary because notifications can occur from a processor not in the current layout
		std::set<int> pids;
		for(IndexLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
			pids.insert(iter->first);
		for(IndexLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
			pids.insert(iter->first);


		create_mark_map(masterLayout);
		m_vMasterLayouts.resize(maxOverlap+1);
		m_vSlaveLayouts.resize(maxOverlap+1);
		AddLayout(m_vMasterLayouts[0], masterLayout);
		AddLayout(m_vSlaveLayouts[0], slaveLayout);
		m_overlapSize.clear();
		for(size_t current_overlap=0; current_overlap <= maxOverlap; current_overlap++)
		{
			m_overlapSize.push_back(m_newMat.num_rows());

			IF_DEBUG(LIB_ALG_MATRIX, 4)
			{
				UG_DLOG(LIB_ALG_MATRIX, 4, "\n---------------------\ncurrentOL: " << current_overlap << "\n");
				UG_DLOG(LIB_ALG_MATRIX, 4, "---------------------\n\n");
			}

			if(current_overlap <= m_overlapDepthMaster)
			{
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

				if(current_overlap == m_overlapDepthMaster && m_masterDirichletLast)
				{
					std::vector<IndexLayout::Element> vIndex;
					CollectUniqueElements(vIndex,  *receive_layout);
					SetDirichletRow(m_newMat, vIndex);
				}
				else
				{
					bool bCreateNewNodes = (current_overlap == m_overlapDepthMaster ? false : true);
					communicate(*send_layout, *receive_layout, bCreateNewNodes,
						slaveOLLayouts[current_overlap], masterOLLayouts[current_overlap], pids,
						nodeNummerator, false, current_overlap);
				}
			}

			// backwards
			if(current_overlap <= m_overlapDepthSlave)
			{
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

				if(current_overlap == m_overlapDepthSlave && m_slaveDirichletLast)
				{
					std::vector<IndexLayout::Element> vIndex;
					CollectUniqueElements(vIndex,  *backward_receive_layout);
					SetDirichletRow(m_newMat, vIndex);
				}
				else
				{
					bool bCreateNewNodes = (current_overlap == m_overlapDepthSlave ? false : true);
					communicate(*backward_send_layout, *backward_receive_layout, bCreateNewNodes,
						backward_slaveOLLayouts[current_overlap], backward_masterOLLayouts[current_overlap], pids,
									nodeNummerator, true, current_overlap+1);
				}
			}

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
		IndexLayout &totalMasterLayout, IndexLayout &totalSlaveLayout, std::vector<IndexLayout> &vMasterLayouts, std::vector<IndexLayout> &vSlaveLayouts,
		std::vector<size_t> &overlapSize,
		size_t overlapDepth=1)
{
	PROFILE_FUNC();
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat, totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts);
	c.m_overlapDepthMaster = overlapDepth;
	c.m_overlapDepthSlave = overlapDepth;
	c.m_masterDirichletLast = false;
	c.m_slaveDirichletLast = false;
	bool b = c.calculate();
	overlapSize = c.m_overlapSize;
	return b;
}

// TODO: one "bug" remains: dirichlet nodes, which have only connection to themselfs = 1.0, get afterwards 2.0 (because rows are not additive there)
template<typename matrix_type>
bool GenerateOverlap2(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat,
		IndexLayout &totalMasterLayout, IndexLayout &totalSlaveLayout, std::vector<IndexLayout> &vMasterLayouts, std::vector<IndexLayout> &vSlaveLayouts,
		size_t overlapDepthMaster, size_t overlapDepthSlave, bool masterDirichletLast, bool slaveDirichletLast)
{
	PROFILE_FUNC();
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat, totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts);
	c.m_overlapDepthMaster = overlapDepthMaster;
	c.m_overlapDepthSlave = overlapDepthSlave;
	c.m_masterDirichletLast = masterDirichletLast;
	c.m_slaveDirichletLast = slaveDirichletLast;
	return c.calculate();
}

// TODO: one "bug" remains: dirichlet nodes, which have only connection to themselfs = 1.0, get afterwards 2.0 (because rows are not additive there)
template<typename matrix_type>
bool MakeConsistent(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat)
{
	PROFILE_FUNC();
	IndexLayout totalMasterLayout, totalSlaveLayout;
	std::vector<IndexLayout> vMasterLayouts;
	std::vector<IndexLayout> vSlaveLayouts;
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat, totalMasterLayout, totalSlaveLayout, vMasterLayouts, vSlaveLayouts);
	c.m_overlapDepthMaster = 0;
	c.m_overlapDepthSlave = 0;
	c.m_masterDirichletLast = false;
	c.m_slaveDirichletLast = false;
	bool b = c.calculate();
	newMat.set_layouts(mat.get_master_layout(), mat.get_slave_layout());
	return b;
}

template<typename matrix_type>
bool MakeFullRowsMatrix(const ParallelMatrix<matrix_type> &mat, ParallelMatrix<matrix_type> &newMat)
{
	// GenerateOverlap2(mat, newMat, totalMasterLayout, totalSlaveLayout, masterLayouts, slaveLayouts, 1, 0, true, true);
	return true;
}


}
#endif
