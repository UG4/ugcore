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


template<typename matrix_type>
class GenerateOverlapClass
{
private:
	typedef IndexLayout::Interface Interface;

	pcl::ParallelCommunicator<IndexLayout> &m_com;

	std::vector<AlgebraID> m_globalIDs;

	std::map<int, std::vector<bool> > m_mapMark;
	matrix_type &m_mat;
	matrix_type &m_newMat;
	IndexLayout &m_masterOLLayout;
	IndexLayout &m_slaveOLLayout;
	size_t m_overlapDepth;


	// neuer algorithmus:
	// overlap 0
	// 1. slave-knoten verschicken ihre Zeile an den Master.
	// 2. Master nimmt diese entgegen, zeile wird für vorhandene Knoten addiert
	// 3. fertig.

	// overlap 1
	// 1. slave-knoten verschicken ihre Zeile an den Master.
	//    /!\ werden verknüpfungen zu anderen prozessoren verschickt, werden die prozessoren informiert
	//    /!\ unter umständen wird so ein prozessor "von 2 seiten" informiert. dann muss es eine
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

	void collect_matrix_row_data(const matrix_type &mat,
			IndexLayout &layout, bool create_new_nodes,
			std::map<int, std::vector<size_t> > &masterOLLayout, std::set<int> pids)
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

				AlgebraID &global_row_index = m_globalIDs[local_index];

				std::vector<size_t> *p_masterOLInterface = NULL;

				// serialize matrix row
				Serialize(stream, global_row_index);
				size_t num_connections = 0;
				for(typename matrix_type::cRowIterator conn = mat.beginRow(local_index); !conn.isEnd(); ++conn)
				{
					if(conn.value() == 0.0) continue;
					num_connections++;
				}
				Serialize(stream, num_connections);

				UG_DLOG(LIB_ALG_MATRIX, 4, "GID " << global_row_index.first << " | " << global_row_index.second << "\n");

				//UG_DLOG(LIB_ALG_MATRIX, 4, mat.num_connections(local_index) << " connections: \n");
				for(typename matrix_type::cRowIterator conn = mat.beginRow(local_index); !conn.isEnd(); ++conn)
				{
					if(conn.value() == 0.0) continue;
					AlgebraID &global_col_index = m_globalIDs[conn.index()];
					size_t local_col_index = conn.index();

					if(create_new_nodes == false)
					{
						UG_DLOG(LIB_ALG_MATRIX, 4, "not adding node to masterOLInterface because create_new_node = false\n");
					}
					else
					{
						if(global_col_index.first == pcl::GetProcRank())
						{
							UG_ASSERT(global_col_index.second == local_col_index, "");

							// only add nodes once to interface
							// for master nodes: mark[pid][i] = true means: we already send some matrix row with this index, so there is already a interface.
							if(mark[local_col_index] == false)
							{
								mark[local_col_index] = true;
								UG_DLOG(LIB_ALG_MATRIX, 4, " adding node to masterOLInterface .interface(" << pid << ").push_back(" << local_col_index << ")\n");
								if(p_masterOLInterface == NULL)
									p_masterOLInterface = &masterOLLayout[pid];
								p_masterOLInterface->push_back(local_col_index);
							}
							else
								UG_DLOG(LIB_ALG_MATRIX, 4, "  not adding node to masterOLInterface because it is already\n");
							//node_marked[local_col_index] = true;
						}
						else
							// wir informieren hiermit den "Besitzer" des Knotens darüber,
							// dass eine kopie des knotens nun auf dem Prozessor mit prozessor id pid liegt.
						{
							UG_DLOG(LIB_ALG_MATRIX, 4, "  not adding node to masterOLInterface because it is on another processor\n");
							if(global_col_index.first != pid) // obviously the owner knows that his node is on his processor
							{
								UG_DLOG(LIB_ALG_MATRIX, 4, " processor " << global_col_index.first << " gets informed that his local node "
																		<< global_col_index.second << " has a copy on processor " << pid << "\n");
								// for slave nodes: mark[pid][i]=true means we already send notification to the owner that local index i has a copy on processor pid.
								// commented out because we would need to resize m_mapMark[global_col_index.first]
								std::vector<bool> &otherMark = m_mapMark[pid];
								otherMark.resize(m_globalIDs.size(), false);

								if(otherMark[local_col_index] == false)
								{
									otherMark[local_col_index] = true;
									BinaryStream &notifications_stream = *notifications_pack.get_stream(global_col_index.first);

									Serialize(notifications_stream, global_col_index.second);
									Serialize(notifications_stream, pid);
									//UG_ASSERT(find(pids.begin(), pids.end(), global_col_index.first) != pids.end(), "sending notification to a processor " << global_col_index.first << ", not in pid list???");
								}
								else
									UG_DLOG(LIB_ALG_MATRIX, 4, " not sending notification because already send one\n");
							}
							else
								UG_DLOG(LIB_ALG_MATRIX, 4, "  not sending notification because owner knows that his node is on his\n");

						}
					}


					UG_DLOG(LIB_ALG_MATRIX, 4, "  " << conn.index() << " (global id " <<
							global_col_index.first << " | " << global_col_index.second << ") -> " << conn.value() << "\n");
					Serialize(stream, global_col_index);
					Serialize(stream, conn.value());
				}
			}

			// use communicator to send data
			UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": sending matrixrow data to proc " << pid << " (size " << stream.size() << ")\n");
			m_com.send_raw(pid, stream.buffer(), stream.size(), false);
		}

		for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			int pid = *iter;
			BinaryStream &stream = *notifications_pack.get_stream(pid);
			m_com.send_raw(pid, stream.buffer(), stream.size(), false);
			UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": sending notifications data to proc " << pid << " (size " << stream.size() << ") \n");
		}

	}

	void receive_notifications(StreamPack &notifications_pack,
			std::map<int, std::vector<size_t> > &masterOLLayout)
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

	// for each processor we need to have a list which of our master nodes exist on their processor
	// this is important because we sometimes will need to add them to interfaces

	// b = m_mapMark[pid][i] :
	// if i is a master node on this processor, b is true iff i has a copy on processor pid
	// if i is a slave node on this processor, b is true iff we sent a notification to its owner that i has a copy on pid.

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

	void get_new_indices(StreamPack &matrixrow_pack, std::map<int, std::vector<size_t> > &overlapLayout,
			NewNodesNummerator &nodeNummerator)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\nGenerateOverlap_GetNewIndices\n");
		typedef IndexLayout::Interface Interface;

		size_t num_connections;
		AlgebraID global_row_index, global_col_index;
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
				Deserialize(stream, global_row_index);

				size_t local_row_index = nodeNummerator.get_index_if_available(global_row_index, has_index);
				UG_ASSERT(has_index, "global id " << global_row_index.first << " | " << global_row_index.second << " has no associated local id");
				UG_DLOG(LIB_ALG_MATRIX, 4, "GID " << global_row_index.first << " | " << global_row_index.second << " has local ID " << local_row_index << "\n");

				Deserialize(stream, num_connections);

				for(size_t i=0; i<num_connections; i++)
				{
					Deserialize(stream, global_col_index);
					size_t local_col_index = nodeNummerator.get_index_or_create_new(global_col_index, bCreated);
					UG_DLOG(LIB_ALG_MATRIX, 4, " connection GID " << global_col_index.first << " | " << global_col_index.second << " has local ID " << local_col_index <<
							(bCreated ? " and was created\n" : "\n"));
					if(bCreated)
						overlapLayout[global_col_index.first].push_back(local_col_index);
					Deserialize(stream, value);
				}
			}

		}

		UG_DLOG(LIB_ALG_MATRIX, 4, "\n");
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
	 * \sa GenerateOverlap
	 */
	void add_matrix_rows(StreamPack &matrixrow_pack, NewNodesNummerator &nodeNummerator, bool bSet)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n\niterate again over all streams to get the matrix lines\n");

		size_t num_connections;
		AlgebraID global_row_index, global_col_index;
		vector<typename matrix_type::connection> cons;
		bool has_index;

		for(StreamPack::iterator iter = matrixrow_pack.begin(); iter != matrixrow_pack.end(); ++iter)
		{
			int pid = iter->first;
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

				size_t j=0;
				for(size_t i=0; i<num_connections; i++)
				{
					// get global column index, and associated local index
					Deserialize(stream, global_col_index);
					Deserialize(stream, cons[j].dValue);
					cons[j].iIndex = nodeNummerator.get_index_if_available(global_col_index, has_index);

					UG_DLOG(LIB_ALG_MATRIX, 4, cons[j].iIndex << " (global index " << global_col_index.first << " | " << global_col_index.second << ")");
					UG_DLOG(LIB_ALG_MATRIX, 4, " -> " << cons[j].dValue << "\n");

					//UG_ASSERT(has_index, "global id " << global_col_index.first << " | " << global_col_index.second << " has no associated local id");
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
					m_newMat.set_matrix_row(local_row_index, &cons[0], j);
				else
				{
					m_newMat.add_matrix_row(local_row_index, &cons[0], j);
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

	void insert_map_into_layout_sorted(std::map<int, std::vector<size_t> > &m, IndexLayout &layout)
	{
		for(std::map<int, std::vector<size_t> >::iterator iter = m.begin(); iter != m.end(); ++iter)
		{
			int pid = iter->first;
			std::vector<size_t> &v = iter->second;
			for(size_t i=0; i<v.size(); i++)

			IndicesSort(v.begin(), v.end(), m_globalIDs);
			Interface &interface = layout.interface(pid);
			for(size_t i=0; i<v.size(); i++)
				interface.push_back(v[i]);
		}
	}
	void communicate(IndexLayout &sendingNodesLayout, IndexLayout &receivingNodesLayout,
			bool bCreateNewNodes,
			IndexLayout &newSlavesLayout, IndexLayout &newMastersLayout,
			std::set<int> &pids, NewNodesNummerator &nodeNummerator, bool bSet)
	{

		std::map<int, std::vector<size_t> > newSlaves;
		std::map<int, std::vector<size_t> > newMasters;

		collect_matrix_row_data(m_newMat, sendingNodesLayout, bCreateNewNodes, newMasters,
				pids);


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

		for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
		{
			size_t pid = *iter;
			UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": preparing buffer for notifications data from proc " << pid << "\n");
			m_com.receive_raw(pid, *notifications_pack.get_stream(pid));
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

			for(std::set<int>::iterator iter = pids.begin(); iter != pids.end(); ++iter)
			{
				size_t pid = *iter;
				UG_ASSERT(notifications_pack.has_stream(pid), "pid = " << pid);
				UG_DLOG(LIB_ALG_MATRIX, 4, "proc " << pcl::GetProcRank() << ": received " << notifications_pack.get_stream(pid)->size() << " bytes of notifications data from proc " << pid << "\n");
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


			insert_map_into_layout_sorted(newSlaves, newSlavesLayout);
			insert_map_into_layout_sorted(newMasters, newMastersLayout);

			AddLayout(m_masterOLLayout, newMastersLayout);
			AddLayout(m_slaveOLLayout, newSlavesLayout);

			for(IndexLayout::iterator iter = newMastersLayout.begin(); iter != newMastersLayout.end(); ++iter)
				pids.insert(iter->first);
			for(IndexLayout::iterator iter = newSlavesLayout.begin(); iter != newSlavesLayout.end(); ++iter)
				pids.insert(iter->first);
		}

		// create as copy & resize matrix to be now new_indices_size x new_indices_size.
		size_t new_indices_size = nodeNummerator.get_new_indices_size();
		UG_DLOG(LIB_ALG_MATRIX, 4, "resize matrix to be now " << new_indices_size << " x " << new_indices_size << ".\n");

		m_newMat.resize(new_indices_size, new_indices_size);

		// iterate again over all streams to get the matrix lines
		add_matrix_rows(matrixrow_pack, nodeNummerator, bSet);

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
	GenerateOverlapClass(matrix_type &mat, matrix_type &newMat, IndexLayout &masterOLLayout, IndexLayout &slaveOLLayout, size_t overlapDepth=1) :
		m_com(mat.get_communicator()), m_mat(mat), m_newMat(newMat), m_masterOLLayout(masterOLLayout), m_slaveOLLayout(slaveOLLayout), m_overlapDepth(overlapDepth)
	{

	}

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
					nodeNummerator, false);

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
								nodeNummerator, true);

			// done!
		}


		m_newMat.set_layouts(m_masterOLLayout, m_slaveOLLayout);
		m_newMat.set_communicator(m_mat.get_communicator());
		m_newMat.set_process_communicator(m_mat.get_process_communicator());
		m_newMat.copy_storage_type(m_mat);



		IF_DEBUG(LIB_ALG_MATRIX, 4)
		{
			UG_DLOG(LIB_ALG_MATRIX, 4, "new matrix\n\n");
			m_newMat.print();

			UG_DLOG(LIB_ALG_MATRIX, 4, "master OL Layout:\n");
			PrintLayout(m_masterOLLayout);
			UG_DLOG(LIB_ALG_MATRIX, 4, "slave OL Layout:\n");
			PrintLayout(m_slaveOLLayout);

			UG_DLOG(LIB_ALG_MATRIX, 4, "OL Layout:\n");
			PrintLayout(m_com, m_masterOLLayout, m_slaveOLLayout);
		}

		return true;

	}
};

// TODO: one "bug" remains: dirichlet nodes, which have only connection to themselfs = 1.0, get afterwards 2.0 (because rows are not additive there)
template<typename matrix_type>
bool GenerateOverlap(const ParallelMatrix<matrix_type> &_mat, ParallelMatrix<matrix_type> &newMat, IndexLayout &masterOLLayout, IndexLayout &slaveOLLayout, size_t overlapDepth=1)
{
	// pcl does not use const much
	//UG_ASSERT(overlap_depth > 0, "overlap_depth has to be > 0");
	ParallelMatrix<matrix_type> &mat = const_cast<ParallelMatrix<matrix_type> &> (_mat);

	GenerateOverlapClass<ParallelMatrix<matrix_type> > c(mat, newMat, masterOLLayout, slaveOLLayout, overlapDepth);
	return c.calculate();
}

}
#endif
