/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION__ROW_SENDING_SCHEME_H_
#define __H__LIB_ALGEBRA__PARALLELIZATION__ROW_SENDING_SCHEME_H_

#include <map>
#include "parallel_nodes.h"
#include "parallelization_util.h"
#include "pcl/pcl.h"
#include "new_layout_creator.h"

namespace ug
{

template<typename matrix_type>
class RowSendingScheme
{
private:
	bool m_bCreateNewNodes;
	using BufferMap = std::map<int, BinaryBuffer>;
	BufferMap rowsBufferMap;
	const matrix_type &mat;
	using connection = typename matrix_type::connection;
	std::map<size_t, std::vector<connection> > connections;
	ParallelNodes &PN;
	size_t rowMax, colMax;

public:
	RowSendingScheme(matrix_type &_mat, ParallelNodes &_PN)
	: mat(_mat), PN(_PN)
	{

	}

	void set_create_new_nodes(bool bCreateNewNodes)
	{
		m_bCreateNewNodes = bCreateNewNodes;
	}


	/**
	 * for processor pid in sendLayout
	 *     for i in interface
	 *     		send row i to processor pid
	 * for processor pid in receivingLayout
	 * 		issue receive data
	 *  issue
	 * @param communicator
	 * @param sendLayout
	 * @param receiveLayout
	 * @sa issue_send
	 */
	void issue_send(pcl::InterfaceCommunicator<IndexLayout> &communicator,
			const IndexLayout &sendLayout, const IndexLayout &receiveLayout)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "*** RowSendingScheme::issue_send: ***\n");
		for(auto it = sendLayout.begin(); it != sendLayout.end(); ++it)
		{
			const IndexLayout::Interface& interface = sendLayout.interface(it);
			int pid = sendLayout.proc_id(it);
			BinaryBuffer buf;
			for(auto iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			{
				issue_send(buf, pid, interface.get_element(iter2));
			}
			UG_DLOG(LIB_ALG_MATRIX, 4, "Sending " << buf.write_pos() << " of data to processor " << pid << "\n");
			communicator.send_raw(pid, buf.buffer(), buf.write_pos(), false);
		}

		rowsBufferMap.clear();

		for(auto it = receiveLayout.begin(); it != receiveLayout.end(); ++it)
		{
			int pid = receiveLayout.proc_id(it);
			UG_DLOG(LIB_ALG_MATRIX, 4, "issue receive from processor " << pid << "\n");
			communicator.receive_raw(pid, rowsBufferMap[pid]);
		}

		PN.issue(communicator);

	}

	void process(const IndexLayout &receiveLayout)
	{
		rowMax=0; colMax=0;
		UG_DLOG(LIB_ALG_MATRIX, 4, "*** RowSendingScheme::process: ***\n");
		connections.clear();

		for(auto it = receiveLayout.begin(); it != receiveLayout.end(); ++it)
		{
			int pid = receiveLayout.proc_id(it);
			const IndexLayout::Interface& interface = receiveLayout.interface(it);

			BinaryBuffer &buf = rowsBufferMap[pid];
			UG_DLOG(LIB_ALG_MATRIX, 4, "rowsBufferMap: received " << buf.write_pos() << " bytes from processor " << pid << "\n");
			for(auto it = interface.begin(); it != interface.end(); ++it)
				process(buf, pid, interface.get_element(it));
		}

		PN.process();
	}


	void set_rows_in_matrix(matrix_type &mat)
	{
		resize_mat(mat);
		for(typename std::map<size_t, std::vector<connection> >::iterator it = connections.begin();
				it != connections.end(); ++it)
		{
			std::vector<connection> &cons = it->second;
			if(!cons.empty())
			  mat.set_matrix_row(it->first, &cons[0], cons.size());
		}
	}

	void add_rows_to_matrix(matrix_type &mat)
	{
		resize_mat(mat);
		for(typename std::map<size_t, std::vector<connection> >::iterator it = connections.begin();
				it != connections.end(); ++it)
		{
			std::vector<connection> &cons = it->second;
			if(!cons.empty())
			  mat.add_matrix_row(it->first, &cons[0], cons.size());
		}
	}


private:
	void resize_mat(matrix_type &mat)
	{
		size_t cols = std::max(colMax, mat.num_cols());
		size_t rows = std::max(rowMax, mat.num_rows());
		if(rows > mat.num_rows() || cols > mat.num_cols())
			mat.resize_and_keep_values(rows, cols);
	}

	/**
	 *  Put row 'localRowIndex' into the binary buffer 'buf'
	 *  PN keeps track about which indices exist on the recieving pid
	 * @param buf				buf to write it to
	 * @param pid				receiving pid
	 * @param localRowIndex		local row index
	 * \note modifies PN by specifying which nodes we sent to pid
	 */
	void issue_send(BinaryBuffer &buf, int pid, int localRowIndex)
	{
		size_t num_connections = mat.num_connections(localRowIndex);

		// serialize number of connections
		Serialize(buf, num_connections);
		UG_DLOG(LIB_ALG_MATRIX, 4, "sending to pid " << pid << " row " << localRowIndex << " (" << PN.local_to_global(localRowIndex) << "), " << num_connections << " connections \n");

		for(typename matrix_type::const_row_iterator conn = mat.begin_row(localRowIndex);
							conn != mat.end_row(localRowIndex); ++conn)
		{
			size_t localColIndex = conn.index();
			const AlgebraID &globalColIndex = PN.local_to_global(localColIndex);
			if(m_bCreateNewNodes)
				PN.create_node(globalColIndex, localColIndex, pid);

			// serialize connection
			Serialize(buf, globalColIndex);
			Serialize(buf, conn.value());
			UG_DLOG(LIB_ALG_MATRIX, 4, " " << localColIndex  << " (" << globalColIndex << ") -> " << conn.value() << "\n");
		}
		UG_DLOG(LIB_ALG_MATRIX, 4, "\n");
	}

	void process(BinaryBuffer &buf, int pid, size_t localRowIndex)
	{
		size_t num_connections;
		rowMax = std::max(rowMax, localRowIndex+1);

		// serialize number of connections
		Deserialize(buf, num_connections);

		UG_DLOG(LIB_ALG_MATRIX, 4, "processing received row " << localRowIndex << ", " << num_connections << " connections \n");

		size_t distanceToMasterOrInner = PN.distance_to_master_or_inner(localRowIndex);

		std::vector<connection> &cons = connections[localRowIndex];
		size_t i = cons.size();
		num_connections += cons.size();
		cons.resize(num_connections);
		size_t j=i;
		for(; i<num_connections; i++)
		{
			AlgebraID globalColIndex;
			Deserialize(buf, globalColIndex);
			Deserialize(buf, cons[j].dValue);
			bool bHasIndex=true;
			if(m_bCreateNewNodes)
				cons[j].iIndex = PN.create_slave_node(globalColIndex, distanceToMasterOrInner+1);
			else
				cons[j].iIndex = PN.get_local_index_if_available(globalColIndex, bHasIndex);

			UG_DLOG(LIB_ALG_MATRIX, 4, " " << (int)(cons[j].iIndex) << " (" << globalColIndex << ") -> " << cons[j].dValue << "\n");
			if(bHasIndex)
			{
				size_t k;
				for(k=0; k<j; k++)
				{
					if(cons[k].iIndex == cons[j].iIndex)
					{
						cons[k].dValue += cons[j].dValue;
						break;
					}
				}
				if(k==j)
				{
					colMax = std::max(colMax, cons[j].iIndex+1);
					++j;
				}
			}
		}
		cons.resize(j);
	}
};

}
#endif