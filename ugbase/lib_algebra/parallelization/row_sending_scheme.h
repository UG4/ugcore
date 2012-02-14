/*
 * row_sending_scheme.h
 *
 *  Created on: 05.10.2011
 *      Author: mrupp
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
	typedef std::map<int, BinaryBuffer>	BufferMap;
	BufferMap rowsBufferMap;
	const matrix_type &mat;
	typedef typename matrix_type::connection connection;
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

	void issue_send(pcl::ParallelCommunicator<IndexLayout> &communicator,
			IndexLayout &sendLayout, IndexLayout &receiveLayout)
	{
		UG_DLOG(LIB_ALG_MATRIX, 4, "*** RowSendingScheme::issue_send: ***\n");
		for(IndexLayout::iterator it = sendLayout.begin(); it != sendLayout.end(); ++it)
		{
			IndexLayout::Interface interface = sendLayout.interface(it);
			int pid = sendLayout.proc_id(it);
			BinaryBuffer buf;
			for(IndexLayout::Interface::iterator iter2 = interface.begin();
					iter2 != interface.end(); ++iter2)
			{
				issue_send(buf, pid, interface.get_element(iter2));
			}
			UG_DLOG(LIB_ALG_MATRIX, 4, "Sending " << buf.write_pos() << " of data to processor " << pid << "\n");
			communicator.send_raw(pid, buf.buffer(), buf.write_pos(), false);
		}

		rowsBufferMap.clear();

		for(IndexLayout::iterator it = receiveLayout.begin(); it != receiveLayout.end(); ++it)
		{
			int pid = receiveLayout.proc_id(it);
			UG_DLOG(LIB_ALG_MATRIX, 4, "issue receive from processor " << pid << "\n");
			communicator.receive_raw(pid, rowsBufferMap[pid]);
		}

		PN.issue(communicator);

	}

	void process(IndexLayout &receiveLayout)
	{
		rowMax=0; colMax=0;
		UG_DLOG(LIB_ALG_MATRIX, 4, "*** RowSendingScheme::process: ***\n");
		connections.clear();

		for(IndexLayout::iterator it = receiveLayout.begin(); it != receiveLayout.end(); ++it)
		{
			int pid = receiveLayout.proc_id(it);
			IndexLayout::Interface interface = receiveLayout.interface(it);

			BinaryBuffer &buf = rowsBufferMap[pid];
			UG_DLOG(LIB_ALG_MATRIX, 4, "rowsBufferMap: received " << buf.write_pos() << " bytes from processor " << pid << "\n");
			for(IndexLayout::Interface::iterator it = interface.begin(); it != interface.end(); ++it)
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
			mat.add_matrix_row(it->first, &cons[0], cons.size());
		}
	}


private:
	void resize_mat(matrix_type &mat)
	{
		size_t cols = std::max(colMax, mat.num_cols());
		size_t rows = std::max(rowMax, mat.num_rows());
		if(rows > mat.num_rows() || cols > mat.num_cols())
			mat.resize(rows, cols);
	}
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
		int j=i;
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
					j++;
				}
			}
		}
		cons.resize(j);
	}
};

}
#endif /* __H__LIB_ALGEBRA__PARALLELIZATION__ROW_SENDING_SCHEME_H_ */
