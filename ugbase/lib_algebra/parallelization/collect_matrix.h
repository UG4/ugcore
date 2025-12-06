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

#ifndef __H__LIB_ALGEBRA__COLLECT_MATRIX_H_
#define __H__LIB_ALGEBRA__COLLECT_MATRIX_H_

#include "parallel_nodes.h"
#include "serialize_interfaces.h"
#include "common/debug_print.h"

namespace ug{

template<typename matrix_type>
void SerializeRow(BinaryBuffer &stream, const matrix_type &mat, size_t localRowIndex, ParallelNodes &PN)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	const AlgebraID &globalRowIndex = PN.local_to_global(localRowIndex);
	UG_COND_THROW(globalRowIndex.master_proc() > pcl::NumProcs() || globalRowIndex.master_proc() < 0, globalRowIndex);
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
		const AlgebraID &globalColIndex = PN.local_to_global(localColIndex);
		UG_COND_THROW(globalColIndex.master_proc() > pcl::NumProcs() || globalColIndex.master_proc() < 0, globalColIndex);

		UG_DLOG(LIB_ALG_AMG, 4, localColIndex << " (" << globalColIndex << ") -> " << conn.value() << " ");

		// serialize connection
		Serialize(stream, globalColIndex);
		Serialize(stream, conn.value());
	}
	UG_DLOG(LIB_ALG_AMG, 4, "\n");
}


template<typename matrix_type>
void SendMatrix(const matrix_type &A, IndexLayout &verticalSlaveLayout,	int destproc, ParallelNodes &PN)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** SendMatrix ************\n\n");

	pcl::InterfaceCommunicator<IndexLayout> &communicator = A.layouts()->comm();
	BinaryBuffer stream;

	Serialize(stream, A.num_rows());
	for(size_t i=0; i<A.num_rows(); i++)
		SerializeRow(stream, A, i, PN);

	SerializeLayout(stream, A.layouts()->master(), PN);
	SerializeLayout(stream, A.layouts()->slave(), PN);

	IndexLayout::Interface &verticalInterface = verticalSlaveLayout.interface(destproc);
	for(size_t i=0; i<A.num_rows(); i++)
		verticalInterface.push_back(i);

	UG_DLOG(LIB_ALG_AMG, 3, "Srcproc " << pcl::ProcRank() << " is sending " << stream.write_pos() << " bytes of data to destproc " << destproc << "\n");
	communicator.send_raw(destproc, stream.buffer(), stream.write_pos(), false);
	communicator.communicate();
}

template<typename TConnectionType>
size_t DeserializeRow(BinaryBuffer &stream, stdvector<TConnectionType> &cons, ParallelNodes &PN)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	AlgebraID globalRowIndex;

	// serialize global row index
	Deserialize(stream, globalRowIndex);
	size_t localRowIndex = PN.global_to_local(globalRowIndex);

	UG_DLOG(LIB_ALG_AMG, 4, "Got row " << localRowIndex << " (" << globalRowIndex << "), ");
	size_t num_connections;

	// serialize number of connections
	Deserialize(stream, num_connections);

	UG_DLOG(LIB_ALG_AMG, 4, num_connections << " connections: ");

	cons.resize(num_connections);
	for(size_t i =0; i<num_connections; i++)
	{
		AlgebraID globalColIndex;
		Deserialize(stream, globalColIndex);
		cons[i].iIndex = PN.global_to_local(globalColIndex);
		Deserialize(stream, cons[i].dValue);
		UG_DLOG(LIB_ALG_AMG, 4, cons[i].iIndex << " (" << globalColIndex << ") -> " << cons[i].dValue << " ");
	}
	UG_DLOG(LIB_ALG_AMG, 4, "\n");
	return localRowIndex;
}

// ReceiveMatrix
//---------------------------------------------------------------------------
/**
 *	Receives a distributed matrix from several processors
 * \param A				(in) input matrix
 * \param M				(out) collected matrix
 * \param verticalMasterLayout	(out) created master layout to processors in srcprocs
 * \param
 * \param srcprocs		list of source processors
 *
 */
template<typename matrix_type>
void ReceiveMatrix(const matrix_type &A, matrix_type &M, IndexLayout &verticalMasterLayout,	const std::vector<int> &srcprocs,
		ParallelNodes &PN)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** ReceiveMatrix ************\n\n");
	pcl::InterfaceCommunicator<IndexLayout> &communicator = A.layouts()->comm();

	M = A;
	//M.print();
	M.set_layouts(SmartPtr(new AlgebraLayouts));
	using BufferMap = std::map<int, BinaryBuffer>;
	BufferMap streams;

	UG_DLOG(LIB_ALG_AMG, 3, "DestProc " << pcl::ProcRank() << " is waiting on data from ");
	for(size_t i=0; i<srcprocs.size(); i++)
	{
		UG_DLOG(LIB_ALG_AMG, 3, srcprocs[i] << " ");
		communicator.receive_raw(srcprocs[i], streams[srcprocs[i]]);
	}
	UG_DLOG(LIB_ALG_AMG, 3, "\n");
	communicator.communicate();

	AlgebraID globalRowIndex, globalColIndex;
	size_t num_connections, numRows;

	for(size_t i=0; i<srcprocs.size(); i++)
	{
		int pid = srcprocs[i];
		BinaryBuffer &stream = streams[pid];
		stream.set_read_pos(0);

		UG_DLOG(LIB_ALG_AMG, 4, "received " << stream.write_pos() << " bytes of data from process " << pid << "\n");
		IndexLayout::Interface &verticalInterface = verticalMasterLayout.interface(pid);
		typename matrix_type::connection con;

		Deserialize(stream, numRows);
		for(size_t i=0; i<numRows; i++)
		{
			// serialize global row index, number of connections
			Deserialize(stream, globalRowIndex);
			Deserialize(stream, num_connections);
			UG_COND_THROW(globalRowIndex.master_proc() > pcl::NumProcs() || globalRowIndex.master_proc() < 0, i << " " << globalRowIndex << " " << pid);

			size_t localRowIndex = PN.get_local_index_or_create_new(globalRowIndex, 0);
			verticalInterface.push_back(localRowIndex);
			UG_DLOG(LIB_ALG_AMG, 4, "Got row " << localRowIndex << " (" << globalRowIndex << "), ");
			UG_DLOG(LIB_ALG_AMG, 4, num_connections << " connections: ");

			for(size_t pid =0; pid<num_connections; pid++)
			{
				Deserialize(stream, globalColIndex);
				Deserialize(stream, con.dValue);

				con.iIndex = PN.get_local_index_or_create_new(globalColIndex, 0);
				UG_DLOG(LIB_ALG_AMG, 4, con.iIndex << " (" << globalColIndex << ") -> " << con.dValue << " ");
			}
			UG_DLOG(LIB_ALG_AMG, 4, "\n");
		}
	}

	M.resize_and_keep_values(PN.local_size(), PN.local_size());
	//M.print();

	for(size_t i=0; i<srcprocs.size(); i++)
	{
		int pid = srcprocs[i];
		BinaryBuffer &stream = streams[pid];
		stream.set_read_pos(0);
		stdvector<typename matrix_type::connection> cons;

		Deserialize(stream, numRows);
		for(size_t i=0; i<numRows; i++)
		{
			size_t localRowIndex = DeserializeRow(stream, cons, PN);
			if(cons.size())
				M.add_matrix_row(localRowIndex, &cons[0], cons.size());
		}
	}

	//UG_DLOG(LIB_ALG_AMG, 4, "\n** the matrix M: \n\n");
	//M.print();
	//UG_DLOG(LIB_ALG_AMG, 4, "\n");

	//UG_LOG("COLLECTED LAYOUT:\n");
	//PrintLayout(processCommunicator, communicator, masterLayout, slaveLayout);
}

/**
 * 1. constructs global indices
 * 2. for pid != proc_id(0) :
 * 		a) send the whole matrix with global ids to proc_id(0)
 * 		b) create slaveLayout to proc_id(0)
 *
 * 3. for pid = proc_id(0) :
 * 		a) receives the matrices
 * 		b) builds up collectedA
 * 		c) creates masterLayout to all other pids.
 *
 * @param A				(input) the distributed parallel matrix A
 * @param collectedA	(output) the collected matrix A on pid = proc_id(0)
 * @param masterLayout	the agglomeration master layout (only defined on pid = proc_id(0))
 * @param slaveLayout	the agglomeration slave layout (only defined on pid != proc_id(0))
 */
template<typename matrix_type>
void CollectMatrixOnOneProc(const matrix_type &A, matrix_type &collectedA, IndexLayout &masterLayout, IndexLayout &slaveLayout)
{
	try{
	PROFILE_FUNC_GROUP("algebra parallelization");
	UG_DLOG(LIB_ALG_AMG, 1, "\n*********** SendMatrix ************\n\n");
	std::vector<int> srcprocs;
	masterLayout.clear();
	slaveLayout.clear();

	const pcl::ProcessCommunicator &pc = A.layouts()->proc_comm();
	ParallelNodes PN(A.layouts(), A.num_rows());

	if(pcl::ProcRank() == pc.get_proc_id(0))
	{
		srcprocs.resize(pc.size()-1);
		for(size_t i=1; i<pc.size(); i++)
			srcprocs[i-1] = pc.get_proc_id(i);
		ReceiveMatrix(A, collectedA, masterLayout, srcprocs, PN);
	}
	else
		SendMatrix(A, slaveLayout, pc.get_proc_id(0), PN);
	}UG_CATCH_THROW(__FUNCTION__ << " failed");
}

/**
 * gathers the vector vec to collectedVec on one processor
 * @param agglomeratedMaster	master agglomeration layout. only nonempty if Root=true
 * @param agglomeratedSlave		slave agglomeration layout. only nonempty if Root=false
 * @param com
 * @param collectedVec			(output) result on proc with bRoot=true
 * @param vec					(input) the distributed vec
 * @param type can be PST_ADDITIVE or PST_CONSISTENT
 * @param bRoot
 */
template<typename T>
void GatherVectorOnOne(IndexLayout &agglomeratedMaster, IndexLayout &agglomeratedSlave,
		pcl::InterfaceCommunicator<IndexLayout> &com,
		ParallelVector<T> &collectedVec,
		const ParallelVector<T> &vec,
		ParallelStorageType type, bool bRoot)
{
	try{
	PROFILE_FUNC_GROUP("algebra parallelization");
	using vector_type = ParallelVector<T>;
	if(!bRoot)
	{
		ComPol_VecAdd<vector_type > compolAdd(&collectedVec, &vec);
		com.send_data(agglomeratedSlave, compolAdd);
		com.communicate();
	}
	else
	{
		//UG_LOG("gather_vertical: receiving data at level " << level << "\n");
		UG_COND_THROW(&vec == &collectedVec, "vec and collected vec may not be same");
		collectedVec.set(0.0);
		for(size_t i=0; i<vec.size(); i++)
			collectedVec[i] = vec[i];

		UG_COND_THROW(!vec.has_storage_type(type), "storage type is " << vec.get_storage_type() << ", not " << type);
		if(type == PST_ADDITIVE)
		{
			ComPol_VecAdd<vector_type > compolAdd(&collectedVec, &vec);
			com.receive_data(agglomeratedMaster, compolAdd);
			com.communicate();
			collectedVec.set_storage_type(PST_ADDITIVE);
		}
		else if(type == PST_CONSISTENT)
		{
			ComPol_VecCopy<vector_type > compolCopy(&collectedVec, &vec);
			com.receive_data(agglomeratedMaster, compolCopy);
			com.communicate();
			collectedVec.set_storage_type(PST_CONSISTENT);
		}
		else { UG_THROW("storage type " << type << "unsupported."); }
	}
	}UG_CATCH_THROW(__FUNCTION__ << " failed");
}

/**
 * broadcasts the vector collectedVec to the distributed vec
 * @param agglomeratedMaster	master agglomeration layout. only nonempty if bRoot=true
 * @param agglomeratedSlave		slave agglomeration layout. only nonempty if bRoot=false
 * @param com
 * @param vec					(output) the distributed vec
 * @param collectedVec			(input) collectedVec
 * @param type can be PST_ADDITIVE or PST_CONSISTENT
 * @param bRoot
 */
template<typename T>
void BroadcastVectorFromOne(IndexLayout &agglomeratedMaster, IndexLayout &agglomeratedSlave,
		pcl::InterfaceCommunicator<IndexLayout> &com,
		ParallelVector<T> &vec,
		const ParallelVector<T> &collectedVec,
		ParallelStorageType type, bool bRoot)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	try{
		using vector_type = ParallelVector<T>;
	if(!bRoot)
	{
		ComPol_VecCopy<vector_type> compolCopy(&vec, &collectedVec);
		com.receive_data(agglomeratedSlave, compolCopy);
		com.communicate();
		vec.set_storage_type(type);
	}
	else
	{
		UG_COND_THROW(&vec == &collectedVec, "vec and collected vec may not be same");
		for(size_t i=0; i<vec.size(); i++)
			vec[i] = collectedVec[i];

		UG_COND_THROW(!collectedVec.has_storage_type(type), "storage type is " << collectedVec.get_storage_type() << ", not " << type);
		vec.set_storage_type(type);

		ComPol_VecAdd<vector_type > compolCopy(&vec, &collectedVec);
		com.send_data(agglomeratedMaster, compolCopy);
		com.communicate();
	}

	if(type == PST_ADDITIVE)
	{
		UG_THROW("ONLY CONSISTENT!");
		// das problem ist, dass der vektor noch slave-interfaces nach "außen" haben kann,
		// diese werden dann fälschlicherweise auch 0 gesetzt.

		//!!! WRONG !!! SetLayoutValues(&vec, vec.layouts()->slave(), 0.0); //!!!
		//vec.set_storage_type(PST_ADDITIVE);

	}
	else if(type == PST_CONSISTENT) {	}
	else { UG_THROW("storage type " << type << "unsupported."); }

	}UG_CATCH_THROW(__FUNCTION__ << " failed");
}

} // namespace ug

#endif