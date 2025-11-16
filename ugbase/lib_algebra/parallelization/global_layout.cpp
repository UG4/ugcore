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

#ifndef UG_PARALLEL
#error "This only works with a UG_PARALLEL define."
#endif

#include "global_layout.h"

namespace ug
{

void DeserializeAndAddGlobalInterface(BinaryBuffer &stream,	std::vector<AlgebraID> &interface)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	std::vector<AlgebraID> tmp;
	Deserialize(stream, tmp);
	interface.insert(interface.end(), tmp.begin(), tmp.end());
}

void DeserializeAndAddGlobalLayout(BinaryBuffer &stream, GlobalLayout &globalLayout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	size_t size;
	Deserialize(stream, size);
	for(size_t i=0; i<size; i++)
	{
		int pid = Deserialize<int>(stream);
		DeserializeAndAddGlobalInterface(stream, globalLayout[pid]);
	}
}

void PrintGlobalLayout(const GlobalLayout &globalLayout, const char *name)
{
	PROFILE_FUNC_GROUP("algebra parallelization debug");
	UG_LOG("GlobalLayout ");
	if(name) UG_DLOG(LIB_ALG_AMG, 4, name);
	UG_LOG("\n");

	for(GlobalLayout::const_iterator iter = globalLayout.begin(); iter != globalLayout.end(); ++iter)
	{
		int pid = iter->first;
		const std::vector<AlgebraID> &v = iter->second;
		UG_LOG("to processor " << pid << "\n");
		for(size_t i=0; i<v.size(); i++)
		{
			UG_LOG(" " << v[i] << "\n");
		}
	}
}

void ReceiveGlobalLayout(pcl::InterfaceCommunicator<IndexLayout> &comm, const std::vector<int> &srcprocs,
		GlobalLayout &globalMasterLayout, GlobalLayout &globalSlaveLayout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	using BufferMap = std::map<int, BinaryBuffer>;
	BufferMap streams;

	for(size_t i=0; i<srcprocs.size(); i++)
		comm.receive_raw(srcprocs[i], streams[srcprocs[i]]);

	comm.communicate();
	for(size_t i=0; i<streams.size(); i++)
	{
		BinaryBuffer &stream = streams[srcprocs[i]];
		UG_DLOG(LIB_ALG_AMG, 4, "received " << stream.write_pos() << " bytes from " << srcprocs[i] << "\n");
		DeserializeAndAddGlobalLayout(stream, globalMasterLayout);
		DeserializeAndAddGlobalLayout(stream, globalSlaveLayout);
	}
}

void SerializeGlobalLayout(BinaryBuffer &stream, const GlobalLayout &globalLayout)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	size_t size = globalLayout.size();
	Serialize(stream, size);
	UG_DLOG(LIB_ALG_AMG, 4, " size = " << size << "\n");
	for(GlobalLayout::const_iterator iter = globalLayout.begin(); iter != globalLayout.end(); ++iter)
	{
		int pid = iter->first;
		const std::vector<AlgebraID> &v = iter->second;
		//UG_DLOG(LIB_ALG_AMG, 4, " " << pid << " - " << v << "\n");
		//if(v.size() == 0) continue;
		Serialize(stream, pid);
		Serialize(stream, v);
	}
}

void SendGlobalLayout(pcl::InterfaceCommunicator<IndexLayout> &comm,
		const GlobalLayout &globalMasterLayout, const GlobalLayout &globalSlaveLayout, int pid)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	UG_DLOG(LIB_ALG_AMG, 4, "sending to " << pid << "\n");
	BinaryBuffer stream;
	UG_DLOG(LIB_ALG_AMG, 4, "serialize globalMasterLayout\n");
	SerializeGlobalLayout(stream, globalMasterLayout);
	UG_DLOG(LIB_ALG_AMG, 4, "serialize globalSlaveLayout\n");
	SerializeGlobalLayout(stream, globalSlaveLayout);

	UG_DLOG(LIB_ALG_AMG, 4, "size is " << stream.write_pos() << "\n");
	comm.send_raw(pid, stream.buffer(), stream.write_pos(), false);
	comm.communicate();
}


void MergeGlobalLayout(GlobalLayout &globalLayout, const std::map<int, int> &merge)
{
	PROFILE_FUNC_GROUP("algebra parallelization");
	for(std::map<int, int>::const_iterator it = merge.begin(); it != merge.end(); ++it)
	{
		std::vector<AlgebraID> &a = globalLayout[it->first];

		if(it->second != pcl::ProcRank())
		{
			std::vector<AlgebraID> &b = globalLayout[it->second];
			b.insert(b.end(), a.begin(), a.end());
		}
		globalLayout.erase(it->first);
	}
	// we don't want interfaces to ourselfs
	globalLayout.erase(pcl::ProcRank());
}

}
