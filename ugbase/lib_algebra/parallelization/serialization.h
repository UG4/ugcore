/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef SERIALIZATION_H_
#define SERIALIZATION_H_

namespace ug{

template<class TOStream>
void Serialize(TOStream &buf, const IndexLayout::Interface &interface)
{
	Serialize(buf, (size_t)interface.size());
	for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
		Serialize(buf, interface.get_element(iter));
}

template<class TIStream>
void Deserialize(TIStream &buf, IndexLayout::Interface &interface)
{
	size_t s = Deserialize<size_t> (buf);
	IndexLayout::Interface::Element el;
	for(size_t i=0; i<s; i++)
	{
		Deserialize(buf, el);
		interface.push_back(el);
	}
}


template<class TOStream>
void Serialize(TOStream &buf, const IndexLayout &layout)
{
	Serialize(buf, (size_t)layout.num_interfaces());
	for(IndexLayout::const_iterator iter = layout.begin(); iter != layout.end(); ++iter)
	{
		Serialize(buf, (int) layout.proc_id(iter));
		Serialize(buf, layout.interface(iter));
	}
}


template<class TIStream>
void Deserialize(TIStream &buf, IndexLayout &layout)
{
	size_t num_interfaces = Deserialize<size_t>(buf);
	for(size_t i=0; i<num_interfaces; i++)
	{
		int pid = Deserialize<int>(buf);
		Deserialize(buf, layout.interface(pid));
	}
}

template<typename T>
void SetParallelData(T &t,
		IndexLayout &masterLayout, IndexLayout &slaveLayout,
		pcl::InterfaceCommunicator<IndexLayout> &ic,
		pcl::ProcessCommunicator& pc)
{
	t.set_master_layout(masterLayout);
	t.set_slave_layout(slaveLayout);

	t.set_communicator(ic);
	t.set_process_communicator(pc);
}

template<typename T, class TOStream>
void SerializeParallelData(TOStream &buf, T &t)
{
	Serialize(buf, t.layouts()->master());
	Serialize(buf, t.layouts()->slave());

	Serialize(buf, t.layouts()->comm());

	Serialize(buf, t.layouts()->proc_comm());
}

template<typename T, class TIStream>
void DeserializeParallelData(TIStream &buf, T &t,
		IndexLayout &masterLayout, IndexLayout &slaveLayout,
		pcl::InterfaceCommunicator<IndexLayout> &ic,
		pcl::ProcessCommunicator& pc)
{
	Deserialize(buf, masterLayout);
	Deserialize(buf, slaveLayout);

	Deserialize(buf, ic);
	Deserialize(buf, pc);
}

///////////////////////

template<typename T, class TOStream>
void SerializeUniquePart(TOStream &buf, const ParallelMatrix<T> &A)
{
	Serialize(buf, A.get_storage_mask());
}

template<typename T, class TIStream>
void DeserializeUniquePart(TIStream &buf, ParallelMatrix<T> &A)
{
	A.set_storage_type( Deserialize<uint>(buf) );
}

/*template<typename T, class TOStream>
void Serialize(TOStream &buf, const ParallelMatrix<T> &A)
{
	UG_ASSERT(0, "use SerializeUniquePart");
}

template<typename T, class TIStream>
void Deserialize(TIStream &buf, ParallelMatrix<T> &A)
{
	UG_ASSERT(0, "use DeserializeUniquePart");
}*/

//////////////////
template<typename TValueType, class TOStream>
void SerializeUniquePart(TOStream &buf, const ParallelVector<TValueType> &v)
{
	Serialize(buf, v.get_storage_mask());
}

template<typename TValueType, class TIStream>
void DeserializeUniquePart(TIStream &buf, ParallelVector<TValueType> &v)
{
	v.set_storage_type( Deserialize<uint>(buf) );
}

/*template<typename TValueType, class TOStream>
void Serialize(TOStream &buf, const ParallelVector<TValueType> &v)
{
	UG_ASSERT(0, "use SerializeUniquePart and SerializeParallelData");
}

template<typename TValueType, class TIStream>
void Deserialize(TIStream &buf, ParallelVector<TValueType> &v)
{
	UG_ASSERT(0, "use DeserializeUniquePart and DeserializeParallelData");
}*/
/////////////////

template<class TOStream>
void Serialize(TOStream &buf, const pcl::InterfaceCommunicator<IndexLayout> &ic)
{
	//Serialize<bool>(buf, ic.communication_debugging_enabled());
}

template<class TIStream>
void Deserialize(TIStream &buf, pcl::InterfaceCommunicator<IndexLayout> &ic)
{
	//
}

template<class TOStream>
void Serialize(TOStream &buf, const pcl::ProcessCommunicator &ic)
{
	int s;
	if(ic.is_world())
		s = -1;
	else
		s = ic.size();
	Serialize(buf, s);

	if(s > 0)
	{
		for(int i=0; i<s; i++)
			Serialize(buf, ic.get_proc_id(i));
	}
}


template<class TIStream>
void Deserialize(TIStream &buf, pcl::ProcessCommunicator &ic)
{
	int s = Deserialize<int>(buf);

	if(s == -1)
	{
		if(ic.is_world()) return;
		ic = pcl::ProcessCommunicator(pcl::PCD_WORLD);
	}
	else if (s == 0)
	{
		if(ic.empty()) return;
		ic = pcl::ProcessCommunicator(pcl::PCD_EMPTY);
	}
	else
	{
		std::vector<int> procs;
		for(int i=0; i<s; i++)
			procs.push_back(Deserialize<int>(buf));
		ic = pcl::ProcessCommunicator::create_communicator(procs);
	}
}

}
#endif