/*
 * serialization.h
 *
 *  Created on: May 29, 2012
 *      Author: rupp4
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
#endif /* SERIALIZATION_H_ */
