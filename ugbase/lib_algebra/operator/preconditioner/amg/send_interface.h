/*
 * send_interface.h
 *
 *  Created on: 13.09.2011
 *      Author: mrupp
 */

#ifndef SEND_INTERFACE_H_
#define SEND_INTERFACE_H_

namespace ug
{
// send full data on Interface

template<typename T>
class StdArrayCommunicationScheme
{
public:
	typedef typename T::value_type value_type;
	StdArrayCommunicationScheme(T &t) : m_arr(t)
	{

	}

	inline const value_type &send(int pid, size_t index) const
	{
		return m_arr[index];
	}

	inline void receive(int pid, size_t index, value_type &v)
	{
		m_arr[index] = v;
	}

	inline int get_element_size() const
	{
		if(block_traits<value_type>::is_static) return sizeof(value_type);
		else
			return -1;
	}
private:
	T &m_arr;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TSendingScheme, typename TPIDs>
void SendOnInterfaces(pcl::ParallelCommunicator<IndexLayout> &communicator, TPIDs &pids,
		IndexLayout &layout, const TSendingScheme &sender)
{
	AMG_PROFILE_FUNC();

	for(size_t i=0; i<pids.size(); i++)
	{
		int pid = pids[i];
		IndexLayout::Interface &interface = layout.interface(pid);
		SendOnInterface(communicator, pid, interface, sender);
	}
	communicator.communicate();

}

template<typename TPIDs, typename TReceiveScheme>
void ReceiveOnInterfaces(pcl::ParallelCommunicator<IndexLayout> &communicator, TPIDs &pids,
		IndexLayout &layout, TReceiveScheme &receiver)
{
	stdvector< BinaryBuffer > bufs(pids.size());
	for(size_t i=0; i<pids.size(); i++)
		communicator.receive_raw(pids[i], bufs[i], receiver.get_element_size() * layout.interface(pids[i]).size());

	communicator.communicate();

	for(size_t i=0; i<pids.size(); i++)
		ReceiveOnInterface(bufs[i], pids[i], layout.interface(pids[i]), receiver); //, TLocalData::value_type());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TCommunicationScheme>
void CommunicateOnInterfaces(pcl::ParallelCommunicator<IndexLayout> &communicator, IndexLayout &receivingLayout, IndexLayout &sendingLayout, TCommunicationScheme &scheme)
{
	AMG_PROFILE_FUNC();
	int i=0;
	for(IndexLayout::iterator it = sendingLayout.begin(); it != sendingLayout.end(); ++it, ++i)
		SendOnInterface(communicator, sendingLayout.proc_id(it), sendingLayout.interface(it), scheme);

	int receivingLayoutSize=0;
	for(IndexLayout::iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it) receivingLayoutSize++;

	stdvector< BinaryBuffer > bufs(receivingLayoutSize);
	i=0;
	for(IndexLayout::iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it, ++i)
		communicator.receive_raw(receivingLayout.proc_id(it), bufs[i], scheme.get_element_size());

	communicator.communicate();
	i=0;
	for(IndexLayout::iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it, ++i)
		ReceiveOnInterface(bufs[i], receivingLayout.proc_id(it), receivingLayout.interface(it), scheme); //, TLocalData::value_type());


}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TSendingScheme> //, typename TValueType>
void SendOnInterface(pcl::ParallelCommunicator<IndexLayout> &communicator, int pid, IndexLayout::Interface &interface,
		const TSendingScheme &sender) //, TValueType dummy)
{
	BinaryBuffer s;
	for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		Serialize(s, sender.send(pid, interface.get_element(iter))); //, TLocalData::value_type());

	communicator.send_raw(interface.get_target_proc(), s.buffer(), s.write_pos(), sender.get_element_size() != -1);
}


template<typename TReceivingScheme> //, typename TValueType>
void ReceiveOnInterface(BinaryBuffer &s, int pid,
		IndexLayout::Interface &interface, TReceivingScheme &receiver) //, TValueType dummy)
{
	for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
	{
		typename TReceivingScheme::value_type val;
		Deserialize(s, val);
		receiver.receive(pid, interface.get_element(iter), val);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void SendOnInterface(pcl::ParallelCommunicator<IndexLayout> &communicator, int pid, IndexLayout::Interface &interface,
		T &localData, bool dummy)
{
	BinaryBuffer s;
	int j=0;
	char a=0;
	for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
	{
		bool b = localData.send(pid, interface.get_element(iter));
		a &= (b << j);
		if(++j == 8)
		{
			Serialize(s, a);
			a = 0;
		}
	}
	if(j)
		Serialize(s, a);
	communicator.send_raw(interface.get_target_proc(), s.buffer(), s.write_pos(), true);
}
template<typename T>
void ReceiveOnInterface(BinaryBuffer &s, int pid,
		IndexLayout::Interface &interface, T &localData, bool dummy)
{
	int j=8;
	char a=0;
	for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
	{
		size_t index = interface.get_element(iter);
		if(j==8)
		{
			Deserialize(s, a);
			j=0;
		}
		bool b = a & (1 << j);
		j++;
		localData.receive(pid, index, b);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// send sparse data on interface?

} // namespace ug
#endif /* SEND_INTERFACE_H_ */
