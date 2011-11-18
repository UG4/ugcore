/*
 * send_interface.h
 *
 *  Created on: 13.09.2011
 *      Author: mrupp
 */

#ifndef SEND_INTERFACE_H_
#define SEND_INTERFACE_H_

#ifdef UG_PARALLEL
namespace ug
{

template<typename TDerived>
class BoolCommunicationScheme
{
public:
	TDerived& derived() { return static_cast<TDerived&>(*this); }

	void send_on_interface(pcl::ParallelCommunicator<IndexLayout> &communicator, int pid, IndexLayout::Interface &interface)
	{
		BinaryBuffer s;
		int j=0;
		char a=0;
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			bool b = derived().send(pid, interface.get_element(iter));
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
	void receive_on_interface(BinaryBuffer &s, int pid, IndexLayout::Interface &interface)
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
			derived().receive(pid, index, b);
		}
	}

	int get_size(IndexLayout::Interface &interface)
	{
		return sizeof(char) * (interface.size()/8 + (interface.size()%8==0?0:1) );
	}

};

template<typename TDerived>
class CommunicationScheme
{
public:
	TDerived& derived() {
	       return static_cast<TDerived&>(*this);
	}
	void send_on_interface(pcl::ParallelCommunicator<IndexLayout> &communicator, int pid, IndexLayout::Interface &interface)
	{
		BinaryBuffer s;
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			Serialize(s, derived().send(pid, interface.get_element(iter)));

		communicator.send_raw(interface.get_target_proc(), s.buffer(), s.write_pos(), get_size(interface) != -1);
	}


	void receive_on_interface(BinaryBuffer &s, int pid,	IndexLayout::Interface &interface)
	{
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			typename TDerived::value_type val;
			Deserialize(s, val);
			derived().receive(pid, interface.get_element(iter), val);
		}
	}

	int get_size(IndexLayout::Interface &interface)
	{
		int s = derived().get_element_size();
		if(s == -1)
			return -1;
		else
			return s * interface.size();
	}
};

template<typename T>
class StdArrayCommunicationScheme : public CommunicationScheme<StdArrayCommunicationScheme<T> >
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

template<typename T>
class BoolArrayCommunicationScheme : public BoolCommunicationScheme<BoolArrayCommunicationScheme<T> >
{
public:
	BoolArrayCommunicationScheme(T &t) : m_arr(t)
	{

	}

	inline bool send(int pid, size_t index) const
	{
		return m_arr[index];
	}

	inline void receive(int pid, size_t index, bool v)
	{
		m_arr[index] = v;
	}

private:
	T &m_arr;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TSendingScheme, typename TPIDs>
void SendOnInterfaces(pcl::ParallelCommunicator<IndexLayout> &communicator, TPIDs &pids,
		IndexLayout &layout, TSendingScheme &sender)
{
	AMG_PROFILE_FUNC();

	for(size_t i=0; i<pids.size(); i++)
	{
		int pid = pids[i];
		IndexLayout::Interface &interface = layout.interface(pid);
		sender.send_on_interface(communicator, pid, interface);
	}
	communicator.communicate();

}

template<typename TPIDs, typename TReceiveScheme>
void ReceiveOnInterfaces(pcl::ParallelCommunicator<IndexLayout> &communicator, TPIDs &pids,
		IndexLayout &layout, TReceiveScheme &receiver)
{
	stdvector< BinaryBuffer > bufs(pids.size());
	for(size_t i=0; i<pids.size(); i++)
		communicator.receive_raw(pids[i], bufs[i], receiver.get_size(layout.interface(pids[i])));

	communicator.communicate();

	for(size_t i=0; i<pids.size(); i++)
		receiver.receive_on_interface(bufs[i], pids[i], layout.interface(pids[i]));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename TCommunicationScheme>
void CommunicateOnInterfaces(pcl::ParallelCommunicator<IndexLayout> &communicator, IndexLayout &receivingLayout, IndexLayout &sendingLayout, TCommunicationScheme &scheme)
{
	AMG_PROFILE_FUNC();
	int i=0;
	for(IndexLayout::iterator it = sendingLayout.begin(); it != sendingLayout.end(); ++it, ++i)
		scheme.send_on_interface(communicator, sendingLayout.proc_id(it), sendingLayout.interface(it));

	int receivingLayoutSize=0;
	for(IndexLayout::iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it) receivingLayoutSize++;

	stdvector< BinaryBuffer > bufs(receivingLayoutSize);
	i=0;
	for(IndexLayout::iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it, ++i)
		communicator.receive_raw(receivingLayout.proc_id(it), bufs[i], scheme.get_size(receivingLayout.interface(it)));

	communicator.communicate();
	i=0;
	for(IndexLayout::iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it, ++i)
		scheme.receive_on_interface(bufs[i], receivingLayout.proc_id(it), receivingLayout.interface(it));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// send sparse data on interface?

} // namespace ug
#endif /* SEND_INTERFACE_H_ */

#endif /* UG_PARALLEL */
