/*
 * communication_scheme.h
 *
 *  Created on: 13.09.2011
 *      Author: mrupp
 */

#ifndef COMMUNICATION_SCHEME_H_
#define COMMUNICATION_SCHEME_H_

/**
 * CommunicationScheme is a base class for a recurring programming task:
 * You want to send data over a layout/interface. Now you just need to
 * specify what to to at each index on the receiver side, and on each index on the sender side

 * example with own CommunicationScheme
 class AddCommunicationScheme : public CommunicationScheme<AddCommunicationScheme, double>
 {
 public:
	AddCommunicationScheme(const std::vector<double> &v) : vec(v) {}
	double send(int pid, size_t index) const	{ return vec[index]; }
	void receive(int pid, size_t index, double val) { vec[index] += val; }
	inline int get_element_size() const { return sizeof(double); }
private:
	const std::vector<double> &vec;
};

* to use it, you just need to write
 AddCommunicationScheme scheme(vec)
 CommunicateOnInterfaces(parallelCommunicator, masterLayout, slaveLayout, scheme);

* and all slave nodes have been added the master nodes.
* (here slaveLayout is the receiving, and masterLayout is the sending Scheme)


* example with StdArrayCommunicationScheme:
 * std::vector<int> vec;
 * ... do sth with the vector on both processors
 * transfer data from master to slave
 * StdArrayCommunicationScheme<std::vector<int> > arrScheme(vec)
 * CommunicateOnIntefaces(parallelCommunicator, masterLayout, slaveLayout, arrScheme);
 * ... done. Now all master nodes have overwritten the data on the slave nodes.
 * note that instead of int you can use any type for there is Serialize/Deserialize,
 * so also: string, map, int, char, double and so on.
 * there is also a specialization for bool which uses only 1 bit for each bool.

 *
 *
 *
 *
 */

#include "pcl/pcl.h"

#ifdef UG_PARALLEL
namespace ug
{

// CommunicationScheme
//----------------------
/**
 * \brief CRTP Base class for communications on layout/interfaces
 * \tparam TDerived Derived class
 * \tparam TValue
 *
 * The derived class has to have the functions.
 * \sa StdArrayCommunicationScheme
 * const TValue &send(int pid, size_t index) const
 * void receive(int pid, size_t index, value_type &v)
 * inline int get_element_size() const
 *
 */
template<typename TDerived, typename TValue>
class CommunicationScheme
{
public:
	//! get the derived class
	TDerived& derived() {
	       return static_cast<TDerived&>(*this);
	}

	/** send_on_interface
	 * \brief send data on the interface based on TDerived::send
	 * \param	communicator	ParallelCommunicator used to send
	 * \param	pid				where to send data
	 * \param	interface		interface to pid over which to send
	 * This function does the following
	 * - loop through the interface. serialize each item we get from TDerived::send(pid, index) to a buffer
	 * - send the buffer to process pid over ParallelCommunicator com.
	*/
	void send_on_interface(pcl::ParallelCommunicator<IndexLayout> &communicator, int pid, IndexLayout::Interface &interface)
	{
		BinaryBuffer s;
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
			Serialize(s, derived().send(pid, interface.get_element(iter)));

		communicator.send_raw(interface.get_target_proc(), s.buffer(), s.write_pos(), get_size(interface) != -1);
	}


	/** receive_on_interface
	 * \brief receive data on the interface and hand it over to TDerived::receive
	 * \param	s			BinaryBuffer with the data we received from processor pid.
	 * \param	pid			processor we got s from
	 * \param	interface	interface to processor pid
	 * This function does the following
	 * - loop through the interface. deserialize each item from s, and hand it to TDerived::receive(pid, index, item)
	*/
	void receive_on_interface(BinaryBuffer &s, int pid,	IndexLayout::Interface &interface)
	{
		for(IndexLayout::Interface::iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			TValue val;
			Deserialize(s, val);
			derived().receive(pid, interface.get_element(iter), val);
		}
	}

	/** get_size
	 * \brief get the size of the required buffer for an interface
	 * \param	interface	The interface so we can calculate the size
	 * If TDerived can specify an exact size for each item it serializes, we can calculate the size of the buffer
	 * otherwise, it returns -1, and we return -1 -> Buffer size is sent.
	*/
	int get_size(IndexLayout::Interface &interface)
	{
		int s = derived().get_element_size();
		if(s == -1)
			return -1;
		else
			return s * interface.size();
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \brief CRTP Base class for communications on layout/interfaces
 * \tparam TDerived Derived class
 * Template specialization for bool values
 * \sa CommunicationScheme<TDerived, TValue>
 */
template<typename TDerived>
class CommunicationScheme<TDerived, bool>
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
			if(b) a |= (1 << j);
			j++;
			if(j == 8)
			{
				Serialize(s, a);
				a = 0;
				j = 0;
			}
		}
		if(j) Serialize(s, a);
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// StdArrayCommunicationScheme
//------------------------------
/**
 * \brief Communication Scheme for synchronization of arrays or similar datastructures
 * \tparam TArray	array to sync. has to support TArray::value_type operator [] (size_t i) const.
 *
 * example:
 * std::vector<int> vec;
 * ... do sth with the vector on both processors
 * transfer data from master to slave
 * StdArrayCommunicationScheme<std::vector<int> > arrScheme(vec)
 * CommunicateOnIntefaces(parallelCommunicator, slaveLayout, masterLayout, arrScheme);
 * ... done. Now all master nodes have overwritten the data on the slave nodes.
 */
template<typename TArray>
class StdArrayCommunicationScheme : public CommunicationScheme<StdArrayCommunicationScheme<TArray>,
	typename TArray::value_type >
{
	typedef typename TArray::value_type value_type;
public:
	StdArrayCommunicationScheme(TArray &t) : m_arr(t)
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
	TArray &m_arr;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** CommunicateOnInterfaces
 * \brief sends data over a CommunicationScheme from a sendingLayout to a receivingLayout
 * \tparam 	TCommunicationScheme	type of the sending/receiving CommunicationScheme
 * \param	communicator			ParallelCommunicator used to send data
 * \param	sendingLayout			layout to send data from
 * \param	receivingLayout			layout to receive data from
 * \param	scheme					CommunicationScheme to use
 */
template<typename TCommunicationScheme>
void CommunicateOnInterfaces(pcl::ParallelCommunicator<IndexLayout> &communicator,
		IndexLayout &sendingLayout, IndexLayout &receivingLayout, TCommunicationScheme &scheme)
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

/** SendOnInterfaces
 * \brief sends data over a interface based on a CommunicationScheme to a subgroup of processes
 * \tparam 	TSendingScheme	type of the sending CommunicationScheme
 * \tparam 	TPIDs			type of the container of the pids
 * \param	communicator	ParallelCommunicator used to send data
 * \param	pids			pids to send to
 * \param	layout			layout to use to send
 * \param	sender			sending CommunicationScheme
 */
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

/** ReceiveOnInterfaces
 * \brief receives data over a interface based on a CommunicationScheme on a subgroup of processes
 * \tparam 	TReceiveScheme	type of the receiving CommunicationScheme
 * \tparam 	TPIDs			type of the container of the pids
 * \param	communicator	ParallelCommunicator used to send data
 * \param	pids			pids to receive from
 * \param	layout			layout to use to receive
 * \param	receiver		receiving CommunicationScheme
 */
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



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// send sparse data on interface?

} // namespace ug
#endif /* SEND_INTERFACE_H_ */

#endif /* UG_PARALLEL */
