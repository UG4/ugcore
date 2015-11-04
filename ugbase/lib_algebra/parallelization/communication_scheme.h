
#ifndef COMMUNICATION_SCHEME_H_
#define COMMUNICATION_SCHEME_H_

#include "pcl/pcl.h"

#ifdef UG_PARALLEL
namespace ug
{

/**
 * \defgroup lib_algebra_parallelization_scheme Parallel Algebra Communication Scheme
 * \brief Communication Scheme for parallel Algebra
 * \ingroup lib_algebra_parallelization
 * \{
 */

// CommunicationScheme
//----------------------
/**
 * \brief CRTP Base class for communications on layout/interfaces
 *
 * \details CommunicationScheme is a base class for a recurring programming task:
 *   You want to send data over a layout/interface. Now you just need to
 *   specify what to to at each index on the receiver side, and on each index 
 *   on the sender side
 *   Example with own CommunicationScheme
 * \code{.cpp}
 * class AddCommunicationScheme : public CommunicationScheme<AddCommunicationScheme, double>
 * {
 * public:
 * 	AddCommunicationScheme(const std::vector<double> &v) : vec(v) {}
 * 	double send(int pid, size_t index) const	{ return vec[index]; }
 * 	void receive(int pid, size_t index, double val) { vec[index] += val; }
 * 	inline int get_element_size() const { return sizeof(double); }
 * private:
 * 	const std::vector<double> &vec;
 * };
 * \endcode
 * to use it, you just need to write
 * \code{.cpp}
 * AddCommunicationScheme scheme(vec);
 * CommunicateOnInterfaces(parallelCommunicator, masterLayout, slaveLayout, scheme);
 * \endcode
 * and all slave nodes have been added the master nodes.
 * (here slaveLayout is the receiving, and masterLayout is the sending Scheme)
 * 
 * Example with <tt>StdArrayCommunicationScheme</tt>:
 * \code{.cpp}
 * std::vector<int> vec;
 * // ... do sth with the vector on both processors
 * // transfer data from master to slave
 * StdArrayCommunicationScheme<std::vector<int> > arrScheme(vec)
 * CommunicateOnIntefaces(parallelCommunicator, masterLayout, slaveLayout, arrScheme);
 * // ... done.
 * \endcode
 * Now all master nodes have overwritten the data on the slave nodes.
 * note that instead of int you can use any type for there is Serialize/Deserialize,
 * so also: \c string, \c map, \c int, \c char, \c double and so on.
 * There is also a specialization for \c bool which uses only 1 bit for each 
 * \c bool.
 * 
 * \note The derived class has to have the functions:
 *  - <tt>%StdArrayCommunicationScheme</tt>
 *  - <tt>const TValue &send(int pid, size_t index) const</tt>
 *  - <tt>void receive(int pid, size_t index, value_type &v)</tt>
 *  - <tt>inline int get_element_size() const</tt>
 *
 * \tparam TDerived Derived class
 * \tparam TValue
 */
template<typename TDerived, typename TValue>
class CommunicationScheme : public pcl::ICommunicationPolicy<IndexLayout>
{
public:
	//! get the derived class
	TDerived& derived() {
	       return static_cast<TDerived&>(*this);
	}

	/** collect
	 * \brief send data on the interface based on TDerived::send
	 * \param	buff			BinaryBuffer to write data to
	 * \param	interface		interface to pid over which to send
	 * This function does the following
	 * - loop through the interface. serialize each item we get from TDerived::send(pid, index) to a buffer
	 * - send the buffer to process pid over InterfaceCommunicator com.
	*/
	virtual bool
	collect(ug::BinaryBuffer& buff, const Interface& interface)

	{
		int pid = interface.get_target_proc();
		for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
			Serialize(buff, derived().send(pid, interface.get_element(iter)));
		return true;
	}


	/** extract
	 * \brief receive data on the interface and hand it over to TDerived::receive
	 * \param	buff		BinaryBuffer with the data we received from processor pid.
	 * \param	interface	interface to processor pid
	 * This function does the following	 *
	 * - loop through the interface. deserialize each item from s, and hand it to TDerived::receive(pid, index, item)
	*/
	virtual bool
	extract(ug::BinaryBuffer& buff, const Interface& interface)
	{
		int pid = interface.get_target_proc();
		TValue val;
		for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			Deserialize(buff, val);
			derived().receive(pid, interface.get_element(iter), val);
		}
		return true;
	}

	/** get_required_buffer_size
	 * \brief get the size of the required buffer for an interface
	 * \param	interface	The interface so we can calculate the size
	 * If TDerived can specify an exact size for each item it serializes, we can calculate the size of the buffer
	 * otherwise, it returns -1, and we return -1 -> Buffer size is sent.
	*/
	virtual int
	get_required_buffer_size(const Interface& interface)
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
class CommunicationScheme<TDerived, bool> : public pcl::ICommunicationPolicy<IndexLayout>
{
public:
	TDerived& derived() { return static_cast<TDerived&>(*this); }

	virtual bool
	collect(ug::BinaryBuffer& buff, const Interface& interface)
	{
		int pid = interface.get_target_proc();
		int j=0;
		char a=0;
		for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			bool b = derived().send(pid, interface.get_element(iter));
			if(b) a |= (1 << j);
			j++;
			if(j == 8)
			{
				Serialize(buff, a);
				a = 0;
				j = 0;
			}
		}
		if(j) Serialize(buff, a);
		return true;
	}

	virtual bool
	extract(ug::BinaryBuffer& buff, const Interface& interface)
	{
		int pid = interface.get_target_proc();
		int j=8;
		char a=0;
		for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			size_t index = interface.get_element(iter);
			if(j==8)
			{
				Deserialize(buff, a);
				j=0;
			}
			bool b = a & (1 << j);
			j++;
			derived().receive(pid, index, b);
		}
		return true;
	}

	virtual int
	get_required_buffer_size(const Interface& interface)
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
 * \param	communicator			InterfaceCommunicator used to send data
 * \param	sendingLayout			layout to send data from
 * \param	receivingLayout			layout to receive data from
 * \param	scheme					CommunicationScheme to use
 */
template<typename TCommunicationScheme>
void CommunicateOnInterfaces(pcl::InterfaceCommunicator<IndexLayout> &communicator,
		const IndexLayout &sendingLayout, const IndexLayout &receivingLayout, TCommunicationScheme &scheme)
{
	AMG_PROFILE_FUNC();
	int i=0;

	for(IndexLayout::const_iterator it = sendingLayout.begin(); it != sendingLayout.end(); ++it, ++i)
	{
		BinaryBuffer buff;
		const IndexLayout::Interface &interface = sendingLayout.interface(it);
		scheme.collect(buff, interface);
		// todo: don't use parallelcommunicator send_raw
		// todo: reserve and reuse buff
		communicator.send_raw(interface.get_target_proc(), buff.buffer(), buff.write_pos(),
				scheme.get_required_buffer_size(interface) != -1);
	}

	int receivingLayoutSize=0;
	for(IndexLayout::const_iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it)
		receivingLayoutSize++;

	stdvector< BinaryBuffer > bufs(receivingLayoutSize);
	i=0;
	for(IndexLayout::const_iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it, ++i)
		communicator.receive_raw(receivingLayout.proc_id(it), bufs[i],
				scheme.get_required_buffer_size(receivingLayout.interface(it)));

	communicator.communicate();
	i=0;
	for(IndexLayout::const_iterator it = receivingLayout.begin(); it != receivingLayout.end(); ++it, ++i)
		scheme.extract(bufs[i], receivingLayout.interface(it));
}

template<typename TCommunicationScheme>
void CommunicateFromSlaveToMaster(HorizontalAlgebraLayouts &layouts, TCommunicationScheme &scheme)
{
	CommunicateOnInterfaces(layouts.proc_comm(), layouts.slave(), layouts.master(), scheme);
}
template<typename TCommunicationScheme>
void CommunicateFromMasterToSlave(HorizontalAlgebraLayouts &layouts, TCommunicationScheme &scheme)
{
	CommunicateOnInterfaces(layouts.proc_comm(), layouts.master(), layouts.slave(), scheme);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** SendOnInterfaces
 * \brief sends data over a interface based on a CommunicationScheme to a subgroup of processes
 * \tparam 	TSendingScheme	type of the sending CommunicationScheme
 * \tparam 	TPIDs			type of the container of the pids
 * \param	communicator	InterfaceCommunicator used to send data
 * \param	pids			pids to send to
 * \param	layout			layout to use to send
 * \param	sender			sending CommunicationScheme
 */
template<typename TSendingScheme, typename TPIDs>
void SendOnInterfaces(pcl::InterfaceCommunicator<IndexLayout> &communicator, TPIDs &pids,
		IndexLayout &layout, TSendingScheme &sender)
{
	AMG_PROFILE_FUNC();

	for(size_t i=0; i<pids.size(); i++)
	{
		int pid = pids[i];
		IndexLayout::Interface &interface = layout.interface(pid);

		BinaryBuffer buff;
		sender.collect(buff, interface);
		// todo: don't use parallelcommunicator send_raw
		// todo: reserve and reuse buff
		communicator.send_raw(pid, buff.buffer(), buff.write_pos(),
				sender.get_required_buffer_size(interface) != -1);
	}
	communicator.communicate();

}

/** ReceiveOnInterfaces
 * \brief receives data over a interface based on a CommunicationScheme on a subgroup of processes
 * \tparam 	TReceiveScheme	type of the receiving CommunicationScheme
 * \tparam 	TPIDs			type of the container of the pids
 * \param	communicator	InterfaceCommunicator used to send data
 * \param	pids			pids to receive from
 * \param	layout			layout to use to receive
 * \param	receiver		receiving CommunicationScheme
 */
template<typename TPIDs, typename TReceiveScheme>
void ReceiveOnInterfaces(pcl::InterfaceCommunicator<IndexLayout> &communicator, TPIDs &pids,
		IndexLayout &layout, TReceiveScheme &receiver)
{
	stdvector< BinaryBuffer > bufs(pids.size());
	for(size_t i=0; i<pids.size(); i++)
		communicator.receive_raw(pids[i], bufs[i], receiver.get_required_buffer_size(layout.interface(pids[i])));

	communicator.communicate();

	for(size_t i=0; i<pids.size(); i++)
		receiver.extract(bufs[i], layout.interface(pids[i]));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0
template<typename TDerived, typename TValue>
class SparseCommunicationScheme : public pcl::ICommunicationPolicy<IndexLayout>
{
public:
	//! get the derived class
	TDerived& derived() {
	       return static_cast<TDerived&>(*this);
	}

	/** collect
	 * \brief send data on the interface based on TDerived::send
	 * \param	buff			BinaryBuffer to write data to
	 * \param	interface		interface to pid over which to send
	 * This function does the following
	 * - loop through the interface. serialize each item we get from TDerived::send(pid, index) to a buffer
	 * - send the buffer to process pid over InterfaceCommunicator com.
	*/
	virtual bool
	collect(ug::BinaryBuffer& buff, const Interface& interface)

	{
		int pid = interface.get_target_proc();
		for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			int index;
			char a = 0;
			if(derived().need_send(pid, index))
				a = 1;
			Serialize(buff, a);
			if(a) Serialize(buff, derived().send(pid, index);
		}
	}


	/** extract
	 * \brief receive data on the interface and hand it over to TDerived::receive
	 * \param	buff		BinaryBuffer with the data we received from processor pid.
	 * \param	interface	interface to processor pid
	 * This function does the following	 *
	 * - loop through the interface. deserialize each item from s, and hand it to TDerived::receive(pid, index, item)
	*/
	virtual bool
	extract(ug::BinaryBuffer& buff, const Interface& interface)
	{
		int pid = interface.get_target_proc();
		TValue val;
		for(IndexLayout::Interface::const_iterator iter = interface.begin(); iter != interface.end(); ++iter)
		{
			char a;
			Deserialize(buff, a);
			if(a)
			{
				Deserialize(buff, val);
				derived().receive(pid, interface.get_element(iter), val);
			}
		}
	}

	/** get_required_buffer_size
	 * \brief get the size of the required buffer for an interface
	 * \param	interface	The interface so we can calculate the size
	 * If TDerived can specify an exact size for each item it serializes, we can calculate the size of the buffer
	 * otherwise, it returns -1, and we return -1 -> Buffer size is sent.
	*/
	virtual int
	get_required_buffer_size(const Interface& interface)
	{
		int s = derived().get_element_size();
		if(s == -1)
			return -1;
		else
			return s * interface.size();
	}
};
#endif

// end group lib_algebra_parallelization_scheme
/// \}

} // namespace ug
#endif /* SEND_INTERFACE_H_ */

#endif /* UG_PARALLEL */
