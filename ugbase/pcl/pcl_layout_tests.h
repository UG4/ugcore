// created by Martin Rupp, 31.01.2011

#ifndef __H__PCL_TEST_LAYOUTS__
#define __H__PCL_TEST_LAYOUTS__

#include "pcl_base.h"
#include "pcl_methods.h"
#include "pcl_communication_structs.h"
#include "pcl_communicator.h"
#include "pcl_process_communicator.h"

#include <vector>
#include "common/util/stream_pack.h"
#include "common/log.h"
#include "common/assert.h"
#include "common/serialization.h"
#include <iomanip> // for 'std::setw()' etc.

#include <boost/function.hpp>

namespace pcl
{

///	Trivial implementation of a to-value callback.
/**	TValue has to be constructable from TElem.*/
template <class TElem>
TElem TrivialToValue(TElem e)
{
	return e;
}

////////////////////////////////////////////////////////////////////////
//  TestLayoutIsDoubleEnded
/// tests if masterLayouts proc id's find a match in correspoding slaveLayouts proc ids.
/**
 * that is, iff processor P1 has processor P2 in his masterLayout, then processor P2 needs to have P1 in his slaveLayout
 */
template<typename TLayout>
bool TestLayoutIsDoubleEnded(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout, TLayout &slaveLayout)
{
	// check if connections are double-ended
	std::vector<char> bMasterToProcess; bMasterToProcess.resize(pcl::GetNumProcesses(), 0x00);
	std::vector<char> bSlaveToProcess; bSlaveToProcess.resize(pcl::GetNumProcesses(), 0x00);

	for(typename TLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
		bMasterToProcess[masterLayout.proc_id(iter)] = true;
	for(typename TLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
		bSlaveToProcess[slaveLayout.proc_id(iter)] = true;

	for(int i=0; i<pcl::GetNumProcesses(); i++)
	{
		com.send_raw(i, &bMasterToProcess[i], sizeof(char));
		com.send_raw(i, &bSlaveToProcess[i], sizeof(char));
	}

	ug::StreamPack masterToThisProcessPack, slaveToThisProcessPack;
	for(int i=0; i<pcl::GetNumProcesses(); i++)
	{
		com.receive_raw(i, *masterToThisProcessPack.get_stream(i));
		com.receive_raw(i, *slaveToThisProcessPack.get_stream(i));
	}

	com.communicate();

	for(int i=0; i<pcl::GetNumProcesses(); i++)
	{
		char bMasterToThisProcess, bSlaveToThisProcess;
		ug::Deserialize(*masterToThisProcessPack.get_stream(i), bMasterToThisProcess);
		ug::Deserialize(*slaveToThisProcessPack.get_stream(i), bSlaveToThisProcess);

		if(bMasterToThisProcess != bSlaveToProcess[i])
		{
			UG_LOG("Process " << std::setw(4) << i << " has " << (bMasterToThisProcess ? "a" : "no") << " master connection to this process (" << std::setw(4) << pcl::GetProcRank()
				<< "), but we have " << (bSlaveToProcess[i] ? "a" : "no") << " slave connection to " << std::setw(4) << i << std::endl);
			return false;
		}
		if(bSlaveToThisProcess != bMasterToProcess[i]){
			UG_LOG("Process " << std::setw(4) << i << " has " << (bSlaveToThisProcess ? "a" : "no") << " slave connection to this process (" << pcl::GetProcRank()
				<< "), but we have " << (bMasterToProcess[i] ? "a" : "no") << " master connection to " << std::setw(4) << i << std::endl);
			return false;
		}
	}

	return true;
}

/// if processor P1 has a interface to P2, then the size of the interface P1->P2 has to be the same as the size of interface P2->P1
/**	You may specify a callback, which allows to print not the elements in the
 * interfaces directly, but to print associated values. This is useful, if e.g.
 * pointers are stored in the layouts.
 *
 * The callback-method takes a parameter of type TLayout::Element and has to
 * return a value of type TValue (TValue = TLayout::Element by default).
 */
template<typename TLayout, typename TValue>
bool TestSizeOfInterfacesInLayoutsMatch(pcl::ParallelCommunicator<TLayout> &com,
					TLayout &masterLayout, TLayout &slaveLayout, bool bPrint=false,
					boost::function<TValue (typename TLayout::Element)> cbToValue
						= TrivialToValue<typename TLayout::Element>)
{
	typedef typename TLayout::Interface Interface;
	ug::StreamPack sendpack, receivepack;
	for(typename TLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		Interface &interface = slaveLayout.interface(iter);
		int pid = slaveLayout.proc_id(iter);
		ug::BinaryStream &stream = *sendpack.get_stream(pid);

		for(typename TLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			typename Interface::Element &element = interface.get_element(iter2);
			Serialize(stream, cbToValue(element));
		}
		com.send_raw(pid, stream.buffer(), stream.size(), false);
	}


	for(typename TLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		int pid = masterLayout.proc_id(iter);
		com.receive_raw(pid, *receivepack.get_stream(pid));
	}

	com.communicate();

	bool layout_broken=false;
	for(typename TLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		int pid = masterLayout.proc_id(iter);
		ug::BinaryStream &stream = *receivepack.get_stream(pid);
		bool broken=false;
		if(bPrint) UG_LOG("Interface processor " << pcl::GetProcRank() << " <-> processor " << pid << " (Master <-> Slave):\n");
		typename TLayout::Interface &interface = masterLayout.interface(iter);
		for(typename TLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			if(stream.can_read_more() == false)
			{
				broken =true;
				if(bPrint){
					TValue value = cbToValue(interface.get_element(iter2));
					UG_LOG(" " << std::setw(9) << value << " <-> " << "BROKEN!" << std::endl);
				}
			}
			else
			{
				TValue val1 = cbToValue(interface.get_element(iter2));
				TValue val2; Deserialize(stream, val2);
				if(bPrint){
					UG_LOG(" " << std::setw(9) << val1 << " <-> " << val2 << std::endl);
				}
			}
		}
		while(stream.can_read_more())
		{
			broken = true;
			if(!bPrint) break;
			TValue value; Deserialize(stream, value);

			UG_LOG(" BROKEN! -> " << std::setw(9) << value << std::endl);
		}

		if(broken)
		{
			if(!bPrint) break;
			UG_LOG("Interface from processor " << std::setw(4) << pcl::GetProcRank() << " to processor " << std::setw(4) << pid << " is BROKEN!\n");
			layout_broken=true;
		}
	}

	if(!bPrint && layout_broken)
	{
		UG_LOG("\n\nOne or more interfaces are broken, printing interfaces:\n")
		TestSizeOfInterfacesInLayoutsMatch<TLayout, TValue>(com, masterLayout, slaveLayout,
															true, cbToValue);
		return false;
	}
	else if(layout_broken == true){
		UG_LOG("One or more interfaces are broken\n");
		return false;
	}
	return true;
}

///	Checks whether the given layouts are consistent.
/**	Checks whether the given master and slave layout have the same number of elements
 * and prints associated elements if bPrint = true. If the callback cbToValue is
 * specified, then it is used during printing to print values associated with the
 * elements in the layouts instead of printing the elements directly.
 *
 * The methods returns true if all interfaces match and false if not. The return
 * values are consistent among all processes.
 */
template<typename TLayout, typename TValue>
bool TestLayout(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout,
				TLayout &slaveLayout, bool bPrint=false,
				boost::function<TValue (typename TLayout::Element)> cbToValue
					= TrivialToValue<typename TLayout::Element>)
{
	bool bDoubleEnded = TestLayoutIsDoubleEnded(com, masterLayout, slaveLayout);

	if(!pcl::AllProcsTrue(bDoubleEnded))
		return false;

	if(bPrint)
	{
		UG_LOG("MasterLayout is to processes ");
		for(typename TLayout::iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)	{
			UG_LOG(" " << std::setw(4) << masterLayout.proc_id(iter) << " ");
		}
		UG_LOG("\nSlave Layout is to processes ");
		for(typename TLayout::iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter) {
			UG_LOG(" " << std::setw(4) << slaveLayout.proc_id(iter) << " ");
		}
		UG_LOG("\n");
	}

	bool bSuccess = TestSizeOfInterfacesInLayoutsMatch<TLayout, TValue>(com, masterLayout,
															slaveLayout, bPrint, cbToValue);
	return pcl::AllProcsTrue(bSuccess);
}

template<typename TLayout>
bool TestLayout(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout,
				TLayout &slaveLayout, bool bPrint=false)
{
	return TestLayout<TLayout, typename TLayout::Element>(com, masterLayout, slaveLayout, bPrint);
}

template<typename TLayout, typename TValue>
bool PrintLayout(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout,
				TLayout &slaveLayout,
				boost::function<TValue (typename TLayout::Element)> cbToValue
					= TrivialToValue<typename TLayout::Element>)
{
	return TestLayout(com, masterLayout, slaveLayout, true, cbToValue);
}

template<typename TLayout>
bool PrintLayout(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout,
				TLayout &slaveLayout)
{
	return TestLayout(com, masterLayout, slaveLayout, true);
}

}

#endif
