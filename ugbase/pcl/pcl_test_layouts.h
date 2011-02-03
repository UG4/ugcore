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

namespace pcl
{

////////////////////////////////////////////////////////////////////////
//  TestLayoutIsDoubleEnded
/// tests if masterLayouts proc id's find a match in correspoding slaveLayouts proc ids.
/**
 * that is, iff processor P1 has processor P2 in his masterLayout, then processor P2 needs to have P1 in his slaveLayout
 */
template<typename TLayout>
void TestLayoutIsDoubleEnded(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout, TLayout &slaveLayout)
{
	//using namespace std;

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

		UG_ASSERT(bMasterToThisProcess == bSlaveToProcess[i],
			"Process " << std::setw(4) << i << " has " << (bMasterToThisProcess ? "a" : "no") << " master connection to this process (" << std::setw(4) << pcl::GetProcRank()
			<< "), but we have " << (bSlaveToProcess[i] ? "a" : "no") << " slave connection to " << std::setw(4) << i);
		UG_ASSERT(bSlaveToThisProcess == bMasterToProcess[i],
			"Process " << std::setw(4) << i << " has " << (bSlaveToThisProcess ? "a" : "no") << " slave connection to this process (" << pcl::GetProcRank()
			<< "), but we have " << (bMasterToProcess[i] ? "a" : "no") << " master connection to " << std::setw(4) << i);
	}
}

/// if processor P1 has a interface to P2, then the size of the interface P1->P2 has to be the same as the size of interface P2->P1
template<typename TLayout>
void TestSizeOfInterfacesInLayoutsMatch(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout, TLayout &slaveLayout, bool bPrint=false)
{
	//using namespace std;

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
			Serialize(stream, element);
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
				if(bPrint) UG_LOG(" " << std::setw(9) << interface.get_element(iter2) << " <-> BROKEN!\n");
			}
			else
			{
				typename Interface::Element element; Deserialize(stream, element);
				if(bPrint) UG_LOG(" " << std::setw(9) << interface.get_element(iter2) << " <-> " << std::setw(9) << element << "\n");
			}
		}
		while(stream.can_read_more())
		{
			broken = true;
			if(!bPrint) break;
			typename Interface::Element element; Deserialize(stream, element);
			UG_LOG(" BROKEN! -> " << element << "\n");
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
		TestSizeOfInterfacesInLayoutsMatch(com, masterLayout, slaveLayout, true);
	}
	else
	{
		UG_ASSERT(layout_broken == false, "One or more interfaces are broken");
	}
}

template<typename TLayout>
void TestLayout(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout, TLayout &slaveLayout, bool bPrint=false)
{
	TestLayoutIsDoubleEnded(com, masterLayout, slaveLayout);

	//using namespace std;

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

	TestSizeOfInterfacesInLayoutsMatch(com, masterLayout, slaveLayout, bPrint);
}


template<typename TLayout>
void PrintLayout(pcl::ParallelCommunicator<TLayout> &com, TLayout &masterLayout, TLayout &slaveLayout)
{
	TestLayout(com, masterLayout, slaveLayout, true);
}

}

#endif
