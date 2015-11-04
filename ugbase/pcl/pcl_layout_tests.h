#ifndef __H__PCL_TEST_LAYOUTS__
#define __H__PCL_TEST_LAYOUTS__

#include <vector>
#include <iomanip> // for 'std::setw()' etc.
#include <map>

#include "pcl_base.h"
#include "pcl_methods.h"
#include "pcl_communication_structs.h"
#include "pcl_interface_communicator.h"
#include "pcl_process_communicator.h"

#include "common/util/binary_buffer.h"
#include "common/log.h"
#include "common/assert.h"
#include "common/serialization.h"
#include "common/profiler/profiler.h"

#include <boost/function.hpp>

namespace pcl
{

/// \addtogroup pcl
/// \{

///	Trivial implementation of a to-value callback.
/**	TValue has to be constructable from TElem.*/
template <class TElem>
TElem TrivialToValue(TElem e)
{
	return e;
}
inline void PrintPC(const pcl::ProcessCommunicator &processCommunicator)
{
	UG_LOG(processCommunicator.size() << " involved procs: ");
	for(size_t i=0; i<processCommunicator.size(); i++)
	{	UG_LOG(processCommunicator.get_proc_id(i) << " "); }
	UG_LOG("\n");

}

#define PRINTPC(com) UG_LOG(__FUNCTION__ << " : " << __LINE__ << "\n"); PrintPC(com)
////////////////////////////////////////////////////////////////////////
//  TestLayoutIsDoubleEnded
/// tests if masterLayouts proc id's find a match in corresponding slaveLayouts proc ids.
/**
 * that is, iff processor P1 has processor P2 in his masterLayout, then processor P2 needs to have P1 in his slaveLayout
 */
template<typename TLayout>
bool TestLayoutIsDoubleEnded(const pcl::ProcessCommunicator processCommunicator,
		pcl::InterfaceCommunicator<TLayout> &com,
		const TLayout &masterLayout, const TLayout &slaveLayout)
{
	PROFILE_FUNC_GROUP("debug");
	//PRINTPC(processCommunicator);
//	check if connections are double-ended
	std::vector<char> bMasterToProcess; bMasterToProcess.resize(processCommunicator.size(), 0x00);
	std::vector<char> bSlaveToProcess; bSlaveToProcess.resize(processCommunicator.size(), 0x00);

	for(typename TLayout::const_iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		int id = processCommunicator.get_local_proc_id(masterLayout.proc_id(iter));
		if(id == -1)
			UG_LOG("Processor " << masterLayout.proc_id(iter) << " not in processCommunicator, but in MasterLayout\n")
		else bMasterToProcess[id] = true;
	}
	for(typename TLayout::const_iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		int id = processCommunicator.get_local_proc_id(slaveLayout.proc_id(iter));
		if(id == -1)
			UG_LOG("Processor " << slaveLayout.proc_id(iter) << " not in processCommunicator, but in SlaveLayout\n")
		else bSlaveToProcess[id] = true;
	}

	for(size_t i=0; i<processCommunicator.size(); i++)
	{
		com.send_raw(processCommunicator.get_proc_id(i), &bMasterToProcess[i], sizeof(char));
		com.send_raw(processCommunicator.get_proc_id(i), &bSlaveToProcess[i], sizeof(char));
	}

	std::vector<ug::BinaryBuffer> masterToThisProcessMap(processCommunicator.size());
	std::vector<ug::BinaryBuffer> slaveToThisProcessMap(processCommunicator.size());
	for(size_t i=0; i < processCommunicator.size(); i++)
	{
		com.receive_raw(processCommunicator.get_proc_id(i), masterToThisProcessMap[i]);
		com.receive_raw(processCommunicator.get_proc_id(i), slaveToThisProcessMap[i]);
	}

	com.communicate();

	for(size_t i=0; i<processCommunicator.size(); i++)
	{
		char bMasterToThisProcess, bSlaveToThisProcess;
		ug::Deserialize(masterToThisProcessMap[i], bMasterToThisProcess);
		ug::Deserialize(slaveToThisProcessMap[i], bSlaveToThisProcess);

		int pid = processCommunicator.get_proc_id(i);
		if(bMasterToThisProcess != bSlaveToProcess[i])
		{
			UG_LOG("Process " << std::setw(4) << pid << " has " << (bMasterToThisProcess ? "a" : "no") << " master connection to this process (" << std::setw(4) << pcl::ProcRank()
				<< "), but we have " << (bSlaveToProcess[i] ? "a" : "no") << " slave connection to " << std::setw(4) << pid << std::endl);
			return false;
		}
		if(bSlaveToThisProcess != bMasterToProcess[i]){
			UG_LOG("Process " << std::setw(4) << pid << " has " << (bSlaveToThisProcess ? "a" : "no") << " slave connection to this process (" << pcl::ProcRank()
				<< "), but we have " << (bMasterToProcess[i] ? "a" : "no") << " master connection to " << std::setw(4) << pid << std::endl);
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
 * If compareValues = true is specified, then all entries supplied by cbToValue
 * are compared between connected interface entries and mismatches are reported.
 *
 * The callback-method takes a parameter of type TLayout::Element and has to
 * return a value of type TValue (TValue = TLayout::Element by default).
 */
template<typename TLayout, typename TValue>
bool TestSizeOfInterfacesInLayoutsMatch(pcl::InterfaceCommunicator<TLayout> &com,
                                        const TLayout &masterLayout,
                                        const TLayout &slaveLayout, bool bPrint=false,
					boost::function<TValue (typename TLayout::Element)> cbToValue
						= TrivialToValue<typename TLayout::Element>,
					bool compareValues = false)
{
	PROFILE_FUNC_GROUP("debug");
	typedef std::map<int, ug::BinaryBuffer>	BufferMap;
	typedef typename TLayout::Interface Interface;

	BufferMap sendMap, receiveMap;

	for(typename TLayout::const_iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter)
	{
		const Interface &interface = slaveLayout.interface(iter);
		int pid = slaveLayout.proc_id(iter);
		ug::BinaryBuffer& buffer = sendMap[pid];

		for(typename TLayout::Interface::const_iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			const typename Interface::Element &element = interface.get_element(iter2);
			Serialize(buffer, cbToValue(element));
		}
		com.send_raw(pid, buffer.buffer(), buffer.write_pos(), false);
	}


	for(typename TLayout::const_iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		int pid = masterLayout.proc_id(iter);
		com.receive_raw(pid, receiveMap[pid]);
	}

	com.communicate();

	bool layoutBroken=false;
	bool valueMismatch = false;

	for(typename TLayout::const_iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)
	{
		int pid = masterLayout.proc_id(iter);
		ug::BinaryBuffer &buffer = receiveMap[pid];
		bool broken=false;
		if(bPrint) { UG_LOG("      Interface processor " << pcl::ProcRank() << " <-> processor " << pid << " (Master <-> Slave):\n"); }
		const typename TLayout::Interface &interface = masterLayout.interface(iter);
		std::stringstream brokenInformation;
		size_t onMaster=0, onSlave=0;
		for(typename TLayout::Interface::const_iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
		{
			TValue val1 = cbToValue(interface.get_element(iter2));

			if(buffer.eof())
			{
				 broken =true;
				 if(bPrint){
					 UG_LOG("        " << onMaster << ": " << std::setw(9) << val1 << " <-> " << "-?-" << "\n");
				 }
			}
			else
			{
				TValue val2;
				Deserialize(buffer, val2);
				onSlave++;

				bool mismatch = false;
				if(compareValues){
					mismatch = (val1 != val2);
				}

				if(mismatch){
					UG_LOG("        " << std::setw(9) << val1 << " <-> " << val2 << "  --- MISMATCH! ---"
							<< " (interface to proc " << pid << ")" << std::endl);
				}
				valueMismatch |= mismatch;
				if(bPrint)
				{
					UG_LOG("        " << onMaster << ": " << std::setw(9) << val1 << " <-> " << val2 << "\n");
				}
			}
			onMaster++;
		}
                            
		if(!buffer.eof())
		{
			while(!buffer.eof())
			{
				TValue val2;
				Deserialize(buffer, val2);
				if(bPrint)
				{UG_LOG("        " << onSlave << ": " << std::setw(9) << "-?-" << " <-> " << val2 << "\n");}
				onSlave++;
			}
			broken =true;
		}
		if(bPrint)
		{
			if(onMaster!=onSlave)
			{
				UG_LOG("   Interface sizes do not match:\n");
				UG_LOG("   - Slave Interface " << pid << " -> " << pcl::ProcRank() << ", size = " << onSlave << "\n");
				UG_LOG("   - Master Interface " << pcl::ProcRank() << " -> " << pid << ", size = " << onMaster << "\n");
			}
			else
			{	UG_LOG("      In total " << std::setw(9) << onMaster << " entries in this interface." << std::endl);	}
		}

		if(broken)
		{
			layoutBroken=true;
			UG_LOG("      Interface from processor " << std::setw(4) << pcl::ProcRank() << " to processor " << std::setw(4) << pid << " is BROKEN!\n");

		}

	}

	if(layoutBroken == true)
	{
		UG_LOG("One or more interfaces are broken\n");
		return false;
	}

	if(valueMismatch){
		UG_LOG("MISMATCH! Not all values at connected interface elements did match!\n");
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
 * If compareValues = true is specified, then all entries supplied by cbToValue
 * are compared between connected interface entries and mismatches are reported.
 *
 * The methods returns true if all interfaces match and false if not. The return
 * values are consistent among all processes.
 */
template<typename TLayout, typename TValue>
bool TestLayout(const pcl::ProcessCommunicator &processCommunicator,
                pcl::InterfaceCommunicator<TLayout> &com,
                const TLayout &masterLayout,
                const TLayout &slaveLayout, bool bPrint=false,
				boost::function<TValue (typename TLayout::Element)> cbToValue
					= TrivialToValue<typename TLayout::Element>,
				bool compareValues = false)
{
	PROFILE_FUNC_GROUP("debug");
	if(bPrint)
	{
		UG_LOG("proc " << std::setw(4) << pcl::ProcRank() << ":\n");
		UG_LOG("   MasterLayout is to processes ");
		for(typename TLayout::const_iterator iter = masterLayout.begin(); iter != masterLayout.end(); ++iter)	{
			UG_LOG(" " << std::setw(4) << masterLayout.proc_id(iter) << " ");
		}
		UG_LOG("\n   Slave Layout is to processes ");
		for(typename TLayout::const_iterator iter = slaveLayout.begin(); iter != slaveLayout.end(); ++iter) {
			UG_LOG(" " << std::setw(4) << slaveLayout.proc_id(iter) << " ");
		}
		UG_LOG("\n");
	}
	bool bDoubleEnded = TestLayoutIsDoubleEnded(processCommunicator,
			com, masterLayout, slaveLayout);
	if(!pcl::AllProcsTrue(bDoubleEnded, processCommunicator))
		return false;

	bool bSuccess = TestSizeOfInterfacesInLayoutsMatch<TLayout, TValue>(com, masterLayout,
											slaveLayout, bPrint, cbToValue, compareValues);
	return pcl::AllProcsTrue(bSuccess, processCommunicator);
}

template<typename TLayout>
bool TestLayout(const pcl::ProcessCommunicator &processCommunicator,
				pcl::InterfaceCommunicator<TLayout> &com, const TLayout &masterLayout,
				const TLayout &slaveLayout, bool bPrint=false,
				bool compareValues = false)
{
	return TestLayout<TLayout, typename TLayout::Element>(processCommunicator, com,
			masterLayout, slaveLayout, bPrint,
			TrivialToValue<typename TLayout::Element>, compareValues);
}

template<typename TLayout, typename TValue>
bool PrintLayout(const pcl::ProcessCommunicator &processCommunicator,
                 pcl::InterfaceCommunicator<TLayout> &com, const TLayout &masterLayout,
                 const TLayout &slaveLayout,
                 boost::function<TValue (typename TLayout::Element)> cbToValue
					= TrivialToValue<typename TLayout::Element>)
{
	return TestLayout(processCommunicator, com, masterLayout, slaveLayout, true, cbToValue);
}


template<typename TLayout>
inline bool PrintLayout(const pcl::ProcessCommunicator &processCommunicator,
                 pcl::InterfaceCommunicator<TLayout> &com,
                 const TLayout &masterLayout,
                 const TLayout &slaveLayout)
{
	return TestLayout(processCommunicator, com, masterLayout, slaveLayout, true);
}



template<typename TLayout>
inline void PrintLayout(const TLayout &layout)
{
	for(typename TLayout::const_iterator iter = layout.begin(); iter != layout.end(); ++iter)
	{
		size_t pid = layout.proc_id(iter);
		UG_LOG("to processor " << pid << ": ");
		const typename TLayout::Interface &interface = layout.interface(iter);
		for(typename TLayout::Interface::const_iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
			UG_LOG(interface.get_element(iter2) << "  ");
		UG_LOG("\n");
	}
}

#define TESTLAYOUT(processCommunicator, com, master, slave) { if(TestLayout(processCommunicator,\
		com, master, slave) == false) { PrintLayout(processCommunicator, com, master, slave); UG_COND_THROW(true, "layout broken"); } }

#ifndef NDEBUG
#define DEBUG_TESTLAYOUT(processCommunicator, com, master, slave) TESTLAYOUT(processCommunicator, com, master, slave)
#else
	#define DEBUG_TESTLAYOUT(processCommunicator, com, master, slave)
#endif

#define DEBUG_TESTLAYOUTS(layout) DEBUG_TESTLAYOUT(layout->proc_comm(), layout->comm(), layout->master(), layout->slave())
#define TESTLAYOUTS(layout) TESTLAYOUT(layout->proc_comm(), layout->comm(), layout->master(), layout->slave())

// end group pcl
/// \}

} // end namespace pcl

#endif
