/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__PCL_UTIL__
#define __H__PCL_UTIL__

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <functional>  // for std::less

#include "common/util/binary_buffer.h"

#include "pcl_layout_util.h"
#include "pcl_base.h"
#include "pcl_communication_structs.h"
#include "pcl_process_communicator.h"

#ifdef PCL_DEBUG_BARRIER_ENABLED
///	A debug barrier. Halts program execution until all processes of the communicator
///	have called the barrier function.
/**	Note that this barrier is only has effect, if the defined PCL_DEBUG_BARRIER_ENABLED
 * is enabled.*/
	#define PCL_DEBUG_BARRIER(communicator) communicator.barrier()
	#define PCL_DEBUG_BARRIER_ALL()	{pcl::ProcessCommunicator com; com.barrier();}

#else
	#define PCL_DEBUG_BARRIER(communicator)
	#define PCL_DEBUG_BARRIER_ALL()
#endif

namespace pcl
{

/// \addtogroup pcl
/// \{

/// synchronizes all processes.
///	This is just a shortcut if all processes have to be synchronized.
/**	This is equivalent to calling pcl::ProcessCommunicator()::barrier().
 * \sa pcl::ProcessCommunicator::barrier
 */
void SynchronizeProcesses();

/// performs an allreduce and returns true if and only if all procs called the
/// function with bFlag = true
bool AllProcsTrue(bool bFlag, ProcessCommunicator comm = ProcessCommunicator());

/// performs an allreduce and returns true at least one process called the
/// function with bFlag = true
bool OneProcTrue(bool bFlag, ProcessCommunicator comm = ProcessCommunicator());

////////////////////////////////////////////////////////////////////////
///	exchanges information about which process wants to communicate with which other process.
/**	If processes want to send data to other processes, but the target
 * processes do not know from where to receive their data, then this
 * method can be used to notify each process from where he has to
 * receive data.
 *
 * Each process has to call this method, passing an array of process-ids
 * to where it will send data (vSendToRanks).
 *
 * The processes from which it has to receive data are then collected in
 * vReceiveFromRanksOut.
 *
 * Note that this method involves global communication (with respect to
 * the used process-communicator) and should thus be rarely used - at least
 * if the used process communicator is big.
 */
void CommunicateInvolvedProcesses(std::vector<int>& vReceiveFromRanksOut,
								  const std::vector<int>& vSendToRanks,
								  const ProcessCommunicator& procComm
								  	= ProcessCommunicator());

////////////////////////////////////////////////////////////////////////
/**
 * Removes unselected entries from the interfaces in the given layout.
 * Empty interfaces are removed from the layout, too.
 *
 * TLayout has to be compatible with pcl::Layout or pcl::MultiLevelLayout.
 *
 * \code
 * bool TSelector::is_selected(TType* t);
 * \endcode
 *
 * \returns:	true, if the layout has changed, false if not.
 */
template <typename TLayout, typename TSelector>
bool RemoveUnselectedInterfaceEntries(TLayout& layout, TSelector& sel);



////////////////////////////////////////////////////////////////////////
/**
 * Removes interface-entries, empty interfaces and empty layouts from
 * the given layoutMap for the given type.
 *
 * TLayoutMap has to be compatible with pcl::LayoutMap.
 *
 * TSelector has to feature a method
 * \code
 * bool TSelector::is_selected(TType* t);
 * \endcode
 *
 * \returns:	true, if the layout-map has changed, false if not.
 */
template <typename TType, typename TLayoutMap, typename TSelector>
bool RemoveUnselectedInterfaceEntries(TLayoutMap& lm,
										TSelector& sel);


////////////////////////////////////////////////////////////////////////
///	communicates selection-status of interface elements
/**
 *	TLayout has to be compatible with the pcl::layout_tags.
 *
 *	TSelectorIn has to feature a method
 *	\code
 *	bool TSelectorIn::is_selected(TLayout::Element e);
 *	\endcode
 *
 *	TSelectorOut has to feature methods
 *	\code
 *	void TSelectorOut::select(TLayout::Element e);
 *	void TSelectorOut::deselect(TLayout::Element e);
 *	\endcode
 */
template <typename TLayout, typename TSelectorIn, typename TSelectorOut>
class SelectionCommPol : public ICommunicationPolicy<TLayout>
{
	public:
		using Interface = typename ICommunicationPolicy<TLayout>::Interface;
		
	public:
		SelectionCommPol(TSelectorIn& selIn, TSelectorOut& selOut) :
			m_selIn(selIn), m_selOut(selOut)	{}
			
	///	iterates over the interface entries. Writes 1 for selected, 0 for unselected.
		bool collect(ug::BinaryBuffer& buff, Interface& interface) override;

		///	iterates over the interface entries. selects for 1, deselects for 0.
		bool extract(ug::BinaryBuffer& buff, Interface& interface) override;

	protected:
		TSelectorIn& m_selIn;
		TSelectorOut& m_selOut;
};


///	checks whether proc-entries in send- and recv-lists on participating processes match
/** The return value is the same for all participating processes.
 * \sa pcl::SendRecvBuffersMatch
 */
bool SendRecvListsMatch(const std::vector<int>& recvFrom,
						const std::vector<int>& sendTo,
						const ProcessCommunicator& involvedProcs = ProcessCommunicator());

///	checks whether matching buffers in send- and recv-lists have the same size

/**	Note that this method does not check whether matching buffers exist. This is assumed.
 * Checks are only performed on the sizes of associated buffers.
 *
 * Make sure that recvBufSizes specifies the buffer size for each entry in recvFrom.
 * Same for sendBufSizes and sendTo.
 *
 * The return value is the same for all participating processes.
 * \sa pcl::SendRecvMapsMatch
 *
 * @param recvFrom
 * @param recvBufSizes
 * @param sendTo
 * @param sendBufSizes
 * @param involvedProcs
 * @return
 */
bool SendRecvBuffersMatch(const std::vector<int>& recvFrom,
							const std::vector<int>& recvBufSizes,
							const std::vector<int>& sendTo,
							const std::vector<int>& sendBufSizes,
							const ProcessCommunicator& involvedProcs = ProcessCommunicator());


/**
 *
 * @tparam TLayout
 * @param destLayout
 * @param sourceLayout
 */
template<typename TLayout>
void AddLayout(TLayout &destLayout, const TLayout &sourceLayout);


/// util function to read a file in parallel.
/** \sa ReadFile in serial
 * \param filename filename (in/out!)
 * \param file the file
 * \param bText if true, use r instead of rb for file functions
 * \param bDistributedLoad if false, just use ReadFile.
 * \param pc for rank which should open files.
 * core with rank pc.get_proc_id(0) reads this filename. status, filename, and file are transmitted to other cores
 * \return false on all cores if
 */
bool ParallelReadFile(std::string &filename,
						std::vector<char> &file,
						bool bText,
						bool bDistributedLoad,
						const ProcessCommunicator& pc = ProcessCommunicator());


/**
 * @brief Find minimal key/value pair across processes
 * This function will receive one key/value pair from each process.
 * They will be gathered on proc 0, where the minimal key (w.r.t. given Compare object, e.g., std::less<TKey> )
 * and corresponding value will be determined.
 * The minimal pair will be made known to all procs and returned.
 *
 * This is a procedure that is required fairly often, e.g., when position-attached
 * data is to be located by a nearest-neighbor search.
 *
 * Key and value types must be serializable in a binary buffer.
 *
 * @param keyInOut  each proc's key (used as input and as output)
 * @param valInOut  each proc's value (used as input and as output)
 * @param cmp       object to use for key comparison (typically std::less<Key>)
 */
template <typename TKey, typename TValue, typename Compare>
void MinimalKeyValuePairAcrossAllProcs(TKey& keyInOut,
										TValue& valInOut,
										const Compare& cmp = Compare());


// end group pcl
/// \}

}//	end of namespace


#include "pcl_util_impl.h"

#endif
