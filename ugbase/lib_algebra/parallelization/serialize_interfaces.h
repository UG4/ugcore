/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_ALGEBRA__SERIALIZE_INTERFACES_H_
#define __H__LIB_ALGEBRA__SERIALIZE_INTERFACES_H_


#ifndef UG_PARALLEL
#error "This only works with a UG_PARALLEL define."
#endif

#include "pcl/pcl.h"
#include "parallelization.h"
#include "common/serialization.h"

namespace ug
{


template<typename TLocalToGlobal>
void SerializeInterface(BinaryBuffer &stream,
                        const IndexLayout::Interface &interface, const TLocalToGlobal &localToGlobal)
{
	Serialize(stream, interface.size());
	for(IndexLayout::Interface::const_iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
	{
		size_t localIndex = interface.get_element(iter2);
		AlgebraID &globalIndex = localToGlobal[localIndex];
		Serialize(stream, globalIndex);
	}
}

template<typename TLocalToGlobal>
void SerializeLayout(BinaryBuffer &stream,
                     const IndexLayout &layout, const TLocalToGlobal &localToGlobal)
{
	size_t size = 0;
	for(IndexLayout::const_iterator iter = layout.begin(); iter != layout.end(); ++iter)
		size++;
	Serialize(stream, size);
	for(IndexLayout::const_iterator iter = layout.begin(); iter != layout.end(); ++iter)
	{
		const IndexLayout::Interface &interface = layout.interface(iter);
		int pid = layout.proc_id(iter);
		Serialize(stream, pid);
		Serialize(stream, interface);
	}
}


template<typename TGlobalToLocal>
void DeserializeInterface(BinaryBuffer &stream,
		IndexLayout::Interface &interface, const TGlobalToLocal &globalToLocal)
{

	size_t size;
	Deserialize(stream, size);
	for(size_t i=0; i<size; i++)
	{
		AlgebraID globalIndex;
		Deserialize(stream, globalIndex);
		size_t localIndex = globalToLocal[globalIndex];
		interface.push_back(localIndex);
	}
}

template<typename TGlobalToLocal>
void DeserializeLayout(BinaryBuffer &stream,
		IndexLayout &layout, const TGlobalToLocal &globalToLocal)
{
	size_t size;
	Deserialize(stream, size);
	for(size_t i=0; i<size; i++)
	{
		int pid = Deserialize<int>(stream);
		DeserializeInterface(stream, layout.interface(pid), globalToLocal);
	}
}


template<typename TLayout>
bool RemoveInterface(TLayout &layout, int pid)
{
	typename TLayout::iterator it = find(layout, pid);
	if(it == layout.end()) return false;
	layout.erase(it);
	return true;
}
// removes all interfaces which are within group
template<typename TLayout>
void RemoveInterfaces(TLayout &layout, std::vector<int> group)
{
	for(size_t i=0; i<group.size(); i++)
		RemoveInterface(layout, group[i]);
}

template<typename TLayout>
bool AppendInterface(TLayout &layout, int pidSource, int pidAppendTo)
{
	typename TLayout::iterator it = find(layout, pidSource);
	if(it == layout.end()) return false;
	typename TLayout::Interface interfaceSource = layout.interface(it);
	if(interfaceSource.empty()) return false;
	typename TLayout::Interface interfaceAppendTo = layout.interface(pidAppendTo);

	for(typename TLayout::Interface::iterator itSource = interfaceSource.begin(); itSource != interfaceSource.end(); ++itSource)
		interfaceAppendTo.push_back(interfaceSource.get_element(itSource));
	return true;
}

template<typename TLayout>
void MergeInterfaces(TLayout &layout, const std::vector<int> pidSources, int pidAppendTo)
{
	for(size_t i=0; i<pidSources.size(); i++)
	{
		AppendInterface(layout, pidSources[i], pidAppendTo);
		RemoveInterface(layout, pidSources[i]);
	}
}


template<typename TLayout>
void MergeInterfaces(TLayout &layout, const std::map<int, int> merge)
{
	for(std::map<int, int>::const_iterator it = merge.begin(); it != merge.end(); ++it)
	{
		AppendInterface(layout, it->first, it->second);
		RemoveInterface(layout, it->first);
	}
}

}
#endif /* SERIALIZE_INTERFACES.H */
