/*
 * serialize_interfaces.h
 *
 *  Created on: 02.08.2011
 *      Author: mrupp
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
