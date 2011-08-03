/*
 * parallelization.h
 *
 *  Created on: 02.08.2011
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__PARALLELIZATION_H_
#define __H__LIB_ALGEBRA__PARALLELIZATION_H_

#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallelization.h"
namespace ug
{

template<typename T>
T Deserialize(BinaryBuffer &stream)
{
	T t;
	Deserialize(stream, t);
	return t;
}
template<typename TLocalToGlobal>
void SerializeInterface(BinaryBuffer &stream,
		IndexLayout::Interface &interface, const TLocalToGlobal &localToGlobal)
{
	Serialize(stream, interface.size());
	for(IndexLayout::Interface::iterator iter2 = interface.begin(); iter2 != interface.end(); ++iter2)
	{
		size_t localIndex = interface.get_element(iter2);
		AlgebraID &globalIndex = localToGlobal[localIndex];
		Serialize(stream, globalIndex);
	}
}

template<typename TLocalToGlobal>
void SerializeLayout(BinaryBuffer &stream,
		IndexLayout &layout, const TLocalToGlobal &localToGlobal)
{
	size_t size = 0;
	for(IndexLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
		size++;
	Serialize(stream, size);
	for(IndexLayout::iterator iter = layout.begin(); iter != layout.end(); ++iter)
	{
		IndexLayout::Interface &interface = layout.interface(iter);
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
typename TLayout::iterator find(TLayout &layout, int pid)
{
	for(typename TLayout::iterator it = layout.begin(); it != layout.end(); ++it)
		if(layout.proc_id(it) == pid) return it;
	return layout.end();
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
#endif /* PARALLELIZATION_H_ */
