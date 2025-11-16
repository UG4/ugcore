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
#ifndef PCL_UTIL_IMPL_H
#define PCL_UTIL_IMPL_H

#include "common/util/binary_buffer.h"
#include "common/serialization.h"

#include "pcl_util.h"

namespace pcl {
	template<typename TLayout, typename TSelector>
	bool RemoveUnselectedInterfaceEntries(TLayout &layout, TSelector &sel) {
		//	iterate over all interfaces of the layout.
		//	for each we'll create a new one, into which elements selected
		//	elements will be inserted.
		//	Finally we'll swap the content of the those interfaces.
		//	if the interface is empty at the end of the operation, it will be
		//	removed from the layout.
		bool retVal = false;

		using Interface = typename TLayout::Interface;
		//ø using InterfaceIter = typename TLayout::iterator;
		using Elem = typename Interface::Element;
		//ø using ElemIter = typename Interface::iterator;

		//	iterate over all levels
		for(size_t level = 0; level < layout.num_levels(); ++level)
		{
			//	iterate over all interfaces
			for(auto iiter = layout.begin(level); iiter != layout.end(level);)
			{
				bool interfaceChanged = false;
				Interface& interface = layout.interface(iiter);

				//	create a temporary interface and fill it with the selected entries
				Interface tInterface;

				for(auto iter = interface.begin(); iter != interface.end(); ++iter)
				{
					Elem& e = interface.get_element(iter);
					if(sel.is_selected(e))
						tInterface.push_back(e);
					else
						interfaceChanged = true;
				}

				//	now swap the interface contents.
				if(interfaceChanged){
					interface.swap(tInterface);

					//	if the interface is empty, erase it.
					//	if not, simply increase the iterator
					if(interface.size() == 0){
						iiter = layout.erase(iiter, level);
					}
					else{
						++iiter;
					}

					retVal = true;
				}
				else{
					++iiter;
				}
			}
		}

		return retVal;
	}

	template<typename TType, typename TLayoutMap, typename TSelector>
	bool RemoveUnselectedInterfaceEntries(TLayoutMap &lm, TSelector &sel) {
		//ø using iterator = typename TLayoutMap::template Types<TType>::Map::iterator;
		using Layout = typename TLayoutMap::template Types<TType>::Layout;

		bool retVal = false;
		for(auto iter = lm.template layouts_begin<TType>(); iter != lm.template layouts_end<TType>();)
		{
			//	get the layout
			Layout& layout = iter->second;
			//	remove unnecessary interface entries and interfaces
			retVal |= RemoveUnselectedInterfaceEntries(layout, sel);
			//	if the layout is empty, it can be removed from the map
			//	if not we'll simply increase the iterator
			if(layout.empty()){
				iter = lm.template erase_layout<TType>(iter);
			}
			else{
				++iter;
			}
		}
		return retVal;
	}

	template<typename TLayout, typename TSelectorIn, typename TSelectorOut>
	bool SelectionCommPol<TLayout, TSelectorIn, TSelectorOut>::collect(ug::BinaryBuffer &buff, Interface &interface) {
		char zero = 0;
		char one = 1;

		for(auto iter = interface.begin(); iter != interface.end(); ++iter)
		{
			if(m_selIn.is_selected(interface.get_element(iter)))
				buff.write(&one, sizeof(char));
			else
				buff.write(&zero, sizeof(char));
		}

		return true;
	}

	template<typename TLayout, typename TSelectorIn, typename TSelectorOut>
	bool SelectionCommPol<TLayout, TSelectorIn, TSelectorOut>::extract(ug::BinaryBuffer &buff, Interface &interface) {
		char tmp;
		for(auto iter = interface.begin(); iter != interface.end(); ++iter)
		{
			buff.read(&tmp, sizeof(char));
			if(tmp == 0)
				m_selOut.deselect(interface.get_element(iter));
			else
				m_selOut.select(interface.get_element(iter));
		}

		return true;
	}

template<typename TLayout>
void AddLayout(TLayout &destLayout, const TLayout &sourceLayout) {
	for(typename TLayout::const_iterator iter = sourceLayout.begin(); iter != sourceLayout.end(); ++iter)
	{
		const typename TLayout::Interface &source_interface = sourceLayout.interface(iter);
		typename TLayout::Interface &dest_interface = destLayout.interface(sourceLayout.proc_id(iter));
		for(auto iter2 = source_interface.begin(); iter2 != source_interface.end(); ++iter2)
			dest_interface.push_back(source_interface.get_element(iter2));
	}
}

template <typename TKey, typename TValue, typename Compare>
void MinimalKeyValuePairAcrossAllProcs(TKey& keyInOut, TValue& valInOut, const Compare& cmp)
{
	// in the serial case, input key and value are already output key and value
#ifndef UG_PARALLEL
	return;
#endif
	size_t nProc = NumProcs();
	if (nProc == 1)
		return;

	// write key/value pair into binary buffer
	ug::BinaryBuffer buf;
	ug::Serialize(buf, keyInOut);
	ug::Serialize(buf, valInOut);

	// gather buffers on root
	ProcessCommunicator pc;
	pc.gather(buf, 0);

	// find minimum on root
	size_t rk = ProcRank();
	if (rk == 0)
	{
		TKey key;
		TValue val;

		for (size_t i = 0; i < nProc; ++i)
		{
			ug::Deserialize(buf, key);
			ug::Deserialize(buf, val);

			if (cmp(key, keyInOut))
			{
				keyInOut = key;
				valInOut = val;
			}
		}

		buf.clear();
		ug::Serialize(buf, keyInOut);
		ug::Serialize(buf, valInOut);
	}

	// communicate minimum to all procs
	pc.broadcast(buf, 0);

	// copy to return values
	ug::Deserialize(buf, keyInOut);
	ug::Deserialize(buf, valInOut);
}


}
#endif