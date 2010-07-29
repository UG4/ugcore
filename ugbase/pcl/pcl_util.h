//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d29

#ifndef __H__PCL_UTIL__
#define __H__PCL_UTIL__

#include <iostream>
#include "pcl_base.h"
#include "pcl_communication_structs.h"
#include "pcl_communicator.h"
#include "pcl_process_communicator.h"

namespace pcl
{
////////////////////////////////////////////////////////////////////////
/**
 * Removes unselected entries from the interfaces in the given layout.
 * Empty interfaces are removed from the layout, too.
 *
 * TLayout has to be compatible with pcl::Layout or pcl::MultiLevelLayout.
 *
 * TSelector has to feature a method bool TSelector::is_selected(TType* t);
 */
template <class TLayout, class TSelector>
void RemoveUnselectedInterfaceEntries(TLayout& layout, TSelector& sel)
{
	using namespace pcl;
//	iterate over all interfaces of the layout.
//	for each we'll create a new one, into which elements selected
//	elements will be inserted.
//	Finally we'll swap the content of the those interfaces.
//	if the interface is empty at the end of the operation, it will be
//	removed from the layout.

//	some typedefs first
	typedef typename TLayout::Interface Interface;
	typedef typename TLayout::iterator InterfaceIter;
	typedef typename Interface::Element	Elem;
	typedef typename Interface::iterator ElemIter;
	
//	iterate over all levels
	for(size_t level = 0; level < layout.num_levels(); ++level)
	{
	//	iterate over all interfaces
		for(InterfaceIter iiter = layout.begin(level);
			iiter != layout.end(level);)
		{
			Interface& interface = layout.interface(iiter);
		
		//	create a temporary interface and fill it with the selected entries
			Interface tInterface;
			
			for(ElemIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Elem& e = interface.get_element(iter);
				if(sel.is_selected(e))
					tInterface.push_back(e);
			}
			
		//	now swap the interface contents.
			interface.swap(tInterface);
		
		//	if the interface is empty, erase it.
		//	if not, simply increase the iterator
			if(interface.size() == 0)
				iiter = layout.erase(iiter, level);
			else
				++iiter;
		}
	}
}

////////////////////////////////////////////////////////////////////////
/**
 * Removes interface-entries, empty interfaces and empty layouts from
 * the given layoutMap for the given type.
 *
 * TLayoutMap has to be compatible with pcl::LayoutMap.
 *
 * TSelector has to feature a method bool TSelector::is_selected(TType* t);
 */
template <class TType, class TLayoutMap, class TSelector>
void RemoveUnselectedInterfaceEntries(TLayoutMap& lm, TSelector& sel)
{
	typedef typename TLayoutMap::template Types<TType>::Map::iterator iterator;
	typedef typename TLayoutMap::template Types<TType>::Layout		Layout;

	for(iterator iter = lm.template layouts_begin<TType>();
		iter != lm.template layouts_end<TType>();)
	{
	//	get the layout
		Layout& layout = iter->second;
	//	remove unnecessary interface entries and interfaces
		RemoveUnselectedInterfaceEntries(layout, sel);
	//	if the layout is empty, it can be removed from the map
	//	if not we'll simply increase the iterator
		if(layout.empty())
			iter = lm.template erase_layout<TType>(iter);
		else
			++iter;
	}
}

}//	end of namespace

#endif
