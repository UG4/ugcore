//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m07 d29

#ifndef __H__PCL_UTIL__
#define __H__PCL_UTIL__

#include <iostream>
#include <vector>
#include <algorithm>
#include "common/util/hash.h"
#include "pcl_base.h"
#include "pcl_communication_structs.h"
#include "pcl_communicator.h"
#include "pcl_process_communicator.h"

namespace pcl
{

////////////////////////////////////////////////////////////////////////
///	holds an id
/**	Default constructor creates an invalid id.*/
/*
template <class TLocalID>
struct ID
{
	ID(int procID = -1, TLocalID localID = 0)
		:m_procID(procID), m_localID(localID)	{}

///	returns true if m_procID >= 0.
	bool is_valid()		{return m_procID >= 0;}

///	sets m_procID to -1.
	void invalidate()	{m_procID = -1;}


	int m_procID;
	TLocalID m_localID;
};


////////////////////////////////////////////////////////////////////////
///	generates global ids for the elements in the given layout.
template <class TLayout, class TLocalID>
void DistributeGlobalIDs(TLayout& masterLayout, TLayout& slaveLayout,
						boost::function<void (typename TLayout::Element elem,
							ID<TLocalID>)> cbSetID,
						boost::function<ID<TLocalID>
							(typename TLayout::Element elem)> cbGetID)
{
	typedef ID<TLocalID> TID;
	TID invalidID;

//	iterate over all elements in the layouts and invalidate their ids.

}
*/

////////////////////////////////////////////////////////////////////////
///	collects the ids of all processes to which interfaces exist.
/**
 * Fills a vector with the process-ids, to which interfaces exist in
 * the given layout.
 *
 * TLayout has to be compatible with pcl::Layout or pcl::MultiLevelLayout.
 *
 * \returns the number of associated processes.
 */
template <class TLayout>
size_t CollectAssociatedProcesses(std::vector<int>& procIDsOut,
								  TLayout& layout)
{
	procIDsOut.clear();
	
//	iterate through the levels of the layout
	for(size_t i = 0; i < layout.num_levels(); ++i){
	//	iterate through the interfaces on that level
		for(typename TLayout::iterator iIter = layout.begin(i);
			iIter != layout.end(i); ++iIter)
		{
			int procID = layout.proc_id(iIter);
		//	check whether the process is already contained in procIDsOut
			if(i > 0){
				if(find(procIDsOut.begin(), procIDsOut.end(), procID)
				   == procIDsOut.end())
				{
				//	the entry has not yet been added
					procIDsOut.push_back(procID);
				}
			}
			else{
			//	on level 0 each process exists only once
				procIDsOut.push_back(procID);
			}
		}
	}
	
	return procIDsOut.size();
}

///	writes all elements in the interfaces into the vector.
/**
 * This function extracts all elements from a layout and stores them into
 * a std::vector. Doubles may occur and are not removed. The container is
 * clear as default, before extracting.
 */
template <class TLayout>
void CollectElements(std::vector<typename TLayout::Element>& elemsOut,
					TLayout& layout,
					bool clearContainer = true)
{
	typedef typename TLayout::Interface Interface;
	typedef typename TLayout::Element TElem;

//	clear the return value
	if(clearContainer) elemsOut.clear();

//	iterate over all interfaces
	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
		for(typename TLayout::iterator interfaceIter = layout.begin(lvl);
			interfaceIter != layout.end(lvl); ++interfaceIter)
		{
		//	iterate over the entries of the interface
			Interface& interface = layout.interface(interfaceIter);
			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	add elem to vector
				elemsOut.push_back(interface.get_element(iter));
			}
		}
	}
}

///	writes all elements in the interfaces into the resulting vector. avoids doubles.
template <class TLayout>
void CollectUniqueElements(std::vector<typename TLayout::Element>& elemsOut,
							TLayout& layout)
{
	typedef typename TLayout::Interface Interface;
	typedef typename TLayout::Element TElem;

//	clear the return value
	elemsOut.clear();

//	we'll use a hash to make sure that each element only exists once
	ug::Hash<int, TElem> hash(layout.num_interface_elements());

//	iterate over all interfaces
	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
		for(typename TLayout::iterator interfaceIter = layout.begin(lvl);
			interfaceIter != layout.end(lvl); ++interfaceIter)
		{
		//	iterate over the entries of the interface
			Interface& interface = layout.interface(interfaceIter);
			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	check whether the entry already exists in the hash
				if(!hash.has_entries(interface.get_element(iter))){
				//	we don't care about the value
					hash.add(0, interface.get_element(iter));
					elemsOut.push_back(interface.get_element(iter));
				}
			}
		}
	}
}

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
template <class TLayout, class TSelector>
bool RemoveUnselectedInterfaceEntries(TLayout& layout, TSelector& sel)
{
//	iterate over all interfaces of the layout.
//	for each we'll create a new one, into which elements selected
//	elements will be inserted.
//	Finally we'll swap the content of the those interfaces.
//	if the interface is empty at the end of the operation, it will be
//	removed from the layout.
	bool retVal = false;
	
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
			bool interfaceChanged = false;
			Interface& interface = layout.interface(iiter);
		
		//	create a temporary interface and fill it with the selected entries
			Interface tInterface;
			
			for(ElemIter iter = interface.begin();
				iter != interface.end(); ++iter)
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
template <class TType, class TLayoutMap, class TSelector>
bool RemoveUnselectedInterfaceEntries(TLayoutMap& lm, TSelector& sel)
{
	typedef typename TLayoutMap::template Types<TType>::Map::iterator iterator;
	typedef typename TLayoutMap::template Types<TType>::Layout		Layout;

	bool retVal = false;
	for(iterator iter = lm.template layouts_begin<TType>();
		iter != lm.template layouts_end<TType>();)
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
template <class TLayout, class TSelectorIn, class TSelectorOut>
class SelectionCommPol : public ICommunicationPolicy<TLayout>
{
	public:
		typedef typename ICommunicationPolicy<TLayout>::Interface Interface;
		
	public:
		SelectionCommPol(TSelectorIn& selIn, TSelectorOut& selOut) :
			m_selIn(selIn), m_selOut(selOut)	{}
			
	///	iterates over the interface entries. Writes 1 for selected, 0 for unselected.
		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
			char zero = 0;
			char one = 1;
			
			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				if(m_selIn.is_selected(interface.get_element(iter)))
					buff.write(&one, sizeof(char));
				else
					buff.write(&zero, sizeof(char));
			}
			
			return true;
		}
		
	///	iterates over the interface entries. selects for 1, deselects for 0.
		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
			char tmp;
			for(typename Interface::iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				buff.read(&tmp, sizeof(char));
				if(tmp == 0)
					m_selOut.deselect(interface.get_element(iter));
				else
					m_selOut.select(interface.get_element(iter));
			}
			
			return true;
		}
		
	protected:
		TSelectorIn&	m_selIn;
		TSelectorOut&	m_selOut;
};

/// performs an allreduce and returns true if and only if all procs called the
/// function with bFlag = true
inline bool AllProcsTrue(bool bFlag,
						ProcessCommunicator comm = ProcessCommunicator())
{
//	local int bool flag
	int boolFlag = (bFlag) ? 1 : 0;

// 	local return flag
	int retBoolFlag;

//	all reduce
	comm.allreduce(&boolFlag, &retBoolFlag, 1, PCL_DT_INT, PCL_RO_LAND);

//	return global flag
	if(retBoolFlag != 0)
		return true;
	return false;
}

}//	end of namespace

#endif
