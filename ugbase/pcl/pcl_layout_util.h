#ifndef __H__PCL__pcl_layout_util__
#define __H__PCL__pcl_layout_util__

#include <vector>
#include "common/util/hash.h"
#include "pcl_communication_structs.h"


namespace pcl{

/// \addtogroup pcl
/// \{

///	removes all empty interfaces from the given layout.
template <class TLayout>
void RemoveEmptyInterfaces(TLayout& layout)
{
	typedef typename TLayout::iterator TInterfaceIter;
	typedef typename TLayout::Interface TInterface;

	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
		for(TInterfaceIter iter = layout.begin(lvl); iter != layout.end(lvl);)
		{
			TInterface& intfc = layout.interface(iter);
			if(intfc.empty())
				iter = layout.erase(iter, lvl);
			else
				++iter;
		}
	}
}

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
						   const TLayout& layout)
{
	typedef typename TLayout::Interface Interface;
	typedef typename TLayout::Element TElem;

//	clear the return value
	elemsOut.clear();

//	we'll use a hash to make sure that each element only exists once
	ug::Hash<TElem, int> hash(layout.num_interface_elements());
	hash.reserve(layout.num_interface_elements());

//	iterate over all interfaces
	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
		for(typename TLayout::const_iterator interfaceIter = layout.begin(lvl);
			interfaceIter != layout.end(lvl); ++interfaceIter)
		{
		//	iterate over the entries of the interface
			const Interface& interface = layout.interface(interfaceIter);
			for(typename Interface::const_iterator iter = interface.begin();
				iter != interface.end(); ++iter)
			{
			//	check whether the entry already exists in the hash
				if(!hash.has_entry(interface.get_element(iter))){
				//	we don't care about the value
					hash.insert(interface.get_element(iter), 0);
					elemsOut.push_back(interface.get_element(iter));
				}
			}
		}
	}
}

// end group pcl
/// \}

}// end of namespace

#endif
