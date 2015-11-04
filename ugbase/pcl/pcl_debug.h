#ifndef __H__PCL_DEBUG__
#define __H__PCL_DEBUG__

#include <iostream>
#include "pcl_base.h"
#include "pcl_methods.h"
#include "pcl_communication_structs.h"
#include "pcl_interface_communicator.h"
#include "pcl_process_communicator.h"

////////////////////////////////////////////////////////////////////////
///	this allows us to print messages to the users terminal
/** Wrapper for UG_LOG.*/
#define PCLLOG(msg) UG_LOG(msg)


namespace pcl
{

/// \addtogroup pcl
/// \{

////////////////////////////////////////////////////////////////////////
//	LogLayoutStructure
///	Logs the internals of a layout.
/**
 * Supported layouts are pcl::SingleLevelLayout and pcl::MultiLevelLayout.
 */
template <class TLayout>
void LogLayoutStructure(TLayout& layout, const char* prefix = "")
{
	using namespace std;
	
	typedef typename TLayout::Interface Interface;
	typedef typename TLayout::iterator InterfaceIter;

	PCLLOG(prefix << "-- PCL_DEBUG: layout-structure --\n");
	PCLLOG(prefix << "---- num_levels: " << layout.num_levels() << endl);
	
	for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl)
	{
		PCLLOG(prefix << "---- interfaces on level " << lvl << " (proc id, size): ");
				
		for(InterfaceIter iiter = layout.begin(lvl);
			iiter != layout.end(lvl); ++iiter)
		{
			Interface& interface = layout.interface(iiter);
			PCLLOG("(" << layout.proc_id(iiter) << ", " << interface.size() << "), ");
		}
		PCLLOG(endl);
	}
}

////////////////////////////////////////////////////////////////////////
//	LogLayoutMapStrucure
///	Logs the internals of a layout-map for a given type.
/**
 * Supported layouts are pcl::SingleLevelLayout and pcl::MultiLevelLayout.
 */
template <class TType, class TLayoutMap>
void LogLayoutMapStructure(TLayoutMap& lm)
{
	using namespace std;
	
	typedef typename TLayoutMap::template Types<TType>::Map::iterator iterator;
	typedef typename TLayoutMap::template Types<TType>::Layout		Layout;

	PCLLOG("-- PCL_DEBUG: layout-map-structure --\n");
	
	for(iterator iter = lm.template layouts_begin<TType>();
		iter != lm.template layouts_end<TType>(); ++iter)
	{
	//	get the key
		PCLLOG("---- has key: " << iter->first << endl);
		
	//	log the layout structure
		Layout& layout = iter->second;
		LogLayoutStructure(layout, "----");
	}
}

// end group pcl
/// \}

}//	end of namespace

#endif
