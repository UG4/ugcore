//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d18

#ifndef __H__LIB_GRID__DISTRIBUTION_UTIL_IMPL__
#define __H__LIB_GRID__DISTRIBUTION_UTIL_IMPL__

#include <iostream>
#include <vector>

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerializeLayoutInterfaces(std::ostream& out, TLayout& layout,
								std::vector<int>* pProcessMap)
{
//	write the number of levels
//	then for each level the number of interfaces
//	then for each interface the number of nodes and the local-ids of the nodes.
	int tmp;
	
	tmp = (int) layout.num_levels();
	out.write((char*)&tmp, sizeof(int));
	
//	iterate through the levels of the layout
	for(size_t level = 0; level < layout.num_levels(); ++level)
	{
		pcl::InterfaceMap& imap = layout.interface_map(level);
	//	write the number of interfaces for this level
		tmp = (int)imap.size();
		out.write((char*)&tmp, sizeof(int));
		
	//	iterate through the interfaces
		for(pcl::InterfaceMap::iterator iter = imap.begin();
			iter != imap.end(); ++iter)
		{
		//	write the connected process-id
		//	if a process map is supplied perform a lookup
			int procID = iter->first;
			if(pProcessMap)
			{
				assert(pProcessMap->size() > procID && "process-map to small.");
				if(pProcessMap->size() > procID)
					procID = (*pProcessMap)[procID];
			}
			out.write((char*)&procID, sizeof(int));
			
		//	write the number of entries that are contained in the interface
			pcl::Interface& interface = iter->second;
			tmp = (int)interface.size();
			out.write((char*)&tmp, sizeof(int));
			
		//	write the interface-entries
			for(size_t i = 0; i < interface.size(); ++i)
				out.write((char*)&interface[i], sizeof(pcl::InterfaceEntry));
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	DeserializeLayoutInterfaces
template <class TLayout>
void DeserializeLayoutInterfaces(TLayout& layoutOut,
								std::istream& in)
{
//	for conveniance
	TLayout& layout = layoutOut;
	
//	read the number of levels
//	then for each level the number of interfaces
//	then for each interface the number of nodes and the local-ids of the nodes.
	int numLevels;
	in.read((char*)&numLevels, sizeof(int));
	layout.set_num_levels((size_t)numLevels);
	
//	iterate through the levels of the layout
	for(size_t level = 0; level < layout.num_levels(); ++level)
	{
		pcl::InterfaceMap& imap = layout.interface_map(level);
	//	read the number of interfaces for this level
		int numInterfaces;
		in.read((char*)&numInterfaces, sizeof(int));
		
	//	iterate through the interfaces
		for(int i = 0; i < numInterfaces; ++i)
		{
		//	read the connected process-id
			int procID;
			in.read((char*)&procID, sizeof(int));
		
		//	get the interface
			pcl::Interface& interface = layout.interface(procID, level);
			
		//	read the number of entries that are contained in the interface
			int numEntries;
			in.read((char*)&numEntries, sizeof(int));
			
		//	resize the interface
			interface.resize(numEntries);
			
		//	read the interface-entries
			for(size_t j = 0; j < interface.size(); ++j)
				in.read((char*)&interface[i], sizeof(pcl::InterfaceEntry));
		}
	}
}

}//	end of namespace

#endif
