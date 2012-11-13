//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d18

#ifndef __H__LIB_GRID__DISTRIBUTION_UTIL_IMPL__
#define __H__LIB_GRID__DISTRIBUTION_UTIL_IMPL__

#include <iostream>
#include <vector>
#include <cassert>

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void SerializeDistributionLayoutInterfaces(BinaryBuffer& out, TLayout& layout,
											Grid& grid, AInt& aLocalInd)
{
	typedef typename TLayout::NodeType			Node;
	typedef typename TLayout::NodeVec			NodeVec;
	typedef typename TLayout::InterfaceEntry		InterfaceEntry;
	typedef typename PtrTypeToGeomObjBaseType<Node>::base_type	Elem;

	assert(grid.has_attachment<Elem>(aLocalInd));
	Grid::AttachmentAccessor<Elem, AInt> aaLocalInd(grid, aLocalInd);

	NodeVec nodes = layout.node_vec();

//	write source-proc and the number of levels
//	then for each level the number of interfaces
//	then for each interface the number of nodes and the local-ids of the nodes.
	int tmp;
	tmp = layout.get_source_proc();
	out.write((char*)&tmp, sizeof(int));

	tmp = (int) layout.num_levels();
	out.write((char*)&tmp, sizeof(int));
	
//	iterate through the levels of the layout
	for(size_t level = 0; level < layout.num_levels(); ++level)
	{
		typename TLayout::InterfaceMap& imap = layout.interface_map(level);
	//	write the number of interfaces for this level
		tmp = (int)imap.size();
		out.write((char*)&tmp, sizeof(int));
		
	//	iterate through the interfaces
		for(typename TLayout::InterfaceMap::iterator iter = imap.begin();
			iter != imap.end(); ++iter)
		{
		//	write the connected process-id
		//	if a process map is supplied perform a lookup
			int procID = iter->first.first;
			int oldTargetProc = iter->first.second;
			out.write((char*)&procID, sizeof(int));
			out.write((char*)&oldTargetProc, sizeof(int));
			
		//	write the number of entries that are contained in the interface
			typename TLayout::Interface& interface = iter->second;
			tmp = (int)interface.size();
			out.write((char*)&tmp, sizeof(int));
			
		//	write the interface-entries
		//	redirect the local index to the index where the associated
		//	element was serialized.
			InterfaceEntry entry;
			for(size_t i = 0; i < interface.size(); ++i){
				entry = DistributionInterfaceEntry(
									aaLocalInd[nodes[interface[i].local_id()]],
									interface[i].type());
				out.write((char*)&entry, sizeof(InterfaceEntry));
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
template <class TLayout>
void DeserializeDistributionLayoutInterfaces(TLayout& layout,
											BinaryBuffer& in)
{
	typedef typename TLayout::ProcPair	ProcPair;
//	read the source-proc and the number of levels
//	then for each level the number of interfaces
//	then for each interface the number of nodes and the local-ids of the nodes.
	int sourceProc;
	in.read((char*)&sourceProc, sizeof(int));
	layout.set_source_proc(sourceProc);

	int numLevels;
	in.read((char*)&numLevels, sizeof(int));
	layout.set_num_levels(numLevels);

//	iterate through the levels of the layout
	for(int level = 0; level < numLevels; ++level)
	{
	//	read the number of interfaces for this level
		int numInterfaces;
		in.read((char*)&numInterfaces, sizeof(int));

	//	iterate through the interfaces
		for(int i = 0; i < numInterfaces; ++i)
		{
		//	read the connected process-id
			int procID, oldTargetProc;
			in.read((char*)&procID, sizeof(int));
			in.read((char*)&oldTargetProc, sizeof(int));

		//	access the interface
			typename TLayout::Interface& interface =
						layout.interface(ProcPair(procID, oldTargetProc), level);

		//	read the number of entries that are contained in the interface
			int num;
			in.read((char*)&num, sizeof(int));
			interface.resize(num);

		//	write the interface-entries
			for(int j = 0; j < num; ++j){
				in.read((char*)&interface[j], sizeof(typename TLayout::InterfaceEntry));
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	DeserializeLayoutInterfaces
template <class TGeomObj, class TLayoutMap>
void DeserializeDistributionLayoutInterfaces(
								TLayoutMap& layoutMapOut,
								std::vector<TGeomObj*>& vGeomObjs,
								BinaryBuffer& in)
{
//	for conveniance
	//typedef typename TLayoutMap::mapped_type::Layout	TLayout;
	//typedef typename TLayoutMap::mapped_type::Interface	TInterface;
	typedef typename TLayoutMap::template
			Types<TGeomObj>::Layout::LevelLayout 	TLayout;
	typedef typename TLayout::Interface				TInterface;
	TLayout* pLayout = NULL;
	unsigned int lastLayoutKey = 0;
	TInterface* pInterface = NULL;
	DistributionInterfaceEntry entry;
	
//	read the source proc and the number of levels
//	then for each level the number of interfaces
//	then for each interface the number of nodes and the local-ids of the nodes.
	int sourceProc;
	in.read((char*)&sourceProc, sizeof(int));

	int numLevels;
	in.read((char*)&numLevels, sizeof(int));
	
//	iterate through the levels of the layout
	for(int level = 0; level < numLevels; ++level)
	{
	//	read the number of interfaces for this level
		int numInterfaces;
		in.read((char*)&numInterfaces, sizeof(int));
		
	//	iterate through the interfaces
		for(int i = 0; i < numInterfaces; ++i)
		{
		//	read the connected process-id
			int procID, oldTargetProc;
			in.read((char*)&procID, sizeof(int));
			in.read((char*)&oldTargetProc, sizeof(int));
					
		//	read the number of entries that are contained in the interface
			int numEntries;
			in.read((char*)&numEntries, sizeof(int));
			
		//	we'll set pLayout to NULL before we read the interface.
		//	This is important since this will cause the actualization of pInterface.
			pLayout = NULL;
			
		//	read the interface-entries
			for(int j = 0; j < numEntries; ++j)
			{
				in.read((char*)&entry, sizeof(DistributionInterfaceEntry));
			//	we're caching the last used layout to avoid too much lookups.
				if((!pLayout) || (lastLayoutKey != entry.type()))
				{
				//	get the matching layout
					pLayout = &layoutMapOut.template get_layout<TGeomObj>(entry.type()).
																layout_on_level(level);
					lastLayoutKey = entry.type();
				//	the interface has changed too
					pInterface = &pLayout->interface(procID);
				}
			//	copy the element into the interface
				pInterface->push_back(vGeomObjs[entry.local_id()]);
			}
		}
	}
}

}//	end of namespace

#endif
