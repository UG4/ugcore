//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d18

#ifndef __H__PCL__DISTRIBUTION_NODE_LAYOUT__
#define __H__PCL__DISTRIBUTION_NODE_LAYOUT__

#include <vector>
#include <map>
#include "parallel_grid_layout.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	The interface entry holds a local id and the type of the entry.
struct DistributionInterfaceEntry
{
	DistributionInterfaceEntry()	{}
	DistributionInterfaceEntry(int nLocalID, int nType) : localID(nLocalID), type(nType)	{}
	
	int localID : 28;
	int type  	: 4;
};

////////////////////////////////////////////////////////////////////////
///	Holds nodes and interfaces. Used during distribution only.
/**
 * This class is used in the process of grid-distribution to assemble
 * the nodes and interfaces that belong to the different processes.
 * It is however not used during parallel communication.
 */
template <class TNode>
struct DistributionNodeLayout
{
// typedefs
///	the interface entry type
	typedef DistributionInterfaceEntry		InterfaceEntry;
///	an interface consists of a list of local ids.
	typedef std::vector<InterfaceEntry>		Interface;
///	an interface-map is a list of interfaces, each associated with a process id.
	typedef std::map<int, Interface>		InterfaceMap;
///	a list of interface-maps. Required for multilevel / hierarchical approaches.
	typedef std::vector<InterfaceMap>		InterfaceMapVec;

///	the type of the nodes
	typedef TNode	NodeType;
///	a vector that holds nodes.
	typedef std::vector<TNode>				NodeVec;	
	
//	some methods
/*
///	set process id. Use with care to not invalidate other ProcessLayouts.
	inline void set_proc_id(int procID)						{m_procID = procID;}
	
///	returns the process id.
	inline int get_proc_id() const							{return m_procID;}
*/	
///	returns a reference to the vector that holds the nodes.
	inline NodeVec& node_vec()								{return m_vNodes;}
	
///	returns the interface to the given process on the given level.
	/**	if you don't specify a level, level = 0 will be used.*/
	inline Interface& interface(int procID, int level = 0)	{return interface_map(level)[procID];}
	
///	returns the interface-map for the given level.
	/**	if you don't specify a level, level = 0 will be used.*/
	inline InterfaceMap& interface_map(int level = 0)		{if(level >= m_vInterfaceMaps.size()) m_vInterfaceMaps.resize(level + 1); return m_vInterfaceMaps[level];}
	
///	sets the number of levels.
	/**	Setting the number of levels is optional. Increases performance for #levels > 1.*/
	void set_num_levels(size_t num)							{m_vInterfaceMaps.resize(num);}
	
///	returns the number of levels.
	inline size_t num_levels() const						{return m_vInterfaceMaps.size();}
	
protected:
	//int				m_procID;
	NodeVec			m_vNodes;
	InterfaceMapVec	m_vInterfaceMaps;
};

}//	end of namespace pcl

#endif
