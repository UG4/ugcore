//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#ifndef __H__LIB_GRID__DISTRIBUTION_UTIL__
#define __H__LIB_GRID__DISTRIBUTION_UTIL__

#include <vector>
#include <map>
#include "../parallel_grid_layout.h"
#include "../distributed_grid.h"
#include "lib_grid/tools/selector_multi_grid.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_util
///	@{

////////////////////////////////////////////////////////////////////////
///	The interface entry holds a local id and the type of the entry.
/** Note that the type is restricted to 4 bytes currently.
 *
 * Note that this is pointed out in the documentation for InterfaceNodeTypes,
 * too.*/
struct DistributionInterfaceEntry
{
	DistributionInterfaceEntry() :
		m_localID(0),
		m_type(0)	{}
		
	DistributionInterfaceEntry(int nLocalID, int nType) :
		m_localID(nLocalID),
		m_type(nType)	{}

	unsigned int local_id()	{return m_localID;}
	unsigned int type()		{return m_type;}

	private:
		unsigned int m_localID : 28;
		unsigned int m_type  	: 4;
};

////////////////////////////////////////////////////////////////////////
///	This small helper gives us information about a node move during redistribution
struct RedistributionNodeTransferInfo{
	RedistributionNodeTransferInfo() :
		srcProc(-1), targetProc(-1), bMove(false)	{}
	RedistributionNodeTransferInfo(int nSrcProc, int nTargetProc, bool nMove) :
		srcProc(nSrcProc), targetProc(nTargetProc), bMove(nMove)	{}

	int srcProc;
	int targetProc;
	bool bMove;
};

///	Attachment for a vector of RedistributionNodeTransferInfo
typedef Attachment<std::vector<RedistributionNodeTransferInfo> >
		ARedistributionNodeTransferInfoVec;

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
///	a proc-pair combines the target-proc (first) and the old target proc(second).
	typedef std::pair<int, int>				ProcPair;
///	an interface-map is a list of interfaces, each associated with a new and old process id.
/**	The old process-id is only used during redistribution.*/
	typedef std::map<ProcPair, Interface>	InterfaceMap;
///	a list of interface-maps. Required for multilevel / hierarchical approaches.
	typedef std::vector<InterfaceMap>		InterfaceMapVec;

///	the type of the nodes
	typedef TNode	NodeType;
///	a vector that holds nodes.
	typedef std::vector<TNode>				NodeVec;	
	
///	sets m_srcProcID to pcl::GetProcRank.
	DistributionNodeLayout() : m_srcProcID(pcl::GetProcRank())	{}
	virtual ~DistributionNodeLayout()	{}

//	some methods
/*
///	set process id. Use with care to not invalidate other ProcessLayouts.
	inline void set_proc_id(int procID)						{m_procID = procID;}
	
///	returns the process id.
	inline int get_proc_id() const							{return m_procID;}
*/	
///	returns a reference to the vector that holds the nodes.
	inline NodeVec& node_vec()								{return m_vNodes;}

///	returns a const reference to the vector that holds the nodes.
	inline const NodeVec& node_vec() const					{return m_vNodes;}
		
///	returns the interface to the given process on the given level.
/**	if you don't specify a level, level = 0 will be used.
 * Since the interface is internally associated with a ProcPair,
 * -1 is used for the old target proc.*/
	inline Interface& interface(int procID, size_t level = 0)	{return interface_map(level)[ProcPair(procID, -1)];}

///	returns the interface to the given process and the given old process on the given level.
	inline Interface& interface(ProcPair procPair, size_t level = 0)	{return interface_map(level)[procPair];}

///	returns true if the interface already exists
/** Since the interface is internally associated with a ProcPair,
 * -1 is used for the old target proc.*/
	inline bool has_interface(int procID, size_t level = 0)
	{
		InterfaceMap& m = interface_map(level);
		return (m.find(ProcPair(procID, -1)) != m.end());
	}

///	returns true if the interface already exists
	inline bool has_interface(ProcPair procPair, size_t level = 0)
	{
		InterfaceMap& m = interface_map(level);
		return (m.find(procPair) != m.end());
	}
	
///	returns the interface-map for the given level.
/**	if you don't specify a level, level = 0 will be used.*/
	inline InterfaceMap& interface_map(size_t level = 0)		{if(level >= m_vInterfaceMaps.size()) m_vInterfaceMaps.resize(level + 1); return m_vInterfaceMaps[level];}
	
///	sets the number of levels.
/**	Setting the number of levels is optional. Increases performance for \#levels > 1.*/
	void set_num_levels(size_t num)							{m_vInterfaceMaps.resize(num);}
	
///	returns the number of levels.
	inline size_t num_levels() const						{return m_vInterfaceMaps.size();}
	
///	the source-proc defines which process created the layout
/**	\{	*/
	inline void set_source_proc(int procID)		{m_srcProcID = procID;}
	inline int get_source_proc()				{return m_srcProcID;}
/**	\}	*/

protected:
	int				m_srcProcID;
	NodeVec			m_vNodes;
	InterfaceMapVec	m_vInterfaceMaps;
};


////////////////////////////////////////////////////////////////////////
///	Used during grid redistribution.
/**	It holds an instance of DistributionNodeLayout<TNode> and a vector
 * of global ids, which store a globalID for each node in the
 * distribution layout.
 */
template <class TNode>
struct RedistributionNodeLayout : public DistributionNodeLayout<TNode>
{
///	this array is of the same size as m_distLayout.node_vec()
	std::vector<GeomObjID>	m_globalIDs;
};


////////////////////////////////////////////////////////////////////////
//	typedefs
typedef DistributionNodeLayout<VertexBase*>	DistributionVertexLayout;
typedef DistributionNodeLayout<EdgeBase*>	DistributionEdgeLayout;
typedef DistributionNodeLayout<Face*>		DistributionFaceLayout;
typedef DistributionNodeLayout<Volume*>		DistributionVolumeLayout;

typedef RedistributionNodeLayout<VertexBase*>	RedistributionVertexLayout;
typedef RedistributionNodeLayout<EdgeBase*>		RedistributionEdgeLayout;
typedef RedistributionNodeLayout<Face*>			RedistributionFaceLayout;
typedef RedistributionNodeLayout<Volume*>		RedistributionVolumeLayout;


////////////////////////////////////////////////////////////////////////
//	CreateDistributionLayouts
///	Creates distribution layouts for vertices, edges, faces and volumes
/**
 * Given a MultiGrid and a SubsetHandler, this method creates distribution
 * layouts for vertices, edges, ...
 * Those layouts can then be used to distribute a grid onto different
 * processes. Please note that those layouts are not used to perform
 * communication later on. Their sole purpose is to help to distribute
 * a grid.
 *
 * For each subset a separate distribution-layout is created. That means 
 * that i.e. vertexLayoutsOut[k] holds the vertex-layout for the k-th subset.
 *
 * If you pass a pointer to a valid selector (which is registered at mg),
 * the selector will be used for internal calculations. The only reason
 * for this parameter is a speed increase.
 * You shouldn't assume anything about the content of pSel after the
 * method finished.
 *
 * TVertexDistributionLayout has to be either DistributionNodeLayout<VertexBase>
 * or RedistributionLayout<VertexBase>. Analogous for the other types.
 *
 * If CreateDistributionLayouts is used during redistribution, a
 * DistributedGridManager should be specified, since existing interface
 * entries will have to be ignored during assignment of distribution-
 * interface entries. Only specify pDistGridMgr for redistribution.
 */
template <class TVertexDistributionLayout, class TEdgeDistributionLayout,
		  class TFaceDistributionLayout, class TVolumeDistributionLayout>
void CreateDistributionLayouts(
						std::vector<TVertexDistributionLayout>& vertexLayoutsOut,
						std::vector<TEdgeDistributionLayout>& edgeLayoutsOut,
						std::vector<TFaceDistributionLayout>& faceLayoutsOut,
						std::vector<TVolumeDistributionLayout>& volumeLayoutsOut,
						MultiGrid& mg, SubsetHandler& sh,
						bool distributeGenealogy,
						bool createVerticalInterfaces = false,
						MGSelector* pSel = NULL,
						DistributedGridManager* pDistGridMgr = NULL,
						std::vector<int>* processMap = NULL);

////////////////////////////////////////////////////////////////////////
//	CreateRedistributionLayouts
///	Creates distribution layouts for vertices, edges, faces and volumes
/**
 * Given a DistributedGridManager and a SubsetHandler, this method
 * creates redistribution layouts for vertices, edges, ...
 * Those layouts can then be used to redistribute a distributed grid onto
 * different processes.
 * Please note that those layouts are not used to perform
 * communication later on. Their sole purpose is to help to redistribute
 * a grid.
 *
 * For each subset a separate redistribution-layout is created. That means
 * that i.e. vertexLayoutsOut[k] holds the vertex-layout for the k-th subset.
 *
 * If nodes that lie in an existing interface are to be distributed to
 * another process, interfaces are automatically generated, which will
 * connect the nodes on the new process to associated nodes on the old
 * neighbor or - if those nodes move too, to their new destinations.
 * Note that if a processMap was specified, the interfaces are build
 * with respect to the processMap.
 *
 * If you pass a pointer to a valid selector (which is registered at mg),
 * the selector will be used for internal calculations. This increases
 * performance and avoids unnecessary data allocation.
 *
 * Make sure that global ids are generated and that aGlobalID points to
 * the attachment which holds them (aGeomObjID by default). Those are
 * only required at this point, since they are written into the
 * redistribution layouts.
 */
void CreateRedistributionLayouts(
					std::vector<RedistributionVertexLayout>& vertexLayoutsOut,
					std::vector<RedistributionEdgeLayout>& edgeLayoutsOut,
					std::vector<RedistributionFaceLayout>& faceLayoutsOut,
					std::vector<RedistributionVolumeLayout>& volumeLayoutsOut,
					DistributedGridManager& distGridMgr,
					SubsetHandler& shPartition, bool distributeGenealogy,
					bool createVerticalInterfaces,
					MGSelector* pSel = NULL, std::vector<int>* processMap = NULL,
					AGeomObjID& aGlobalID = aGeomObjID);
						
////////////////////////////////////////////////////////////////////////
//	SerializeGridAndDistributionLayouts
/**
 * Writes the elements of a grid, which are referenced by the given
 * gridLayout to a binary stream.
 * You may pass a selector via pSel, which increases performance of this
 * method. After this method finished pSel will contain all elements
 * that have been written to the stream (in the same order as they were
 * written).
 *
 * This method first uses SerializeMultiGridElements to write the multi-grid
 * to the stream and calls SerializeDistributionLayoutInterfaces for the
 * vertex, edge, face and volume layouts afterwards.
 *
 * Make sure to only pass a process-map if your layouts do not already
 * reference the final processes.
 *
 * \todo	layouts should already reference the real target procs -
 * 			pProcessMap should thus be removed.
 */
void SerializeGridAndDistributionLayouts(
							std::ostream& out, MultiGrid& mg,
							DistributionVertexLayout& vrtLayout,
							DistributionEdgeLayout& edgeLayout,
							DistributionFaceLayout& faceLayout,
							DistributionVolumeLayout& volLayout,
							AInt& aLocalIndVRT, AInt& aLocalIndEDGE,
							AInt& aLocalIndFACE, AInt& aLocalIndVOL,
							MGSelector* pSel = NULL);

////////////////////////////////////////////////////////////////////////
//	SerializeDistributionLayoutInterfaces
/**	Make sure to only pass a process-map if your layouts do not already
 * reference the final processes.
 *
 * aLocalInd has to be attached to the elements in the distrbution-layout.
 * They have to store the index for each element, where it was written to
 * the stream during serialization of the grid.
 * This is important, since the nodes in the node-vec of layout are not
 * necessarily ordered by their mg-level. This however is the order in which
 * elements are serialized. We thus have to redirect the indices stored in
 * the distribution-interfaces to the real indices of the elements.
 *
 * \todo	layouts should already reference the real target procs -
 * 			pProcessMap should thus be removed.
 */
template <class TLayout>
void SerializeDistributionLayoutInterfaces(std::ostream& out, TLayout& layout,
											Grid& grid, AInt& aLocalInd);

////////////////////////////////////////////////////////////////////////
///	deserializes a DistributionLayoutInterface.
template <class TLayout>
void DeserializeDistributionLayoutInterfaces(TLayout& layout,
											std::istream& in);

////////////////////////////////////////////////////////////////////////
//	DeserializeDistributionLayoutInterfaces
///	Reads serialized DistributionLayout data and creates GridLayoutMaps on the fly.
/**
 * TLayoutMap has to be compatible with an
 * std::map<int, ParallelELEMENTLayout>, where ELEMENT can be either
 * Vertex, Edge, Face or Volume.
 */
template <class TGeomObj, class TLayoutMap>
void DeserializeDistributionLayoutInterfaces(
									TLayoutMap& layoutMapOut,
									std::vector<TGeomObj*>& vGeomObjs,
									std::istream& in);

////////////////////////////////////////////////////////////////////////
//	DeserializeGridAndDistributionLayouts
void DeserializeGridAndDistributionLayouts(
									MultiGrid& mgOut,
									GridLayoutMap& gridLayoutOut,
									std::istream& in);

////////////////////////////////////////////////////////////////////////
///	selects all elements in a DistributionLayout into a selector
/**	TSelector should be either Selector or MGSelector.
 *	TDistributionLayout has to be a specialization of
 *	DistributionNodeLayout for either VertexBase, EdgeBase,
 *	Face or Volume.
 */
template <class TSelector, class TDistributionLayout>
void SelectNodesInLayout(TSelector& sel, TDistributionLayout& layout)
{
	typename TDistributionLayout::NodeVec& nodes = layout.node_vec();
	for(size_t i = 0; i < nodes.size(); ++i)
		sel.select(nodes[i]);
}

////////////////////////////////////////////////////////////////////////
///	marks all elements in a DistributionLayout
/**	TDistributionLayout has to be a specialization of
 * DistributionNodeLayout for either VertexBase, EdgeBase,
 * Face or Volume.
 * Make sure that this method is called between calls to g.begin_marking()
 * and g.end_marking().
 */
template <class TDistributionLayout>
void MarkNodesInLayout(Grid& g, TDistributionLayout& layout)
{
	typename TDistributionLayout::NodeVec& nodes = layout.node_vec();
	for(size_t i = 0; i < nodes.size(); ++i)
		g.mark(nodes[i]);
}

////////////////////////////////////////////////////////////////////////
///	marks all elements in a series of DistributionLayouts
/**	TIterator has to be a stl-compatible iterator type and
 * TIterator::value_type has to be an instance of a specialization of
 * DistributionNodeLayout for either VertexBase, EdgeBase, Face or Volume.
 *
 * Make sure that this method is called between calls to g.begin_marking()
 * and g.end_marking().
 */
template <class TIterator>
void MarkNodesInLayouts(Grid& g, TIterator layoutsBegin, TIterator layoutsEnd)
{
	for(TIterator iter = layoutsBegin; iter != layoutsEnd; iter++)
		MarkNodesInLayout(g, *iter);
}

///	Counts how many entries with the given type are contained in the given interface.
size_t NumEntriesOfTypeInDistributionInterface(unsigned int type,
			std::vector<DistributionInterfaceEntry>& interface);

///	checks whether the interconnections between the layouts are consistent.
/**	Note that the i-th entry of distLayouts is considered to
 * be the layout of process i (eventually redirected by procMap).
 */
template <class TDistLayout>
bool TestDistributionLayouts(std::vector<TDistLayout>& distLayouts,
							 int* procMap = NULL);

///	Currently simply outputs the connections in the layouts.
/**	Note that the i-th entry of distLayouts is considered to
 * be the layout of process i (eventually redirected by procMap).
 */
template <class TDistLayout>
bool TestRedistributionLayouts(std::vector<TDistLayout>& distLayouts,
								int* procMap = NULL);

///	Simply prints the contents of the given layouts.
template <class TDistLayout>
bool PrintRedistributionLayouts(std::vector<TDistLayout>& distLayouts);
/// @}
}//	end of namespace

////////////////////////////////
//	include implementation
#include "distribution_util_impl.hpp"

#endif
