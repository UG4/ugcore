//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#include "lib_grid/lib_grid.h"
#include "distribution_util.h"
#include "common/util/stream_pack.h"

namespace pcl
{

enum InterfaceNodeTypes
{
	INT_UNKNOWN = 0,
	INT_MASTER = 1,
	INT_SLAVE = 3,
	INT_LINK = 7
};

struct InterfaceEntry
{
	InterfaceEntry()	{}
	InterfaceEntry(int nLocalID, int nType) : localID(nLocalID), type(nType)	{}
	
	int localID : 28;
	int type  	: 4;
};

///	an interface consists of a list of local ids.
typedef std::vector<InterfaceEntry>		Interface;
///	an interface-map is a list of interfaces, each associated with a process id.
typedef std::map<int, Interface>		InterfaceMap;
///	a list of interface-maps. Required for multilevel / hierarchical approaches.
typedef std::vector<InterfaceMap>		InterfaceMapVec;
///	allows iteration over (procID, Interface)-pairs.
typedef InterfaceMap::iterator			InterfaceIterator;

//	move this class to pcl!
template <class TNode>
struct ParallelNodeLayout
{
//	some typedefs
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

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	specialization for lib_grid.
typedef pcl::ParallelNodeLayout<VertexBase*>	ParallelVertexLayout;
typedef pcl::ParallelNodeLayout<EdgeBase*>		ParallelEdgeLayout;
typedef pcl::ParallelNodeLayout<Face*>			ParallelFaceLayout;
typedef pcl::ParallelNodeLayout<Volume*>		ParallelVolumeLayout;

template <class TNodeLayout, class TIterator, class TAIntAccessor>
void AddElementsToLayout(std::vector<TNodeLayout>& layouts,
							int layoutIndex, 
							TIterator nodesBegin, TIterator nodesEnd,
							TAIntAccessor& aaFirstLayout,
							TAIntAccessor& aaFirstProcLocalInd,
							int level = 0)
{
	TNodeLayout& layout = layouts[layoutIndex];	
	
	for(TIterator iter = nodesBegin; iter != nodesEnd; ++iter)
	{
		typename TNodeLayout::NodeType node = *iter;
		int masterLayoutIndex = aaFirstLayout[node];
		if(masterLayoutIndex == -1)
		{
		//	the node has been encountered for the first time.
		//	add it to the layout and set aaFirstLayout and
		//	aaLocalIndex accordingly.
			aaFirstLayout[node] = layoutIndex;
			aaFirstProcLocalInd[node] = (int)layout.node_vec().size();
			layout.node_vec().push_back(node);
		}
		else
		{
		//	this helps debugging: if you assume that no interfaces will be build
		//	during the execution of this method, you may pass a level of -1.
			assert(level != -1 && "bad level index.");
			
		//	the node has already been added to another layout.
		//	add it to the new layout and create interface entries
		//	on both sides.
			int localMasterID = aaFirstProcLocalInd[node];
			int localID = (int)layout.node_vec().size();
			TNodeLayout& masterLayout = layouts[masterLayoutIndex];
			
		//	access the interfaces
			pcl::Interface& masterInterface = masterLayout.interface(layoutIndex, level);
			pcl::Interface& slaveInterface = layout.interface(masterLayoutIndex, level);
			masterInterface.push_back(pcl::InterfaceEntry(localMasterID, pcl::INT_MASTER));
			slaveInterface.push_back(pcl::InterfaceEntry(localID, pcl::INT_SLAVE));
		}
	}
}


void CreateGridLayouts(	std::vector<ParallelVertexLayout>& vertexLayoutsOut,
						std::vector<ParallelEdgeLayout>& edgeLayoutsOut,
						std::vector<ParallelFaceLayout>& faceLayoutsOut,
						std::vector<ParallelVolumeLayout>& volumeLayoutsOut,
						MultiGrid& mg, SubsetHandler& sh,
						MGSelector* pSel = NULL)
{
//	initialize a selector.
	MGSelector tmpSel;
	if(!pSel)
	{
		tmpSel.assign_grid(mg);
		pSel = &tmpSel;
	}
	MGSelector& msel = *pSel;
	
//	resize and clear the layouts
	vertexLayoutsOut = std::vector<ParallelVertexLayout>(sh.num_subsets());
	edgeLayoutsOut = std::vector<ParallelEdgeLayout>(sh.num_subsets());
	faceLayoutsOut = std::vector<ParallelFaceLayout>(sh.num_subsets());
	volumeLayoutsOut = std::vector<ParallelVolumeLayout>(sh.num_subsets());
	
//	attach first-proc-indices and local-ids to the elements of the grid.
	AInt aFirstProc;
	AInt aFirstProcLocalInd;
	mg.attach_to_vertices(aFirstProc);
	mg.attach_to_edges(aFirstProc);
	mg.attach_to_faces(aFirstProc);
	mg.attach_to_volumes(aFirstProc);
	mg.attach_to_vertices(aFirstProcLocalInd);
	mg.attach_to_edges(aFirstProcLocalInd);
	mg.attach_to_faces(aFirstProcLocalInd);
	mg.attach_to_volumes(aFirstProcLocalInd);

//	the attachment accessors
	Grid::VertexAttachmentAccessor<AInt> aaFirstProcVRT(mg, aFirstProc);
	Grid::EdgeAttachmentAccessor<AInt> aaFirstProcEDGE(mg, aFirstProc);
	Grid::FaceAttachmentAccessor<AInt> aaFirstProcFACE(mg, aFirstProc);
	Grid::VolumeAttachmentAccessor<AInt> aaFirstProcVOL(mg, aFirstProc);
	Grid::VertexAttachmentAccessor<AInt> aaFirstProcLocalIndVRT(mg, aFirstProcLocalInd);
	Grid::EdgeAttachmentAccessor<AInt> aaFirstProcLocalIndEDGE(mg, aFirstProcLocalInd);
	Grid::FaceAttachmentAccessor<AInt> aaFirstProcLocalIndFACE(mg, aFirstProcLocalInd);
	Grid::VolumeAttachmentAccessor<AInt> aaFirstProcLocalIndVOL(mg, aFirstProcLocalInd);

//	initialise first-proc attachments
	SetAttachmentValues(aaFirstProcVRT, mg.vertices_begin(), mg.vertices_end(), -1);
	SetAttachmentValues(aaFirstProcEDGE, mg.edges_begin(), mg.edges_end(), -1);
	SetAttachmentValues(aaFirstProcFACE, mg.faces_begin(), mg.faces_end(), -1);
	SetAttachmentValues(aaFirstProcVOL, mg.volumes_begin(), mg.volumes_end(), -1);
	
//	iterate through the subsets and and create the packs.
//	we have to do this in two steps to make sure that all
//	elements are masters on the processes that they are assigned to
//	in the subsethandler.

//	step 1: add the elements to the groups to which they were assigned.
	for(uint i = 0; i < sh.num_subsets(); ++i)
	{
	//	the level is ignored since it won't be used in this phase.
	//	by passing -1 we can assert that no interface is accessed.
		AddElementsToLayout(vertexLayoutsOut, i,
							sh.begin<VertexBase>(i), sh.end<VertexBase>(i),
							aaFirstProcVRT, aaFirstProcLocalIndVRT, -1);
		AddElementsToLayout(edgeLayoutsOut, i,
							sh.begin<EdgeBase>(i), sh.end<EdgeBase>(i),
							aaFirstProcEDGE, aaFirstProcLocalIndEDGE, -1);
		AddElementsToLayout(faceLayoutsOut, i,
							sh.begin<Face>(i), sh.end<Face>(i),
							aaFirstProcFACE, aaFirstProcLocalIndFACE, -1);
		AddElementsToLayout(volumeLayoutsOut, i,
							sh.begin<Volume>(i), sh.end<Volume>(i),
							aaFirstProcVOL, aaFirstProcLocalIndVOL, -1);
	}
	
//	step 2: add all the associated elements to the distribution groups, which
//			have not already been assigned.
	for(uint i = 0; i < sh.num_subsets(); ++i)
	{
		msel.clear_selection();		
//TODO: overlap can be easily handled here! simply increase the selection.
//		eventually we first would have to select all associated elements.

	//	the hierarchy has to be complete. make sure the whole geneology
	//	is selected. By passing true, all associated elements of lower
	//	dimension will be selected, too.
		SelectAssociatedGenealogy(msel, true);
		
	//	make sure that we won't add elements twice.
		msel.deselect(sh.begin<VertexBase>(i), sh.end<VertexBase>(i));
		msel.deselect(sh.begin<EdgeBase>(i), sh.end<EdgeBase>(i));
		msel.deselect(sh.begin<Face>(i), sh.end<Face>(i));
		msel.deselect(sh.begin<Volume>(i), sh.end<Volume>(i));
		
	//	add the elements to the groups
	//	since interfaces are generated during this step, we have to take
	//	care of the levels.
		for(uint level = 0; level < msel.num_levels(); ++level)
		{
			AddElementsToLayout(vertexLayoutsOut, i,
								msel.begin<VertexBase>(level), msel.end<VertexBase>(level),
								aaFirstProcVRT, aaFirstProcLocalIndVRT, level);
			AddElementsToLayout(edgeLayoutsOut, i,
								msel.begin<EdgeBase>(level), msel.end<EdgeBase>(level),
								aaFirstProcEDGE, aaFirstProcLocalIndEDGE, level);
			AddElementsToLayout(faceLayoutsOut, i,
								msel.begin<Face>(level), msel.end<Face>(level),
								aaFirstProcFACE, aaFirstProcLocalIndFACE, level);
			AddElementsToLayout(volumeLayoutsOut, i,
								msel.begin<Volume>(level), msel.end<Volume>(level),
								aaFirstProcVOL, aaFirstProcLocalIndVOL, level);
		}
	}
	
//	The layouts are now complete.
//	we're done in here.

//	clean up
	mg.detach_from_vertices(aFirstProc);
	mg.detach_from_edges(aFirstProc);
	mg.detach_from_faces(aFirstProc);
	mg.detach_from_volumes(aFirstProc);
	mg.detach_from_vertices(aFirstProcLocalInd);
	mg.detach_from_edges(aFirstProcLocalInd);
	mg.detach_from_faces(aFirstProcLocalInd);
	mg.detach_from_volumes(aFirstProcLocalInd);
}


////////////////////////////////////////////////////////////////////////
template <class TSelector, class TLayout>
static
void SelectNodesInLayout(TSelector& sel, TLayout& layout)
{
	typename TLayout::NodeVec& nodes = layout.node_vec();
	for(size_t i = 0; i < nodes.size(); ++i)
		sel.select(nodes[i]);
}

////////////////////////////////////////////////////////////////////////
template <class TLayout, class TAIntAccessor>
void SerializeLayout(std::ostream& out, TLayout& layout,
					TAIntAccessor& aaInt, std::vector<int>& processMap)
{
//...
}

////////////////////////////////////////////////////////////////////////
void SerializeGridAndLayouts(std::ostream& out, MultiGrid& mg,
						std::vector<int>& processMap,
						ParallelVertexLayout& vrtLayout,
						ParallelEdgeLayout& edgeLayout,
						ParallelFaceLayout& faceLayout,
						ParallelVolumeLayout& volLayout,
						AInt& aLocalIndVRT, AInt& aLocalIndEDGE,
						AInt& aLocalIndFACE, AInt& aLocalIndVOL,
						MGSelector* pSel = NULL)
{
//	initialize a selector.
	MGSelector tmpSel;
	if(!pSel)
	{
		tmpSel.assign_grid(mg);
		pSel = &tmpSel;
	}
	MGSelector& msel = *pSel;
	
	
//	select all elements in the layouts so that we can serialize
//	that part of the grid.
	SelectNodesInLayout(msel, vrtLayout);
	SelectNodesInLayout(msel, edgeLayout);
	SelectNodesInLayout(msel, faceLayout);
	SelectNodesInLayout(msel, volLayout);
	
//	write the grid.
	SerializeMultiGridElements(mg,
						msel.get_multi_level_geometric_object_collection(),
						aLocalIndVRT, aLocalIndEDGE,
						aLocalIndFACE, aLocalIndVOL, out);

//	write the layouts
}
						
}//	end of namespace
