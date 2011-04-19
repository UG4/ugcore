//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#include "lib_grid/lib_grid.h"
#include "pcl/pcl.h"
#include "distribution_util.h"
#include "common/util/stream_pack.h"
#include "common/util/binary_stream.h"
#include "common/serialization.h"
#include "common/assert.h"

using namespace std;

namespace ug
{
/*
void PrintData(int* data, int size)
{
	for(int i = 0; i < size; ++i)
		cout << data[i] << ", ";

	cout << endl;
}
*/

////////////////////////////////////////////////////////////////////////
//	AddNodesToLayout
///	adds nodes to a layout and to interfaces if required.
/**
 * Please note that this method not only alters the layout referenced
 * by layoutIndex, but all layouts that share a node with this
 * layout. For each node that is referenced by multiple layouts,
 * corresponding interface entries are automatically generated.
 */
template <class TNodeLayout, class TIterator, class TAIntAccessor>
static
void AddNodesToLayout(std::vector<TNodeLayout>& layouts,
						int layoutIndex,
						TIterator nodesBegin, TIterator nodesEnd,
						TAIntAccessor& aaFirstLayout,
						TAIntAccessor& aaFirstProcLocalInd,
						std::vector<int>* processMap = NULL,
						int level = 0,
						int interfacesOnLevelOnly = -1,
						DistributedGridManager* pDistGridMgr = NULL)
{
	typedef typename TNodeLayout::Interface			Interface;
	typedef typename TNodeLayout::InterfaceEntry	InterfaceEntry;

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
			layout.node_vec().push_back(node);

		//	access the interfaces
		//	if the node already is in a 'real' horizontal interface, we'll ignore it
		//todo: check the type of interface
			if(pDistGridMgr){
				if(pDistGridMgr->contains_status(*iter, INT_H_MASTER)
				  || pDistGridMgr->contains_status(*iter, INT_H_SLAVE))
					continue;
			}

			if((interfacesOnLevelOnly == -1) ||
				(interfacesOnLevelOnly == level))
			{
			//	get the master and slave proc-id
				int masterProc = masterLayoutIndex;
				int slaveProc = layoutIndex;
				if(processMap){
					masterProc = (*processMap)[masterProc];
					slaveProc = (*processMap)[slaveProc];
				}

				Interface& masterInterface = masterLayout.interface(slaveProc, level);
				Interface& slaveInterface = layout.interface(masterProc, level);
				masterInterface.push_back(InterfaceEntry(localMasterID, INT_H_MASTER));
				slaveInterface.push_back(InterfaceEntry(localID, INT_H_SLAVE));
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
/**	Creates vertical interfaces for all elements which are moved away
 * from this process and which do not have a parent or whose parent
 * is not moved along with them.
 *
 * Note that this method may add a new distLayout, a new process-map entry
 * and a new empty subset to shPart.
 */
template <class TDistLayout, class TAAInt>
void CreateVerticalInterfaces(std::vector<TDistLayout>& distLayouts,
								MultiGrid& mg, MGSelector& msel,
								SubsetHandler& shPart, TAAInt& aaInt,
								DistributedGridManager* pDistGridMgr = NULL,
								std::vector<int>* processMap = NULL)
{
//	all elements which do not have a parent on their associated proc
//	have to be added to a vertical interface to localProc. An associated
//	vertical master has to be added in the associated interface on
//	localProc.

	typedef typename TDistLayout::NodeVec			NodeVec;
	typedef typename TDistLayout::NodeType			Node;
	typedef typename TDistLayout::Interface			Interface;
	typedef typename TDistLayout::InterfaceEntry	InterfaceEntry;

	int locProcRank = pcl::GetProcRank();

//	we have to find the local LayoutIndex
	int locLayoutInd = -1;
	if(processMap){
		vector<int>& procMap = *processMap;
		for(size_t i = 0; i < distLayouts.size(); ++i){
			if(procMap[i] == locProcRank){
				locLayoutInd = i;
				break;
			}
		}
	//	if no layout was found for the local proc, we have to add a new one
		if(locLayoutInd == -1){
			locLayoutInd = (int)distLayouts.size();
			shPart.subset_required(shPart.num_subsets());
			distLayouts.push_back(TDistLayout());
			procMap.push_back(locProcRank);
		}
	}
	else{
		locLayoutInd = locProcRank;
	//	if no layout exists for the local proc, we have to add new ones
		if(locLayoutInd >= (int)distLayouts.size()){
			distLayouts.resize(locLayoutInd + 1);
			shPart.subset_required(locLayoutInd);
		}
	}

	TDistLayout& locLayout = distLayouts[locLayoutInd];
	NodeVec& locNodes = locLayout.node_vec();

//	msel will always hold all nodes in localLayout
//	aaInt will be used to store the local index of each selected node
	msel.clear();
	for(size_t i = 0; i < locNodes.size(); ++i){
		msel.select(locNodes[i]);
		aaInt[locNodes[i]] = i;
	}

//	iterate over all distLayouts but the local one
	for(size_t i_layout = 0; i_layout < distLayouts.size(); ++i_layout)
	{
		if((int)i_layout == locLayoutInd)
			continue;

		int curProcRank = i_layout;
		if(processMap)
			curProcRank = processMap->at(i_layout);

		TDistLayout& curLayout = distLayouts[i_layout];

	//	mark all nodes in curLayout
		mg.begin_marking();
		MarkNodesInLayout(mg, curLayout);

	//	we wont access the interfaces directly. Instead we'll cache them
	//	for efficient reuse
		Interface* locInterface = NULL;
		Interface* curInterface = NULL;
		int interfaceLevel = -1;

	//	check for each node whether its parent is not marked or if it has
	//	no parent at all
		NodeVec& curNodes = curLayout.node_vec();
		for(size_t i_node = 0; i_node < curNodes.size(); ++i_node)
		{
			Node node = curNodes[i_node];

		//	if pDistGridMgr is supplied, we first check however whether
		//	the element is already contained in a vertical interface.
		//	If so, we wont add it to another one.
			if(!pDistGridMgr || pDistGridMgr->contains_status(node, INT_V_MASTER)
			  || pDistGridMgr->contains_status(node, INT_V_SLAVE))
				continue;

			GeometricObject* parent = mg.get_parent(node);
		//	!mg.is_marked(parent) is only executed if a parent exists.
			if((!parent) || !mg.is_marked(parent)){
			//	the element has to be put into a vertical interface
			//	get the level and check whether we have to access the interfaces again
				int lvl = mg.get_level(node);
				if(lvl != interfaceLevel){
				//	we have to access the interfaces
					interfaceLevel = lvl;
					locInterface = &locLayout.interface(curProcRank, lvl);
					curInterface = &curLayout.interface(locProcRank, lvl);
				}

			//	before we can insert node into locInterface, we first have to
			//	make sure, that it is contained in locLayout.
				if(!msel.is_selected(node)){
					msel.select(node);
					aaInt[node] = (int)locNodes.size();
					locNodes.push_back(node);
				}

			//	finally insert the node into the interfaces
				locInterface->push_back(InterfaceEntry(aaInt[node], INT_V_MASTER));
				curInterface->push_back(InterfaceEntry(i_node, INT_V_SLAVE));
			}
		}

		mg.end_marking();
	}
}

////////////////////////////////////////////////////////////////////////
template <class TVertexDistributionLayout, class TEdgeDistributionLayout,
		  class TFaceDistributionLayout, class TVolumeDistributionLayout>
void CreateDistributionLayouts(
						std::vector<TVertexDistributionLayout>& vertexLayoutsOut,
						std::vector<TEdgeDistributionLayout>& edgeLayoutsOut,
						std::vector<TFaceDistributionLayout>& faceLayoutsOut,
						std::vector<TVolumeDistributionLayout>& volumeLayoutsOut,
						MultiGrid& mg, SubsetHandler& sh,
						bool distributeGenealogy,
						bool createVerticalInterfaces,
						MGSelector* pSel,
						DistributedGridManager* pDistGridMgr,
						std::vector<int>* processMap)
{
//	initialize a selector.
	MGSelector tmpSel;
	if(!pSel)
	{
		tmpSel.assign_grid(mg);
		pSel = &tmpSel;
	}
	MGSelector& msel = *pSel;

//	make sure that the processMap has the right size
	if(processMap){
		UG_ASSERT((int)processMap->size() == sh.num_subsets(),
				  "ProcessMap has to have as many entries as there are partitions");
	}
//	resize and clear the layouts
	vertexLayoutsOut = std::vector<TVertexDistributionLayout>(sh.num_subsets());
	edgeLayoutsOut = std::vector<TEdgeDistributionLayout>(sh.num_subsets());
	faceLayoutsOut = std::vector<TFaceDistributionLayout>(sh.num_subsets());
	volumeLayoutsOut = std::vector<TVolumeDistributionLayout>(sh.num_subsets());

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
	for(int i = 0; i < sh.num_subsets(); ++i)
	{
	//	the level is ignored since it won't be used in this phase.
	//	by passing -1 we can assert that no interface is accessed.
		AddNodesToLayout(vertexLayoutsOut, i,
							sh.begin<VertexBase>(i), sh.end<VertexBase>(i),
							aaFirstProcVRT, aaFirstProcLocalIndVRT,
							processMap, -1);
		AddNodesToLayout(edgeLayoutsOut, i,
							sh.begin<EdgeBase>(i), sh.end<EdgeBase>(i),
							aaFirstProcEDGE, aaFirstProcLocalIndEDGE,
							processMap, -1);
		AddNodesToLayout(faceLayoutsOut, i,
							sh.begin<Face>(i), sh.end<Face>(i),
							aaFirstProcFACE, aaFirstProcLocalIndFACE,
							processMap, -1);
		AddNodesToLayout(volumeLayoutsOut, i,
							sh.begin<Volume>(i), sh.end<Volume>(i),
							aaFirstProcVOL, aaFirstProcLocalIndVOL,
							processMap, -1);
	}

//	step 2: add all the associated elements to the distribution groups, which
//			have not already been assigned.
	for(int i = 0; i < sh.num_subsets(); ++i)
	{
		msel.clear();
		msel.select(sh.begin<VertexBase>(i), sh.end<VertexBase>(i));
		msel.select(sh.begin<EdgeBase>(i), sh.end<EdgeBase>(i));
		msel.select(sh.begin<Face>(i), sh.end<Face>(i));
		msel.select(sh.begin<Volume>(i), sh.end<Volume>(i));
//TODO: overlap can be easily handled here! simply increase the selection.
//		eventually we first would have to select all associated elements.
	//	the hierarchy has to be complete. make sure the whole genealogy
	//	is selected. By passing true, all associated elements of lower
	//	dimension will be selected, too.

	//	if the whole genealogy shall be distributed, then select it here.
	//	associated elements will automatically be selected.
	//	If however vertical interfaces shall be created, the genealogy
	//	shouldn't be distributed. In this case only associated geometric
	//	objects have to be selected.
		if(distributeGenealogy)
			SelectAssociatedGenealogy(msel, true);
		else
			SelectAssociatedGeometricObjects(msel);

		int interfacesOnLevelOnly = -1;
		if(distributeGenealogy)
			interfacesOnLevelOnly = mg.num_levels() - 1;

	//	now add the missing horizontal interfaces
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
			AddNodesToLayout(vertexLayoutsOut, i,
								msel.begin<VertexBase>(level), msel.end<VertexBase>(level),
								aaFirstProcVRT, aaFirstProcLocalIndVRT,
								processMap, level,
								interfacesOnLevelOnly, pDistGridMgr);
			AddNodesToLayout(edgeLayoutsOut, i,
								msel.begin<EdgeBase>(level), msel.end<EdgeBase>(level),
								aaFirstProcEDGE, aaFirstProcLocalIndEDGE,
								processMap, level,
								interfacesOnLevelOnly, pDistGridMgr);
			AddNodesToLayout(faceLayoutsOut, i,
								msel.begin<Face>(level), msel.end<Face>(level),
								aaFirstProcFACE, aaFirstProcLocalIndFACE,
								processMap, level,
								interfacesOnLevelOnly, pDistGridMgr);
			AddNodesToLayout(volumeLayoutsOut, i,
								msel.begin<Volume>(level), msel.end<Volume>(level),
								aaFirstProcVOL, aaFirstProcLocalIndVOL,
								processMap, level,
								interfacesOnLevelOnly, pDistGridMgr);
		}
	}

//	horizontal layouts are complete by now. All nodes that go to process i
//	are contained in layout i at this point.
//	we can now use this information to add vertical interfaces
	if(createVerticalInterfaces){
	//	we'll reuse aaFirstProc... for a different purpose here.
		CreateVerticalInterfaces(vertexLayoutsOut, mg, msel, sh, aaFirstProcVRT,
								pDistGridMgr, processMap);
		CreateVerticalInterfaces(edgeLayoutsOut, mg, msel, sh, aaFirstProcEDGE,
								pDistGridMgr, processMap);
		CreateVerticalInterfaces(faceLayoutsOut, mg, msel, sh, aaFirstProcFACE,
								pDistGridMgr, processMap);
		CreateVerticalInterfaces(volumeLayoutsOut, mg, msel, sh, aaFirstProcVOL,
								pDistGridMgr, processMap);
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

//	explicit template instantiation
template void CreateDistributionLayouts<DistributionVertexLayout, DistributionEdgeLayout,
										DistributionFaceLayout, DistributionVolumeLayout>(
										std::vector<DistributionVertexLayout>&,
										std::vector<DistributionEdgeLayout>&,
										std::vector<DistributionFaceLayout>&,
										std::vector<DistributionVolumeLayout>&,
										MultiGrid&, SubsetHandler&, bool, bool, MGSelector*,
										DistributedGridManager*, std::vector<int>*);

template void CreateDistributionLayouts<RedistributionVertexLayout, RedistributionEdgeLayout,
										RedistributionFaceLayout, RedistributionVolumeLayout>(
										std::vector<RedistributionVertexLayout>&,
										std::vector<RedistributionEdgeLayout>&,
										std::vector<RedistributionFaceLayout>&,
										std::vector<RedistributionVolumeLayout>&,
										MultiGrid&, SubsetHandler&, bool, bool, MGSelector*,
										DistributedGridManager*, std::vector<int>*);

/*
////////////////////////////////////////////////////////////////////////
void CreateDistributionLayouts_SplitBaseGrid(
						std::vector<DistributionVertexLayout>& vertexLayoutsOut,
						std::vector<DistributionEdgeLayout>& edgeLayoutsOut,
						std::vector<DistributionFaceLayout>& faceLayoutsOut,
						std::vector<DistributionVolumeLayout>& volumeLayoutsOut,
						MultiGrid& mg, SubsetHandler& sh,
						IDomainDecompositionInfo& ddinfo,
						MGSelector* pSel)
{
//	initialize a selector.
	MGSelector tmpSel;
	if(!pSel){
		tmpSel.assign_grid(mg);
		pSel = &tmpSel;
	}
	MGSelector& msel = *pSel;

//	call normal CreateDistributionLayouts first.
	CreateDistributionLayouts(vertexLayoutsOut, edgeLayoutsOut,
							  faceLayoutsOut, volumeLayoutsOut,
							  mg, sh, false, &msel);

//	now we have to create a base grid for each domain partition
//	to do so we'll iterate over all subdomains in ddinfo and
//	collect the base grid of each. On the fly we'll create the
//	horizontal interfaces.

	std::vector<int> subdomProcs;
	for(int i_subdom = 0; i_subdom < ddinfo.num_subdomains(); ++i_subdom)
	{
		ddinfo.get_subdomain_procs(subdomProcs, i_subdom);
	//	the first proc in each subdomain will hold the base grid.

	}

}
*/

////////////////////////////////////////////////////////////////////////
void SerializeGridAndDistributionLayouts(
								std::ostream& out, MultiGrid& mg,
								DistributionVertexLayout& vrtLayout,
								DistributionEdgeLayout& edgeLayout,
								DistributionFaceLayout& faceLayout,
								DistributionVolumeLayout& volLayout,
								AInt& aLocalIndVRT, AInt& aLocalIndEDGE,
								AInt& aLocalIndFACE, AInt& aLocalIndVOL,
								MGSelector* pSel)
{
//	initialize a selector.
	MGSelector tmpSel;
	if(!pSel)
	{
		tmpSel.assign_grid(mg);
		pSel = &tmpSel;
	}
	MGSelector& msel = *pSel;

	msel.clear();

//	select all elements in the layouts so that we can serialize
//	that part of the grid.
	SelectNodesInLayout(msel, vrtLayout);
	SelectNodesInLayout(msel, edgeLayout);
	SelectNodesInLayout(msel, faceLayout);
	SelectNodesInLayout(msel, volLayout);

//	write the grid.
//	during serialization the local indices are automatically generated
//	and written to the aLocalInd... attachments.
	SerializeMultiGridElements(mg,
						msel.get_geometric_object_collection(),
						aLocalIndVRT, aLocalIndEDGE,
						aLocalIndFACE, aLocalIndVOL, out);

//	write the layouts
	SerializeDistributionLayoutInterfaces(out, vrtLayout);
	SerializeDistributionLayoutInterfaces(out, edgeLayout);
	SerializeDistributionLayoutInterfaces(out, faceLayout);
	SerializeDistributionLayoutInterfaces(out, volLayout);

//	done. Please note that no attachments have been serialized in this method.
}
/*
template <class TGeomObj, class TLayout>
static
void
FillLayoutWithNodes(TLayout& layout, Grid& grid)
{
	typedef typename TLayout::NodeType TNode;
	typedef typename geometry_traits<TGeomObj>::iterator iterator;
	typename TLayout::NodeVec& nodes = layout.node_vec();

	for(iterator iter = grid.begin<TGeomObj>();
		iter != grid.end<TGeomObj>(); ++iter)
		nodes.push_back(*iter);
}
*/

////////////////////////////////////////////////////////////////////////
//	DeserializeGridAndLayouts
void DeserializeGridAndDistributionLayouts(MultiGrid& mgOut,
											GridLayoutMap& gridLayoutOut,
											std::istream& in)
{
//	read the grid.
//	we'll need vectors which contain the elements of the grid later on.
//	This is handled by the deserialization routine automatically, if
//	we pass pointers to those vectors to the method.
	vector<VertexBase*>	vVrts;
	vector<EdgeBase*>	vEdges;
	vector<Face*>		vFaces;
	vector<Volume*>		vVols;

	DeserializeMultiGridElements(mgOut, in, &vVrts, &vEdges, &vFaces, &vVols);

//	read the layouts
/*
	DeserializeLayoutInterfaces<VertexBase>(
					gridLayoutOut.vertex_layout_hierarchy_map(), vVrts, in);
	DeserializeLayoutInterfaces<EdgeBase>(
					gridLayoutOut.edge_layout_hierarchy_map(), vEdges, in);
	DeserializeLayoutInterfaces<Face>(
					gridLayoutOut.face_layout_hierarchy_map(), vFaces, in);
	DeserializeLayoutInterfaces<Volume>(
					gridLayoutOut.volume_layout_hierarchy_map(), vVols, in);
*/

	DeserializeDistributionLayoutInterfaces<VertexBase>(gridLayoutOut,
														vVrts, in);
	DeserializeDistributionLayoutInterfaces<EdgeBase>(gridLayoutOut,
														vEdges, in);
	DeserializeDistributionLayoutInterfaces<Face>(gridLayoutOut,
													vFaces, in);
	DeserializeDistributionLayoutInterfaces<Volume>(gridLayoutOut,
													vVols, in);

//DEBUG
/*
	PCLLOG("deserialization done.\n");
	if(gridLayoutOut.has_vertex_layout(INT_H_MASTER))
	{
		ParallelVertexLayout& pvl = gridLayoutOut.vertex_layout(INT_H_MASTER);
		PCLLOG("process has vertex-master-layout with " << pvl.num_levels() << " levels\n");
		ParallelVertexLayout::Layout& layout = pvl.layout(0);
		ParallelVertexLayout::Layout::iterator iter;
		for(iter = layout.begin(); iter != layout.end(); ++iter)
		{
			PCLLOG("master-interface to process " << iter->first);
			PCLLOG(" contains " << iter->second.size() << " elements.\n");
		}
	}

	if(gridLayoutOut.has_vertex_layout(INT_H_SLAVE))
	{
		ParallelVertexLayout& pvl = gridLayoutOut.vertex_layout(INT_H_SLAVE);
		PCLLOG("process has vertex-slave-layout with " << pvl.num_levels() << " levels\n");
		ParallelVertexLayout::Layout& layout = pvl.layout(0);
		ParallelVertexLayout::Layout::iterator iter;
		for(iter = layout.begin(); iter != layout.end(); ++iter)
		{
			PCLLOG("slave-interface to process " << iter->first);
			PCLLOG(" contains " << iter->second.size() << " elements.\n");
		}
	}
*/
//	done. Please note that no attachments have been serialized in this method.
}


size_t NumEntriesOfTypeInDistributionInterface(unsigned int type,
			std::vector<DistributionInterfaceEntry>& interface)
{
	size_t counter = 0;
	for(size_t i = 0; i < interface.size(); ++i){
		if(interface[i].type == type)
			++counter;
	}
	return counter;
}

//todo: copy implementation to ..._impl.hpp
template <class TDistLayout>
bool TestDistributionLayouts(std::vector<TDistLayout>& distLayouts,
							int* procMap)
{
	bool bSuccess = true;

	UG_LOG("Performing DistributionLayout Tests: ...\n")
//	first check whether corresponding interfaces exist
	typedef typename TDistLayout::InterfaceMap 	InterfaceMap;
	typedef typename TDistLayout::Interface		Interface;

	for(int i_curLayout = 0; i_curLayout < (int)distLayouts.size(); ++i_curLayout)
	{
		TDistLayout& curLayout = distLayouts[i_curLayout];

		int curProcID = i_curLayout;
		if(procMap)
			curProcID = procMap[i_curLayout];

		for(size_t lvl = 0; lvl < curLayout.num_levels(); ++lvl)
		{
			InterfaceMap& curMap = curLayout.interface_map(lvl);
			for(typename InterfaceMap::iterator mapIter = curMap.begin();
				mapIter != curMap.end(); ++mapIter)
			{
			//	we'll only compare with connected processes with a higher rank.
			//	All others have already been checked.
				int conProcID = mapIter->first.first;
				if(conProcID <= curProcID)
					continue;

				Interface& curIntf = mapIter->second;
				TDistLayout& conLayout = distLayouts[conProcID];
				Interface& conIntf = conLayout.interface(curProcID, lvl);

			//	make sure that both interfaces have the same number of entries.
				if(curIntf.size() != conIntf.size()){
					bSuccess = false;
					UG_LOG("  WARNING: Sizes do not match between interfaces of procs "
							<< curProcID << " and " << conProcID << " on level " << lvl << endl);
				}

			//	make sure that the different interfaces match each other in size
				size_t numCurMasters = NumEntriesOfTypeInDistributionInterface(
															INT_H_MASTER, curIntf);
				size_t numCurSlaves = NumEntriesOfTypeInDistributionInterface(
															INT_H_SLAVE, curIntf);
				size_t numConMasters = NumEntriesOfTypeInDistributionInterface(
															INT_H_MASTER, conIntf);
				size_t numConSlaves = NumEntriesOfTypeInDistributionInterface(
															INT_H_SLAVE, conIntf);

				size_t numCurVrtMasters = NumEntriesOfTypeInDistributionInterface(
													INT_V_MASTER, curIntf);
				size_t numCurVrtSlaves = NumEntriesOfTypeInDistributionInterface(
													INT_V_SLAVE, curIntf);
				size_t numConVrtMasters = NumEntriesOfTypeInDistributionInterface(
													INT_V_MASTER, conIntf);
				size_t numConVrtSlaves = NumEntriesOfTypeInDistributionInterface(
													INT_V_SLAVE, conIntf);

				if(numCurMasters != numConSlaves){
					UG_LOG("  Master -> Slave Interface mismatch on level " << lvl << ":\n");
					UG_LOG("\t" << numCurMasters << " masters on process " << curProcID << endl);
					UG_LOG("\t" << numConSlaves << " slaves on process " << conProcID << endl);
				}

				if(numCurSlaves != numConMasters){
					UG_LOG("  Slave -> Master Interface mismatch on level " << lvl << ":\n");
					UG_LOG("\t" << numCurSlaves << " slaves on process " << curProcID << endl);
					UG_LOG("\t" << numConMasters << " masters on process " << conProcID << endl);
				}

				if(numCurVrtMasters != numConVrtSlaves){
					UG_LOG("  VerticalMaster -> VerticalSlave Interface mismatch on level " << lvl << ":\n");
					UG_LOG("\t" << numCurVrtMasters << " vertical masters on process " << curProcID << endl);
					UG_LOG("\t" << numConVrtSlaves << " vertical slaves on process " << conProcID << endl);
				}

				if(numCurVrtSlaves != numConVrtMasters){
					UG_LOG("  VerticalSlave -> VerticalMaster Interface mismatch on level " << lvl << ":\n");
					UG_LOG("\t" << numCurVrtSlaves << " vertical slaves on process " << curProcID << endl);
					UG_LOG("\t" << numConVrtMasters << " vertical masters on process " << conProcID << endl);
				}
			}
		}
	}
	UG_LOG("  ... done\n");
	return bSuccess;
}


template bool TestDistributionLayouts<DistributionVertexLayout>(std::vector<DistributionVertexLayout>&, int*);
template bool TestDistributionLayouts<DistributionEdgeLayout>(std::vector<DistributionEdgeLayout>&, int*);
template bool TestDistributionLayouts<DistributionFaceLayout>(std::vector<DistributionFaceLayout>&, int*);
template bool TestDistributionLayouts<DistributionVolumeLayout>(std::vector<DistributionVolumeLayout>&, int*);




template <class TDistLayout>
bool TestRedistributionLayouts(std::vector<TDistLayout>& distLayouts,
								int* procMap)
{
	bool bSuccess = true;

	UG_LOG("Performing RedistributionLayout Tests: ...\n")
	UG_LOG("Layouts: " << distLayouts.size() << endl);

//	first check whether corresponding interfaces exist
	typedef typename TDistLayout::InterfaceMap 	InterfaceMap;
	typedef typename TDistLayout::Interface		Interface;

	for(int i_curLayout = 0; i_curLayout < (int)distLayouts.size(); ++i_curLayout)
	{
		TDistLayout& curLayout = distLayouts[i_curLayout];

		int curProcID = i_curLayout;
		if(procMap)
			curProcID = procMap[i_curLayout];

		for(size_t lvl = 0; lvl < curLayout.num_levels(); ++lvl)
		{
			InterfaceMap& curMap = curLayout.interface_map(lvl);
			for(typename InterfaceMap::iterator mapIter = curMap.begin();
				mapIter != curMap.end(); ++mapIter)
			{
				int conProcID = mapIter->first.first;
				if(conProcID == curProcID)
					continue;

				UG_LOG("  connections " << curProcID << " - " << conProcID << ":");

				Interface& curIntf = mapIter->second;

			//	make sure that the different interfaces match each other in size
				size_t numCurMasters = NumEntriesOfTypeInDistributionInterface(
															INT_H_MASTER, curIntf);
				size_t numCurSlaves = NumEntriesOfTypeInDistributionInterface(
															INT_H_SLAVE, curIntf);

				size_t numCurVrtMasters = NumEntriesOfTypeInDistributionInterface(
													INT_V_MASTER, curIntf);
				size_t numCurVrtSlaves = NumEntriesOfTypeInDistributionInterface(
													INT_V_SLAVE, curIntf);

				if(numCurMasters){
					UG_LOG("    h-masters: " << numCurMasters);
				}

				if(numCurSlaves){
					UG_LOG("    h-slaves: " << numCurSlaves);
				}

				if(numCurVrtMasters){
					UG_LOG("    v-masters: " << numCurVrtMasters);
				}

				if(numCurVrtSlaves){
					UG_LOG("    v-slaves: " << numCurVrtSlaves);
				}

				UG_LOG(endl);
			}
		}
	}
	UG_LOG("  ... done\n");
	return bSuccess;
}


template bool TestRedistributionLayouts<RedistributionVertexLayout>(std::vector<RedistributionVertexLayout>&, int*);
template bool TestRedistributionLayouts<RedistributionEdgeLayout>(std::vector<RedistributionEdgeLayout>&, int*);
template bool TestRedistributionLayouts<RedistributionFaceLayout>(std::vector<RedistributionFaceLayout>&, int*);
template bool TestRedistributionLayouts<RedistributionVolumeLayout>(std::vector<RedistributionVolumeLayout>&, int*);

template <class TDistLayout>
bool PrintRedistributionLayouts(std::vector<TDistLayout>& distLayouts)
{
	bool bSuccess = true;

	UG_LOG("Printing RedistributionLayouts: ...\n")
	UG_LOG("Layouts: " << distLayouts.size() << endl);

//	first check whether corresponding interfaces exist
	typedef typename TDistLayout::InterfaceMap 	InterfaceMap;
	typedef typename TDistLayout::Interface		Interface;

	for(int i_curLayout = 0; i_curLayout < (int)distLayouts.size(); ++i_curLayout)
	{
		TDistLayout& curLayout = distLayouts[i_curLayout];
		UG_LOG("layout with source proc: " << curLayout.get_source_proc() << endl);

		for(size_t lvl = 0; lvl < curLayout.num_levels(); ++lvl)
		{
			InterfaceMap& curMap = curLayout.interface_map(lvl);
			for(typename InterfaceMap::iterator mapIter = curMap.begin();
				mapIter != curMap.end(); ++mapIter)
			{
				int conProcID = mapIter->first.first;
				int oldConProcID = mapIter->first.second;

				UG_LOG("  interface to " << conProcID << ":\n");
				UG_LOG("  old connected proc: " << oldConProcID << endl);

				Interface& curIntf = mapIter->second;

			//	make sure that the different interfaces match each other in size
				size_t numCurMasters = NumEntriesOfTypeInDistributionInterface(
															INT_H_MASTER, curIntf);
				size_t numCurSlaves = NumEntriesOfTypeInDistributionInterface(
															INT_H_SLAVE, curIntf);

				size_t numCurVrtMasters = NumEntriesOfTypeInDistributionInterface(
													INT_V_MASTER, curIntf);
				size_t numCurVrtSlaves = NumEntriesOfTypeInDistributionInterface(
													INT_V_SLAVE, curIntf);

				if(numCurMasters){
					UG_LOG("    h-masters: " << numCurMasters);
				}

				if(numCurSlaves){
					UG_LOG("    h-slaves: " << numCurSlaves);
				}

				if(numCurVrtMasters){
					UG_LOG("    v-masters: " << numCurVrtMasters);
				}

				if(numCurVrtSlaves){
					UG_LOG("    v-slaves: " << numCurVrtSlaves);
				}

				UG_LOG(endl);
			}
		}
	}
	UG_LOG("  ... done\n");
	return bSuccess;
}


template bool PrintRedistributionLayouts<RedistributionVertexLayout>(std::vector<RedistributionVertexLayout>&);
template bool PrintRedistributionLayouts<RedistributionEdgeLayout>(std::vector<RedistributionEdgeLayout>&);
template bool PrintRedistributionLayouts<RedistributionFaceLayout>(std::vector<RedistributionFaceLayout>&);
template bool PrintRedistributionLayouts<RedistributionVolumeLayout>(std::vector<RedistributionVolumeLayout>&);

}//	end of namespace
