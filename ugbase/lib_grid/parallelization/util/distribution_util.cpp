//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#include "lib_grid/lib_grid.h"
#include "pcl/pcl.h"
#include "distribution_util.h"
#include "common/util/stream_pack.h"
#include "common/util/binary_stream.h"

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
						int level = 0,
						int interfacesOnLevelOnly = -1)
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
			layout.node_vec().push_back(node);

		//	access the interfaces
			if((interfacesOnLevelOnly == -1) ||
				(interfacesOnLevelOnly == level))
			{
				typename TNodeLayout::Interface& masterInterface = masterLayout.interface(layoutIndex, level);
				typename TNodeLayout::Interface& slaveInterface = layout.interface(masterLayoutIndex, level);
				masterInterface.push_back(typename TNodeLayout::InterfaceEntry(localMasterID, INT_MASTER));
				slaveInterface.push_back(typename TNodeLayout::InterfaceEntry(localID, INT_SLAVE));
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
void CreateDistributionLayouts(
						std::vector<DistributionVertexLayout>& vertexLayoutsOut,
						std::vector<DistributionEdgeLayout>& edgeLayoutsOut,
						std::vector<DistributionFaceLayout>& faceLayoutsOut,
						std::vector<DistributionVolumeLayout>& volumeLayoutsOut,
						MultiGrid& mg, SubsetHandler& sh,
						bool distributeGenealogy,
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

//	resize and clear the layouts
	vertexLayoutsOut = std::vector<DistributionVertexLayout>(sh.num_subsets());
	edgeLayoutsOut = std::vector<DistributionEdgeLayout>(sh.num_subsets());
	faceLayoutsOut = std::vector<DistributionFaceLayout>(sh.num_subsets());
	volumeLayoutsOut = std::vector<DistributionVolumeLayout>(sh.num_subsets());

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
							aaFirstProcVRT, aaFirstProcLocalIndVRT, -1);
		AddNodesToLayout(edgeLayoutsOut, i,
							sh.begin<EdgeBase>(i), sh.end<EdgeBase>(i),
							aaFirstProcEDGE, aaFirstProcLocalIndEDGE, -1);
		AddNodesToLayout(faceLayoutsOut, i,
							sh.begin<Face>(i), sh.end<Face>(i),
							aaFirstProcFACE, aaFirstProcLocalIndFACE, -1);
		AddNodesToLayout(volumeLayoutsOut, i,
							sh.begin<Volume>(i), sh.end<Volume>(i),
							aaFirstProcVOL, aaFirstProcLocalIndVOL, -1);
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
	//	the hierarchy has to be complete. make sure the whole geneology
	//	is selected. By passing true, all associated elements of lower
	//	dimension will be selected, too.

	//	if the whole genealogy shall be distributed, then select it here.
	//	associated elements will automatically be selected.
	//	If howerver vertical interfaces shall be created, the genealogy
	//	shouldn't be distributed. In this case only associated geometric
	//	objects have to be selected.
//TODO: do it as commented above and use the uncommented code below.
		SelectAssociatedGenealogy(msel, true);//remove this
	/*
		if(distributeGenealogy)
			SelectAssociatedGenealogy(msel, true);
		else
			SelectAssociatedGeometricObjects(msel);
	*/

		int interfacesOnLevelOnly = -1;
		if(!distributeGenealogy)
			interfacesOnLevelOnly = mg.num_levels() - 1;

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
								aaFirstProcVRT, aaFirstProcLocalIndVRT, level,
								interfacesOnLevelOnly);
			AddNodesToLayout(edgeLayoutsOut, i,
								msel.begin<EdgeBase>(level), msel.end<EdgeBase>(level),
								aaFirstProcEDGE, aaFirstProcLocalIndEDGE, level,
								interfacesOnLevelOnly);
			AddNodesToLayout(faceLayoutsOut, i,
								msel.begin<Face>(level), msel.end<Face>(level),
								aaFirstProcFACE, aaFirstProcLocalIndFACE, level,
								interfacesOnLevelOnly);
			AddNodesToLayout(volumeLayoutsOut, i,
								msel.begin<Volume>(level), msel.end<Volume>(level),
								aaFirstProcVOL, aaFirstProcLocalIndVOL, level,
								interfacesOnLevelOnly);
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
template <class TDistLayout>
void AddExistingInterfacesForRedistribution(
					DistributedGridManager& distGridMgr,
					std::vector<TDistLayout>& distLayoutVec)
{
//	access the associated multi-grid
	if(!distGridMgr.get_assigned_grid())
		return;

	MultiGrid& mg = *distGridMgr.get_assigned_grid();

//	we'll use this vector to gather existing interface-entries for a node
	vector<pair<int, size_t> > interfaceEntries;

//	iterate through the processes to which interfaces have to be build.
	for(size_t iDistLayout = 0; iDistLayout < distLayoutVec.size();
		++iDistLayout)
	{
		TDistLayout& distLayout = distLayoutVec[iDistLayout];

	//	iterate through the nodes of the dist-layout.
	//	depending on the status, we'll have to add a new interface entry.
		typename TDistLayout::NodeVec& nodes = distLayout.node_vec();

		for(size_t iNode = 0; iNode < nodes.size(); ++iNode)
		{
			typename TDistLayout::NodeType& node = nodes[iNode];
			if(distGridMgr.is_interface_element(node))
			{
			//	the node lies in at least one interface to an old neighbour.
			//	We have to add the entry to corresponding distribution-interfaces
				byte nodeStatus = distGridMgr.get_status(node);
				int iType = INT_UNKNOWN;
				if(nodeStatus & ES_MASTER){
					iType = INT_MASTER;
				}
				else if(nodeStatus & ES_SLAVE){
					iType = INT_SLAVE;
				}
				else{
					UG_LOG("Only ES_MASTER and ES_SLAVE are supported during redistribution in the moment!\n");
					throw(int(0));
				}

				distGridMgr.collect_interface_entries(interfaceEntries, node);

				if(iType != INT_UNKNOWN){
					for(size_t i = 0; i < interfaceEntries.size(); ++i){
						typename TDistLayout::Interface& interface =
											distLayout.interface(interfaceEntries[i].first,
																 mg.get_level(node));
						interface.push_back(DistributionInterfaceEntry(iNode, iType,
																	   interfaceEntries[i].second));
					}
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
void CreateRedistributionLayouts(
						std::vector<DistributionVertexLayout>& vertexLayoutsOut,
						std::vector<DistributionEdgeLayout>& edgeLayoutsOut,
						std::vector<DistributionFaceLayout>& faceLayoutsOut,
						std::vector<DistributionVolumeLayout>& volumeLayoutsOut,
						DistributedGridManager& distGridMgr, SubsetHandler& sh,
						bool distributeGenealogy,
						MGSelector* pSel)
{
//	access the associated multi-grid
	if(!distGridMgr.get_assigned_grid())
		return;

	MultiGrid& mg = *distGridMgr.get_assigned_grid();

//	first we'll create normal distribution layouts
	CreateDistributionLayouts(vertexLayoutsOut, edgeLayoutsOut,
							  faceLayoutsOut, volumeLayoutsOut,
							  mg, sh, distributeGenealogy, pSel);

//	now we can create distribution-interfaces from existing ones
	AddExistingInterfacesForRedistribution(distGridMgr, vertexLayoutsOut);
	AddExistingInterfacesForRedistribution(distGridMgr, edgeLayoutsOut);
	AddExistingInterfacesForRedistribution(distGridMgr, faceLayoutsOut);
	AddExistingInterfacesForRedistribution(distGridMgr, volumeLayoutsOut);

}

////////////////////////////////////////////////////////////////////////
void SerializeGridAndDistributionLayouts(
								std::ostream& out, MultiGrid& mg,
								DistributionVertexLayout& vrtLayout,
								DistributionEdgeLayout& edgeLayout,
								DistributionFaceLayout& faceLayout,
								DistributionVolumeLayout& volLayout,
								AInt& aLocalIndVRT, AInt& aLocalIndEDGE,
								AInt& aLocalIndFACE, AInt& aLocalIndVOL,
								MGSelector* pSel,
								std::vector<int>* pProcessMap)
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
	SerializeDistributionLayoutInterfaces(out, vrtLayout, pProcessMap);
	SerializeDistributionLayoutInterfaces(out, edgeLayout, pProcessMap);
	SerializeDistributionLayoutInterfaces(out, faceLayout, pProcessMap);
	SerializeDistributionLayoutInterfaces(out, volLayout, pProcessMap);

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
	if(gridLayoutOut.has_vertex_layout(INT_MASTER))
	{
		ParallelVertexLayout& pvl = gridLayoutOut.vertex_layout(INT_MASTER);
		PCLLOG("process has vertex-master-layout with " << pvl.num_levels() << " levels\n");
		ParallelVertexLayout::Layout& layout = pvl.layout(0);
		ParallelVertexLayout::Layout::iterator iter;
		for(iter = layout.begin(); iter != layout.end(); ++iter)
		{
			PCLLOG("master-interface to process " << iter->first);
			PCLLOG(" contains " << iter->second.size() << " elements.\n");
		}
	}

	if(gridLayoutOut.has_vertex_layout(INT_SLAVE))
	{
		ParallelVertexLayout& pvl = gridLayoutOut.vertex_layout(INT_SLAVE);
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


size_t NumEntriesOfTypeInDistributionInterface(int type,
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
bool TestDistributionLayouts(std::vector<TDistLayout>& distLayouts)
{
	bool bSuccess = true;

	UG_LOG("Performing DistributionLayout Tests: ...\n")
//	first check whether corresponding interfaces exist
	typedef typename TDistLayout::InterfaceMap 	InterfaceMap;
	typedef typename TDistLayout::Interface		Interface;

	for(int curProcID = 0; curProcID < (int)distLayouts.size(); ++curProcID){
		TDistLayout& curLayout = distLayouts[curProcID];
		for(size_t lvl = 0; lvl < curLayout.num_levels(); ++lvl){
			InterfaceMap& curMap = curLayout.interface_map(lvl);
			for(typename InterfaceMap::iterator mapIter = curMap.begin();
				mapIter != curMap.end(); ++mapIter)
			{
			//	we'll only compare with connected processes with a higher rank.
			//	All others have already been checked.
				int conProcID = mapIter->first;
				if(conProcID < curProcID)
					continue;

				Interface& curIntf = mapIter->second;
				TDistLayout& conLayout = distLayouts[conProcID];
				Interface& conIntf = conLayout.interface(curProcID, lvl);

			//	make sure that both interfaces have the same number of entries.
				if(curIntf.size() != conIntf.size()){
					bSuccess = false;
					UG_LOG("  WARNING: Sizes do not match between interfaces of procs "
							<< curProcID << " and " << conProcID << "on level " << lvl << endl);
				}

			//	make sure that the different interfaces match each other in size
				size_t numCurMasters = NumEntriesOfTypeInDistributionInterface(
															INT_MASTER, curIntf);
				size_t numCurSlaves = NumEntriesOfTypeInDistributionInterface(
															INT_SLAVE, curIntf);
				size_t numConMasters = NumEntriesOfTypeInDistributionInterface(
															INT_MASTER, conIntf);
				size_t numConSlaves = NumEntriesOfTypeInDistributionInterface(
															INT_SLAVE, conIntf);

				size_t numCurVrtMasters = NumEntriesOfTypeInDistributionInterface(
													INT_VERTICAL_MASTER, curIntf);
				size_t numCurVrtSlaves = NumEntriesOfTypeInDistributionInterface(
													INT_VERTICAL_SLAVE, curIntf);
				size_t numConVrtMasters = NumEntriesOfTypeInDistributionInterface(
													INT_VERTICAL_MASTER, conIntf);
				size_t numConVrtSlaves = NumEntriesOfTypeInDistributionInterface(
													INT_VERTICAL_SLAVE, conIntf);

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


template bool TestDistributionLayouts<DistributionVertexLayout>(std::vector<DistributionVertexLayout>&);
template bool TestDistributionLayouts<DistributionEdgeLayout>(std::vector<DistributionEdgeLayout>&);
template bool TestDistributionLayouts<DistributionFaceLayout>(std::vector<DistributionFaceLayout>&);
template bool TestDistributionLayouts<DistributionVolumeLayout>(std::vector<DistributionVolumeLayout>&);

}//	end of namespace
