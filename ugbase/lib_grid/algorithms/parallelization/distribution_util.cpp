//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#include "lib_grid/lib_grid.h"
#include "pcl/parallel_node_layout.h"
#include "distribution_util.h"
#include "common/util/stream_pack.h"


namespace ug
{

////////////////////////////////////////////////////////////////////////
void DistributeGrid(MultiGrid& grid, SubsetHandler& sh,
					int localProcID, MultiGrid* pLocalGridOut,
					ParallelGridLayout* pLocalGridLayoutOut,
					std::vector<int>* pProcessMap)
{
/*
	vector<GridDistributionPack> distPacks;
	vector<int> fallbackProcessMap;
	vector<int>& processMap = fallbackProcessMap;
	int localProcIndex = -1;//the index at which localProcID is in the process-map.
	
//	if a process-map was provided, we'll use it.
//	if not, we'll create our own.
	if(pProcessMap)
		processMap = *pProcessMap;
	else
		for(int i = 0; i < sh.num_subsets(); ++i)
			processMap.push_back(i);
					
	CreateGridDistributionPacks(distPacks, grid, sh, processMap);

//	set up the receiver-process-map that is used to distribute the binaryStream to
//	the processes. This one differs from the original process-map in that it does
//	not contain the localProcID. If the original process-map did not contain the
//	localProcID, then both maps are equal.
	vector<int> receiverProcMap;
	for(int i = 0; i < processMap.size(); ++i)
	{
		if(processMap[i] != localProcID)
			receiverProcMap.push_back(processMap[i]);
		else
			localProcIndex = i;
	}
	
	int numRecProcs = (int)receiverProcMap.size();
	
//	send to each receiver-process the size of the grid-distribution-pack
//	it will receive.
//	The binary-stream to which we will pack all the data.
	BinaryStream binaryStream;
	vector<int>	blockSize(numRecProcs);
	int blockInd = 0;
	for(int i = 0; i < (int)distPacks.size(); ++i)
	{
		if(i != localProcIndex)
		{
			int oldSize = binaryStream.size();
			WriteGridDistributionPackToBinaryStream(distPacks[i], binaryStream);
			blockSize[blockInd++] = (int)(binaryStream.size() - oldSize);
		}
	}
	
	int tag = 0;
	
	vector<int> bufferSizes(numRecProcs, sizeof(int));
	if(numRecProcs > 0)
	{		
	//	distribute the block-sizes to the different processes
		pcl::DistributeData(localProcID, &receiverProcMap.front(), numRecProcs,
							&blockSize.front(), &bufferSizes.front(), 38);

	//	distribute the grids-distribution-packs
		pcl::DistributeData(localProcID, &receiverProcMap.front(), numRecProcs,
							binaryStream.buffer(), &blockSize.front(), 39);
	}
	
//	create the grid for process localProcID
	if(pLocalGridOut && pLocalGridCommSetOut && (localProcIndex != -1))
		UnpackGridDistributionPack(*pLocalGridOut, *pLocalGridCommSetOut,
									distPacks[localProcIndex]);
*/
}

////////////////////////////////////////////////////////////////////////
void ReceiveGrid(MultiGrid& gridOut, ParallelGridLayout& gridLayoutOut,
					int srcProcID)
{
/*
//	receive the stream-size
	int streamSize;
	pcl::ReceiveData(&streamSize, srcProcID, sizeof(int), 38);

//	receive the buffer
	BinaryStream binaryStream(streamSize);
	pcl::ReceiveData(binaryStream.buffer(), srcProcID, streamSize, 39);

//	create the grid and the communication set.
	GridDistributionPack distPack;
	ReadGridDistributionPackFromBinaryStream(distPack, binaryStream);
	UnpackGridDistributionPack(gridOut, gridCommSetOut, distPack);
*/
}


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

////////////////////////////////////////////////////////////////////////
void CreateGridLayouts(	std::vector<ParallelVertexLayout>& vertexLayoutsOut,
						std::vector<ParallelEdgeLayout>& edgeLayoutsOut,
						std::vector<ParallelFaceLayout>& faceLayoutsOut,
						std::vector<ParallelVolumeLayout>& volumeLayoutsOut,
						MultiGrid& mg, SubsetHandler& sh,
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
			AddNodesToLayout(vertexLayoutsOut, i,
								msel.begin<VertexBase>(level), msel.end<VertexBase>(level),
								aaFirstProcVRT, aaFirstProcLocalIndVRT, level);
			AddNodesToLayout(edgeLayoutsOut, i,
								msel.begin<EdgeBase>(level), msel.end<EdgeBase>(level),
								aaFirstProcEDGE, aaFirstProcLocalIndEDGE, level);
			AddNodesToLayout(faceLayoutsOut, i,
								msel.begin<Face>(level), msel.end<Face>(level),
								aaFirstProcFACE, aaFirstProcLocalIndFACE, level);
			AddNodesToLayout(volumeLayoutsOut, i,
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
void SerializeGridAndLayouts(std::ostream& out, MultiGrid& mg,
						ParallelVertexLayout& vrtLayout,
						ParallelEdgeLayout& edgeLayout,
						ParallelFaceLayout& faceLayout,
						ParallelVolumeLayout& volLayout,
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
	
//	init the attachment accessors
	Grid::VertexAttachmentAccessor<AInt> aaLocalIndVRT(mg, aLocalIndVRT);
	Grid::EdgeAttachmentAccessor<AInt> aaLocalIndEDGE(mg, aLocalIndEDGE);
	Grid::FaceAttachmentAccessor<AInt> aaLocalIndFACE(mg, aLocalIndFACE);
	Grid::VolumeAttachmentAccessor<AInt> aaLocalIndVOL(mg, aLocalIndVOL);
	
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
						msel.get_multi_level_geometric_object_collection(),
						aLocalIndVRT, aLocalIndEDGE,
						aLocalIndFACE, aLocalIndVOL, out);

//	write the layouts
	SerializeLayoutInterfaces(out, vrtLayout, aaLocalIndVRT, pProcessMap);
	SerializeLayoutInterfaces(out, edgeLayout, aaLocalIndEDGE, pProcessMap);
	SerializeLayoutInterfaces(out, faceLayout, aaLocalIndFACE, pProcessMap);
	SerializeLayoutInterfaces(out, volLayout, aaLocalIndVOL, pProcessMap);
	
//	done. Please note that no attachments have been serialized in this method.
}
						
}//	end of namespace
