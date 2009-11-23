//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#include "lib_grid/lib_grid.h"
#include "pcl/parallel_node_layout.h"
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
void DistributeGrid(MultiGrid& mg, SubsetHandler& sh,
					int localProcID, MultiGrid* pLocalGridOut,
					ParallelGridLayout* pLocalGridLayoutOut,
					std::vector<int>* pProcessMap)
{

//	we have to store the layouts for all the processes.
	vector<ParallelVertexLayout> vVertexLayouts;
	vector<ParallelEdgeLayout> vEdgeLayouts;
	vector<ParallelFaceLayout> vFaceLayouts;
	vector<ParallelVolumeLayout> vVolumeLayouts;
	
//	we need some attachments that will speed up the called processes.
	AInt aInt;
	mg.attach_to_vertices(aInt);
	mg.attach_to_edges(aInt);
	mg.attach_to_faces(aInt);
	mg.attach_to_volumes(aInt);
	
//	the selector will help to speed things up a little.
	MGSelector msel(mg);
	
	CreateGridLayouts(vVertexLayouts, vEdgeLayouts, vFaceLayouts,
						vVolumeLayouts, mg, sh, &msel);

//	we will now fill a binary stream with all the grids.
//	this stream will receive the data that has to be copied to the local grid.
	BinaryStream localStream;
//	this stream will receive all the data that is to be sent to other processes.
	BinaryStream globalStream;
//	this vector is required so that we can use distribute-data later on.
	vector<int>	vBlockSizes;
//	here we'll store the ids of the receiving processes.
	vector<int> vReceiverIDs;

	int numProcs = (int)sh.num_subsets();
	if(pProcessMap)
		numProcs = std::min((int)pProcessMap->size(), numProcs);
		
	for(int i = 0; i < numProcs; ++i)
	{
/*
cout << "proc " << i << ":\n";
cout << "  layouts:\n";
cout << "    vrts: " << vVertexLayouts[i].node_vec().size() << endl;
cout << "    edges: " << vEdgeLayouts[i].node_vec().size() << endl;
cout << "    faces: " << vFaceLayouts[i].node_vec().size() << endl;
cout << "    vols: " << vVolumeLayouts[i].node_vec().size() << endl;
*/
		int proc = i;
		if(pProcessMap)
			proc = (*pProcessMap)[i];
		
		if(proc == localProcID)
		{
			SerializeGridAndLayouts(localStream, mg, vVertexLayouts[i],
									vEdgeLayouts[i], vFaceLayouts[i], vVolumeLayouts[i],
									aInt, aInt, aInt, aInt, &msel, pProcessMap);
									
		//	serialize position attachment
			for(uint iLevel = 0; iLevel < mg.num_levels(); ++iLevel)
			{
				SerializeAttachment<VertexBase>(mg, aPosition,
												msel.begin<VertexBase>(iLevel),
												msel.end<VertexBase>(iLevel),
												localStream);
			}

//TODO:		the user should be able to add personal data to those buffers.
		}
		else
		{
			int oldSize = globalStream.size();
			SerializeGridAndLayouts(globalStream, mg, vVertexLayouts[i],
									vEdgeLayouts[i], vFaceLayouts[i], vVolumeLayouts[i],
									aInt, aInt, aInt, aInt, &msel, pProcessMap);

		//	serialize position attachment
			for(uint iLevel = 0; iLevel < mg.num_levels(); ++iLevel)
			{
				SerializeAttachment<VertexBase>(mg, aPosition,
												msel.begin<VertexBase>(iLevel),
												msel.end<VertexBase>(iLevel),
												globalStream);
			}									
//TODO:		the user should be able to add personal data to those buffers.

			vBlockSizes.push_back((int)(globalStream.size() - oldSize));
			vReceiverIDs.push_back(proc);
		}
	}
	
//	send the grids to their target processes
	int numReceivers = (int)vReceiverIDs.size();
	if(numReceivers > 0)
	{
	//	every process receives the size of the data-buffer first.
		vector<int> bufferSizes(numReceivers, sizeof(int));
		
	//	distribute the block-sizes to the different processes
		pcl::DistributeData(localProcID, &vReceiverIDs.front(), numReceivers,
							&vBlockSizes.front(), &bufferSizes.front(), 38);

	//	distribute the grids-distribution-packs
		pcl::DistributeData(localProcID, &vReceiverIDs.front(), numReceivers,
							globalStream.buffer(), &vBlockSizes.front(), 39);
	}

//	fill the local grid and subset-handler
	if(pLocalGridOut && pLocalGridLayoutOut && (localStream.size() > 0))
	{		
		if(!pLocalGridOut->has_vertex_attachment(aPosition))
			pLocalGridOut->attach_to_vertices(aPosition);

		DeserializeGridAndLayouts(*pLocalGridOut, *pLocalGridLayoutOut,
									localStream);

		DeserializeAttachment<VertexBase>(*pLocalGridOut, aPosition,
										pLocalGridOut->begin<VertexBase>(),
										pLocalGridOut->end<VertexBase>(),
										localStream);		

	}
	
//	clean up
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);
	
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
void ReceiveGrid(MultiGrid& mgOut, ParallelGridLayout& gridLayoutOut,
					int srcProcID)
{
//	receive the stream-size
	int streamSize;
	pcl::ReceiveData(&streamSize, srcProcID, sizeof(int), 38);

//	receive the buffer
	BinaryStream binaryStream(streamSize);
	pcl::ReceiveData(binaryStream.buffer(), srcProcID, streamSize, 39);

//	fill the grid and the layout
	DeserializeGridAndLayouts(mgOut, gridLayoutOut, binaryStream);
	
//	read the attached data
	if(!mgOut.has_vertex_attachment(aPosition))
		mgOut.attach_to_vertices(aPosition);

	DeserializeAttachment<VertexBase>(mgOut, aPosition, binaryStream);

//TODO:	allow the user to read his data.
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
			layout.node_vec().push_back(node);
			
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
		msel.select(sh.begin<VertexBase>(i), sh.end<VertexBase>(i));
		msel.select(sh.begin<EdgeBase>(i), sh.end<EdgeBase>(i));
		msel.select(sh.begin<Face>(i), sh.end<Face>(i));
		msel.select(sh.begin<Volume>(i), sh.end<Volume>(i));
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

	msel.clear_selection();
	
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
	SerializeLayoutInterfaces(out, vrtLayout, pProcessMap);
	SerializeLayoutInterfaces(out, edgeLayout, pProcessMap);
	SerializeLayoutInterfaces(out, faceLayout, pProcessMap);
	SerializeLayoutInterfaces(out, volLayout, pProcessMap);
	
//	done. Please note that no attachments have been serialized in this method.
}

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

////////////////////////////////////////////////////////////////////////
//	DeserializeGridAndLayouts
void DeserializeGridAndLayouts(MultiGrid& mgOut,
							ParallelGridLayout& gridLayoutOut,
							std::istream& in)
{	
//	read the grid.
//	during deserialization the local indices are automatically generated
//	and written to the aLocalInd... attachments.
	DeserializeMultiGridElements(mgOut, in);

//	read the layouts
	DeserializeLayoutInterfaces(gridLayoutOut.vertexLayout, in);
	DeserializeLayoutInterfaces(gridLayoutOut.edgeLayout, in);
	DeserializeLayoutInterfaces(gridLayoutOut.faceLayout, in);
	DeserializeLayoutInterfaces(gridLayoutOut.volumeLayout, in);
	
//TODO: this is not optimal. Probably it is not even required...
//	add the nodes to the layouts
	FillLayoutWithNodes<VertexBase>(gridLayoutOut.vertexLayout, mgOut);
	FillLayoutWithNodes<EdgeBase>(gridLayoutOut.edgeLayout, mgOut);
	FillLayoutWithNodes<Face>(gridLayoutOut.faceLayout, mgOut);
	FillLayoutWithNodes<Volume>(gridLayoutOut.volumeLayout, mgOut);
	
//	done. Please note that no attachments have been serialized in this method.
}

}//	end of namespace
