//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d08

#include "grid_distribution.h"
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
/**	All nodes in the given distribution layout will be added to
 *	the INT_VERTICAL_MASTER interface of layoutMapOut.*/
template <class TGeomObj>
static void AddVerticalMasterInterfaces(GridLayoutMap& layoutMapOut,
										MultiGrid& mg,
								  		const DistributionNodeLayout<TGeomObj*>& distLayout,
								  		int targetProc)
{
//	some typedefs
	typedef DistributionNodeLayout<TGeomObj*> 			DistLayout;
	typedef typename GridLayoutMap::Types<TGeomObj>::
								Layout::LevelLayout 	TLayout;
	typedef typename TLayout::Interface					TInterface;
	
//	access the nodes of the distribution layout
	const typename DistLayout::NodeVec& nodeVec = distLayout.node_vec();

//	we'll cache the layout and interface.
	TLayout* pLayout = NULL;
	TInterface* pInterface = NULL;
	int currentLevel = -1;

//	iterate over the nodes
	for(size_t i = 0; i < nodeVec.size(); ++i){
		TGeomObj* node = nodeVec[i];
	//	get the level and the matching layout and interface
		int level = (int)mg.get_level(node);
		if(level != currentLevel){
			currentLevel = level;
		//	get the layout and the interface to targetProc
			pLayout = &layoutMapOut.template get_layout<TGeomObj>(INT_VERTICAL_MASTER).
														layout_on_level(level);
			pInterface = &pLayout->interface(targetProc);
		}
		
	//	add node to the interface
		pInterface->push_back(node);
	}
}

template <class TGeomObj>
static void AddHorizontalInterfaces(GridLayoutMap& layoutMapOut,
								DistributionNodeLayout<TGeomObj*>& distLayout,
								std::vector<int>* pProcessMap)
{
//	some typedefs
	typedef DistributionNodeLayout<TGeomObj*> 			DistLayout;
	typedef typename GridLayoutMap::Types<TGeomObj>::
								Layout::LevelLayout 	TLayout;
	typedef typename TLayout::Interface					TInterface;
	
//	access the nodes of the distribution layout
	const typename DistLayout::NodeVec& nodeVec = distLayout.node_vec();

//	we'll cache the layout and interface.
	TLayout* pLayout = NULL;
	TInterface* pInterface = NULL;
	int elemType = -1;

//	iterate over the nodes
	for(size_t level = 0; level < distLayout.num_levels(); ++level)
	{
		typename DistLayout::InterfaceMap& imap = distLayout.interface_map(level);
		for(typename DistLayout::InterfaceMap::iterator iter = imap.begin();
			iter != imap.end(); ++iter)
		{
		//	get the proc-id of this interface
			int procID = iter->first;
			if(pProcessMap)
			{
				assert((int)pProcessMap->size() > procID && "process-map to small.");
				if((int)pProcessMap->size() > procID)
					procID = (*pProcessMap)[procID];
			}
			
		//	copy the interface-elements to the layoutMap
			typename DistLayout::Interface& distInterface = iter->second;

		//	iterate over the nodes in the interface
			for(size_t i = 0; i < distInterface.size(); ++i){
				int newElemType = distInterface[i].type;
				TGeomObj* node = nodeVec[distInterface[i].localID];
			//	get the layout and the matching interface
				if(elemType != newElemType){
					elemType = newElemType;
				//	get the layout and the interface to targetProc
					pLayout = &layoutMapOut.template get_layout<TGeomObj>(elemType).
																layout_on_level(level);
					pInterface = &pLayout->interface(procID);
				}
				
			//	add the node to the local interface.
				pInterface->push_back(node);
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////
//	DistributeGrid_KeepSrcGrid
bool DistributeGrid_KeepSrcGrid(MultiGrid& mg, ISubsetHandler& sh,
								GridLayoutMap& layoutMap,
								SubsetHandler& shPartition,
								int localProcID,
								std::vector<int>* pProcessMap)
{

//	we have to store the layouts for all the processes.
	vector<DistributionVertexLayout> vVertexLayouts;
	vector<DistributionEdgeLayout> vEdgeLayouts;
	vector<DistributionFaceLayout> vFaceLayouts;
	vector<DistributionVolumeLayout> vVolumeLayouts;
	
//	we need some attachments that will speed up the called processes.
	AInt aInt;
	mg.attach_to_vertices(aInt);
	mg.attach_to_edges(aInt);
	mg.attach_to_faces(aInt);
	mg.attach_to_volumes(aInt);
	
//	the selector will help to speed things up a little.
	MGSelector msel(mg);
	
//	if vertical interfaces shall be created then we won't distribute
//	the whole genealogy
	CreateDistributionLayouts(vVertexLayouts, vEdgeLayouts, vFaceLayouts,
							  vVolumeLayouts, mg, shPartition,
							  false, &msel);
	
//	we will now fill a binary stream with all the grids.
//	this stream will receive all the data that is to be sent to other processes.
	BinaryStream globalStream;
//	this vector is required so that we can use distribute-data later on.
	vector<int>	vBlockSizes;
//	here we'll store the ids of the receiving processes.
	vector<int> vReceiverIDs;

	int numProcs = (int)shPartition.num_subsets();
	if(pProcessMap)
		numProcs = std::min((int)pProcessMap->size(), numProcs);
		
	for(int i = 0; i < numProcs; ++i)
	{
		int proc = i;
		if(pProcessMap)
			proc = (*pProcessMap)[i];
		
		if(proc == localProcID)
		{
		//	create local horizontal interfaces
			AddHorizontalInterfaces<VertexBase>(layoutMap, vVertexLayouts[i], pProcessMap);
			AddHorizontalInterfaces<EdgeBase>(layoutMap, vEdgeLayouts[i], pProcessMap);
			AddHorizontalInterfaces<Face>(layoutMap, vFaceLayouts[i], pProcessMap);
			AddHorizontalInterfaces<Volume>(layoutMap, vVolumeLayouts[i], pProcessMap);
		}
		else
		{
		//	since the original grid is kept, we have to add vertical interfaces.
		//	all nodes in the layouts will be vertical interface members.
		//	Here we create the local vertical interfaces.
			AddVerticalMasterInterfaces<VertexBase>(layoutMap, mg,
													vVertexLayouts[i], proc);
			AddVerticalMasterInterfaces<EdgeBase>(layoutMap, mg,
												  vEdgeLayouts[i], proc);
			AddVerticalMasterInterfaces<Face>(layoutMap, mg,
											  vFaceLayouts[i], proc);
			AddVerticalMasterInterfaces<Volume>(layoutMap, mg,
											  vVolumeLayouts[i], proc);

		//	serialize the parts of the grid that will be send to other processes.
			int oldSize = globalStream.size();
			SerializeGridAndDistributionLayouts(
									globalStream, mg, vVertexLayouts[i],
									vEdgeLayouts[i], vFaceLayouts[i], vVolumeLayouts[i],
									aInt, aInt, aInt, aInt, &msel, pProcessMap);

		//	serialize subset indices
			SerializeSubsetHandler(mg, sh,
								   msel.get_multi_level_geometric_object_collection(),
								   globalStream);
								   
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
	
//	clean up
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);
	
	return true;
}

////////////////////////////////////////////////////////////////////////
//	DistributeGrid
bool DistributeGrid(MultiGrid& mg, ISubsetHandler& sh,
					SubsetHandler& shPartition,
					int localProcID, MultiGrid* pLocalGridOut,
					ISubsetHandler* pLocalSHOut,
					GridLayoutMap* pLocalGridLayoutMapOut,
					std::vector<int>* pProcessMap)
{

//	we have to store the layouts for all the processes.
	vector<DistributionVertexLayout> vVertexLayouts;
	vector<DistributionEdgeLayout> vEdgeLayouts;
	vector<DistributionFaceLayout> vFaceLayouts;
	vector<DistributionVolumeLayout> vVolumeLayouts;
	
//	we need some attachments that will speed up the called processes.
	AInt aInt;
	mg.attach_to_vertices(aInt);
	mg.attach_to_edges(aInt);
	mg.attach_to_faces(aInt);
	mg.attach_to_volumes(aInt);
	
//	the selector will help to speed things up a little.
	MGSelector msel(mg);

	CreateDistributionLayouts(vVertexLayouts, vEdgeLayouts, vFaceLayouts,
							  vVolumeLayouts, mg, shPartition,
							  true, &msel);
	
//	we will now fill a binary stream with all the grids.
//	this stream will receive the data that has to be copied to the local grid.
	BinaryStream localStream;
//	this stream will receive all the data that is to be sent to other processes.
	BinaryStream globalStream;
//	this vector is required so that we can use distribute-data later on.
	vector<int>	vBlockSizes;
//	here we'll store the ids of the receiving processes.
	vector<int> vReceiverIDs;

	int numProcs = (int)shPartition.num_subsets();
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
			SerializeGridAndDistributionLayouts(
								localStream, mg, vVertexLayouts[i],
								vEdgeLayouts[i], vFaceLayouts[i], vVolumeLayouts[i],
								aInt, aInt, aInt, aInt, &msel, pProcessMap);
		
		//	serialize subset indices
			SerializeSubsetHandler(mg, sh,
								   msel.get_multi_level_geometric_object_collection(),
								   localStream);
								   
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
			SerializeGridAndDistributionLayouts(
									globalStream, mg, vVertexLayouts[i],
									vEdgeLayouts[i], vFaceLayouts[i], vVolumeLayouts[i],
									aInt, aInt, aInt, aInt, &msel, pProcessMap);

		//	serialize subset indices
			SerializeSubsetHandler(mg, sh,
								   msel.get_multi_level_geometric_object_collection(),
								   globalStream);
								   
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
	if(pLocalGridOut && pLocalSHOut && 
		pLocalGridLayoutMapOut && (localStream.size() > 0))
	{		
		if(!pLocalGridOut->has_vertex_attachment(aPosition))
			pLocalGridOut->attach_to_vertices(aPosition);

		DeserializeGridAndDistributionLayouts(
									*pLocalGridOut, *pLocalGridLayoutMapOut,
									localStream);

	//	deserialize subset handler
		if(!DeserializeSubsetHandler(*pLocalGridOut, *pLocalSHOut,
									localStream))
		{
			goto bailout_false;
		}
					
		if(!DeserializeAttachment<VertexBase>(*pLocalGridOut, aPosition,
										pLocalGridOut->begin<VertexBase>(),
										pLocalGridOut->end<VertexBase>(),
										localStream))
		{
			goto bailout_false;
		}
	}
	
//	clean up
//bailout_true:
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);
	
	return true;

bailout_false:
	mg.detach_from_vertices(aInt);
	mg.detach_from_edges(aInt);
	mg.detach_from_faces(aInt);
	mg.detach_from_volumes(aInt);
	
	return false;

}


////////////////////////////////////////////////////////////////////////
/**	All nodes in the given distribution layout will be added to
 *	the INT_VERTICAL_MASTER interface of layoutMapOut.
 *	This method is a little unsafe...
 *	It is absolutly crucial, that all elements have been created
 *	in the same order as they were received through the in-stream
 *	in ReceiveInitialGrid. Especially auto-create options can be
 *	problematic here. Make sure to avoid auto-element-creation.*/
template <class TGeomObj>
static void AddVerticalSlaveInterfaces(GridLayoutMap& layoutMapOut,
									   MultiGrid& mg,
									   int srcProc)
{
//	some typedefs
	typedef typename GridLayoutMap::Types<TGeomObj>::
								Layout::LevelLayout 	TLayout;
	typedef typename TLayout::Interface					TInterface;
	typedef typename geometry_traits<TGeomObj>::iterator	ObjIter;
	
//	iterate over the levels of the multi-grid
	for(size_t level = 0; level < mg.num_levels(); ++level){
	//	if there are no nodes in this level, we can return immediatly
		if(mg.num<TGeomObj>(level) == 0)
			break;
	//	get the appropriate layout and interface
		TLayout& layout = layoutMapOut.template
							get_layout<TGeomObj>(INT_VERTICAL_SLAVE).
								layout_on_level(level);
								
		TInterface& interface = layout.interface(srcProc);
		
	//	iterate over the elements of the level
		for(ObjIter iter = mg.begin<TGeomObj>(level);
			iter != mg.end<TGeomObj>(level); ++iter)
		{
			interface.push_back(*iter);
		}
	}
}

////////////////////////////////////////////////////////////////////////
//	ReceiveGrid
bool ReceiveGrid(MultiGrid& mgOut, ISubsetHandler& shOut,
				 GridLayoutMap& gridLayoutMapOut,
				 int srcProcID, bool createVerticalLayouts)
{
//	receive the stream-size
	int streamSize;
	pcl::ReceiveData(&streamSize, srcProcID, sizeof(int), 38);

//	receive the buffer
	BinaryStream binaryStream(streamSize);
	pcl::ReceiveData(binaryStream.buffer(), srcProcID, streamSize, 39);

//	fill the grid and the layout
	DeserializeGridAndDistributionLayouts(mgOut, gridLayoutMapOut,
											binaryStream);
	
//	if vertical layouts shall be created, do it now.
	if(createVerticalLayouts){
		AddVerticalSlaveInterfaces<VertexBase>(gridLayoutMapOut, mgOut, srcProcID);
		AddVerticalSlaveInterfaces<EdgeBase>(gridLayoutMapOut, mgOut, srcProcID);
		AddVerticalSlaveInterfaces<Face>(gridLayoutMapOut, mgOut, srcProcID);
		AddVerticalSlaveInterfaces<Volume>(gridLayoutMapOut, mgOut, srcProcID);
	}
	
//	deserialize subset handler
	if(!DeserializeSubsetHandler(mgOut, shOut, binaryStream))
		return false;
	
//	read the attached data
	if(!mgOut.has_vertex_attachment(aPosition))
		mgOut.attach_to_vertices(aPosition);

	if(!DeserializeAttachment<VertexBase>(mgOut, aPosition, binaryStream))
		return false;

//TODO:	allow the user to read his data.

	return true;
}

}//	end of namespace
