//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d08

#include "grid_distribution.h"
#include "lib_grid/lib_grid.h"
#include "distribution_node_layout.h"
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
//	DistributeGrid
void DistributeGrid(MultiGrid& mg, SubsetHandler& sh,
					int localProcID, MultiGrid* pLocalGridOut,
					ParallelGridLayout* pLocalGridLayoutOut,
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
}


////////////////////////////////////////////////////////////////////////
//	ReceiveGrid
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

}//	end of namespace
