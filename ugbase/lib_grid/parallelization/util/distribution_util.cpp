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
			typename TNodeLayout::Interface& masterInterface = masterLayout.interface(layoutIndex, level);
			typename TNodeLayout::Interface& slaveInterface = layout.interface(masterLayoutIndex, level);
			masterInterface.push_back(typename TNodeLayout::InterfaceEntry(localMasterID, INT_MASTER));
			slaveInterface.push_back(typename TNodeLayout::InterfaceEntry(localID, INT_SLAVE));
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

}//	end of namespace
