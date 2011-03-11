//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#include "lib_grid/lib_grid.h"
#include "pcl/pcl.h"
#include "distribution_util.h"
#include "common/util/stream_pack.h"
#include "common/util/binary_stream.h"
#include "common/serialization.h"

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
template <class TVertexDistributionLayout, class TEdgeDistributionLayout,
		  class TFaceDistributionLayout, class TVolumeDistributionLayout>
void CreateDistributionLayouts(
						std::vector<TVertexDistributionLayout>& vertexLayoutsOut,
						std::vector<TEdgeDistributionLayout>& edgeLayoutsOut,
						std::vector<TFaceDistributionLayout>& faceLayoutsOut,
						std::vector<TVolumeDistributionLayout>& volumeLayoutsOut,
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

//	explicit template instantiation
template void CreateDistributionLayouts<DistributionVertexLayout, DistributionEdgeLayout,
										DistributionFaceLayout, DistributionVolumeLayout>(
										std::vector<DistributionVertexLayout>&,
										std::vector<DistributionEdgeLayout>&,
										std::vector<DistributionFaceLayout>&,
										std::vector<DistributionVolumeLayout>&,
										MultiGrid&, SubsetHandler&, bool, MGSelector*);

template void CreateDistributionLayouts<RedistributionVertexLayout, RedistributionEdgeLayout,
										RedistributionFaceLayout, RedistributionVolumeLayout>(
										std::vector<RedistributionVertexLayout>&,
										std::vector<RedistributionEdgeLayout>&,
										std::vector<RedistributionFaceLayout>&,
										std::vector<RedistributionVolumeLayout>&,
										MultiGrid&, SubsetHandler&, bool, MGSelector*);

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
///	Copies the associated global id from each interface entry into the redist-layouts.
/**	Make sure that aGlobalID contains valid global ids.*/
template <class TGeomObj>
void AssociateGlobalIDs(Grid& g,
						std::vector<RedistributionNodeLayout<TGeomObj*> >& layouts,
						AGeomObjID& aGlobalID)
{
	if(!g.has_attachment<TGeomObj>(aGlobalID)){
		g.attach_to<TGeomObj>(aGlobalID);
	}

	Grid::AttachmentAccessor<TGeomObj, AGeomObjID> aaID(g, aGlobalID);

	for(size_t i_layout = 0; i_layout < layouts.size(); ++i_layout){
		RedistributionNodeLayout<TGeomObj*>& layout = layouts[i_layout];
		const typename DistributionNodeLayout<TGeomObj*>::NodeVec& nodes =
												layout.node_vec();
		layout.m_globalIDs.resize(nodes.size());

		for(size_t i_node = 0; i_node < nodes.size(); ++i_node){
			layout.m_globalIDs[i_node] = aaID[nodes[i_node]];
		}
	}
}

////////////////////////////////////////////////////////////////////////
///	Collects all target processes for each entry in the grid.
/**	For each element of the grid of type TGeomObj, the target processes can
 * be found in the attachment to which TAAIntVec points, after the algorithm
 * is done.
 *
 * Note that the current process itself is not regarded as a target proc.
 *
 * Make sure that TGeomObj is either VertexBase, EdgeBase, Face or Volume.
 * Make also sure that TAAIntVec is accessing a valid Attachment<vector<int> >
 * attachment at grid.
 */
template <class TGeomObj, class TAAIntVec>
static void
CollectRedistributionTargetProcs(
						vector<RedistributionNodeLayout<TGeomObj> >& layouts,
						TAAIntVec& aaIntVec, int localLayoutIndex, int* processMap)
{
	typedef typename DistributionNodeLayout<TGeomObj>::NodeVec		NodeVec;
	typedef typename DistributionNodeLayout<TGeomObj>::Interface	Interface;
	typedef typename DistributionNodeLayout<TGeomObj>::InterfaceMap	InterfaceMap;

	for(size_t i_layout = 0; i_layout < layouts.size(); ++i_layout){
		if((int)i_layout == localLayoutIndex)
			continue;

		int targetProc = i_layout;
		if(processMap)
			targetProc = processMap[i_layout];

	//	get the nodes of the current layout
		NodeVec& nodes = layouts[i_layout].node_vec();

	//	now iterate over the nodes and add the target-proc
		for(size_t i = 0; i < nodes.size(); ++i)
			aaIntVec[nodes[i]].push_back(targetProc);
	}
}

////////////////////////////////////////////////////////////////////////
///	A highly specialized communication policy used during grid redistribution
/**	All nodes that will reside on the local process have to be selected
 * in the specified selector.
 *
 * TAAIntVec has to be an attachment accessor to a vector<int> attachment and
 * TAATransferInfoVec has to be an attachment accessor to a
 * vector<RedistributionNodeTransferInfo> attachment.
 *
 * The fist accessor is the source-accessor from which data will be taken,
 * the second one is the target to which data will be written.
 *
 * Make sure that both accessors work with the elements in the given layout.
 */
template <class TLayout, class TAAIntVec, class TAATransferInfoVec>
class ComPol_SynchronizeNodeTransfer : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename Interface::iterator	InterfaceIter;

	///	Construct the communication policy with a ug::Selector.
	/**	Through the parameters select and deselect one may specify whether
	 * a process selects and/or deselects elements based on the received
	 * selection status.*/
		ComPol_SynchronizeNodeTransfer(ISelector& sel, TAAIntVec& aaIntVec,
										TAATransferInfoVec& aaTransInfoVec)
			 :	m_sel(sel), m_aaIntVec(aaIntVec), m_aaTransInfoVec(aaTransInfoVec)
		{}

		virtual int
		get_required_buffer_size(Interface& interface)		{return -1;}

	///	write target processes and move-flag
		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
			byte bTrue = 1; byte bFalse = 0;

		//	write the entry indices of selected elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
			//	first write whether the elem is moved away from this proc
				if(m_sel.is_selected(elem))
					buff.write((char*)&bFalse, sizeof(byte));
				else
					buff.write((char*)&bTrue, sizeof(byte));

			//	now write the target processes
				Serialize(buff, m_aaIntVec[elem]);
			}

			return true;
		}

	///	read target processes and move-flag
		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
			byte bMove;
			int srcProc = interface.get_target_proc();
			vector<int> targetProcs;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&bMove, sizeof(byte));
				targetProcs.clear();
				Deserialize(buff, targetProcs);

				for(size_t i = 0; i < targetProcs.size(); ++i){
					m_aaTransInfoVec[elem].push_back(
						RedistributionNodeTransferInfo(srcProc, targetProcs[i],
														bMove != 0));
				}
			}
			return true;
		}

	protected:
		ISelector&			m_sel;
		TAAIntVec&			m_aaIntVec;
		TAATransferInfoVec&	m_aaTransInfoVec;
};


////////////////////////////////////////////////////////////////////////
///	This method distributes information between neighbor procs on where which node is send.
/**	All nodes that will reside on the local process have to be selected
 * in the specified selector.
 *
 * Synchronization of the node movement is required, since
 * neighbored processes may redistribute connected parts of a grid at the
 * same time. They thus have to inform associated processes about those
 * changes. This information can be used during creation of the distribution
 * interfaces.
 */
template <class TGeomObj, class TAAIntVec, class TAATransferInfoVec>
static void
SynchronizeNodeTransfer(Grid& g, GridLayoutMap& glm,
						pcl::ParallelCommunicator<typename GridLayoutMap::Types<TGeomObj>::Layout> & comm,
						TAAIntVec& aaIntVec,
						TAATransferInfoVec& aaTransInfoVec,
						int localLayoutIndex,
						ISelector& sel,
						int* processMap)
{
	typedef typename GridLayoutMap::Types<TGeomObj>		Types;
	typedef typename Types::Layout		Layout;
	typedef typename Types::Interface	Interface;

//	this communication policy fills aaTransInfoVec
	ComPol_SynchronizeNodeTransfer<Layout, TAAIntVec, TAATransferInfoVec>
		compol(sel, aaIntVec, aaTransInfoVec);

//	exchange data, both from masters to slaves and vice versa
	comm.exchange_data(glm, INT_MASTER, INT_SLAVE, compol);
	comm.exchange_data(glm, INT_SLAVE, INT_MASTER, compol);
	comm.communicate();
}

////////////////////////////////////////////////////////////////////////
/**	\brief iterates through existing interfaces and adds interface entries to
 *			the redistribution layouts.
 *
 * All nodes that will reside on the local process have to be selected
 * in the specified selector.
 *
 * Note that the interfaces of the resulting redistribution layouts will not
 * match the original interfaces in the distGridMgr. This is because node-
 * movement on other processes is already considered. The resulting
 * interfaces are final, which means that they exactly describe the connections
 * which will be created on other processes.
 *
 * Note that the redistribution layouts already contain interfaces to the
 * target processes to which parts of the grid are sent. Those interfaces
 * may be enhanced, but existing entries won't be changed nor removed.
 */
template <class TDistLayout, class TAAProcTargets,
		  class TAATransferInfos>
static
void FinalizeRedistributionLayoutInterfaces(
					DistributedGridManager& distGridMgr,
					std::vector<TDistLayout>& distLayoutVec,
					TAAProcTargets& aaTargetProcs,
					TAATransferInfos& aaTransferInfos,
					ISelector& sel,
					int* processMap = NULL)
{
	typedef typename TDistLayout::NodeType 	Node;
	typedef typename TDistLayout::Interface	Interface;

//	access the associated multi-grid
	if(!distGridMgr.get_assigned_grid())
		return;

	MultiGrid& mg = *distGridMgr.get_assigned_grid();

//	we'll use this vector to gather existing interface-entries for a node
	vector<pair<int, size_t> > interfaceEntries;

//	iterate over all nodes in the distribution layouts
	for(size_t i_layout = 0; i_layout < distLayoutVec.size(); ++i_layout)
	{
		int targetProcID = i_layout;
		if(processMap)
			targetProcID = processMap[i_layout];

		TDistLayout& layout = distLayoutVec[i_layout];
		const typename TDistLayout::NodeVec& nodes = layout.node_vec();

	//	iterate over all nodes
		for(size_t i_node = 0; i_node < nodes.size(); ++i_node){
		//	Node is a VertexBase*, EdgeBase*, Face* or Volume*
			const Node& node = nodes[i_node];

		//	only coninue for interface nodes
			if(!distGridMgr.is_interface_element(node))
				continue;

		//	the node lies in at least one interface to an old neighbor.
		//	We have to add the entry to corresponding distribution-interfaces
			byte nodeStatus = distGridMgr.get_status(node);
			distGridMgr.collect_interface_entries(interfaceEntries, node);

		//todo: cache interfaces
		//	target procs of the current node
			const vector<int>& targetProcs = aaTargetProcs[node];
			const vector<RedistributionNodeTransferInfo>& transferInfos =
													aaTransferInfos[node];

		//	copyAll means that slave interfaces to all associated nodes
		//	are created in the target redistribution layout.
			bool copyAll = false;

		//	if the node is a master node and a copy is created on the target
		//	process, then we only have to create a slave interface on target
		//	process to this master.
		//todo: consider horizontal and virtual interfaces


//todo:	FIX IT!


			if((nodeStatus & ES_MASTER)){
				if(sel.is_selected(node)){
					UG_LOG("typ1 -");
				//	a copy is created. This is already handled by CreateDistributionLayouts.
				/*
					Interface& interface = layout.interface(localProcID,
														 mg.get_level(node));
					interface.push_back(DistributionInterfaceEntry(i_node,
																INT_SLAVE));
				*/
				}
				else{
				//	check whether the target-process is the new master.
				//	The new master is always the target proc with the lowest rank.
					int newMaster = pcl::GetNumProcesses();
					for(size_t i = 0; i < targetProcs.size(); ++i){
						if(targetProcs[i] < newMaster){
							newMaster = targetProcs[i];
						}
					}

					if(newMaster == targetProcID){
						for(size_t i_int = 0; i_int < interfaceEntries.size(); ++i_int){
							int assProc = interfaceEntries[i_int].first;
							bool createDefault = true;

						//	check whether the associated node on the associated process
						//	is also sent somewhere else.
							for(size_t i = 0; i < transferInfos.size(); ++i){
								const RedistributionNodeTransferInfo& info = transferInfos[i];
								if(info.srcProc == assProc){
									UG_LOG("typ2_1 -");
								//	we have to create a new interface to the infos targetProc
									Interface& interface = layout.interface(info.targetProc,
																			 	 mg.get_level(node));
									interface.push_back(DistributionInterfaceEntry(i_node, INT_MASTER));
									if(info.bMove)
										createDefault = false;
								}
							}

						//	create the connection to assProc, if the associated node resides there
							if(createDefault){
								UG_LOG("typ2_2 -");
								Interface& interface = layout.interface(assProc,
																			mg.get_level(node));
								interface.push_back(DistributionInterfaceEntry(i_node, INT_MASTER));
							}
						}
					}
					else
						copyAll = true;
				}
			}
			else if(nodeStatus & ES_SLAVE)
				copyAll = true;
			else{
				throw(UGError("Only ES_MASTER and ES_SLAVE are supported during "
							  "redistribution in the current implementation!"));
			}

			if(copyAll){
				for(size_t i_int = 0; i_int < interfaceEntries.size(); ++i_int){
					int assProc = interfaceEntries[i_int].first;
					bool createDefault = true;
					int newMasterRank = pcl::GetProcRank();

				//	check whether the associated node on the associated process
				//	is also sent somewhere else.
					for(size_t i = 0; i < transferInfos.size(); ++i){
						const RedistributionNodeTransferInfo& info = transferInfos[i];
						if(info.srcProc == assProc){
						//	we have to create a new interface to the infos targetProc only
						//	if the node moves. If not, the old master resides where it was.
							if(info.bMove){
							//	the master node moves. The new master is the process with the
							//	lowest rank
								newMasterRank = min(newMasterRank, info.targetProc);
								createDefault = false;
							}
						}
					}

					if(createDefault){
						UG_LOG("typ3_1 -");
					//	create the connection to assProc, if the associated master node resides there
						Interface& interface = layout.interface(assProc,
																	mg.get_level(node));
						interface.push_back(DistributionInterfaceEntry(i_node, INT_SLAVE));
					}
					else{
						UG_LOG("typ3_2 -");
					//	The master moves. Create an interface to the new master
						Interface& interface = layout.interface(newMasterRank,
																 mg.get_level(node));
						interface.push_back(DistributionInterfaceEntry(i_node, INT_SLAVE));
					}
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
void CreateRedistributionLayouts(
					std::vector<RedistributionVertexLayout>& vertexLayoutsOut,
					std::vector<RedistributionEdgeLayout>& edgeLayoutsOut,
					std::vector<RedistributionFaceLayout>& faceLayoutsOut,
					std::vector<RedistributionVolumeLayout>& volumeLayoutsOut,
					DistributedGridManager& distGridMgr, SubsetHandler& sh,
					bool distributeGenealogy, MGSelector* pSel, int* processMap,
					Attachment<std::vector<int> >* paTargetProcs,
					ARedistributionNodeTransferInfoVec* paTransferInfoVec,
					AGeomObjID& aGlobalID)
{
//	access the associated multi-grid
	if(!distGridMgr.get_assigned_grid())
		return;

	MultiGrid& mg = *distGridMgr.get_assigned_grid();

	MGSelector tmpSel;
	if(!pSel){
		tmpSel.assign_grid(mg);
		pSel = &tmpSel;
	}
	MGSelector& msel = *pSel;

//	first we'll create normal distribution layouts
	CreateDistributionLayouts(vertexLayoutsOut, edgeLayoutsOut,
							  faceLayoutsOut, volumeLayoutsOut,
							  mg, sh, distributeGenealogy, &msel);

//	gather global ids for each node
	AssociateGlobalIDs<VertexBase>(mg, vertexLayoutsOut, aGlobalID);
	AssociateGlobalIDs<EdgeBase>(mg, edgeLayoutsOut, aGlobalID);
	AssociateGlobalIDs<Face>(mg, faceLayoutsOut, aGlobalID);
	AssociateGlobalIDs<Volume>(mg, volumeLayoutsOut, aGlobalID);

//	we have to know which redistribution layout corresponds with this process
	int localLayoutInd = -1;
	if(processMap){
	//	find the local process in the processMap
		int myRank = pcl::GetProcRank();
		for(int i = 0; i < sh.num_subsets(); ++i){
			if(processMap[i] == myRank){
				localLayoutInd = i;
				break;
			}
		}
	}
	else{
		localLayoutInd = pcl::GetProcRank();
	}

	if(localLayoutInd >= sh.num_subsets())
		localLayoutInd = -1;

	msel.clear();
	if(localLayoutInd != -1){
		SelectNodesInLayout(msel, vertexLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, edgeLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, faceLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, volumeLayoutsOut[localLayoutInd]);
	}

//	to each node we will now attach a vector<int>, into which we'll write all
//	processes to which the node will be sent.
//	we'll use the autoattach mechanism of the AttachmentAccessor.
	typedef Attachment<vector<int> > AIntVec;
	AIntVec aTargetProcs;
	if(paTargetProcs)
		aTargetProcs = *paTargetProcs;

	Grid::AttachmentAccessor<VertexBase, AIntVec>
		aaTargetProcsVRT(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<EdgeBase, AIntVec>
		aaTargetProcsEDGE(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<Face, AIntVec>
		aaTargetProcsFACE(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<Volume, AIntVec>
		aaTargetProcsVOL(mg, aTargetProcs, true);

	CollectRedistributionTargetProcs(vertexLayoutsOut, aaTargetProcsVRT,
									 localLayoutInd, processMap);
	CollectRedistributionTargetProcs(edgeLayoutsOut, aaTargetProcsEDGE,
									 localLayoutInd, processMap);
	CollectRedistributionTargetProcs(faceLayoutsOut, aaTargetProcsFACE,
									 localLayoutInd, processMap);
	CollectRedistributionTargetProcs(volumeLayoutsOut, aaTargetProcsVOL,
									 localLayoutInd, processMap);

//	we need a second attachment, which tells us for each node, to where
//	associated nodes on other processes go. We store this through a
//	vector<RedistributionNodeMoveInfo>, which tells us on which
//	process the copy/move operation was scheduled (srcProc) and to which
//	process the element will be transfered (targetProc) and which operation shall be
//	performed (bMove).
	typedef ARedistributionNodeTransferInfoVec ATransferInfoVec;
	ATransferInfoVec aTransferInfoVec;
	if(paTransferInfoVec)
		aTransferInfoVec = *paTransferInfoVec;

	Grid::AttachmentAccessor<VertexBase, ATransferInfoVec>
		aaTransInfoVecVRT(mg, aTransferInfoVec, true);
	Grid::AttachmentAccessor<EdgeBase, ATransferInfoVec>
		aaTransInfoVecEDGE(mg, aTransferInfoVec, true);
	Grid::AttachmentAccessor<Face, ATransferInfoVec>
		aaTransInfoVecFACE(mg, aTransferInfoVec, true);
	Grid::AttachmentAccessor<Volume, ATransferInfoVec>
		aaTransInfoVecVOL(mg, aTransferInfoVec, true);

//	the synchronization is performed in SynchronizeNodeTransfer.
//	It fills data associated with aTransferInfoVec.
//todo: the communicators could also be declared in SynchronizeNodeTransfer.
	pcl::ParallelCommunicator<VertexLayout> commVRT;
	pcl::ParallelCommunicator<EdgeLayout> commEDGE;
	pcl::ParallelCommunicator<FaceLayout> commFACE;
	pcl::ParallelCommunicator<VolumeLayout> commVOL;

	GridLayoutMap& glm = distGridMgr.grid_layout_map();
	SynchronizeNodeTransfer<VertexBase>(mg, glm, commVRT, aaTargetProcsVRT,
							aaTransInfoVecVRT, localLayoutInd, msel, processMap);
	SynchronizeNodeTransfer<EdgeBase>(mg, glm, commEDGE, aaTargetProcsEDGE,
							aaTransInfoVecEDGE, localLayoutInd, msel, processMap);
	SynchronizeNodeTransfer<Face>(mg, glm, commFACE, aaTargetProcsFACE,
							aaTransInfoVecFACE, localLayoutInd, msel, processMap);
	SynchronizeNodeTransfer<Volume>(mg, glm, commVOL, aaTargetProcsVOL,
							aaTransInfoVecVOL, localLayoutInd, msel, processMap);


//	finalize the redistribution layouts by constructing all interfaces
	FinalizeRedistributionLayoutInterfaces(distGridMgr, vertexLayoutsOut,
										aaTargetProcsVRT, aaTransInfoVecVRT,
										msel, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, edgeLayoutsOut,
										aaTargetProcsEDGE, aaTransInfoVecEDGE,
										msel, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, faceLayoutsOut,
										aaTargetProcsFACE, aaTransInfoVecFACE,
										msel, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, volumeLayoutsOut,
										aaTargetProcsVOL, aaTransInfoVecVOL,
										msel, processMap);

// detach temporary attachments, only if they were not specified from outside.
	if(!paTargetProcs)
		mg.detach_from_all(aTargetProcs);
	if(!paTransferInfoVec)
		mg.detach_from_all(aTransferInfoVec);
}
/*
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
				int conProcID = mapIter->first;
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
				int conProcID = mapIter->first;
				if(conProcID == curProcID)
					continue;

				UG_LOG("  connections " << curProcID << " - " << conProcID << ":");

				Interface& curIntf = mapIter->second;

			//	make sure that the different interfaces match each other in size
				size_t numCurMasters = NumEntriesOfTypeInDistributionInterface(
															INT_MASTER, curIntf);
				size_t numCurSlaves = NumEntriesOfTypeInDistributionInterface(
															INT_SLAVE, curIntf);
/*
				size_t numCurVrtMasters = NumEntriesOfTypeInDistributionInterface(
													INT_VERTICAL_MASTER, curIntf);
				size_t numCurVrtSlaves = NumEntriesOfTypeInDistributionInterface(
													INT_VERTICAL_SLAVE, curIntf);
*/
				if(numCurMasters){
					UG_LOG("    masters: " << numCurMasters);
				}

				if(numCurSlaves){
					UG_LOG("    slaves: " << numCurSlaves);
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

}//	end of namespace
