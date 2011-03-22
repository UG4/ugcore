// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 15.03.2011 (m,d,y)
 
#include "lib_grid/lib_grid.h"
#include "pcl/pcl.h"
#include "distribution_util.h"
#include "common/util/stream_pack.h"
#include "common/util/binary_stream.h"
#include "common/serialization.h"

using namespace std;

namespace ug{

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
		//if((int)i_layout == localLayoutIndex)
		//	continue;

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
 * This method uses Grid::mark.
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
template <class TGeomObj,
		  class TAAProcTargets, class TAATransferInfos>
static
void FinalizeRedistributionLayoutInterfaces(
					DistributedGridManager& distGridMgr,
					std::vector<RedistributionNodeLayout<TGeomObj*> >& distLayoutVec,
					TAAProcTargets& aaTargetProcs,
					TAATransferInfos& aaTransferInfos,
					ISelector& sel,
					int* processMap = NULL)
{
	typedef RedistributionNodeLayout<TGeomObj*> TDistLayout;
	typedef typename TDistLayout::NodeType 	Node;
	typedef typename TDistLayout::Interface	DistInterface;

	typedef typename GridLayoutMap::Types<TGeomObj>	Types;
	typedef typename Types::Layout 		Layout;
	typedef typename Layout::Interface 	Interface;
	typedef typename Layout::iterator	IntfIter;
	typedef typename Interface::iterator IntfElemIter;

//	access the associated multi-grid
	if(!distGridMgr.get_assigned_grid())
		return;

	int localProcRank = pcl::GetProcRank();

	MultiGrid& mg = *distGridMgr.get_assigned_grid();

	Layout* masters = NULL;
	if(distGridMgr.grid_layout_map().template has_layout<TGeomObj>(INT_MASTER))
		masters = &distGridMgr.grid_layout_map().template get_layout<TGeomObj>(INT_MASTER);

	Layout* slaves = NULL;
	if(distGridMgr.grid_layout_map().template has_layout<TGeomObj>(INT_SLAVE))
		slaves = &distGridMgr.grid_layout_map().template get_layout<TGeomObj>(INT_SLAVE);

//	store the localLayoutIndex on the fly (the layout index of the current proc).
//	int localLayoutIndex = -1;

//	we have to attach integers to the elements, in which we'll store the
//	redistribution layouts node indices.
	AInt aNodeInd;
	Grid::AttachmentAccessor<TGeomObj, AInt> aaNodeInd(mg, aNodeInd, true);

//	we'll mark the nodes of the currently processed redistribution layout
	mg.begin_marking();

//	iterate over all nodes in the distribution layouts
	for(size_t i_layout = 0; i_layout < distLayoutVec.size(); ++i_layout)
	{
		int targetProc = i_layout;
		if(processMap)
			targetProc = processMap[i_layout];

		TDistLayout& distLayout = distLayoutVec[i_layout];
		const typename TDistLayout::NodeVec& nodes = distLayout.node_vec();

	//	mark all nodes in the current layout and assign their redist-layout node index
		mg.clear_marks();
		for(size_t i = 0; i < nodes.size(); ++i){
			aaNodeInd[nodes[i]] = (int)i;
			mg.mark(nodes[i]);
		}

	//	iterate over all master nodes of this process. If the node is marked we'll
	//	have to further process it.
		if(masters){
			for(size_t lvl = 0; lvl < masters->num_levels(); ++lvl){
				for(IntfIter intfIter = masters->begin(lvl);
					intfIter != masters->end(lvl); ++intfIter)
				{
					Interface& intf = masters->interface(intfIter);
					int connProc = intf.get_target_proc();

					//UG_LOG("checking master interface: " << localProcRank << "->" << connProc);
					//UG_LOG(" for target proc " << targetProc << endl);

					for(IntfElemIter elemIter = intf.begin();
						elemIter != intf.end(); ++elemIter)
					{
					//	check for each interface node whether it is marked
						TGeomObj* elem = intf.get_element(elemIter);
						if(!mg.is_marked(elem))
							continue;

					//	access target and trasferInfo attachments
						const vector<int>& targets = aaTargetProcs[elem];
						const vector<RedistributionNodeTransferInfo>& transferInfos =
																aaTransferInfos[elem];

					//	we have to find out which process will hold this master
					//	when redistribution is done.
						int newMasterProc = -1;
					//	if the elem stays on local-proc, the master stays where it was.
						if(sel.is_selected(elem)){
						//	If we're not on the master proc himself, then we can ignore
						//	such an entry.
//							if(targetProc != localProcRank)
//								continue;
							newMasterProc = localProcRank;
						}
						else{
						//	we have to find the target process with the smallest rank
							newMasterProc = pcl::GetNumProcesses();
							for(size_t i = 0; i < targets.size(); ++i){
								if(targets[i] < newMasterProc){
									newMasterProc = targets[i];
								}
							}
						}

						assert(newMasterProc != -1 && "A new master process has to exist!");

					//	we have to add interface entries for all target-procs
					//	between all target processes of the node (those were
					//	ignored during CreateDistributionLayouts).
						for(size_t i = 0; i < targets.size(); ++i){
							int target = targets[i];
							if(target == targetProc){
								if(target == newMasterProc){
								//	build master interfaces to all other targets
									for(size_t j = 0; j < targets.size(); ++j){
										if(targets[j] == target)
											continue;

										DistInterface& interface = distLayout.interface(
																	targets[j], lvl);
										interface.push_back(DistributionInterfaceEntry(
																aaNodeInd[elem], INT_MASTER));
									}
									break;
								}
								else{
								//	create interface to master
									DistInterface& interface = distLayout.interface(
																	newMasterProc, lvl);
									interface.push_back(DistributionInterfaceEntry(
															aaNodeInd[elem], INT_SLAVE));
									break;
								}
							}
						}

					//	if newMasterProc and targetProc match, we'll build a new master entry.
					//	however, we have to check whether the associated slave is sent
					//	to a new proc, too.
						if(newMasterProc == targetProc){
							for(size_t i = 0; i < transferInfos.size(); ++i){
								const RedistributionNodeTransferInfo& info = transferInfos[i];
								if(info.srcProc == connProc){
								//	the node is moved. only create an interface
								//	to the target proc.
									DistInterface& interface = distLayout.interface(
																info.targetProc, lvl);
									interface.push_back(DistributionInterfaceEntry(
															aaNodeInd[elem], INT_MASTER));
								}
							}
						}
					}
				}
			}
		}

	//	iterate over all slave nodes of this process. If the node is marked we'll
	//	have to further process it.
		if(slaves){
			for(size_t lvl = 0; lvl < slaves->num_levels(); ++lvl){
				for(IntfIter intfIter = slaves->begin(lvl);
					intfIter != slaves->end(lvl); ++intfIter)
				{
					Interface& intf = slaves->interface(intfIter);
					int connProc = intf.get_target_proc();

					//UG_LOG("checking slave interface: " << localProcRank << "->" << connProc);
					//UG_LOG(" for target proc " << targetProc << endl);

					for(IntfElemIter elemIter = intf.begin();
						elemIter != intf.end(); ++elemIter)
					{
					//	check for each interface node whether it is marked
						TGeomObj* elem = intf.get_element(elemIter);
						if(!mg.is_marked(elem))
							continue;

					//	access target and trasferInfo attachments
						const vector<int>& targets = aaTargetProcs[elem];
						const vector<RedistributionNodeTransferInfo>& transferInfos =
																aaTransferInfos[elem];

					//	we have to find out which process will hold the
					//	associated master when redistribution is done.
						int newMasterProc = connProc;
						for(size_t i = 0; i < transferInfos.size(); ++i){
							const RedistributionNodeTransferInfo& info = transferInfos[i];
							if(info.srcProc == connProc){
							//	we found the old master. if it moves, we have to
							//	find the new master location.
								if(info.bMove){
								//	iterate through the transfer infos again and
								//	find the smallest target.
									newMasterProc = pcl::GetNumProcesses();
									for(size_t j = 0; j < transferInfos.size(); ++j){
										if(transferInfos[j].targetProc < newMasterProc)
											newMasterProc = transferInfos[j].targetProc;
									}
								}
								break;
							}
						}

						assert(newMasterProc != -1 && "A new master process has to exist!");

					//	we have to add interface entries for all target-procs
					//	between all target processes of the node (those were
					//	ignored during CreateDistributionLayouts).
						for(size_t i = 0; i < targets.size(); ++i){
							int target = targets[i];
							if(target == targetProc){
								if(target == newMasterProc){
								//	build master interfaces to all other targets
									for(size_t j = 0; j < targets.size(); ++j){
										if(targets[j] == target)
											continue;

										DistInterface& interface = distLayout.interface(
																	targets[j], lvl);
										interface.push_back(DistributionInterfaceEntry(
																aaNodeInd[elem], INT_MASTER));
									}
									break;
								}
								else{
								//	create interface to master
									DistInterface& interface = distLayout.interface(
																	newMasterProc, lvl);
									interface.push_back(DistributionInterfaceEntry(
															aaNodeInd[elem], INT_SLAVE));
									break;
								}
							}
						}

					//	if newMasterProc and targetProc match, we'll build a new master entry.
					//	however, we have to check whether the associated slave is sent
					//	to a new proc, too.
						if(newMasterProc == targetProc){
							for(size_t i = 0; i < transferInfos.size(); ++i){
								const RedistributionNodeTransferInfo& info = transferInfos[i];
								if(info.srcProc == connProc){
								//	the node is moved. only create an interface
								//	to the target proc.
									DistInterface& interface = distLayout.interface(
																info.targetProc, lvl);
									interface.push_back(DistributionInterfaceEntry(
															aaNodeInd[elem], INT_MASTER));
								}
							}
						}
					}
				}
			}
		}
	}

	mg.end_marking();
	mg.detach_from<TGeomObj>(aNodeInd);
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
//todo:	One has to consider existing interfaces during creation of
//		the distribution layouts. Otherwise existing slave nodes
//		could accidently be marked as masters.
	CreateDistributionLayouts(vertexLayoutsOut, edgeLayoutsOut,
							  faceLayoutsOut, volumeLayoutsOut,
							  mg, sh, distributeGenealogy, &msel,
							  &distGridMgr);

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

}// end of namespace
