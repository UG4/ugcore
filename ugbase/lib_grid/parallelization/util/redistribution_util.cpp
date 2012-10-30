// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 15.03.2011 (m,d,y)
 
#include "lib_grid/lg_base.h"
#include "pcl/pcl.h"
#include "distribution_util.h"
#include "common/util/binary_stream.h"
#include "common/serialization.h"
#include "lib_grid/algorithms/selection_util.h"

using namespace std;

namespace ug{

////////////////////////////////////////////////////////////////////////////////
///	A Debug method. Not used during normal program runs.
template <class TGeomObj, class TAATargetProcs, class TAAGlobalIDs>
static
void LogTargetProcs(MultiGrid& mg, TAATargetProcs& aaTargetProcs,
					  TAAGlobalIDs& aaIDs)
{
	typedef typename geometry_traits<TGeomObj>::iterator GeomObjIter;
	UG_LOG("  TARGET-PROCS:\n");
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		UG_LOG("    level " << lvl << ":\n");
		for(GeomObjIter iter = mg.begin<VertexBase>(lvl);
			iter != mg.end<VertexBase>(lvl); ++iter)
		{
		//	log the id
			TGeomObj* o = *iter;
			UG_LOG("      " << aaIDs[o].first << "_" << aaIDs[o].second << ":\t");

		//	log the target procs
			vector<int>& v = aaTargetProcs[o];
			for(size_t i = 0; i < v.size(); ++i){
				UG_LOG(" " << v[i]);
			}

			UG_LOG(endl);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
///	A Debug method. Not used during normal program runs.
template <class TGeomObj, class TAATransferInfos, class TAAGlobalIDs>
static
void LogTransferInfos(MultiGrid& mg, TAATransferInfos& aaTransferInfos,
					  TAAGlobalIDs& aaIDs)
{
	typedef typename geometry_traits<TGeomObj>::iterator GeomObjIter;
	UG_LOG("  TRANSFER-INFOS:\n");
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		UG_LOG("    level " << lvl << ":\n");
		for(GeomObjIter iter = mg.begin<VertexBase>(lvl);
			iter != mg.end<VertexBase>(lvl); ++iter)
		{
		//	log the id
			TGeomObj* o = *iter;
			UG_LOG("      " << aaIDs[o].first << "_" << aaIDs[o].second << ":\t");

		//	log the target procs
			vector<RedistributionNodeTransferInfo>& infos = aaTransferInfos[o];
			for(size_t i = 0; i < infos.size(); ++i){
				UG_LOG("src: " << infos[i].srcProc);
				UG_LOG(", target: " << infos[i].targetProc);
				UG_LOG(", bMove: " << infos[i].bMove << "  ||  ");
			}

			UG_LOG(endl);
		}
	}
}

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

////////////////////////////////////////////////////////////////////////////////
///	Collects all target processes for each entry in the grid.
/** Collects the target processes for each node, as specified through the
 * given redistribution layouts. The target processes are pushed to the
 * vector<int> attachment in each node.
 *
 * A selector may optionally be specified through psel. For the layout of the
 * current process, values are then only assigned to selected elements. This can
 * be used to exclude vertical masters.
 */
/* THIS IS THE ORIGINAL AND A LITTLE MISLEADING COMMENT...
 * For each element of the grid of type TGeomObj, the target processes can
 * be found in the attachment to which TAAIntVec points, after the algorithm
 * is done.
 *
 * Note that the current process itself is not regarded as a target proc.
 *
 * Make sure that TGeomObj is either VertexBase, EdgeBase, Face or Volume.
 * Make also sure that TAAIntVec is accessing a valid Attachment<vector<int> >
 * attachment at grid.
 *
 * if psel is specified, then unselected entries of the layout with
 * localLayoutIndex will only recieve receive the current-proc as their target.
 */
template <class TGeomObj, class TAAIntVec>
static void
CollectRedistributionTargetProcs(
						vector<RedistributionNodeLayout<TGeomObj> >& layouts,
						TAAIntVec& aaIntVec, int localLayoutIndex,
						std::vector<int>* processMap,
						ISelector* psel = NULL)
{
	typedef typename DistributionNodeLayout<TGeomObj>::NodeVec		NodeVec;
	typedef typename DistributionNodeLayout<TGeomObj>::Interface	Interface;
	typedef typename DistributionNodeLayout<TGeomObj>::InterfaceMap	InterfaceMap;

	for(size_t i_layout = 0; i_layout < layouts.size(); ++i_layout){
		//if((int)i_layout == localLayoutIndex)
		//	continue;

		int targetProc = i_layout;
		if(processMap)
			targetProc = (*processMap)[i_layout];

	//	get the nodes of the current layout
		NodeVec& nodes = layouts[i_layout].node_vec();

	//	if we're checking the local layout and if psel was specified, then
	//	we have to only add targetProcs to selected entries.
		if(localLayoutIndex == (int)i_layout && psel != NULL){
		//	now iterate over the nodes and add the target-proc
			for(size_t i = 0; i < nodes.size(); ++i){
				if(psel->is_selected(nodes[i]))
					aaIntVec[nodes[i]].push_back(targetProc);
			}
		}
		else{
		//	now iterate over the nodes and add the target-proc
			for(size_t i = 0; i < nodes.size(); ++i)
				aaIntVec[nodes[i]].push_back(targetProc);
		}
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
 *
 * \param pureVertical	if this parameter is set to true, data is only communicated
 * 						for elements which only lie in a vertical interface.
 * 						This is important to avoid double-entries.
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
		ComPol_SynchronizeNodeTransfer(DistributedGridManager& distGridMgr,
										ISelector& sel, TAAIntVec& aaIntVec,
										TAATransferInfoVec& aaTransInfoVec,
										bool pureVertical) :
			m_distGridMgr(distGridMgr),
			m_sel(sel),
			m_aaIntVec(aaIntVec),
			m_aaTransInfoVec(aaTransInfoVec),
			m_pureVertical(pureVertical)
		{}

		virtual int
		get_required_buffer_size(Interface& interface)		{return -1;}

	///	write target processes and move-flag
		virtual bool
		collect(ug::BinaryBuffer& buff, Interface& interface)
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
		extract(ug::BinaryBuffer& buff, Interface& interface)
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

				if(m_pureVertical){
				//	if pure vertical is enabled, we wont use data for
				//	elements which also lie in a horizontal interface
					if(m_distGridMgr.contains_status(elem, INT_H_MASTER)
					  || m_distGridMgr.contains_status(elem, INT_H_SLAVE))
						continue;
				}

				for(size_t i = 0; i < targetProcs.size(); ++i){
					m_aaTransInfoVec[elem].push_back(
						RedistributionNodeTransferInfo(srcProc, targetProcs[i],
														bMove != 0));
				}
			}

			return true;
		}

	protected:
		DistributedGridManager&	m_distGridMgr;
		ISelector&			m_sel;
		TAAIntVec&			m_aaIntVec;
		TAATransferInfoVec&	m_aaTransInfoVec;
		bool				m_pureVertical;
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
SynchronizeHNodeTransfer(DistributedGridManager& distGridMgr,
						pcl::InterfaceCommunicator<typename GridLayoutMap::Types<TGeomObj>::Layout> & comm,
						TAAIntVec& aaIntVec,
						TAATransferInfoVec& aaTransInfoVec,
						int localLayoutIndex,
						ISelector& sel)
{
	typedef typename GridLayoutMap::Types<TGeomObj>::Layout		Layout;

	GridLayoutMap& glm = distGridMgr.grid_layout_map();

//	this communication policy fills aaTransInfoVec
	ComPol_SynchronizeNodeTransfer<Layout, TAAIntVec, TAATransferInfoVec>
		compolH(distGridMgr, sel, aaIntVec, aaTransInfoVec, false);

//	exchange data, both from masters to slaves and vice versa
	comm.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolH);
	comm.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolH);
	comm.communicate();
}

template <class TGeomObj, class TAAIntVec, class TAATransferInfoVec>
static void
SynchronizeVNodeTransfer(DistributedGridManager& distGridMgr,
						pcl::InterfaceCommunicator<typename GridLayoutMap::Types<TGeomObj>::Layout> & comm,
						TAAIntVec& aaIntVec,
						TAATransferInfoVec& aaTransInfoVec,
						int localLayoutIndex,
						ISelector& sel)
{
	typedef typename GridLayoutMap::Types<TGeomObj>::Layout		Layout;

	GridLayoutMap& glm = distGridMgr.grid_layout_map();

//	this communication policy fills aaTransInfoVec
	ComPol_SynchronizeNodeTransfer<Layout, TAAIntVec, TAATransferInfoVec>
		compolV(distGridMgr, sel, aaIntVec, aaTransInfoVec, false);

//	exchange data, both from masters to slaves and vice versa
	comm.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolV);
	comm.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolV);
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
					byte masterKey,
					byte slaveKey,
					vector<int>* processMap = NULL)
{
	typedef RedistributionNodeLayout<TGeomObj*> TDistLayout;
	typedef typename TDistLayout::NodeType 	Node;
	typedef typename TDistLayout::Interface	DistInterface;
	typedef typename TDistLayout::ProcPair	ProcPair;

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

//TEMPORARY POSITION ATTACHMENT
//	Grid::VertexAttachmentAccessor<APosition2> aaPos(mg, aPosition2);

	Layout* masters = NULL;
	if(distGridMgr.grid_layout_map().template has_layout<TGeomObj>(masterKey))
		masters = &distGridMgr.grid_layout_map().template get_layout<TGeomObj>(masterKey);

	Layout* slaves = NULL;
	if(distGridMgr.grid_layout_map().template has_layout<TGeomObj>(slaveKey))
		slaves = &distGridMgr.grid_layout_map().template get_layout<TGeomObj>(slaveKey);

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
			targetProc = (*processMap)[i_layout];

		//UG_LOG("TARGET_PROC: " << targetProc << endl);

		TDistLayout& distLayout = distLayoutVec[i_layout];
		const typename TDistLayout::NodeVec& nodes = distLayout.node_vec();

	//	mark all nodes in the current layout and assign their redist-layout node
	//	index
		mg.clear_marks();
		for(size_t i = 0; i < nodes.size(); ++i){
			aaNodeInd[nodes[i]] = (int)i;
			mg.mark(nodes[i]);
		}

	//	iterate over all existing master nodes of this process. If the node is marked we'll
	//	have to further process it.
		if(masters){
			for(size_t lvl = 0; lvl < masters->num_levels(); ++lvl){
			//	Note that interfaces are sorted by associated proc-id.
				for(IntfIter intfIter = masters->begin(lvl);
					intfIter != masters->end(lvl); ++intfIter)
				{
					Interface& intf = masters->interface(intfIter);
					int connProc = intf.get_target_proc();
					//UG_LOG("CONN_PROC: " << connProc << endl);
					//UG_LOG("INTF_SIZE: " << intf.size() << endl);
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
															ProcPair(targets[j], connProc), lvl);
										interface.push_back(DistributionInterfaceEntry(
																aaNodeInd[elem], masterKey));
									}
									break;
								}
								else{
								//	create interface to master
									DistInterface& interface = distLayout.interface(
															ProcPair(newMasterProc, connProc), lvl);
									interface.push_back(DistributionInterfaceEntry(
															aaNodeInd[elem], slaveKey));
									break;
								}
							}
						}

					//	if newMasterProc and targetProc match, we'll build a new master entry.
					//	however, we have to check whether the associated slave is sent
					//	to a new proc, too.
						if(newMasterProc == targetProc){
						/*
							UG_LOG("TI:");
							if(elem->base_object_id() == VERTEX){
								UG_LOG(" [" << aaPos[elem] << "] -");
							}
						*/
							for(size_t i = 0; i < transferInfos.size(); ++i){
								const RedistributionNodeTransferInfo& info = transferInfos[i];
								if(info.srcProc == connProc){
									//UG_LOG(" (" << info.srcProc << ", " << info.targetProc << ")");
								//	the node is moved. only create an interface
								//	to the target proc.
								//	If targetProc and info.targetProc are the same, then
								//	we don't have to create the entry...
									if(targetProc == info.targetProc)
										continue;

									DistInterface& interface = distLayout.interface(
														ProcPair(info.targetProc, connProc), lvl);

									interface.push_back(DistributionInterfaceEntry(
															aaNodeInd[elem], masterKey));
								}
								else{
									//UG_LOG(" *(" << info.srcProc << ", " << info.targetProc << ")*");
								}
							}
							//UG_LOG(endl);
						}
					}
				}
			}
		}

	//	iterate over all slave nodes of this process. If the node is marked we'll
	//	have to further process it.
		if(slaves){
			for(size_t lvl = 0; lvl < slaves->num_levels(); ++lvl){
			//	Note that interfaces are sorted by associated proc-id.
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

					//	access target and transferInfo attachments
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
															ProcPair(targets[j], connProc), lvl);
										interface.push_back(DistributionInterfaceEntry(
																aaNodeInd[elem], masterKey));
									}
									break;
								}
								else{
								//	create interface to master
									DistInterface& interface = distLayout.interface(
															ProcPair(newMasterProc, connProc), lvl);
									interface.push_back(DistributionInterfaceEntry(
															aaNodeInd[elem], slaveKey));
									break;
								}
							}
						}

					//	if newMasterProc and targetProc match, we'll build a new master entry.
					//	however, we have to check whether the associated slave is sent
					//	to a new proc, too.
					/*
					 * As far as I can see, this piece of code is unnecessary.
					 * One should remove it as soon as the code proves to run stable.
						if(newMasterProc == targetProc){
							for(size_t i = 0; i < transferInfos.size(); ++i){
								const RedistributionNodeTransferInfo& info = transferInfos[i];
								if(info.srcProc == connProc){
								//	the node is moved. only create an interface
								//	to the target proc.
								//	If targetProc and info.targetProc are the same, then
								//	we don't have to create the entry...
									if(targetProc == info.targetProc)
										continue;

									DistInterface& interface = distLayout.interface(
																info.targetProc, lvl);
									interface.push_back(DistributionInterfaceEntry(
															aaNodeInd[elem], INT_H_MASTER));
								}
							}
						}
					*/
					}
				}
			}
		}
	}

	mg.end_marking();
	mg.detach_from<TGeomObj>(aNodeInd);
}


///	Selects elements in layout, which are not vertical-masters
template <class TGeomObj, class TRedistLayout>
void SelectAllButVerticalMasters(DistributedGridManager& distGridMgr,
								 ISelector& sel, TRedistLayout& layout)
{

	typename TRedistLayout::NodeVec& nodes = layout.node_vec();
	for(size_t i = 0; i < nodes.size(); ++i)
	{
		if(!(distGridMgr.contains_status(nodes[i], INT_V_MASTER)))
		{
			sel.select(nodes[i]);
		}
	}
}

///	Deselects all elements in layouts other than localLayoutInd.
template <class TGeomObj>
void DeselectElementsInOtherLayouts(ISelector& sel, int localLayoutInd,
				std::vector<RedistributionNodeLayout<TGeomObj*> >& distLayoutVec)
{
//	iterate through all layouts.
//	deselect all elements in layouts other than loaclLayoutInd.
	for(size_t i_layout = 0; i_layout < distLayoutVec.size(); ++i_layout)
	{
		if((int)i_layout == localLayoutInd)
			continue;

		typename RedistributionNodeLayout<TGeomObj*>::NodeVec&
			nodes = distLayoutVec[i_layout].node_vec();

		for(size_t i = 0; i < nodes.size(); ++i)
			sel.deselect(nodes[i]);
	}
}

///	Calls clear on all vectors in the given attachments.
template <class TGeomObj, class TAATargetProcs, class TAATransferInfos>
void ClearTargetAndTransferInfos(Grid& g, TAATargetProcs& aaTP,
								 TAATransferInfos& aaTF)
{
//	iterate over all geom objs and clear the attachments
	typedef typename geometry_traits<TGeomObj>::iterator TIter;
	for(TIter iter = g.begin<TGeomObj>(); iter != g.end<TGeomObj>(); ++iter)
	{
		aaTP[*iter].clear();
		aaTF[*iter].clear();
	}
}

////////////////////////////////////////////////////////////////////////
void CreateRedistributionLayouts(
					std::vector<RedistributionVertexLayout>& vertexLayoutsOut,
					std::vector<RedistributionEdgeLayout>& edgeLayoutsOut,
					std::vector<RedistributionFaceLayout>& faceLayoutsOut,
					std::vector<RedistributionVolumeLayout>& volumeLayoutsOut,
					DistributedGridManager& distGridMgr, SubsetHandler& sh,
					bool distributeGenealogy, bool createVerticalInterfaces,
					MGSelector* pSel, std::vector<int>* processMap,
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
//		... I think this has already been taken care of...
	CreateDistributionLayouts(vertexLayoutsOut, edgeLayoutsOut,
							  faceLayoutsOut, volumeLayoutsOut,
							  mg, sh, distributeGenealogy,
							  createVerticalInterfaces,
							  &msel, &distGridMgr, processMap);

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
		for(size_t i = 0; i < processMap->size(); ++i){
			if(processMap->at(i) == myRank){
				localLayoutInd = (int)i;
				break;
			}
		}
	}
	else{
		localLayoutInd = pcl::GetProcRank();
	}

	if(localLayoutInd >= sh.num_subsets())
		localLayoutInd = -1;

//	to each node we will now attach a vector<int>, into which we'll write all
//	processes to which the node will be sent.
//	we'll use the autoattach mechanism of the AttachmentAccessor.
	typedef Attachment<vector<int> > AIntVec;
	AIntVec aTargetProcs;

	Grid::AttachmentAccessor<VertexBase, AIntVec>
		aaTargetProcsVRT(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<EdgeBase, AIntVec>
		aaTargetProcsEDGE(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<Face, AIntVec>
		aaTargetProcsFACE(mg, aTargetProcs, true);
	Grid::AttachmentAccessor<Volume, AIntVec>
		aaTargetProcsVOL(mg, aTargetProcs, true);

//	select all elements, which are not contained in a vertical master interface.
//	Include all associated elements of lower dimension in that selection, too.
//	However - don't select elements from other layouts. This is important,
//	since the selection should afterwards only hold the elements it would
//	hold, if vertical interfaces were not present.
//	This is important, since otherwise we would accidentally build h-interfaces,
//	where only v-interfaces are required.
	if(localLayoutInd != -1){
		SelectAllButVerticalMasters<VertexBase>(distGridMgr, msel,
											vertexLayoutsOut[localLayoutInd]);
		SelectAllButVerticalMasters<EdgeBase>(distGridMgr, msel,
											edgeLayoutsOut[localLayoutInd]);
		SelectAllButVerticalMasters<Face>(distGridMgr, msel,
											faceLayoutsOut[localLayoutInd]);
		SelectAllButVerticalMasters<Volume>(distGridMgr, msel,
											volumeLayoutsOut[localLayoutInd]);

		DeselectElementsInOtherLayouts(msel, localLayoutInd, vertexLayoutsOut);
		DeselectElementsInOtherLayouts(msel, localLayoutInd, edgeLayoutsOut);
		DeselectElementsInOtherLayouts(msel, localLayoutInd, faceLayoutsOut);
		DeselectElementsInOtherLayouts(msel, localLayoutInd, volumeLayoutsOut);

		SelectAssociatedGeometricObjects(msel);
	}
/*
UG_LOG("DEBUG: Positions of selected vertices:\n");
if(mg.num_levels() > 0){
	Grid::AttachmentAccessor<VertexBase, APosition2> aaPos(mg, aPosition2);
	for(VertexBaseIterator iter = msel.begin<VertexBase>(mg.num_levels()-1);
		iter != msel.end<VertexBase>(mg.num_levels()-1); ++iter)
	{
		UG_LOG(aaPos[*iter] << endl);
	}
}
*/
	CollectRedistributionTargetProcs(vertexLayoutsOut, aaTargetProcsVRT,
									 localLayoutInd, processMap, &msel);
	CollectRedistributionTargetProcs(edgeLayoutsOut, aaTargetProcsEDGE,
									 localLayoutInd, processMap, &msel);
	CollectRedistributionTargetProcs(faceLayoutsOut, aaTargetProcsFACE,
									 localLayoutInd, processMap, &msel);
	CollectRedistributionTargetProcs(volumeLayoutsOut, aaTargetProcsVOL,
									 localLayoutInd, processMap, &msel);

//	we need a second attachment, which tells us for each node, to where
//	associated nodes on other processes go. We store this through a
//	vector<RedistributionNodeMoveInfo>, which tells us on which
//	process the copy/move operation was scheduled (srcProc) and to which
//	process the element will be transfered (targetProc) and which operation shall be
//	performed (bMove).
	typedef ARedistributionNodeTransferInfoVec ATransferInfoVec;
	ATransferInfoVec aTransferInfoVec;

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
	pcl::InterfaceCommunicator<VertexLayout> commVRT;
	pcl::InterfaceCommunicator<EdgeLayout> commEDGE;
	pcl::InterfaceCommunicator<FaceLayout> commFACE;
	pcl::InterfaceCommunicator<VolumeLayout> commVOL;

////////////////////////////////
//	BUILD HORIZONTAL INTERFACES
//	NOTE on the following comment: Shortening may be intended.
//	before synchronizing horizontal node movement, we have to make sure, that
//	all elements of the current layout are selected. Otherwise old h-interfaces
//	would be lost or shortened.
/*
	msel.clear();
	if(localLayoutInd != -1){
		SelectNodesInLayout(msel, vertexLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, edgeLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, faceLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, volumeLayoutsOut[localLayoutInd]);
	}
*/

	SynchronizeHNodeTransfer<VertexBase>(distGridMgr, commVRT, aaTargetProcsVRT,
							aaTransInfoVecVRT, localLayoutInd, msel);
	SynchronizeHNodeTransfer<EdgeBase>(distGridMgr, commEDGE, aaTargetProcsEDGE,
							aaTransInfoVecEDGE, localLayoutInd, msel);
	SynchronizeHNodeTransfer<Face>(distGridMgr, commFACE, aaTargetProcsFACE,
							aaTransInfoVecFACE, localLayoutInd, msel);
	SynchronizeHNodeTransfer<Volume>(distGridMgr, commVOL, aaTargetProcsVOL,
							aaTransInfoVecVOL, localLayoutInd, msel);

//	finalize the redistribution layouts by constructing all interfaces
	FinalizeRedistributionLayoutInterfaces(distGridMgr, vertexLayoutsOut,
										aaTargetProcsVRT, aaTransInfoVecVRT,
										msel, INT_H_MASTER, INT_H_SLAVE, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, edgeLayoutsOut,
										aaTargetProcsEDGE, aaTransInfoVecEDGE,
										msel, INT_H_MASTER, INT_H_SLAVE, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, faceLayoutsOut,
										aaTargetProcsFACE, aaTransInfoVecFACE,
										msel, INT_H_MASTER, INT_H_SLAVE, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, volumeLayoutsOut,
										aaTargetProcsVOL, aaTransInfoVecVOL,
										msel, INT_H_MASTER, INT_H_SLAVE, processMap);

////////////////////////////////
//	BUILD VERTICAL INTERFACES
//	targetProcs and transferInfos now have to be cleared. We will then
//	synchronize the targetProcs for all horizontal interfaces.
	msel.clear();
	if(localLayoutInd != -1){
		SelectNodesInLayout(msel, vertexLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, edgeLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, faceLayoutsOut[localLayoutInd]);
		SelectNodesInLayout(msel, volumeLayoutsOut[localLayoutInd]);
	}


	ClearTargetAndTransferInfos<VertexBase>(mg, aaTargetProcsVRT, aaTransInfoVecVRT);
	ClearTargetAndTransferInfos<EdgeBase>(mg, aaTargetProcsEDGE, aaTransInfoVecEDGE);
	ClearTargetAndTransferInfos<Face>(mg, aaTargetProcsFACE, aaTransInfoVecFACE);
	ClearTargetAndTransferInfos<Volume>(mg, aaTargetProcsVOL, aaTransInfoVecVOL);

	CollectRedistributionTargetProcs(vertexLayoutsOut, aaTargetProcsVRT,
									 localLayoutInd, processMap);
	CollectRedistributionTargetProcs(edgeLayoutsOut, aaTargetProcsEDGE,
									 localLayoutInd, processMap);
	CollectRedistributionTargetProcs(faceLayoutsOut, aaTargetProcsFACE,
									 localLayoutInd, processMap);
	CollectRedistributionTargetProcs(volumeLayoutsOut, aaTargetProcsVOL,
									 localLayoutInd, processMap);

	SynchronizeVNodeTransfer<VertexBase>(distGridMgr, commVRT, aaTargetProcsVRT,
							aaTransInfoVecVRT, localLayoutInd, msel);
	SynchronizeVNodeTransfer<EdgeBase>(distGridMgr, commEDGE, aaTargetProcsEDGE,
							aaTransInfoVecEDGE, localLayoutInd, msel);
	SynchronizeVNodeTransfer<Face>(distGridMgr, commFACE, aaTargetProcsFACE,
							aaTransInfoVecFACE, localLayoutInd, msel);
	SynchronizeVNodeTransfer<Volume>(distGridMgr, commVOL, aaTargetProcsVOL,
							aaTransInfoVecVOL, localLayoutInd, msel);

	FinalizeRedistributionLayoutInterfaces(distGridMgr, vertexLayoutsOut,
										aaTargetProcsVRT, aaTransInfoVecVRT,
										msel, INT_V_MASTER, INT_V_SLAVE, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, edgeLayoutsOut,
										aaTargetProcsEDGE, aaTransInfoVecEDGE,
										msel, INT_V_MASTER, INT_V_SLAVE, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, faceLayoutsOut,
										aaTargetProcsFACE, aaTransInfoVecFACE,
										msel, INT_V_MASTER, INT_V_SLAVE, processMap);
	FinalizeRedistributionLayoutInterfaces(distGridMgr, volumeLayoutsOut,
										aaTargetProcsVOL, aaTransInfoVecVOL,
										msel, INT_V_MASTER, INT_V_SLAVE, processMap);

// detach temporary attachments.
	mg.detach_from_all(aTargetProcs);
	mg.detach_from_all(aTransferInfoVec);
}

}// end of namespace
