//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m11 d17

#include <utility>
#include <vector>
#include "lib_grid/lg_base.h"
#include "pcl/pcl.h"
#include "distribution_util.h"
#include "common/util/binary_stream.h"
#include "common/serialization.h"
#include "common/assert.h"
#include "lib_grid/algorithms/selection_util.h"
#include "lib_grid/algorithms/serialization.h"

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

///	The InfoVec is used to store pairs of (layoutIndex, indexInLayout).
/**	The InfoVec will be associated with each element and stores where
 * in which distribution layout the element can be found.
 */
typedef vector<pair<int, int> > InterfaceEntries;
struct DistInfo{
	DistInfo() : vrtMasterInterface(-1)	{}
	InterfaceEntries	interfaceEntries;
	int					vrtMasterInterface;
};
typedef Attachment<DistInfo>	ADistInfo;

///	Automatically attaches ADistInfo to all elements of a grid.
/**	On destruction, the attached dist info attachments are removed again.
 * You may access the dist-info in an element through the get method.
 * The get method returns a reference to the attached DistInfo object.
 *
 * Make sure that the given grid is valid while the DistInfoSupplier exists.
 */
class DistInfoSupplier{
	public:
		DistInfoSupplier(Grid& grid) : m_grid(grid)
		{
			m_grid.attach_to_all(m_aDistInfo);
			m_aaDistInfoVRT.access(grid, m_aDistInfo);
			m_aaDistInfoEDGE.access(grid, m_aDistInfo);
			m_aaDistInfoFACE.access(grid, m_aDistInfo);
			m_aaDistInfoVOL.access(grid, m_aDistInfo);
		}

		~DistInfoSupplier()
		{
			m_grid.detach_from_all(m_aDistInfo);
		}

		DistInfo& get(VertexBase* vrt)	{return m_aaDistInfoVRT[vrt];}
		DistInfo& get(EdgeBase* edge)		{return m_aaDistInfoEDGE[edge];}
		DistInfo& get(Face* face)			{return m_aaDistInfoFACE[face];}
		DistInfo& get(Volume* vol)		{return m_aaDistInfoVOL[vol];}
		DistInfo& get(GeometricObject* obj)
		{
			int objType = obj->base_object_id();
			switch(objType){
				case VERTEX:	return get(static_cast<VertexBase*>(obj));
				case EDGE:		return get(static_cast<EdgeBase*>(obj));
				case FACE:		return get(static_cast<Face*>(obj));
				case VOLUME:	return get(static_cast<Volume*>(obj));
				default:	UG_THROW("Unknown geometric object base type."); break;
			}
		}

	private:
		Grid& 		m_grid;
		ADistInfo	m_aDistInfo;
		Grid::AttachmentAccessor<VertexBase, ADistInfo> m_aaDistInfoVRT;
		Grid::AttachmentAccessor<EdgeBase, ADistInfo> m_aaDistInfoEDGE;
		Grid::AttachmentAccessor<Face, ADistInfo> m_aaDistInfoFACE;
		Grid::AttachmentAccessor<Volume, ADistInfo> m_aaDistInfoVOL;
};


///	manages the given layouts.
/**	Note that internally only references to the layouts specified in the constructor
 * are stored. Those vectors thus have to stay valid while the DistLayoutSupplier
 * is in use.
 */
template <template <class T> class TDistNodeLayout>
class DistLayoutSupplier{
	public:
		template <class TElem>
		struct traits{
			typedef TDistNodeLayout<TElem*>	layout;
			typedef typename layout::NodeVec			node_vec;
			typedef typename layout::Interface			interface;
			typedef typename layout::InterfaceEntry		interface_entry;
		};

		DistLayoutSupplier(
							vector<TDistNodeLayout<VertexBase*> >& vertexLayouts,
							vector<TDistNodeLayout<EdgeBase*> >& edgeLayouts,
							vector<TDistNodeLayout<Face*> >& faceLayouts,
							vector<TDistNodeLayout<Volume*> >& volumeLayouts) :
			m_vertexLayouts(vertexLayouts), m_edgeLayouts(edgeLayouts),
			m_faceLayouts(faceLayouts), m_volumeLayouts(volumeLayouts)
		{}

		void clear_and_resize(size_t newSize)
		{
			m_vertexLayouts.clear();
			m_vertexLayouts.resize(newSize);
			m_edgeLayouts.clear();
			m_edgeLayouts.resize(newSize);
			m_faceLayouts.clear();
			m_faceLayouts.resize(newSize);
			m_volumeLayouts.clear();
			m_volumeLayouts.resize(newSize);
		}

		template <class TElem>
		typename traits<TElem>::layout&
		layout(size_t layoutInd)
		{
			TElem dummy;
			return layout(layoutInd, dummy);
		}

		size_t size() const	{return m_vertexLayouts.size();}

	private:
		TDistNodeLayout<VertexBase*>& layout(size_t layoutInd, const VertexBase&)
		{return m_vertexLayouts[layoutInd];}

		TDistNodeLayout<EdgeBase*>& layout(size_t layoutInd, const EdgeBase&)
		{return m_edgeLayouts[layoutInd];}

		TDistNodeLayout<Face*>& layout(size_t layoutInd, const Face&)
		{return m_faceLayouts[layoutInd];}

		TDistNodeLayout<Volume*>& layout(size_t layoutInd, const Volume&)
		{return m_volumeLayouts[layoutInd];}

	private:
		vector<TDistNodeLayout<VertexBase*> >& m_vertexLayouts;
		vector<TDistNodeLayout<EdgeBase*> >& m_edgeLayouts;
		vector<TDistNodeLayout<Face*> >& m_faceLayouts;
		vector<TDistNodeLayout<Volume*> >& m_volumeLayouts;
};


////////////////////////////////////////////////////////////////////////
//	AddNodesToLayout
///	adds nodes to a layout and to interfaces if required.
/**
 * Please note that this method not only alters the layout referenced
 * by layoutIndex, but all layouts that share a node with this
 * layout. For each node that is referenced by multiple layouts,
 * corresponding interface entries are automatically generated.
 *
 * aDistInfo has to be attachmed at all elements.
 */
//	Implementation note: Instead of aaInfoVec we previously used
//	aaFirstProc and aaFirstProcLocalIndex attachments. This resulted in less
//	memory usage and still worked perfectly fine. However, during the creation
//	of vertical interfaces the InfoVecs are required. Since vertical
//	interfaces occur in all ug-multi-grid applications, we consider this
//	to be the standard case and try to optimize the code for this.
template <class TNodeLayout, class TIterator>
static
void AddNodesToLayout(std::vector<TNodeLayout>& layouts,
						int layoutIndex,
						TIterator nodesBegin, TIterator nodesEnd,
						Grid& grid,
						DistInfoSupplier& distInfos,
						std::vector<int>* processMap = NULL,
						int level = 0,
						int interfacesOnLevelOnly = -1,
						DistributedGridManager* pDistGridMgr = NULL)
{
	typedef typename TNodeLayout::Interface			Interface;
	typedef typename TNodeLayout::InterfaceEntry	InterfaceEntry;
	typedef typename TNodeLayout::NodeType			Node;
	typedef typename PtrTypeToGeomObjBaseType<Node>::base_type	Elem;

	TNodeLayout& layout = layouts[layoutIndex];

	for(TIterator iter = nodesBegin; iter != nodesEnd; ++iter)
	{
		typename TNodeLayout::NodeType node = *iter;
		InterfaceEntries& intfcEntries = distInfos.get(node).interfaceEntries;

		if(intfcEntries.empty())
		{
		//	the node has been encountered for the first time.
		//	add it to the layout and the info-vec
			intfcEntries.push_back(make_pair(layoutIndex,
										(int)layout.node_vec().size()));
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
			int masterLayoutIndex = intfcEntries.front().first;
			TNodeLayout& masterLayout = layouts[masterLayoutIndex];

			int localMasterID = intfcEntries.front().second;
			int localID = (int)layout.node_vec().size();

			intfcEntries.push_back(make_pair(layoutIndex, localID));
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
///	Selects all elements whose parents are sent to other processes only.
/**	This method uses Grid::mark.
 *
 * The method is useful to find all elements in a layout to which vertical interfaces shall
 * be build. One should use SelectAssociatedGeometricObjects to select all associated
 * lower dimensional elements. Only then all elements for vertical interfaces
 * are selected.
 */
template <class TDistLayout>
void SelectNewGhosts(TDistLayout& distLayout, MultiGrid& mg,
					 MGSelector& msel, DistInfoSupplier& distInfos)
{
	typedef typename TDistLayout::NodeVec			NodeVec;
	typedef typename TDistLayout::NodeType			Node;
	typedef typename TDistLayout::Interface			Interface;
	typedef typename TDistLayout::InterfaceEntry	InterfaceEntry;
	typedef typename PtrTypeToGeomObjBaseType<Node>::base_type	Elem;

	msel.clear<Elem>();

//	mark all nodes in curLayout
	mg.begin_marking();
	MarkNodesInLayout(mg, distLayout);

//	select all nodes, which have no parent in the same layout
	NodeVec& curNodes = distLayout.node_vec();
	for(size_t i_node = 0; i_node < curNodes.size(); ++i_node)
	{
		Node node = curNodes[i_node];
		if(distInfos.get(node).vrtMasterInterface != -1)
			continue;

		GeometricObject* parent = mg.get_parent(node);

	//	only consider parents of the same base-type.
	//	Other parents can be ignored, since we will adjust the selection
	//	after having executed this method for all object-types.
		if(parent && (parent->base_object_id() != Elem::BASE_OBJECT_ID))
			continue;

		if(!parent || !mg.is_marked(parent))
			msel.select(node);
	}

	mg.end_marking();
}

///	adds the given element to vertical layouts based on the given distInfos.
/**	Calls MakeVerticalMaster recursively for sides of elem.*/
template <template <class T> class TDistNodeLayout, class TElem>
void MakeVerticalMaster(TElem* elem, MultiGrid& mg,
						int masterLayoutInd,
						DistInfoSupplier& distInfos,
						DistLayoutSupplier<TDistNodeLayout>& distLayouts,
						std::vector<int>* processMap = NULL)
{
	typedef typename DistLayoutSupplier<TDistNodeLayout>::template traits<TElem>::node_vec NodeVec;
	typedef typename DistLayoutSupplier<TDistNodeLayout>::template traits<TElem>::interface Interface;
	typedef typename DistLayoutSupplier<TDistNodeLayout>::template traits<TElem>::interface_entry InterfaceEntry;

//	get the interface entries of the current element
	DistInfo& elemInfo = distInfos.get(elem);
	InterfaceEntries& intfcEntries = elemInfo.interfaceEntries;
	NodeVec& masterNodes = distLayouts.template layout<TElem>(masterLayoutInd).node_vec();

//	get the masters process and entry index
	int masterEntryInd = -1;
	for(size_t i = 0; i < intfcEntries.size(); ++i){
		if(intfcEntries[i].first == masterLayoutInd){
			masterEntryInd = intfcEntries[i].second;
			break;
		}
	}

	bool copyElemToMasterLayout = false;

	if(masterEntryInd == -1){
		copyElemToMasterLayout = true;
		masterEntryInd = (int)masterNodes.size();
	}

	int masterProcInd = masterLayoutInd;
	if(processMap)
		masterProcInd = (*processMap)[masterLayoutInd];

	size_t lvl = mg.get_level(elem);
	bool interfaceCreated = false;
	for(size_t i = 0; i < intfcEntries.size(); ++i){
		int slaveLayoutInd = intfcEntries[i].first;

		if(slaveLayoutInd == masterLayoutInd)
			continue;

		int slaveEntryInd = intfcEntries[i].second;

		int slaveProcInd = slaveLayoutInd;
		if(processMap)
			slaveProcInd = (*processMap)[slaveLayoutInd];

		Interface& masterInterface = distLayouts.template layout<TElem>(masterLayoutInd).
										interface(slaveProcInd, lvl);

		Interface& slaveInterface = distLayouts.template layout<TElem>(slaveLayoutInd).
										interface(masterProcInd, lvl);

		masterInterface.push_back(InterfaceEntry(masterEntryInd, INT_V_MASTER));
		slaveInterface.push_back(InterfaceEntry(slaveEntryInd, INT_V_SLAVE));
		interfaceCreated = true;
	}

	if(interfaceCreated){
		if(copyElemToMasterLayout){
		//	we have to create a new entry in the master-layout
		//	we want the new entry to be the first in the vector to
		//	optimize search performance.
			if(!intfcEntries.empty()){
				intfcEntries.push_back(intfcEntries.front());
				intfcEntries.front() = make_pair((int)masterLayoutInd,
												  (int)masterNodes.size());
			}
			else{
				intfcEntries.push_back(make_pair((int)masterLayoutInd,
												(int)masterNodes.size()));
			}

			masterNodes.push_back(elem);
		}
	//	now make sure that master interfaces for children are built accordingly
		elemInfo.vrtMasterInterface = masterLayoutInd;

		if(TElem::HAS_SIDES){
			typename Grid::traits<typename TElem::side>::secure_container sides;
			mg.associated_elements(sides, elem);
			for(size_t i = 0; i < sides.size(); ++i){
			//	Attention: Not much is currently done to really assert this. Indeed it probably
			//	can't even be asserted...
			//todo:	THINK ABOUT IT...
				DistInfo& sideInfo = distInfos.get(sides[i]);
				UG_ASSERT((sideInfo.vrtMasterInterface == -1)
						|| (sideInfo.vrtMasterInterface == masterLayoutInd),
						"An element may only be made 'vertical master' once.");
				if(sideInfo.vrtMasterInterface == -1){
					MakeVerticalMaster(sides[i], mg, masterLayoutInd, distInfos,
									   distLayouts, processMap);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
/**	Creates vertical interfaces for all selected elements.
 * Make sure to call this method first for volumes, then for faces,
 * edges and vertices.
 *
 * Elements which already are in vertical interfaces won't be added to new ones.
 *
 * Make sure that aDistInfo is attached to all elements of the grid and that it
 * stores (layout / local-index) pairs for each layout in which an element is
 * located.
 *
 * pMSel can optionally be passed to the method to avoid temporary allocation of
 * a new MGSelector.
 */
template <class TElem, template <class T> class TDistNodeLayout>
bool CreateVerticalInterfaces(DistLayoutSupplier<TDistNodeLayout>& distLayouts,
								MultiGrid& mg,
								SubsetHandler& shPart,
								DistInfoSupplier& distInfos,
								std::vector<int>* processMap = NULL,
								MGSelector* pMSel = NULL)
{
//	all elements which do not have a parent on their associated proc
//	have to be added to a vertical interface. An associated
//	vertical master has to be added in the associated interface somewhere.
	typedef TDistNodeLayout<TElem*>					TDistLayout;
	typedef typename TDistLayout::NodeVec			NodeVec;
	typedef typename TDistLayout::NodeType			Node;
	typedef typename TDistLayout::Interface			Interface;
	typedef typename TDistLayout::InterfaceEntry	InterfaceEntry;
	typedef typename PtrTypeToGeomObjBaseType<Node>::base_type	Elem;
	typedef typename Grid::traits<Elem>::iterator	ElemIter;

//	we need a multigrid selector to select ghosts of each layout
	MGSelector localMSel;
	if(!pMSel){
		localMSel.assign_grid(mg);
		pMSel = &localMSel;
	}
	MGSelector& msel = *pMSel;

	DistributedGridManager* pDistGridMgr = mg.distributed_grid_manager();

//	the base-layout index is used to create ghosts for elements which do not
//	have a parent at all.
//	Currently we use the first layout to which elements of highest dimension are
//	sent. This could surely be improved. If the local process keeps elements of
//	highest dimension, it might be better to simply choose the local process, since
//	this would reduce communication during redistribution.
	int baseLayoutInd = -1;

	if(mg.num<Volume>() > 0){
		for(int i = 0; i < shPart.num_subsets(); ++i){
			if(shPart.num<Volume>() > 0){
				baseLayoutInd = i;
				break;
			}
		}
	}
	if((baseLayoutInd == -1) && (mg.num<Face>() > 0)){
		for(int i = 0; i < shPart.num_subsets(); ++i){
			if(shPart.num<Face>() > 0){
				baseLayoutInd = i;
				break;
			}
		}
	}
	if((baseLayoutInd == -1) && (mg.num<EdgeBase>() > 0)){
		for(int i = 0; i < shPart.num_subsets(); ++i){
			if(shPart.num<EdgeBase>() > 0){
				baseLayoutInd = i;
				break;
			}
		}
	}
	if((baseLayoutInd == -1) && (mg.num<VertexBase>() > 0)){
		for(int i = 0; i < shPart.num_subsets(); ++i){
			if(shPart.num<VertexBase>() > 0){
				baseLayoutInd = i;
				break;
			}
		}
	}

	if(baseLayoutInd == -1)
		return false;

//	Now create the interfaces. Iterate over all distLayouts
	for(size_t i_layout = 0; i_layout < distLayouts.size(); ++i_layout)
	{
	//	select all new ghosts of the current layout
		SelectNewGhosts(distLayouts.template layout<Elem>(i_layout), mg, msel, distInfos);

	//	find for each selected element the master proc and create interfaces to
	//	all other procs (given in aaInfoVec).
		for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
			for(ElemIter iter = msel.begin<Elem>(lvl);
				iter != msel.end<Elem>(lvl); ++iter)
			{
				Elem* node = *iter;

			//	if pDistGridMgr is supplied, we first check however whether
			//	the element is already contained in a vertical interface.
			//	If so, we wont add it to another one.
				if(!pDistGridMgr || pDistGridMgr->is_in_vertical_interface(node))
					continue;

			//	the element has to be put into a vertical interface. If it
			//	hasn't got a parent, we'll create an interface to baseLayout.
			//	If it has a parent, we'll have to create an interface to the
			//	first copy, which has a parent on the same proc.
				int masterLayoutInd = -1;
				int masterLocalInd = -1;

			//	info vec of node
				InterfaceEntries& nivec = distInfos.get(node).interfaceEntries;
			//	pointer to info-vec of parent
				InterfaceEntries* ppivec = NULL;

				GeometricObject* parent = mg.get_parent(node);

			//	access the parents info vec.
				if(parent)
					ppivec = &distInfos.get(parent).interfaceEntries;

			//	now find the master-layout-index for the element. If required
			//	create copies of the elements on processes where a parent lies.
				if(parent && !ppivec->empty()){
				//	we have to find a copy of elem whose parent is in the same
				//	layout as the copy.
				//	If none can be found, we'll create a copy of elem on the
				//	first proc on which parent lies and make it a vertical master.
					InterfaceEntries& pivec = *ppivec;
					masterLayoutInd = pivec[0].first;

					for(size_t i_ni = 0; i_ni < nivec.size(); ++i_ni){
						int tLayoutInd = nivec[i_ni].first;
						for(size_t i_pi = 0; i_pi < pivec.size(); ++i_pi){
							if(tLayoutInd == pivec[i_pi].first){
								masterLayoutInd = tLayoutInd;
								masterLocalInd = nivec[i_ni].second;
								break;
							}
						}

						if(masterLocalInd != -1)
							break;
					}
				}
				else
					masterLayoutInd = baseLayoutInd;

				if(masterLayoutInd == (int)i_layout)
					continue;

				MakeVerticalMaster<TDistNodeLayout>(node, mg, masterLayoutInd, distInfos,
								   	   	   	   	    distLayouts, processMap);
			}
		}
	}

	mg.end_marking();
	return true;
}

////////////////////////////////////////////////////////////////////////
template <template <class T> class TDistNodeLayout>
void CreateDistributionLayouts(
						std::vector<TDistNodeLayout<VertexBase*> >& vertexLayoutsOut,
						std::vector<TDistNodeLayout<EdgeBase*> >& edgeLayoutsOut,
						std::vector<TDistNodeLayout<Face*> >& faceLayoutsOut,
						std::vector<TDistNodeLayout<Volume*> >& volumeLayoutsOut,
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
		UG_ASSERT((int)processMap->size() >= sh.num_subsets(),
				  "ProcessMap has to have enough entries for the number of partitions");
	}

//	resize and clear the layouts
	DistLayoutSupplier<TDistNodeLayout> distLayouts(vertexLayoutsOut, edgeLayoutsOut,
								   	   	   	   	    faceLayoutsOut, volumeLayoutsOut);
	distLayouts.clear_and_resize(sh.num_subsets());

	DistInfoSupplier distInfos(mg);

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
							mg, distInfos, processMap, -1);
		AddNodesToLayout(edgeLayoutsOut, i,
							sh.begin<EdgeBase>(i), sh.end<EdgeBase>(i),
							mg, distInfos, processMap, -1);
		AddNodesToLayout(faceLayoutsOut, i,
							sh.begin<Face>(i), sh.end<Face>(i),
							mg, distInfos, processMap, -1);
		AddNodesToLayout(volumeLayoutsOut, i,
							sh.begin<Volume>(i), sh.end<Volume>(i),
							mg, distInfos, processMap, -1);
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

	//todo:	replace the loop with a more efficient structure
		bool foundSomething;
		do{
			if(distributeGenealogy)
				SelectAssociatedGenealogy(msel, true);
			else
				SelectAssociatedGeometricObjects(msel);

			foundSomething = false;
		//	we have to make sure that constraining objects are sent to all
		//	processes to which their constrained objects are sent.
			for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
				for(ConstrainedVertexIterator iter = msel.begin<ConstrainedVertex>(lvl);
					iter != msel.end<ConstrainedVertex>(lvl); ++iter)
				{
					GeometricObject* cobj = (*iter)->get_constraining_object();
					if(cobj && !msel.is_selected(cobj)){
						foundSomething = true;
						msel.select(cobj);
					}
				}

				for(ConstrainedEdgeIterator iter = msel.begin<ConstrainedEdge>(lvl);
					iter != msel.end<ConstrainedEdge>(lvl); ++iter)
				{
					GeometricObject* cobj = (*iter)->get_constraining_object();
					if(cobj && !msel.is_selected(cobj)){
						foundSomething = true;
						msel.select(cobj);
					}
				}

				for(ConstrainedTriangleIterator iter = msel.begin<ConstrainedTriangle>(lvl);
					iter != msel.end<ConstrainedTriangle>(lvl); ++iter)
				{
					GeometricObject* cobj = (*iter)->get_constraining_object();
					if(cobj && !msel.is_selected(cobj)){
						foundSomething = true;
						msel.select(cobj);
					}
				}

				for(ConstrainedQuadrilateralIterator iter = msel.begin<ConstrainedQuadrilateral>(lvl);
					iter != msel.end<ConstrainedQuadrilateral>(lvl); ++iter)
				{
					GeometricObject* cobj = (*iter)->get_constraining_object();
					if(cobj && !msel.is_selected(cobj)){
						foundSomething = true;
						msel.select(cobj);
					}
				}
			}

		//	we also have to make sure that children of constraining elements are
		//	also sent to all the processes to which those constraining elements go.
			for(size_t lvl = 0; lvl < msel.num_levels(); ++lvl){
				for(ConstrainingEdgeIterator iter = msel.begin<ConstrainingEdge>(lvl);
					iter != msel.end<ConstrainingEdge>(lvl); ++iter)
				{
					ConstrainingEdge* ce = *iter;
					for(size_t i = 0; i < ce->num_constrained_vertices(); ++i){
						VertexBase* cde = ce->constrained_vertex(i);
						if(!msel.is_selected(cde)){
							foundSomething = true;
							msel.select(cde);
						}
					}

					for(size_t i = 0; i < ce->num_constrained_edges(); ++i){
						EdgeBase* cde = ce->constrained_edge(i);
						if(!msel.is_selected(cde)){
							foundSomething = true;
							msel.select(cde);
						}
					}
				}

				for(ConstrainingTriangleIterator iter = msel.begin<ConstrainingTriangle>(lvl);
					iter != msel.end<ConstrainingTriangle>(lvl); ++iter)
				{
					ConstrainingTriangle* ce = *iter;

					for(size_t i = 0; i < ce->num_constrained_vertices(); ++i){
						VertexBase* cde = ce->constrained_vertex(i);
						if(!msel.is_selected(cde)){
							foundSomething = true;
							msel.select(cde);
						}
					}

					for(size_t i = 0; i < ce->num_constrained_edges(); ++i){
						EdgeBase* cde = ce->constrained_edge(i);
						if(!msel.is_selected(cde)){
							foundSomething = true;
							msel.select(cde);
						}
					}

					for(size_t i = 0; i < ce->num_constrained_faces(); ++i){
						Face* cde = ce->constrained_face(i);
						if(!msel.is_selected(cde)){
							foundSomething = true;
							msel.select(cde);
						}
					}
				}

				for(ConstrainingQuadrilateralIterator iter = msel.begin<ConstrainingQuadrilateral>(lvl);
					iter != msel.end<ConstrainingQuadrilateral>(lvl); ++iter)
				{
					ConstrainingQuadrilateral* ce = *iter;

					for(size_t i = 0; i < ce->num_constrained_vertices(); ++i){
						VertexBase* cde = ce->constrained_vertex(i);
						if(!msel.is_selected(cde)){
							foundSomething = true;
							msel.select(cde);
						}
					}

					for(size_t i = 0; i < ce->num_constrained_edges(); ++i){
						EdgeBase* cde = ce->constrained_edge(i);
						if(!msel.is_selected(cde)){
							foundSomething = true;
							msel.select(cde);
						}
					}

					for(size_t i = 0; i < ce->num_constrained_faces(); ++i){
						Face* cde = ce->constrained_face(i);
						if(!msel.is_selected(cde)){
							foundSomething = true;
							msel.select(cde);
						}
					}
				}
			}
		}while(foundSomething == true);

		int interfacesOnLevelOnly = -1;
		/*
		if(distributeGenealogy)
			interfacesOnLevelOnly = mg.num_levels() - 1;
		*/
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
								mg, distInfos, processMap, level,
								interfacesOnLevelOnly, pDistGridMgr);
			AddNodesToLayout(edgeLayoutsOut, i,
								msel.begin<EdgeBase>(level), msel.end<EdgeBase>(level),
								mg, distInfos, processMap, level,
								interfacesOnLevelOnly, pDistGridMgr);
			AddNodesToLayout(faceLayoutsOut, i,
								msel.begin<Face>(level), msel.end<Face>(level),
								mg, distInfos, processMap, level,
								interfacesOnLevelOnly, pDistGridMgr);
			AddNodesToLayout(volumeLayoutsOut, i,
								msel.begin<Volume>(level), msel.end<Volume>(level),
								mg, distInfos, processMap, level,
								interfacesOnLevelOnly, pDistGridMgr);
		}
	}

//	horizontal layouts are complete by now. All nodes that go to process i
//	are contained in layout i at this point.
//	we can now use this information to add vertical interfaces
	if(createVerticalInterfaces){
	//	we'll call higher dim elems first, since lower dim-elems
	//	may be added to dist-groups during those calls
		CreateVerticalInterfaces<Volume>(distLayouts, mg, sh, distInfos, processMap, &msel);
		CreateVerticalInterfaces<Face>(distLayouts, mg, sh, distInfos, processMap, &msel);
		CreateVerticalInterfaces<EdgeBase>(distLayouts, mg, sh, distInfos, processMap, &msel);
		CreateVerticalInterfaces<VertexBase>(distLayouts, mg, sh, distInfos, processMap, &msel);
	}

//	The layouts are now complete.
}

//	explicit template instantiation
template void CreateDistributionLayouts<DistributionNodeLayout>(
										std::vector<DistributionNodeLayout<VertexBase*> >&,
										std::vector<DistributionNodeLayout<EdgeBase*> >&,
										std::vector<DistributionNodeLayout<Face*> >&,
										std::vector<DistributionNodeLayout<Volume*> >&,
										MultiGrid&, SubsetHandler&, bool, bool, MGSelector*,
										DistributedGridManager*, std::vector<int>*);

template void CreateDistributionLayouts<RedistributionNodeLayout>(
										std::vector<RedistributionNodeLayout<VertexBase*> >&,
										std::vector<RedistributionNodeLayout<EdgeBase*> >&,
										std::vector<RedistributionNodeLayout<Face*> >&,
										std::vector<RedistributionNodeLayout<Volume*> >&,
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
								BinaryBuffer& out, MultiGrid& mg,
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
						msel.get_geometric_objects(),
						aLocalIndVRT, aLocalIndEDGE,
						aLocalIndFACE, aLocalIndVOL, out);

//	write the layouts
	SerializeDistributionLayoutInterfaces(out, vrtLayout, mg, aLocalIndVRT);
	SerializeDistributionLayoutInterfaces(out, edgeLayout, mg, aLocalIndEDGE);
	SerializeDistributionLayoutInterfaces(out, faceLayout, mg, aLocalIndFACE);
	SerializeDistributionLayoutInterfaces(out, volLayout, mg, aLocalIndVOL);

//	done. Please note that no attachments have been serialized in this method.
}

////////////////////////////////////////////////////////////////////////
void SerializeGridAndRedistributionLayouts(
								BinaryBuffer& out, MultiGrid& mg,
								RedistributionVertexLayout& vrtLayout,
								RedistributionEdgeLayout& edgeLayout,
								RedistributionFaceLayout& faceLayout,
								RedistributionVolumeLayout& volLayout,
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

	SerializeGridAndDistributionLayouts(out, mg, vrtLayout, edgeLayout,
								faceLayout, volLayout, aLocalIndVRT,
								aLocalIndEDGE, aLocalIndFACE, aLocalIndVOL, pSel);

//	now serialize the global ids. Consider the order given by the index attachments
	SerializeGlobalIDs<VertexBase>(out, mg, vrtLayout, aLocalIndVRT);
	SerializeGlobalIDs<EdgeBase>(out, mg, edgeLayout, aLocalIndEDGE);
	SerializeGlobalIDs<Face>(out, mg, faceLayout, aLocalIndFACE);
	SerializeGlobalIDs<Volume>(out, mg, volLayout, aLocalIndVOL);
}


////////////////////////////////////////////////////////////////////////
void DeserializeGridAndDistributionLayouts(
								MultiGrid& mg, BinaryBuffer& in,
								DistributionVertexLayout& vrtLayout,
								DistributionEdgeLayout& edgeLayout,
								DistributionFaceLayout& faceLayout,
								DistributionVolumeLayout& volLayout)
{
	DeserializeMultiGridElements(mg, in, &vrtLayout.node_vec(),
								&edgeLayout.node_vec(), &faceLayout.node_vec(),
								&volLayout.node_vec());

	DeserializeDistributionLayoutInterfaces(vrtLayout, in);
	DeserializeDistributionLayoutInterfaces(edgeLayout, in);
	DeserializeDistributionLayoutInterfaces(faceLayout, in);
	DeserializeDistributionLayoutInterfaces(volLayout, in);
}

////////////////////////////////////////////////////////////////////////
void DeserializeGridAndRedistributionLayouts(
								MultiGrid& mg, BinaryBuffer& in,
								RedistributionVertexLayout& vrtLayout,
								RedistributionEdgeLayout& edgeLayout,
								RedistributionFaceLayout& faceLayout,
								RedistributionVolumeLayout& volLayout)
{
	DeserializeGridAndDistributionLayouts(mg, in, vrtLayout, edgeLayout,
										  faceLayout, volLayout);

	DeserializeGlobalIDs<VertexBase>(vrtLayout, in);
	DeserializeGlobalIDs<EdgeBase>(edgeLayout, in);
	DeserializeGlobalIDs<Face>(faceLayout, in);
	DeserializeGlobalIDs<Volume>(volLayout, in);
}


////////////////////////////////////////////////////////////////////////
//	DeserializeGridAndLayouts
void DeserializeGridAndDistributionLayouts(MultiGrid& mgOut,
											GridLayoutMap& gridLayoutOut,
											BinaryBuffer& in)
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
		if(interface[i].type() == type)
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
