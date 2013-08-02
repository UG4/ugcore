/*
 * domain_impl.h
 *
 *  Created on: 02.03.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOMAIN_IMPL__
#define __H__UG__LIB_DISC__DOMAIN_IMPL__

#include "domain.h"
#include "common/serialization.h"

#ifdef UG_PARALLEL
#include "pcl/pcl_process_communicator.h"
#endif


namespace ug{

////////////////////////////////////////////////////////////////////////////////
// IDomain
////////////////////////////////////////////////////////////////////////////////

template <typename TGrid, typename TSubsetHandler>
IDomain<TGrid,TSubsetHandler>::IDomain(bool isAdaptive)
	:
	m_spGrid(new TGrid(GRIDOPT_NONE)),	// Note: actual options are set by the derived class (dimension dependent).
	m_spSH(new TSubsetHandler(*m_spGrid)),
	m_isAdaptive(isAdaptive),
	m_adaptionIsActive(false)
{
	#ifdef UG_PARALLEL
	//	the grid has to be prepared for parallelism
		m_spGrid->set_parallel(true);
	#endif

//	register function for grid adaption
	m_spGridAdaptionCallbackID =
		message_hub()->register_class_callback(this,
		&ug::IDomain<ug::MultiGrid, ug::MultiGridSubsetHandler>::grid_adaption_callback);

	m_spGridCreationCallbackID =
		message_hub()->register_class_callback(this,
		&ug::IDomain<ug::MultiGrid, ug::MultiGridSubsetHandler>::grid_creation_callback);

	m_spGridDistributionCallbackID =
		message_hub()->register_class_callback(this,
		&ug::IDomain<ug::MultiGrid, ug::MultiGridSubsetHandler>::grid_distribution_callback);
}

///	Destructor
template <typename TGrid, typename TSubsetHandler>
IDomain<TGrid,TSubsetHandler>::~IDomain()
{
}


template <typename TGrid, typename TSubsetHandler>
void IDomain<TGrid,TSubsetHandler>::
update_subset_infos(int rootProc)
{
	TSubsetHandler& sh = *m_spSH;
	for(int i = 0; i < sh.num_subsets(); ++i){
		int dim = -1;
		if(sh.contains_volumes(i))
			dim = 3;
		else if(sh.contains_faces(i))
			dim = 2;
		else if(sh.contains_edges(i))
			dim = 1;
		else if(sh.contains_vertices(i))
			dim = 0;

		sh.subset_info(i).set_property("dim", dim);
	}

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator procCom;

//	prepare the subset-info package, send it to all procs and extract the info again.
	BinaryBuffer buf;
	if(pcl::GetProcRank() == rootProc){
		Serialize(buf, sh.num_subsets());
		for(int i = 0; i < sh.num_subsets(); ++i){
			Serialize(buf, sh.subset_info(i).name);
			Serialize(buf, sh.subset_info(i).get_property("dim").to_int());
		}
	}

	procCom.broadcast(buf, rootProc);

	if(pcl::GetProcRank() != rootProc){
		int numSubsets;
		Deserialize(buf, numSubsets);
		if(numSubsets > 0)
			sh.subset_required(numSubsets - 1);

		for(int i = 0; i < numSubsets; ++i){
			Deserialize(buf, sh.subset_info(i).name);
			int dim;
			Deserialize(buf, dim);
			sh.subset_info(i).set_property("dim", dim);
		}
	}
#endif

}


template <typename TGrid, typename TSubsetHandler>
inline
void IDomain<TGrid,TSubsetHandler>::
grid_adaption_callback(const GridMessage_Adaption& msg)
{
	if(msg.adaption_begins())
		m_adaptionIsActive = true;

	else if(m_adaptionIsActive){
		if(msg.adaptive()){
			if(msg.adaption_ends())
			{
				update_domain_info();
				m_adaptionIsActive = false;
			}
		}
	}

	else{
		UG_THROW("Before any grid-adaption may be performed, the domain "
				"has to be informed that grid-adaption shall begin. "
				"You may use IRefiner::grid_adaption_begins() or schedule "
				"an appropriate message to the associated grids message-hub.");
	}
}

template <typename TGrid, typename TSubsetHandler>
inline
void IDomain<TGrid,TSubsetHandler>::
grid_creation_callback(const GridMessage_Creation& msg)
{
	if(msg.msg() == GMCT_CREATION_STOPS){
		if(msg.proc_id() != -1)
			update_subset_infos(msg.proc_id());
		update_domain_info();
	}
}

template <typename TGrid, typename TSubsetHandler>
inline
void IDomain<TGrid,TSubsetHandler>::
grid_distribution_callback(const GridMessage_Distribution& msg)
{
//	this is already handled in grid_creation_callback
	/*if(msg.msg() == GMDT_DISTRIBUTION_STOPS){
		update_domain_info();
	}*/
}


template <typename TGrid, typename TSubsetHandler>
void IDomain<TGrid,TSubsetHandler>::
update_domain_info()
{
	TGrid& mg = *m_spGrid;
	TSubsetHandler& sh = *m_spSH;

	GeometricBaseObject	locElemType;
	if(mg.template num<Volume>() > 0)
		locElemType = VOLUME;
	else if(mg.template num<Face>() > 0)
		locElemType = FACE;
	else if(mg.template num<EdgeBase>() > 0)
		locElemType = EDGE;
	else
		locElemType = VERTEX;

	GeometricBaseObject	elemType;
	#ifdef UG_PARALLEL
		pcl::ProcessCommunicator commWorld;

	//	all processes should have the same number of grid levels
		int numLevLocal = m_spGrid->num_levels();
		int numLevGlobal;
		commWorld.allreduce(&numLevLocal, &numLevGlobal, 1, PCL_DT_INT, PCL_RO_MAX);

		if(numLevGlobal > 0){
			m_spGrid->level_required(numLevGlobal-1);
			m_spSH->level_required(numLevGlobal-1);
		}

	//	communicate the element type of highest dimension present in the global grid.
		elemType = (GeometricBaseObject) commWorld.allreduce((int)locElemType, PCL_RO_MAX);
	#else
		elemType = locElemType;
	#endif

//	the number of levels of the multi-grid is now equal on all processes
	std::vector<int>	locNumElemsOnLvl;
	locNumElemsOnLvl.reserve(mg.num_levels());

	switch(elemType){
		case VOLUME:
			for(size_t i = 0; i < mg.num_levels(); ++i)
				locNumElemsOnLvl.push_back(mg.template num<Volume>(i));
			break;
		case FACE:
			for(size_t i = 0; i < mg.num_levels(); ++i)
				locNumElemsOnLvl.push_back(mg.template num<Face>(i));
			break;
		case EDGE:
			for(size_t i = 0; i < mg.num_levels(); ++i)
				locNumElemsOnLvl.push_back(mg.template num<EdgeBase>(i));
			break;
		case VERTEX:
			for(size_t i = 0; i < mg.num_levels(); ++i)
				locNumElemsOnLvl.push_back(mg.template num<VertexBase>(i));
			break;
		default:
			UG_THROW("Unknown base object type");
			break;
	}

	std::vector<int>	numElemsOnLvl;
	#ifdef UG_PARALLEL
		numElemsOnLvl.resize(locNumElemsOnLvl.size());
		commWorld.allreduce(locNumElemsOnLvl, numElemsOnLvl, PCL_RO_SUM);
	#else
		numElemsOnLvl = locNumElemsOnLvl;
	#endif

	std::vector<int>	subsetDims;
	subsetDims.reserve(sh.num_subsets());
	for(int i = 0; i < sh.num_subsets(); ++i)
		subsetDims.push_back(sh.subset_info(i).get_property("dim").to_int());


	std::vector<int> numLocalGhosts;
	#ifdef UG_PARALLEL
		switch(elemType){
			case VOLUME:
				count_ghosts<Volume>(numLocalGhosts);
				break;
			case FACE:
				count_ghosts<Face>(numLocalGhosts);
				break;
			case EDGE:
				count_ghosts<EdgeBase>(numLocalGhosts);
				break;
			case VERTEX:
				count_ghosts<VertexBase>(numLocalGhosts);
				break;
			default:
				UG_THROW("Unknown base object type");
				break;
		}
	#else
		numLocalGhosts.resize(mg.num_levels(), 0);
	#endif
	m_domainInfo.set_info(elemType, numElemsOnLvl, locNumElemsOnLvl, numLocalGhosts, subsetDims);
}

template <typename TGrid, typename TSubsetHandler>
bool IDomain<TGrid,TSubsetHandler>::
create_additional_subset_handler(std::string name)
{
	if(m_additionalSH[name].valid())
		return false;

	m_additionalSH[name] = SmartPtr<TSubsetHandler>(new TSubsetHandler(*m_spGrid));
	return true;
}

template <typename TGrid, typename TSubsetHandler>
std::vector<std::string> IDomain<TGrid,TSubsetHandler>::
additional_subset_handler_names() const
{
	typedef typename std::map<std::string, SmartPtr<TSubsetHandler> >::const_iterator iterator_t;
	std::vector<std::string> names;
	for(iterator_t iter = m_additionalSH.begin(); iter != m_additionalSH.end(); ++iter){
		if(iter->second.valid())
			names.push_back(iter->first);
	}
	return names;
}

template <typename TGrid, typename TSubsetHandler>
SmartPtr<TSubsetHandler> IDomain<TGrid,TSubsetHandler>::
additional_subset_handler(std::string name)
{
	SmartPtr<TSubsetHandler> sp = m_additionalSH[name];
	if(!sp.valid()){
		UG_THROW("Requested additional subset handler with name '" << name
				 << "' doesn't exist in the given domain!");
	}
	return sp;
}

template <typename TGrid, typename TSubsetHandler>
const ConstSmartPtr<TSubsetHandler> IDomain<TGrid,TSubsetHandler>::
additional_subset_handler(std::string name) const
{
	SmartPtr<TSubsetHandler> sp = m_additionalSH[name];
	if(!sp.valid()){
		UG_THROW("Requested additional subset handler with name '" << name
				 << "' doesn't exist in the given domain!");
	}
	return sp;
}


#ifdef UG_PARALLEL
template <typename TGrid, typename TSubsetHandler>
template <class TElem>
void IDomain<TGrid,TSubsetHandler>::
count_ghosts(std::vector<int>& numGhostsOnLvlOut)
{
	TGrid& mg = *m_spGrid;
	DistributedGridManager& dgm = *mg.distributed_grid_manager();
	GridLayoutMap& glm = dgm.grid_layout_map();

	numGhostsOnLvlOut.clear();
	numGhostsOnLvlOut.resize(mg.num_levels(), 0);

	if(glm.has_layout<TElem>(INT_V_MASTER)){
		typedef typename GridLayoutMap::Types<TElem>::Layout Layout;
		typedef typename GridLayoutMap::Types<TElem>::Interface Interface;
		Layout& layout = glm.get_layout<TElem>(INT_V_MASTER);
		for(size_t lvl = 0; lvl < layout.num_levels(); ++lvl){
			if(lvl >= mg.num_levels())
				break;

			typename Layout::LevelLayout& lvlLayout = layout.layout_on_level(lvl);
			for(typename Layout::LevelLayout::iterator iiter = lvlLayout.begin();
				iiter != lvlLayout.end(); ++iiter)
			{
				Interface& intfc = lvlLayout.interface(iiter);
				for(typename Interface::iterator eiter = intfc.begin();
					eiter != intfc.end(); ++eiter)
				{
					if(dgm.is_ghost(intfc.get_element(eiter)))
						++numGhostsOnLvlOut[lvl];
				}
			}
		}
	}
}
#endif

////////////////////////////////////////////////////////////////////////////////
// Domain
////////////////////////////////////////////////////////////////////////////////

template <int d, typename TGrid, typename TSubsetHandler>
Domain<d,TGrid,TSubsetHandler>::
Domain(bool isAdaptive) : IDomain<TGrid, TSubsetHandler>(isAdaptive)
{
//	Depending on the dimesion, we'll activeate different options.
//	Otherwise we probably would waste memory...
//	In any case, sides of elements should always be present
	uint gridOpts = GRIDOPT_AUTOGENERATE_SIDES;

//	Furthermore vertices should store associated elements.
//	This option depends on the dimension of the domain
	if(dim > 0)
		gridOpts |= VRTOPT_STORE_ASSOCIATED_EDGES;
	if(dim > 1)
		gridOpts |= VRTOPT_STORE_ASSOCIATED_FACES;
	if(dim > 2)
		gridOpts |= VRTOPT_STORE_ASSOCIATED_VOLUMES;

//	thats it for now. One could think about enabling
//	FACEOPT_STORE_ASSOCIATED_EDGES, VOLOPT_STORE_ASSOCIATED_EDGES
//	and VOLOPT_STORE_ASSOCIATED_FACES. However this costs considerably
//	more memory compared to the performance benefits.
//	Now set the options
	this->grid()->set_options(gridOpts);

//	get position attachment
	m_aPos = GetDefaultPositionAttachment<position_attachment_type>();

// 	let position accessor access Vertex Coordinates
	if(!this->grid()->template has_attachment<VertexBase>(m_aPos))
		this->grid()->template attach_to<VertexBase>(m_aPos);
	m_aaPos.access(*(this->grid()), m_aPos);
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOMAIN_IMPL__ */
