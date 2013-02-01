/*
 * domain_impl.h
 *
 *  Created on: 02.03.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOMAIN_IMPL__
#define __H__UG__LIB_DISC__DOMAIN_IMPL__

#include "domain.h"

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
	int msgID = GridMessageId_Adaption(message_hub());
	m_spGridAdaptionCallbackID =
		message_hub()->register_class_callback(msgID, this,
		&ug::IDomain<ug::MultiGrid, ug::MultiGridSubsetHandler>::grid_changed_callback);

	msgID = GridMessageId_Distribution(message_hub());
	m_spGridDistributionCallbackID =
		message_hub()->register_class_callback(msgID, this,
		&ug::IDomain<ug::MultiGrid, ug::MultiGridSubsetHandler>::grid_distributed_callback);
}

///	Destructor
template <typename TGrid, typename TSubsetHandler>
IDomain<TGrid,TSubsetHandler>::~IDomain()
{
}

template <typename TGrid, typename TSubsetHandler>
inline
void IDomain<TGrid,TSubsetHandler>::
grid_changed_callback(int, const GridMessage_Adaption* msg)
{
	if(msg->adaption_begins())
		m_adaptionIsActive = true;

	else if(m_adaptionIsActive){
		if(msg->adaptive()){
			if(msg->adaption_ends())
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
grid_distributed_callback(int, const GridMessage_Distribution* msg)
{
	if(msg->msg() == GMDT_DISTRIBUTION_STOPS){
//todo	update_domain_info may be unnecessary here and other actions may have
//		to be performed.
		update_domain_info();
	}
}


#ifdef UG_PARALLEL
template <typename TGrid, typename TSubsetHandler>
void IDomain<TGrid,TSubsetHandler>::
update_local_multi_grid()
{
//	proc local number of level
	int numLevLocal = m_spGrid->num_levels();

//	storage for global number of levels
	int numLevGlobal;

	pcl::ProcessCommunicator commWorld;
	commWorld.allreduce(&numLevLocal, &numLevGlobal, 1, PCL_DT_INT, PCL_RO_MAX);

	if(numLevGlobal > 0){
		m_spGrid->level_required(numLevGlobal-1);
		m_spSH->level_required(numLevGlobal-1);
	}
}
#endif


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


	update_local_subset_dim_property();

	GeometricBaseObject	elemType;
	#ifdef UG_PARALLEL
		pcl::ProcessCommunicator commWorld;

		update_local_multi_grid();
		update_global_subset_dim_property();


		elemType = (GeometricBaseObject) commWorld.allreduce((int)locElemType, PCL_RO_MAX);

	#else
		elemType = locElemType;
	#endif

//	the number of levels of the multi-grid and the number of subsets
//	is now equal on all processes
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
		numElemsOnLvl.swap(locNumElemsOnLvl);
	#endif

	std::vector<int>	subsetDims;
	subsetDims.reserve(sh.num_subsets());
	for(int i = 0; i < sh.num_subsets(); ++i)
		subsetDims.push_back(sh.subset_info(i).get_property("dim").to_int());

	m_domainInfo.set_info(elemType, numElemsOnLvl, subsetDims);
}


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
