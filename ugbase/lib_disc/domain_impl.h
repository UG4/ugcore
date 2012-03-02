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
	m_spGrid(new TGrid(GRIDOPT_NONE)),
	m_spSH(new TSubsetHandler(*m_spGrid)),
	m_isAdaptive(isAdaptive)
{
#ifdef UG_PARALLEL
//	create Distributed Grid Manager
	m_distGridMgr = new DistributedGridManager(*m_spGrid);
#else
	m_dim_distGridMgr = NULL;
#endif

//	register function for grid adaption
	int msgID = GridMessageId_Adaption(message_hub());
	m_spGridAdaptionCallbackID =
		message_hub()->register_class_callback(msgID, this,
		&ug::IDomain<ug::MultiGrid, ug::MultiGridSubsetHandler>::grid_changed_callback);
}

///	Destructor
template <typename TGrid, typename TSubsetHandler>
IDomain<TGrid,TSubsetHandler>::~IDomain()
{
#ifdef UG_PARALLEL
	if(m_distGridMgr) delete m_distGridMgr;
#endif
}

template <typename TGrid, typename TSubsetHandler>
void IDomain<TGrid,TSubsetHandler>::
grid_changed_callback(int, const GridMessage_Adaption* msg)
{
	if(msg->adaptive())
		if(msg->adaption_ends())
		{
			update_local_subset_dim_property();
			#ifdef UG_PARALLEL
			update_local_multi_grid();
			update_global_subset_dim_property();
			#endif
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

	m_spGrid->level_required(numLevGlobal);
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
