/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__DOMAIN_IMPL__
#define __H__UG__LIB_DISC__DOMAIN_IMPL__

#include "domain.h"
#include "common/serialization.h"
#include "common/profiler/profiler.h"

#ifdef UG_PARALLEL
	#include "lib_grid/refinement/projectors/projectors.h"
	#include "common/boost_serialization_routines.h"
	#include "common/util/archivar.h"
	#include "common/util/factory.h"
	#include <boost/archive/text_oarchive.hpp>
	#include <boost/archive/text_iarchive.hpp>
#endif

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
	PROFILE_FUNC();

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

	// do not communicate if geom is created on all procs
	if (rootProc == -2) return;

#ifdef UG_PARALLEL
	pcl::ProcessCommunicator procCom;

//	prepare the subset-info package, send it to all procs and extract the info again.
	BinaryBuffer buf;
	if(pcl::ProcRank() == rootProc){
		Serialize(buf, sh.num_subsets());
		for(int i = 0; i < sh.num_subsets(); ++i){
			Serialize(buf, sh.subset_info(i).name);
			Serialize(buf, sh.subset_info(i).get_property("dim").to_int());
		}
	}

	procCom.broadcast(buf, rootProc);

	if(pcl::ProcRank() != rootProc){
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

//todo:	distribute projectors from rootProc to all other processors.
//note:	first check whether source-proc has a ProjectionHandler. If so,
//		create a local projection handler first and perform serialization
//		afterwards.
	set_refinement_projector(
			broadcast_refinement_projector(
						rootProc, procCom, geometry3d(), m_refinementProjector));

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
		//if(msg.adaptive()){
			if(msg.adaption_ends())
			{
				update_domain_info();
				m_adaptionIsActive = false;
			}
		//}
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
	PROFILE_FUNC();

	TGrid& mg = *m_spGrid;
	TSubsetHandler& sh = *m_spSH;

	GridBaseObjectId	locElemType;
	if(mg.template num<Volume>() > 0)
		locElemType = VOLUME;
	else if(mg.template num<Face>() > 0)
		locElemType = FACE;
	else if(mg.template num<Edge>() > 0)
		locElemType = EDGE;
	else
		locElemType = VERTEX;

	GridBaseObjectId	elemType;
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
		elemType = (GridBaseObjectId) commWorld.allreduce((int)locElemType, PCL_RO_MAX);
	#else
		elemType = locElemType;
	#endif

//	the number of levels of the multi-grid is now equal on all processes

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
				count_ghosts<Edge>(numLocalGhosts);
				break;
			case VERTEX:
				count_ghosts<Vertex>(numLocalGhosts);
				break;
			default:
				UG_THROW("Unknown base object type");
				break;
		}
	#else
		numLocalGhosts.resize(mg.num_levels(), 0);
	#endif

	std::vector<int>	numLocalElems;
	numLocalElems.reserve(mg.num_levels());

	switch(elemType){
		case VOLUME:
			for(size_t i = 0; i < mg.num_levels(); ++i)
				numLocalElems.push_back(mg.template num<Volume>(i) - numLocalGhosts[i]);
			break;
		case FACE:
			for(size_t i = 0; i < mg.num_levels(); ++i)
				numLocalElems.push_back(mg.template num<Face>(i) - numLocalGhosts[i]);
			break;
		case EDGE:
			for(size_t i = 0; i < mg.num_levels(); ++i)
				numLocalElems.push_back(mg.template num<Edge>(i) - numLocalGhosts[i]);
			break;
		case VERTEX:
			for(size_t i = 0; i < mg.num_levels(); ++i)
				numLocalElems.push_back(mg.template num<Vertex>(i) - numLocalGhosts[i]);
			break;
		default:
			UG_THROW("Unknown base object type");
			break;
	}

	std::vector<int>	numGlobalElems, minNumLocalElems, maxNumLocalElems;
	#ifdef UG_PARALLEL
	//	we have to sum local element counts excluding ghosts.
		numGlobalElems.resize(numLocalElems.size());
		minNumLocalElems.resize(numLocalElems.size());
		maxNumLocalElems.resize(numLocalElems.size());
		commWorld.allreduce(numLocalElems, numGlobalElems, PCL_RO_SUM);
		commWorld.allreduce(numLocalElems, minNumLocalElems, PCL_RO_MIN);
		commWorld.allreduce(numLocalElems, maxNumLocalElems, PCL_RO_MAX);
	#else
		numGlobalElems = numLocalElems;
		minNumLocalElems = numLocalElems;
		maxNumLocalElems = numLocalElems;
	#endif

	m_domainInfo.set_info(elemType, numGlobalElems, numLocalElems, minNumLocalElems,
						  maxNumLocalElems, numLocalGhosts, subsetDims);
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

template <typename TGrid, typename TSubsetHandler>
void IDomain<TGrid, TSubsetHandler>::
set_refinement_projector(SPRefinementProjector proj)
{
	m_refinementProjector = proj;
	if(proj.valid())
		proj->set_geometry(geometry3d());
}

template <typename TGrid, typename TSubsetHandler>
SPRefinementProjector IDomain<TGrid, TSubsetHandler>::
refinement_projector() const
{
	return m_refinementProjector;
}


#ifdef UG_PARALLEL
template <typename TGrid, typename TSubsetHandler>
SPRefinementProjector IDomain<TGrid, TSubsetHandler>::broadcast_refinement_projector
(
	int rootProc,
	pcl::ProcessCommunicator& procCom,
	SPIGeometry3d geometry,
	SPRefinementProjector projector
)
{
	BinaryBuffer 	buf;
	const int		magicNumber	= 3243578;
	const bool 		isRoot		= (pcl::ProcRank() == rootProc);

	static Factory<RefinementProjector, ProjectorTypes>	projFac;

	if(isRoot){
		Archivar<boost::archive::text_oarchive,
				RefinementProjector,
				ProjectorTypes>
			archivar;

	//	if the specified projector is a projection handler, we'll perform a
	//	special operation.
		ProjectionHandler* ph = NULL;
		int projectorType = -1;// -1: none, 0: normal projector, 1: projection handler
		if(projector.valid()){
			ph = dynamic_cast<ProjectionHandler*>(projector.get());
			if(ph)
				projectorType = 1;
			else
				projectorType = 0;
		}

		Serialize(buf, projectorType);
		if(ph){
			const ISubsetHandler* psh = ph->subset_handler();
			if (psh == subset_handler().get())
				Serialize(buf, std::string(""));
			else
			{
				typedef typename std::map<std::string, SmartPtr<TSubsetHandler> >::const_iterator map_it_t;
				map_it_t it = m_additionalSH.begin();
				map_it_t it_end = m_additionalSH.end();
				for (; it != it_end; ++it)
				{
					if (it->second.get() == psh)
					{
						Serialize(buf, it->first);
						break;
					}
				}
				UG_COND_THROW(it == it_end, "Subset handler for projection handler not found "
											"in list of available subset handlers.");
			}

			size_t numProjectors = ph->num_projectors();
			Serialize(buf, numProjectors);
			for(size_t iproj = 0; iproj < numProjectors; ++iproj){
				SPRefinementProjector	proj		= ph->projector(iproj);
				const std::string&			projName 	= projFac.class_name(*proj);
				Serialize(buf, projName);

				std::stringstream ss;
				boost::archive::text_oarchive ar(ss, boost::archive::no_header);
				archivar.archive(ar, *proj);
				Serialize(buf, ss.str());
			}
		}
		else if(projector.valid()){
			RefinementProjector&	proj		= *projector;
			const std::string&			projName 	= projFac.class_name(proj);
			Serialize(buf, projName);

			std::stringstream ss;
			boost::archive::text_oarchive ar(ss, boost::archive::no_header);
			archivar.archive(ar, proj);
			Serialize(buf, ss.str());
		}

		Serialize(buf, magicNumber);
	}

	procCom.broadcast(buf, rootProc);

	if(!isRoot){
		Archivar<boost::archive::text_iarchive,
				RefinementProjector,
				ProjectorTypes>
			archivar;

		int projectorType;
		Deserialize(buf, projectorType);
		if(projectorType == 1){
			std::string sh_name;
			Deserialize(buf, sh_name);
			ProjectionHandler* ph;
			if (sh_name == std::string(""))
				ph = new ProjectionHandler(geometry, subset_handler());
			else
			{
				typedef typename std::map<std::string, SmartPtr<TSubsetHandler> >::const_iterator map_it_t;
				map_it_t it = m_additionalSH.begin();
				map_it_t it_end = m_additionalSH.end();
				for (; it != it_end; ++it)
				{
					if (it->first == sh_name)
					{
						ph = new ProjectionHandler(geometry, it->second);
						break;
					}
				}
				UG_COND_THROW(it == it_end, "Subset handler name for projection handler not found "
											"in list of available names.");
			}
			SPProjectionHandler projHandler = make_sp(ph);

			size_t numProjectors;
			Deserialize(buf, numProjectors);

			for(size_t iproj = 0; iproj < numProjectors; ++iproj){
				std::string name;
				Deserialize(buf, name);
				SPRefinementProjector proj = projFac.create(name);

				std::string data;
				Deserialize(buf, data);
				std::stringstream ss(data, std::ios_base::in);
				boost::archive::text_iarchive ar(ss, boost::archive::no_header);
				archivar.archive(ar, *proj);

				ph->set_projector(iproj, proj);
			}

			projector = projHandler;
		}
		else if(projectorType == 0){
			std::string name;
			Deserialize(buf, name);
			SPRefinementProjector proj = projFac.create(name);
			proj->set_geometry(geometry);

			std::string data;
			Deserialize(buf, data);
			std::stringstream ss(data, std::ios_base::in);
			boost::archive::text_iarchive ar(ss, boost::archive::no_header);
			archivar.archive(ar, *proj);

			projector = proj;
		}
		else if(projectorType == -1){
			projector = SPNULL;
		}
		else{
			UG_THROW("Invalid projector type in 'BroadcastRefinementProjector': "
					 << projectorType);
		}

		int tmp;
		Deserialize(buf, tmp);
		UG_COND_THROW(tmp != magicNumber, "Magic number mismatch in "
					  "'BroadcastRefinementProjector'. Received "
					  << tmp << ", but expected " << magicNumber);
	}

	return projector;
}

#endif


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
	if(!this->grid()->template has_attachment<Vertex>(m_aPos))
		this->grid()->template attach_to<Vertex>(m_aPos);
	m_aaPos.access(*(this->grid()), m_aPos);

	m_geometry3d = MakeGeometry3d(*(this->grid()), m_aPos);
	this->m_refinementProjector = make_sp(new RefinementProjector(m_geometry3d));
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOMAIN_IMPL__ */
