/*
 * level_dof_distribution.cpp
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#include "level_dof_distribution.h"
#include "mg_dof_distribution_impl.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// LevelMGDoFDistribution
////////////////////////////////////////////////////////////////////////////////

LevelMGDoFDistribution::
LevelMGDoFDistribution(SmartPtr<MultiGrid> spMG,
					   SmartPtr<MGSubsetHandler> spMGSH,
					   const DoFDistributionInfo& rDDInfo,
                       bool bGrouped)
	:	MGDoFDistribution(spMG, spMGSH, rDDInfo, bGrouped)
{
	if(num_levels() > 0) level_required(num_levels()-1);
	init();
}

template <typename TBaseElem>
void LevelMGDoFDistribution::init()
{
	typedef typename geometry_traits<TBaseElem>::iterator iterator;
	static const int dim = TBaseElem::dim;

//	check if indices in the dimension
	if(max_dofs(dim) == 0) return;

	for(int l = 0; l < num_levels(); ++l)
	{
		for(int si = 0; si < num_subsets(); ++si)
		{
		// 	skip if no dofs to be distributed
			if(max_dofs(dim, si) == 0) continue;

		//	get iterators of elems
			iterator iter = m_spMGSH->begin<TBaseElem>(si,l);
			iterator iterEnd = m_spMGSH->end<TBaseElem>(si,l);

		// 	loop elems
			for(; iter != iterEnd; ++iter)
			{
			// 	get vertex
				TBaseElem* elem = *iter;

			//	get roid
				const ReferenceObjectID roid = elem->reference_object_id();

			//	add element
				add(elem, roid, si, m_vLev[l]);
			}
		} // end subset
	} // end level
}

void LevelMGDoFDistribution::init()
{
	level_required(num_levels() - 1);

	if(max_dofs(0) > 0) init<VertexBase>();
	if(max_dofs(1) > 0) init<EdgeBase>();
	if(max_dofs(2) > 0) init<Face>();
	if(max_dofs(3) > 0) init<Volume>();

#ifdef UG_PARALLEL
	m_pDistGridMgr = m_spMG->distributed_grid_manager();
	for(int l = 0; l < num_levels(); ++l)
		create_layouts_and_communicator(l);
#endif
}

void LevelMGDoFDistribution::redistribute_dofs()
{
	for(int l = 0; l < num_levels(); ++l){
		lev_info(l).clear_all();
		lev_info(l).vNumIndexOnSubset.resize(num_subsets(), 0);
	}

	init();

	for(int l = 0; l < num_levels(); ++l){
		if(m_managingDoFDists[l])
			m_managingDoFDists[l]->resize_values(num_indices(l));
	}
}

void LevelMGDoFDistribution::
register_managing_dof_distribution(ManagingDoFDistribution* mdd, int lvl)
{
	level_required(lvl);
	m_managingDoFDists[lvl] = mdd;
}

#ifdef UG_PARALLEL
void LevelMGDoFDistribution::create_layouts_and_communicator(int l)
{
	pcl::ProcessCommunicator commWorld;

//  -----------------------------------
//	CREATE PROCESS COMMUNICATOR
//  -----------------------------------
//	The idea  of local processes is to exclude processes from
//	e.g. norm computation that does not have a grid on a given
//	level. If no DoFs exist on the level, that level is excluded from
//	norm computations. In those cases the process votes false for
//	the subcommunicator.

// 	choose if this process participates
	bool participate = !commWorld.empty() && (num_indices(l) > 0);

	UG_DLOG_ALL_PROCS(LIB_DISC_MULTIGRID, 2,
					  ": Participate = "<< participate <<
					  " for level "<<l<<" (num_indices="<<num_indices(l)<<
					  ",!empty=" << !commWorld.empty() << ").\n");

//	create process communicator for interprocess layouts
	level_required(l);
	lev_info(l).layouts().proc_comm() = commWorld.create_sub_communicator(participate);

//  -----------------------------------
//	CREATE INDEX LAYOUTS ON LEVEL
//  -----------------------------------

	create_index_layout(lev_info(l).layouts().master(), INT_H_MASTER, l);
	create_index_layout(lev_info(l).layouts().slave(), INT_H_SLAVE, l);
	create_index_layout(lev_info(l).layouts().vertical_master(), INT_V_MASTER, l);
	create_index_layout(lev_info(l).layouts().vertical_slave(), INT_V_SLAVE, l);
}

void LevelMGDoFDistribution::create_index_layout(IndexLayout& layout,
												 InterfaceNodeTypes keyType,
                                                 int l)
{
//	clear layout
	layout.clear();

//	add the index from grid layouts
	if(max_dofs(VERTEX)) add_indices_from_layouts<VertexBase>(layout, keyType, l);
	if(max_dofs(EDGE))   add_indices_from_layouts<EdgeBase>(layout, keyType, l);
	if(max_dofs(FACE))   add_indices_from_layouts<Face>(layout, keyType, l);
	if(max_dofs(VOLUME)) add_indices_from_layouts<Volume>(layout, keyType, l);

//	touching an interface means creation. Thus we remove the empty interfaces
//	to avoid storage, communication (should not happen any longer) etc...
	pcl::RemoveEmptyInterfaces(layout);
}

template <typename TBaseElem>
void LevelMGDoFDistribution::add_indices_from_layouts(IndexLayout& indexLayout,
													  InterfaceNodeTypes keyType,
                                                      int l)
{
//	get the grid layout map
	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

//	check if layout present
	if(!layoutMap.has_layout<TBaseElem>(keyType)) return;

//	get element layout
	typedef typename GridLayoutMap::Types<TBaseElem>::Layout::LevelLayout TLayout;
	TLayout& elemLayout = layoutMap.get_layout<TBaseElem>(keyType).layout_on_level(l);

//	iterator for grid element interfaces
	typedef typename TLayout::iterator InterfaceIterator;

//	type of grid element interfaces
	typedef typename TLayout::Interface ElemInterface;

//	iterator for grid elements
	typedef typename ElemInterface::iterator ElemIterator;

//	type of index interfaces
	typedef IndexLayout::Interface IndexInterface;

//	vector for algebra indices
	std::vector<size_t> vIndex;

//	iterate over all grid element interfaces
	for(InterfaceIterator iIter = elemLayout.begin();
		iIter != elemLayout.end(); ++iIter)
	{
	//	get a grid element interface
		ElemInterface& elemInterface = elemLayout.interface(iIter);

	//	get a corresponding index interface
		IndexInterface& indexInterface = indexLayout.interface(
											elemLayout.proc_id(iIter));

	//	iterate over entries in the grid element interface
		for(ElemIterator eIter = elemInterface.begin();
			eIter != elemInterface.end(); ++eIter)
		{
		//	get the grid element
			typename ElemInterface::Element elem = elemInterface.get_element(eIter);

		//	get the algebraic indices on the grid element
			inner_algebra_indices(elem, vIndex);

		//	add the indices to the interface
			for(size_t i = 0; i < vIndex.size(); ++i)
				indexInterface.push_back(vIndex[i]);
		}
	}
}
#endif

template <typename TBaseElem>
void LevelMGDoFDistribution::defragment(std::vector<std::pair<size_t,size_t> >& vReplaced, int l)
{
	typedef typename geometry_traits<TBaseElem>::iterator iterator;
	static const int dim = TBaseElem::dim;

	if(l < 0 || l >= num_levels())
		UG_THROW("Level does not exist.");

//	if nothing to do, continue
	if(!m_vLev[l].free_index_available()) return;

	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no dofs to be distributed
		if(max_dofs(dim, si) == 0) continue;

	//	get iterators of elems
		iterator iter = m_spMGSH->begin<TBaseElem>(si,l);
		iterator iterEnd = m_spMGSH->end<TBaseElem>(si,l);

	// 	loop elems
		for(; iter != iterEnd; ++iter)
		{
		// 	get vertex
			TBaseElem* elem = *iter;

		//	get roid
			const ReferenceObjectID roid = elem->reference_object_id();

		//	check correct index and replace if needed
			MGDoFDistribution::defragment(elem, roid, si, m_vLev[l], vReplaced);
		}
	} // end subset
}

void LevelMGDoFDistribution::defragment(std::vector<std::pair<size_t,size_t> >& vReplaced, int lev)
{
	level_required(lev);

	if(max_dofs(0) > 0) defragment<VertexBase>(vReplaced, lev);
	if(max_dofs(1) > 0) defragment<EdgeBase>(vReplaced, lev);
	if(max_dofs(2) > 0) defragment<Face>(vReplaced, lev);
	if(max_dofs(3) > 0) defragment<Volume>(vReplaced, lev);

//	check that only invalid indices left
	for(LevInfo<std::vector<size_t> >::iterator it = lev_info(lev).begin(); it != lev_info(lev).end(); ++it)
		UG_ASSERT(*it >= lev_info(lev).numIndex, "After defragment still index in "
								  "valid range present in free index container.");

//	clear container
	lev_info(lev).sizeIndexSet -= lev_info(lev).num_free_index();
	lev_info(lev).clear();

	if(lev_info(lev).free_index_available())
		UG_THROW("Internal error: Still free indices available after "
						"defragment: " <<  lev_info(lev).num_free_index());

	if(lev_info(lev).numIndex != lev_info(lev).sizeIndexSet)
		UG_THROW("Internal error: numIndex and sizeIndexSet must be "
						"equal after defragment, since the index set does not "
						"contain holes anymore. But numIndex = "<<lev_info(lev).numIndex
						<<", sizeIndexSet = "<<lev_info(lev).sizeIndexSet);
}

template <typename TBaseElem>
inline void LevelMGDoFDistribution::obj_created(TBaseElem* obj, GeometricObject* pParent,
                        bool replacesParent)
{
	if(is_frozen())
		return;

//	check level
	const int lev = m_spMGSH->get_level(obj);
	level_required(lev);

//	add indices
	add_from_free(obj,
	              obj->reference_object_id(),
	              m_spMGSH->get_subset_index(obj),
	              m_vLev[lev]);
}

template <typename TBaseElem>
inline void LevelMGDoFDistribution::obj_to_be_erased(TBaseElem* obj,TBaseElem* replacedBy)
{
	if(is_frozen())
		return;

//	add indices
	erase(obj,
	      obj->reference_object_id(),
	      m_spMGSH->get_subset_index(obj),
	      m_vLev[m_spMGSH->get_level(obj)]);
}

void LevelMGDoFDistribution::vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent, bool replacesParent) {obj_created(vrt, pParent, replacesParent);}
void LevelMGDoFDistribution::edge_created(Grid* grid, EdgeBase* e, GeometricObject* pParent, bool replacesParent) {obj_created(e, pParent, replacesParent);}
void LevelMGDoFDistribution::face_created(Grid* grid, Face* f, GeometricObject* pParent, bool replacesParent) {obj_created(f, pParent, replacesParent);}
void LevelMGDoFDistribution::volume_created(Grid* grid, Volume* vol, GeometricObject* pParent, bool replacesParent) {obj_created(vol, pParent, replacesParent);}

void LevelMGDoFDistribution::vertex_to_be_erased(Grid* grid, VertexBase* vrt, VertexBase* replacedBy) {obj_to_be_erased(vrt, replacedBy);}
void LevelMGDoFDistribution::edge_to_be_erased(Grid* grid, EdgeBase* e, EdgeBase* replacedBy) {obj_to_be_erased(e, replacedBy);}
void LevelMGDoFDistribution::face_to_be_erased(Grid* grid, Face* f, Face* replacedBy) {obj_to_be_erased(f, replacedBy);}
void LevelMGDoFDistribution::volume_to_be_erased(Grid* grid, Volume* vol, Volume* replacedBy) {obj_to_be_erased(vol, replacedBy);}

////////////////////////////////////////////////////////////////////////////////
// LevelDoFDistribution
////////////////////////////////////////////////////////////////////////////////
LevelDoFDistribution::
LevelDoFDistribution(SmartPtr<LevelMGDoFDistribution> spLevMGDD,
                     SmartPtr<SurfaceView> spSurfView,
                     int level)
	: 	DoFDistributionInfoProvider(spLevMGDD->dof_distribution_info()),
	  	DoFDistribution(spLevMGDD, spSurfView, GridLevel(level, GridLevel::LEVEL, true)),
		m_spMGDD(spLevMGDD), m_spSurfView(spSurfView)
{
	spLevMGDD->register_managing_dof_distribution(this, level);
};

LevelDoFDistribution::
~LevelDoFDistribution()
{
	m_spMGDD->register_managing_dof_distribution(NULL, grid_level().level());
}

template <typename TBaseElem>
void LevelDoFDistribution::
get_connections(std::vector<std::vector<size_t> >& vvConnection) const
{
//	dimension of Base Elem
	static const int dim = TBaseElem::dim;

//	Adjacent geometric objects
	std::vector<VertexBase*> vVrts;
	std::vector<EdgeBase*> vEdges;
	std::vector<Face*> vFaces;
	std::vector<Volume*> vVols;

//	Multigrid
	MultiGrid& rMultiGrid = *const_cast<MultiGrid*>(&(*m_spMGDD->multi_grid()));

// 	Iterators
	typedef typename traits<TBaseElem>::const_iterator const_iterator;
	const_iterator iterEnd = end<TBaseElem>();

//	loop elem
	for(const_iterator iter = begin<TBaseElem>(); iter != iterEnd; ++iter)
	{
	// 	Get elem
		TBaseElem* elem = *iter;

	//	vector of indices
		std::vector<size_t> vIndex;

	//	Get connected elements
		if(dim >= VERTEX && max_dofs(VERTEX) > 0) {
			CollectAssociated(vVrts, rMultiGrid, elem);
			m_spMGDD->changable_indices<VertexBase>(vIndex, vVrts);
		}
		if(dim >= EDGE   && max_dofs(EDGE) > 0)	{
			CollectAssociated(vEdges, rMultiGrid, elem);
			m_spMGDD->changable_indices<EdgeBase>(vIndex, vEdges);
		}
		if(dim >= FACE   && max_dofs(FACE) > 0)	{
			CollectAssociated(vFaces, rMultiGrid, elem);
			m_spMGDD->changable_indices<Face>(vIndex, vFaces);
		}
		if(dim >= VOLUME && max_dofs(VOLUME) > 0) {
			CollectAssociated(vVols, rMultiGrid, elem);
			m_spMGDD->changable_indices<Volume>(vIndex, vVols);
		}

	//	remove doubles
        std::sort(vIndex.begin(), vIndex.end());
        vIndex.erase(std::unique(vIndex.begin(), vIndex.end()), vIndex.end());

    //	add coupling to adjacency graph
		std::vector<size_t>::iterator it;

		for(size_t i = 0; i < vIndex.size(); ++i)
			for(size_t j = i; j < vIndex.size(); ++j)
			{
				const size_t iIndex = vIndex[i];
				const size_t jIndex = vIndex[j];

			//	add connection (iIndex->jIndex)
				it = find(vvConnection[iIndex].begin(), vvConnection[iIndex].end(), jIndex);
				if(it == vvConnection[iIndex].end()) vvConnection[iIndex].push_back(jIndex);


			//	add the opposite direction (adjInd -> index)
				it = find(vvConnection[jIndex].begin(), vvConnection[jIndex].end(), iIndex);
				if(it == vvConnection[jIndex].end()) vvConnection[jIndex].push_back(iIndex);
			}
	}
}

bool LevelDoFDistribution::
get_connections(std::vector<std::vector<size_t> >& vvConnection) const
{
//	if no subset given, we're done
	if(num_subsets() == 0) return true;

//	check that in all subsets same number of functions and at least one
//	if this is not the case for non-grouped DoFs, we cannot allow reordering
	if(!grouped())
	{
		size_t numDoFs = 0;
		for(int si = 0; si < num_subsets(); ++si)
			for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
			{
				const ReferenceObjectID roid = (ReferenceObjectID) i;
				if(num_dofs(roid,si) == 0) continue;
				if(numDoFs == 0) {numDoFs = num_dofs(roid,si); continue;}
				if(num_dofs(roid,si) != numDoFs)
					UG_THROW("LevelDoFDistribution::get_connections: "
							"Currently only implemented iff same number of DoFs on"
							" all geometric objects in all subsets: \n"
							"num_dofs("<<roid<<","<<si<<")="<<num_dofs(roid,si)<<
							", but previously found "<<numDoFs);
			}
	}

//	clear neighbors
	vvConnection.clear(); vvConnection.resize(num_indices());

//	get_connections<VertexBase>(vvConnection);
	get_connections<EdgeBase>(vvConnection);
	get_connections<Face>(vvConnection);
	get_connections<Volume>(vvConnection);

	return true;
}

template <typename TBaseElem>
void LevelMGDoFDistribution::permute_indices(const std::vector<size_t>& vNewInd, int lev)
{
// 	loop Vertices
	typedef typename geometry_traits<TBaseElem>::const_iterator const_iterator;
	const_iterator iterEnd = m_spMGSH->multi_grid()->template end<TBaseElem>(lev);

	for(const_iterator iter = m_spMGSH->multi_grid()->template begin<TBaseElem>(lev); iter != iterEnd; ++iter)
	{
	// 	get vertex
		TBaseElem* elem = *iter;

	// 	get current (old) index
		const size_t oldIndex = obj_index(elem);

	//	replace old index by new one
		obj_index(elem) = vNewInd[oldIndex];
	}
}

void LevelMGDoFDistribution::permute_indices(const std::vector<size_t>& vNewInd, int lev)
{
	if(max_dofs(VERTEX)) permute_indices<VertexBase>(vNewInd, lev);
	if(max_dofs(EDGE))   permute_indices<EdgeBase>(vNewInd, lev);
	if(max_dofs(FACE))   permute_indices<Face>(vNewInd, lev);
	if(max_dofs(VOLUME)) permute_indices<Volume>(vNewInd, lev);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void LevelMGDoFDistribution::level_required(int level)
{
//	check if something to do
	if(level < (int)m_vLev.size()) return;

//	resize info vector
	m_vLev.resize(level+1);
	m_managingDoFDists.resize(level+1, NULL);

//	adjust subsets
	for(size_t l = 0; l < m_vLev.size(); ++l)
	{
		m_vLev[l].vNumIndexOnSubset.resize(num_subsets(), 0);
	}
}

//////////////////////////////////////////////////////////////////////////////
// LevelDoFDistribution
//////////////////////////////////////////////////////////////////////////////

void LevelDoFDistribution::
permute_indices(const std::vector<size_t>& vIndNew)
{
	m_spMGDD->permute_indices(vIndNew, grid_level().level());
#ifdef UG_PARALLEL
	m_spMGDD->create_layouts_and_communicator(grid_level().level());
#endif

//	permute values in managed grid functions
	permute_values(vIndNew);
}

void LevelDoFDistribution::
defragment()
{
//	defragment
	std::vector<std::pair<size_t,size_t> > vReplaced;
	m_spMGDD->defragment(vReplaced, grid_level().level());

//	adapt managed vectors
	if(!vReplaced.empty())
		copy_values(vReplaced, true);

//	num indices may have changed
	resize_values(num_indices());

#ifdef UG_PARALLEL
	m_spMGDD->create_layouts_and_communicator(grid_level().level());
#endif
}


} // end namespace ug
