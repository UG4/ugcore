/*
 * surface_dof_distribution.cpp
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#include "surface_dof_distribution.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// SurfaceDoFDistribution
////////////////////////////////////////////////////////////////////////////////

///	constructor
SurfaceDoFDistribution::
SurfaceDoFDistribution(SmartPtr<MGSubsetHandler> spMGSH, FunctionPattern& fctPatt,
                       SmartPtr<SurfaceLevelView> spSurfLevelView,
                       int level
#ifdef UG_PARALLEL
                       , DistributedGridManager* pDistGridMgr
#endif
			)
		:	MGDoFDistribution(spMGSH, fctPatt), m_spSurfLevelView(spSurfLevelView),
		 	m_level(level)
#ifdef UG_PARALLEL
			, m_pDistGridMgr(pDistGridMgr)
#endif
{
	init();
}

template <typename TBaseElem>
void SurfaceDoFDistribution::init()
{
	typedef typename traits<TBaseElem>::iterator iterator;
	static const int dim = TBaseElem::dim;

	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no dofs to be distributed
		if(max_dofs(dim, si) == 0) continue;

	//	get iterators of elems
		iterator iter = begin<TBaseElem>(si);
		iterator iterEnd = end<TBaseElem>(si);

	// 	loop elems
		for(; iter != iterEnd; ++iter)
		{
		// 	get vertex
			TBaseElem* elem = *iter;

		//	a) add new indices

		//	add element
			const ReferenceObjectID roid = elem->reference_object_id();
			add(elem, roid, si, m_levInfo);

		//	b) if copy exists, copy also to parent and grand-parents and ... etc.
			//\todo: save this execution in non-adaptive case
			TBaseElem* parent = parent_if_copy(elem);
			while(parent){
				copy(parent, elem);
				elem = parent;
				parent = parent_if_copy(elem);
			}
		}
	} // end subset
}

void SurfaceDoFDistribution::init()
{
	m_levInfo.vNumIndexOnSubset.resize(num_subsets());

	if(max_dofs(VERTEX) > 0) init<VertexBase>();
	if(max_dofs(EDGE) > 0)   init<EdgeBase>();
	if(max_dofs(FACE) > 0)   init<Face>();
	if(max_dofs(VOLUME) > 0) init<Volume>();


#ifdef UG_PARALLEL
	create_layouts_and_communicator();
#endif
}

#ifdef UG_PARALLEL
void SurfaceDoFDistribution::create_layouts_and_communicator()
{
	pcl::ProcessCommunicator commWorld;

//	no grid level present, therefore no dof distribution.
//	just vote false in participate
	if(m_level >= num_levels())
	{
		commWorld.create_sub_communicator(false);
		return;
	}

//  -----------------------------------
//	CREATE PROCESS COMMUNICATOR
//  -----------------------------------
//	The idea  of local processes is to exclude processes from
//	e.g. norm computation that does not have a grid on a given
//	level. If no DoFs exist on the level, that level is excluded from
//	norm computations. In those cases the process votes false for
//	the subcommunicator.

// 	choose if this process participates
	bool participate = !commWorld.empty() && (num_indices() > 0);

//	create process communicator for interprocess layouts
	lev_info().processCommunicator	= commWorld.create_sub_communicator(participate);

//  -----------------------------------
//	CREATE INDEX LAYOUTS ON LEVEL
//  -----------------------------------

	create_index_layout(lev_info().masterLayout, INT_H_MASTER);
	create_index_layout(lev_info().slaveLayout, INT_H_SLAVE);

//	no vertical layouts in surface dof distribution
	lev_info().verticalMasterLayout.clear();
	lev_info().verticalSlaveLayout.clear();
}

void SurfaceDoFDistribution::create_index_layout(IndexLayout& layout,
                                                 int keyType)
{
//	clear layout
	layout.clear();

//	add the index from grid layouts
	if(has_indices_on(VERTEX)) add_indices_from_layouts<VertexBase>(layout, keyType);
	if(has_indices_on(EDGE))   add_indices_from_layouts<EdgeBase>(layout, keyType);
	if(has_indices_on(FACE))   add_indices_from_layouts<Face>(layout, keyType);
	if(has_indices_on(VOLUME)) add_indices_from_layouts<Volume>(layout, keyType);

//	touching an interface means creation. Thus we remove the empty interfaces
//	to avoid storage, communication (should not happen any longer) etc...
	pcl::RemoveEmptyInterfaces(layout);
}

template <typename TBaseElem>
void SurfaceDoFDistribution::add_indices_from_layouts(IndexLayout& indexLayout,
                                                      int keyType)
{
//	get the grid layout map
	GridLayoutMap& layoutMap = m_pDistGridMgr->grid_layout_map();

//	check if layout present
	if(!layoutMap.has_layout<TBaseElem>(keyType)) return;

//	choose level, that must be looped
	int fromLev = 0, toLev = layoutMap.get_layout<VertexBase>(keyType).num_levels() - 1;
	if(m_level == GridLevel::TOPLEVEL){
		if(toLev < 0) toLev = 0;
	}
	else toLev = m_level;

//	loop all level
	for(int level = fromLev; level <= toLev; ++level)
	{
	//	get element layout
		typedef typename GridLayoutMap::Types<TBaseElem>::Layout::LevelLayout TLayout;
		TLayout& elemLayout = layoutMap.get_layout<TBaseElem>(keyType).layout_on_level(level);

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

			//	check if element is on surface (i.e. has no children). Shadows are
			//	not taken into account here, since their indices are already added
			//	to the interface by the shadowing objects
				if(multi_grid().has_children(elem)) {continue;}

			//	check if element is a ghost element, i.e. it is a surface element
			//	but only due to a hierarchical cut of the grid in order to
			//	refine it further on another process. These cuts lead to so called
			//	vertical interfaces.
				if(m_spSurfLevelView->surface_view()->is_ghost(elem)) {continue;}

			//	get the algebraic indices on the grid element
				inner_algebra_indices(elem, vIndex);

			//	add the indices to the interface
				for(size_t i = 0; i < vIndex.size(); ++i)
					indexInterface.push_back(vIndex[i]);
			}
		}
	}
}
#endif

template <typename TBaseElem>
TBaseElem* SurfaceDoFDistribution::parent_if_copy(TBaseElem* elem)
{
	GeometricObject* pParent = m_rMultiGrid.get_parent(elem);
	TBaseElem* parent = dynamic_cast<TBaseElem*>(pParent);
	if(parent != NULL &&
		m_rMultiGrid.num_children<TBaseElem>(parent) == 1) return parent;
	else return NULL;
}

template <typename TBaseElem>
TBaseElem* SurfaceDoFDistribution::parent_if_same_type(TBaseElem* elem)
{
	GeometricObject* pParent = m_rMultiGrid.get_parent(elem);
	return dynamic_cast<TBaseElem*>(pParent);
}

template <typename TBaseElem>
TBaseElem* SurfaceDoFDistribution::child_if_copy(TBaseElem* elem)
{
	if(m_rMultiGrid.num_children<TBaseElem>(elem) != 1) return NULL;
	return m_rMultiGrid.get_child<TBaseElem>(elem, 0);
}

template <typename TBaseElem>
void SurfaceDoFDistribution::defragment(std::vector<std::pair<size_t,size_t> >& vReplaced)
{
	typedef typename traits<TBaseElem>::iterator iterator;
	static const int dim = TBaseElem::dim;

//	if nothing to do, return
	if(lev_info().vFreeIndex.empty()) return;

//	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	// 	skip if no DoFs to be distributed
		if(max_dofs(dim, si) == 0) continue;

		iterator iter = begin<TBaseElem>(si);
		iterator iterEnd = end<TBaseElem>(si);

	// 	loop elems
		for(; iter != iterEnd; ++iter)
		{
		// 	get element
			TBaseElem* elem = *iter;

		//	get roid
			const ReferenceObjectID roid = elem->reference_object_id();

		//	check correct index and replace if needed
			MGDoFDistribution::defragment(elem, roid, si, lev_info(), vReplaced);
		}
	}
}

void SurfaceDoFDistribution::defragment()
{
//	check, if holes exist. If not, we're done
	if(lev_info().vFreeIndex.empty())
	{
	//	the index set might have been increased. Thus resize grid functions
		ManagingDoFDistribution::resize_values(num_indices());

	//	we're done
		return;
	}

//	defragment
	std::vector<std::pair<size_t,size_t> > vReplaced;
	if(max_dofs(VERTEX) > 0) defragment<VertexBase>(vReplaced);
	if(max_dofs(EDGE) > 0) defragment<EdgeBase>(vReplaced);
	if(max_dofs(FACE) > 0) defragment<Face>(vReplaced);
	if(max_dofs(VOLUME) > 0) defragment<Volume>(vReplaced);

//	adapt managed vectors
	copy_values(vReplaced, true);

//	num indices may have changed
	resize_values(num_indices());
}

template <typename TBaseElem>
inline void SurfaceDoFDistribution::obj_created(TBaseElem* obj, GeometricObject* pParent,
                        bool replacesParent)
{
//	case 1: A real insertion in the multigrid: add indices
	if(!replacesParent){
		TBaseElem* copyParent = parent_if_copy(obj);

		// a) Is a copy of underlying object, simply copy indices (i.e. reuse them)
		if(copyParent != NULL) {
			copy(obj, copyParent);
		}
		// b) Cannot copy from parent, create new indices.
		else {
			add(obj,
			    obj->reference_object_id(),
			    m_spMGSH->get_subset_index(obj),
			    m_levInfo);
		}
	}

//	case 2: if replacesParent == false, only an obj (e.g. HangingVertex)
//			is replaced by a similar obj (e.g. Vertex). This case is
//			handled in obj_to_be_erased(), that will be called in addition
//			to this callback with replacedBy != NULL.
}

template <typename TBaseElem>
inline void SurfaceDoFDistribution::obj_to_be_erased(TBaseElem* obj,
                             TBaseElem* replacedBy)
{
//	case 1: Only replacement. Just copy indices from one obj to the other
	if(replacedBy != NULL) {copy(replacedBy, obj); return;}

//	case 2: Element disappears, but parent is identical. Thus, the parent
//			has the same indices attached. All indices remain valid on
//			every surface level. No resizement in the index set must be
//			performed.
	if(parent_if_copy(obj) != NULL) return;

//	case 3: The object that will be erased has no identical parent on the
//			coarser grid. In this case we have to remove the index from
//			the index set.
	erase(obj,
	      obj->reference_object_id(),
	      m_spMGSH->get_subset_index(obj),
	      m_levInfo);
}


void SurfaceDoFDistribution::vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent, bool replacesParent) {obj_created<VertexBase>(vrt, pParent, replacesParent);}
void SurfaceDoFDistribution::edge_created(Grid* grid, EdgeBase* e, GeometricObject* pParent, bool replacesParent) {obj_created<EdgeBase>(e, pParent, replacesParent);}
void SurfaceDoFDistribution::face_created(Grid* grid, Face* f, GeometricObject* pParent, bool replacesParent) {obj_created<Face>(f, pParent, replacesParent);}
void SurfaceDoFDistribution::volume_created(Grid* grid, Volume* vol, GeometricObject* pParent, bool replacesParent) {obj_created<Volume>(vol, pParent, replacesParent);}

void SurfaceDoFDistribution::vertex_to_be_erased(Grid* grid, VertexBase* vrt, VertexBase* replacedBy) {obj_to_be_erased<VertexBase>(vrt, replacedBy);}
void SurfaceDoFDistribution::edge_to_be_erased(Grid* grid, EdgeBase* e, EdgeBase* replacedBy) {obj_to_be_erased<EdgeBase>(e, replacedBy);}
void SurfaceDoFDistribution::face_to_be_erased(Grid* grid, Face* f, Face* replacedBy) {obj_to_be_erased<Face>(f, replacedBy);}
void SurfaceDoFDistribution::volume_to_be_erased(Grid* grid, Volume* vol, Volume* replacedBy) {obj_to_be_erased<Volume>(vol, replacedBy);}


template <typename TBaseElem>
void SurfaceDoFDistribution::
get_connections(std::vector<std::vector<size_t> >& vvConnection) const
{
//	dimension of Base Elem
	static const int dim = TBaseElem::dim;

//	Adjacent geometric objects
	std::vector<VertexBase*> vVrts;
	std::vector<EdgeBase*> vEdges;
	std::vector<Face*> vFaces;
	std::vector<Volume*> vVols;

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
			m_spSurfLevelView->collect_associated(vVrts, elem);
			changable_indices<VertexBase>(vIndex, vVrts);
		}
		if(dim >= EDGE   && max_dofs(EDGE) > 0)	{
			m_spSurfLevelView->collect_associated(vEdges, elem);
			changable_indices<EdgeBase>(vIndex, vEdges);
		}
		if(dim >= FACE   && max_dofs(FACE) > 0)	{
			m_spSurfLevelView->collect_associated(vFaces, elem);
			changable_indices<Face>(vIndex, vFaces);
		}
		if(dim >= VOLUME && max_dofs(VOLUME) > 0) {
			m_spSurfLevelView->collect_associated(vVols, elem);
			changable_indices<Volume>(vIndex, vVols);
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

bool SurfaceDoFDistribution::
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
			{
				UG_LOG("ERROR in 'get_connections':"
						" Currently only implemented iff same number of DoFs on"
						" all geometric objects in all subsets.\n");
				return false;
			}
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
void SurfaceDoFDistribution::permute_indices(const std::vector<size_t>& vNewInd)
{
// 	loop Vertices
	typedef typename traits<TBaseElem>::const_iterator const_iterator;

	const_iterator iterEnd = end<TBaseElem>();
	for(const_iterator iter = begin<TBaseElem>(); iter != iterEnd; ++iter)
	{
	// 	get vertex
		TBaseElem* elem = *iter;

	// 	get current (old) index
		const size_t oldIndex = obj_index(elem);

	//	replace old index by new one
		obj_index(elem) = vNewInd[oldIndex];
	}
}

void SurfaceDoFDistribution::permute_indices(const std::vector<size_t>& vNewInd)
{
	if(max_dofs(VERTEX)) permute_indices<VertexBase>(vNewInd);
	if(max_dofs(EDGE))   permute_indices<EdgeBase>(vNewInd);
	if(max_dofs(FACE))   permute_indices<Face>(vNewInd);
	if(max_dofs(VOLUME)) permute_indices<Volume>(vNewInd);

#ifdef UG_PARALLEL
	create_layouts_and_communicator();
#endif

//	permute indices in associated vectors
	permute_values(vNewInd);
}


} // end namespace ug
