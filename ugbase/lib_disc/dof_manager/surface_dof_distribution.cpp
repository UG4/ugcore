/*
 * surface_dof_distribution.cpp
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#include "surface_dof_distribution.h"
#include "mg_dof_distribution_impl.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/algorithms/debug_util.h"

using namespace std;

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// SurfaceDoFDistribution
////////////////////////////////////////////////////////////////////////////////

///	constructor
SurfaceDoFDistribution::
SurfaceDoFDistribution(SmartPtr<MultiGrid> spMG,
                       SmartPtr<MGSubsetHandler> spMGSH,
                       ConstSmartPtr<DoFDistributionInfo> spDDInfo,
                       SmartPtr<SurfaceView> spSurfView,
                       int level, bool bGrouped)
		:	DoFDistributionInfoProvider(spDDInfo),
		 	MGDoFDistribution(spMG, spMGSH, spDDInfo, bGrouped),
		 	DoFDistribution(*this, spSurfView, GridLevel(level, GridLevel::SURFACE, false)),
		 	m_spSurfView(spSurfView),
		 	m_level(level)
{
	reinit();
}


template <typename TBaseElem>
void SurfaceDoFDistribution::reinit()
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

		//	b) if shadow exists, handle also to parent and grand-parents and ... etc.
			TBaseElem* parent = dynamic_cast<TBaseElem*>(this->m_spMG->get_parent(elem));
			while(parent && m_spSurfView->is_shadowed(parent)){

				// if it is a copy element: assign same index
				if(this->m_spMG->template num_children<TBaseElem>(parent) == 1){
					obj_index(parent) = obj_index(elem);
				}
				// if no copy, assign new index
				else{
					add(parent, roid, si, m_levInfo);
				}

				elem = parent;
				parent = dynamic_cast<TBaseElem*>(this->m_spMG->get_parent(elem));
			}

		}
	} // end subset
}

void SurfaceDoFDistribution::reinit()
{
	m_levInfo.clear();
	m_levInfo.vNumIndexOnSubset.resize(num_subsets());

	if(max_dofs(VERTEX)) reinit<VertexBase>();
	if(max_dofs(EDGE))   reinit<EdgeBase>();
	if(max_dofs(FACE))   reinit<Face>();
	if(max_dofs(VOLUME)) reinit<Volume>();

#ifdef UG_PARALLEL
	create_layouts_and_communicator();
#endif
}


#ifdef UG_PARALLEL
void SurfaceDoFDistribution::create_layouts_and_communicator()
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
	bool participate = !commWorld.empty() && (num_indices() > 0);

//	create process communicator for interprocess layouts
	lev_info().layouts()->proc_comm() = commWorld.create_sub_communicator(participate);

//  -----------------------------------
//	CREATE INDEX LAYOUTS ON LEVEL
//  -----------------------------------

	create_index_layout(lev_info().layouts()->master(), INT_H_MASTER);
	create_index_layout(lev_info().layouts()->slave(), INT_H_SLAVE);

//	no vertical layouts in surface dof distribution
	lev_info().layouts()->vertical_master().clear();
	lev_info().layouts()->vertical_slave().clear();
}

void SurfaceDoFDistribution::create_index_layout(IndexLayout& layout,
                                                 int keyType)
{
//	clear layout
	layout.clear();

//	add the index from grid layouts
	if(max_dofs(VERTEX)) add_indices_from_layouts<VertexBase>(layout, keyType);
	if(max_dofs(EDGE))   add_indices_from_layouts<EdgeBase>(layout, keyType);
	if(max_dofs(FACE))   add_indices_from_layouts<Face>(layout, keyType);
	if(max_dofs(VOLUME)) add_indices_from_layouts<Volume>(layout, keyType);

//	touching an interface means creation. Thus we remove the empty interfaces
//	to avoid storage, communication (should not happen any longer) etc...
	pcl::RemoveEmptyInterfaces(layout);
}

template <typename TBaseElem>
void SurfaceDoFDistribution::add_indices_from_layouts(IndexLayout& indexLayout,
                                                      int keyType)
{
//	get the grid layout map
	GridLayoutMap& layoutMap = m_spMG->distributed_grid_manager()->grid_layout_map();

//	check if layout present
	if(!layoutMap.has_layout<TBaseElem>(keyType)) return;

//	choose level, that must be looped
	int fromLev = 0, toLev = layoutMap.get_layout<TBaseElem>(keyType).num_levels() - 1;
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
				//if(multi_grid()->has_children(elem)) {continue;}

			//	check if element is a ghost element, i.e. it is a surface element
			//	but only due to a hierarchical cut of the grid in order to
			//	refine it further on another process. These cuts lead to so called
			//	vertical interfaces.
				//if(m_spSurfView->is_ghost(elem)) {continue;}
				if(!m_spSurfView->is_surface_element(elem)) {continue;}

			//	get the algebraic indices on the grid element
				MGDoFDistribution::inner_algebra_indices(elem, vIndex);

			//	add the indices to the interface
				for(size_t i = 0; i < vIndex.size(); ++i)
					indexInterface.push_back(vIndex[i]);
			}
		}
	}
}
#endif

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
			m_spSurfView->collect_associated(vVrts, elem);
			changable_indices<VertexBase>(vIndex, vVrts);
		}
		if(dim >= EDGE   && max_dofs(EDGE) > 0)	{
			m_spSurfView->collect_associated(vEdges, elem);
			changable_indices<EdgeBase>(vIndex, vEdges);
		}
		if(dim >= FACE   && max_dofs(FACE) > 0)	{
			m_spSurfView->collect_associated(vFaces, elem);
			changable_indices<Face>(vIndex, vFaces);
		}
		if(dim >= VOLUME && max_dofs(VOLUME) > 0) {
			m_spSurfView->collect_associated(vVols, elem);
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
	if(!MGDoFDistribution::grouped())
	{
		size_t numDoFs = 0;
		for(int si = 0; si < num_subsets(); ++si)
			for(int i = 0; i < NUM_REFERENCE_OBJECTS; ++i)
			{
				const ReferenceObjectID roid = (ReferenceObjectID) i;
				if(num_dofs(roid,si) == 0) continue;
				if(numDoFs == 0) {numDoFs = num_dofs(roid,si); continue;}
				if(num_dofs(roid,si) != numDoFs)
					UG_THROW("SurfaceDoFDistribution::get_connections: "
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

	//	if copy exists, copy also to parent and grand-parents and ... etc.
		//\todo: save this execution in non-adaptive case
		TBaseElem* parent = parent_if_shadowed_copy(elem);
		while(parent){
			obj_index(parent) = obj_index(elem);
			elem = parent;
			parent = parent_if_shadowed_copy(elem);
		}
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
