
#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__ADAPTION_SURFACE_GRID_FUNCTION_IMPL__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__ADAPTION_SURFACE_GRID_FUNCTION_IMPL__

#include "adaption_surface_grid_function.h"
#include "common/profiler/profiler.h"
#include "lib_grid/algorithms/debug_util.h"

#define GFUNCADAPT_PROFILE_FUNC()	PROFILE_FUNC_GROUP("gfunc-adapt")


namespace ug{

template <typename TDomain>
template <typename TElem, typename TAlgebra>
void AdaptionSurfaceGridFunction<TDomain>::
copy_from_surface(const GridFunction<TDomain,TAlgebra>& rSurfaceFct, TElem* elem)
{
	std::vector<DoFIndex> vInd;
	const size_t numFct = m_spDDInfo->num_fct();

	std::vector<std::vector<number> >& vvVal = m_aaValue[elem];
	vvVal.resize(numFct);

	for(size_t fct = 0; fct < numFct; ++fct){

		std::vector<number>& vVal = vvVal[fct];

		rSurfaceFct.inner_dof_indices(elem, fct, vInd);
		vVal.resize(vInd.size());
		for(size_t i = 0; i < vInd.size(); ++i)
			vVal[i] = DoFRef(rSurfaceFct, vInd[i]);
	}
}

template <typename TDomain>
template <typename TElem, typename TAlgebra>
void AdaptionSurfaceGridFunction<TDomain>::
copy_from_surface(const GridFunction<TDomain,TAlgebra>& rSurfaceFct)
{
	ConstSmartPtr<SurfaceView> spSurfView = rSurfaceFct.approx_space()->surface_view();
	ConstSmartPtr<MultiGrid> spGrid = m_spDomain->grid();
	typedef typename GridFunction<TDomain,TAlgebra>::template traits<TElem>::const_iterator iter_type;
	iter_type iter = rSurfaceFct.template begin<TElem>(SurfaceView::ALL);
	iter_type iterEnd = rSurfaceFct.template end<TElem>(SurfaceView::ALL);

	for( ; iter != iterEnd; ++iter){

		// get surface element
		TElem* elem = *iter;

		// copy indices to attachment
		copy_from_surface(rSurfaceFct, elem);
	}
}

template <typename TDomain>
template <typename TAlgebra>
void AdaptionSurfaceGridFunction<TDomain>::
copy_from_surface(const GridFunction<TDomain,TAlgebra>& rSurfaceFct)
{
	GFUNCADAPT_PROFILE_FUNC();
	attach_entries(rSurfaceFct.dof_distribution()->dof_distribution_info());

	if(rSurfaceFct.max_dofs(VERTEX))copy_from_surface<Vertex,TAlgebra>(rSurfaceFct);
	if(rSurfaceFct.max_dofs(EDGE))	copy_from_surface<Edge,TAlgebra>(rSurfaceFct);
	if(rSurfaceFct.max_dofs(FACE))	copy_from_surface<Face,TAlgebra>(rSurfaceFct);
	if(rSurfaceFct.max_dofs(VOLUME))copy_from_surface<Volume,TAlgebra>(rSurfaceFct);

	#ifdef UG_PARALLEL
	m_ParallelStorageType = rSurfaceFct.get_storage_mask();
	#endif

//	\todo: At this point the vert. slaves are not provided with the correct values
//		   Thus, the vert. masters should sent their values there.
}


template <typename TDomain>
template <typename TElem, typename TAlgebra>
void AdaptionSurfaceGridFunction<TDomain>::
copy_to_surface(GridFunction<TDomain,TAlgebra>& rSurfaceFct, TElem* elem)
{
	std::vector<DoFIndex> vInd;
	const std::vector<std::vector<number> >& vvVal = m_aaValue[elem];

	UG_ASSERT(vvVal.size() == m_spDDInfo->num_fct(), "Array says numFct: "<<
	          vvVal.size()<<", but should be "<<m_spDDInfo->num_fct()<<" on "
	          << ElementDebugInfo(*m_spDomain->grid(), elem));

	for(size_t fct = 0; fct < vvVal.size(); ++fct){

		const std::vector<number>& vVal = vvVal[fct];
		rSurfaceFct.inner_dof_indices(elem, fct, vInd);

		UG_ASSERT(vVal.size() == vInd.size(), "Stored dofs are "<<vVal.size()<<
		          ", but fct "<<fct<<" has "<<vInd.size()<<" dofs on "<<
		          elem->reference_object_id())

		for(size_t i = 0; i < vInd.size(); ++i){
			 UG_ASSERT(vInd[i][0] != size_t(-1),
			 		   "Bad dof-index encountered on: "
			 		   << ElementDebugInfo(*rSurfaceFct.domain()->grid(), elem));
			 DoFRef(rSurfaceFct, vInd[i]) = vVal[i];
		}
	}
}

template <typename TDomain>
template <typename TElem, typename TAlgebra>
void AdaptionSurfaceGridFunction<TDomain>::
copy_to_surface(GridFunction<TDomain,TAlgebra>& rSurfaceFct)
{
	ConstSmartPtr<SurfaceView> spSurfView = rSurfaceFct.approx_space()->surface_view();
	ConstSmartPtr<MultiGrid> spGrid = m_spDomain->grid();

	typedef typename GridFunction<TDomain,TAlgebra>::template traits<TElem>::const_iterator iter_type;
	iter_type iter = rSurfaceFct.template begin<TElem>(SurfaceView::ALL);
	iter_type iterEnd = rSurfaceFct.template end<TElem>(SurfaceView::ALL);

	for( ; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

		copy_to_surface(rSurfaceFct, elem);
	}
}

template <typename TDomain>
template <typename TAlgebra>
void AdaptionSurfaceGridFunction<TDomain>::
copy_to_surface(GridFunction<TDomain,TAlgebra>& rSurfaceFct)
{
	GFUNCADAPT_PROFILE_FUNC();
	if(rSurfaceFct.max_dofs(VERTEX))copy_to_surface<Vertex,TAlgebra>(rSurfaceFct);
	if(rSurfaceFct.max_dofs(EDGE)) 	copy_to_surface<Edge,TAlgebra>(rSurfaceFct);
	if(rSurfaceFct.max_dofs(FACE))	copy_to_surface<Face,TAlgebra>(rSurfaceFct);
	if(rSurfaceFct.max_dofs(VOLUME))copy_to_surface<Volume,TAlgebra>(rSurfaceFct);

	#ifdef UG_PARALLEL
	rSurfaceFct.set_storage_type(m_ParallelStorageType);
	#endif

	detach_entries();
}


template <typename TDomain>
template <typename TBaseElem>
inline void AdaptionSurfaceGridFunction<TDomain>::obj_created(TBaseElem* elem){

//	 get value attachment
	std::vector<std::vector<number> >& vvVal = m_aaValue[elem];

//	resize to number of functions
	vvVal.resize(m_spDDInfo->num_fct());

//	get geom obj infos
	const ReferenceObjectID roid = elem->reference_object_id();
	const int si = m_spDomain->subset_handler()->get_subset_index(elem);

//	for each fct cmp resize to number of dofs associated with elem
	for(size_t fct = 0; fct < vvVal.size(); ++fct){
		vvVal[fct].resize(m_spDDInfo->num_fct_dofs(fct, roid, si));
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__ADAPTION_SURFACE_GRID_FUNCTION_IMPL__ */
