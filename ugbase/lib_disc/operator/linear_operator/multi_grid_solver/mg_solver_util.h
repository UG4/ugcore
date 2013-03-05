/*
 * mg_solver_util.h
 *
 *  Created on: 01.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER_UTIL__
#define __H__UG__LIB_DISC__MULTI_GRID_SOLVER__MG_SOLVER_UTIL__

// extern header
#include <vector>

// library intern headers
#include "common/common.h"
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
#endif
#include "lib_algebra/operator/interface/operator_iterator.h"
#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/operator.h"

#include "lib_disc/dof_manager/mg_dof_distribution.h"

//only for debugging!!!
#include "lib_grid/algorithms/debug_util.h"


namespace ug{

////////////////////////////////////////////////////////////////////////////////
// SurfaceToToplevelMap
////////////////////////////////////////////////////////////////////////////////

void CreateSurfaceToToplevelMap(std::vector<size_t>& vMap,
                                ConstSmartPtr<DoFDistribution> surfDD,
                                ConstSmartPtr<DoFDistribution> topDD);

////////////////////////////////////////////////////////////////////////////////
// Projections
////////////////////////////////////////////////////////////////////////////////

/// projects surface function to level functions
/**
 * This function copies the values of a vector defined for the surface grid
 * to the vectors corresponding to the level grids. It loops all elements of the
 * surface grid and copies the values for all DoFs on the element.
 *
 * \param[in]	si				Subset index
 * \param[out]	vLevelVector	level vectors
 * \param[in]	vLevelDD		DoF distribution on levels
 * \param[in]	surfaceVector	surface vector
 * \param[in]	surfaceDD		DoF distribution on surface
 * \param[in]	surfaceView		Surface view
 *
 * \tparam		TElem					element type to process
 * \tparam		TVector					vector type
 * \tparam		TDoFDistributionImpl	type of DoF distribution
 */
template <typename TElem, typename TVector>
void ProjectSurfaceToLevel(const std::vector<TVector*>& vLevelVector,
                           std::vector<ConstSmartPtr<DoFDistribution> > vLevelDD,
                           const TVector& surfaceVector,
                           ConstSmartPtr<DoFDistribution> surfaceDD,
                           const SurfaceView& surfaceView,
                           const int baseLvl = 0)
{
	PROFILE_FUNC_GROUP("gmg");
//	type of element iterator
	typedef typename DoFDistribution::traits<TElem>::const_iterator iter_type;

//	vector of indices
	std::vector<size_t> surfaceInd, levelInd;

//	loop all elements of type
	for(int si = 0; si < surfaceDD->num_subsets(); ++si)
	{
	//	iterators for subset
		iter_type iter = surfaceDD->begin<TElem>(si);
		iter_type iterEnd = surfaceDD->end<TElem>(si);

		for( ; iter != iterEnd; ++iter)
		{
		//	get elem
			TElem* elem = *iter;

		//	get level of element in hierarchy
			int level = surfaceView.get_level(elem);

		//	make sure that we're not below baseLvl
			if(level < baseLvl)
				continue;

		//	get corresponding level vector for element
			UG_ASSERT(vLevelVector[level] != NULL, "Vector missing on level " << level);
			TVector& levelVector = *(vLevelVector[level]);

		//	extract all algebra indices for the element on surface
			surfaceDD->inner_algebra_indices(elem, surfaceInd);

		//	check that level is correct
			UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, that is not passed.");

		//	extract all algebra indices for the element on level
			UG_ASSERT(vLevelDD[level].valid(), "DoF Distribution missing");
			vLevelDD[level]->inner_algebra_indices(elem, levelInd);

		//	check that index sets have same cardinality
			UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

		//	copy all elements of the vector
			for(size_t i = 0; i < surfaceInd.size(); ++i)
			{
			//	copy entries into level vector
				levelVector[levelInd[i]] = surfaceVector[surfaceInd[i]];
			}
		}
	}
}

/// projects surface function to level functions
template <typename TVector>
void ProjectSurfaceToLevel(const std::vector<TVector*>& vLevelVector,
                           std::vector<ConstSmartPtr<DoFDistribution> > vLevelDD,
                           const TVector& surfVector,
                           ConstSmartPtr<DoFDistribution> surfDD,
                           const SurfaceView& surfView,
                           const int baseLvl = 0)
{
	PROFILE_FUNC_GROUP("gmg");
//	check, that levelFuntions and level DoFDistributions are the same number
	if(vLevelVector.size() > vLevelDD.size())
		UG_THROW("ProjectSurfaceToLevel: Number of level Vectors ("
				<< vLevelVector.size() <<") greater that level DoF Distributions ("
				<< vLevelDD.size() << "). Aborting.\n");

//	forward for all BaseObject types
	if(surfDD->max_dofs(VERTEX))
		ProjectSurfaceToLevel<VertexBase, TVector>
					(vLevelVector, vLevelDD, surfVector, surfDD, surfView, baseLvl);
	if(surfDD->max_dofs(EDGE))
		ProjectSurfaceToLevel<EdgeBase, TVector>
					(vLevelVector, vLevelDD, surfVector, surfDD, surfView, baseLvl);
	if(surfDD->max_dofs(FACE))
		ProjectSurfaceToLevel<Face, TVector>
					(vLevelVector, vLevelDD, surfVector, surfDD, surfView, baseLvl);
	if(surfDD->max_dofs(VOLUME))
		ProjectSurfaceToLevel<Volume, TVector>
					(vLevelVector, vLevelDD, surfVector, surfDD, surfView, baseLvl);

#ifdef UG_PARALLEL
//	copy storage type into all vectors
	for(size_t lev = 0; lev < vLevelVector.size(); ++lev)
		if(vLevelVector[lev] != NULL)
			vLevelVector[lev]->set_storage_type(surfVector.get_storage_mask());
#endif
}

/// projects surface function to level functions
template <typename TElem, typename TVector>
void ProjectLevelToSurface(TVector& surfaceVector,
                           ConstSmartPtr<DoFDistribution> surfaceDD,
                           const SurfaceView& surfaceView,
						   const std::vector<const TVector*>& vLevelVector,
						   std::vector<ConstSmartPtr<DoFDistribution> > vLevelDD,
                           const int baseLevel = 0)
{
	PROFILE_FUNC_GROUP("gmg");
//	type of element iterator
	typedef typename DoFDistribution::traits<TElem>::const_iterator iter_type;

//	vector of indices
	std::vector<size_t> surfaceInd, levelInd;

//	loop all elements of type
	for(int si = 0; si < surfaceDD->num_subsets(); ++si)
	{
	//	iterators for subset
		iter_type iter = surfaceDD->begin<TElem>(si);
		iter_type iterEnd = surfaceDD->end<TElem>(si);

		for( ; iter != iterEnd; ++iter)
		{
		//	get elem
			TElem* elem = *iter;

		//	skip shadows
			if(surfaceView.is_shadowed(elem)) continue;

		//	get level of element in hierarchy
			int level = surfaceView.get_level(elem);

		//	make sure that we're not below base-level
			if(level < baseLevel)
				continue;

		//	extract all algebra indices for the element on surface
			surfaceDD->inner_algebra_indices(elem, surfaceInd);

		//	get corresponding level vector for element
			const TVector& levelVector = *(vLevelVector[level]);

		//	check that level is correct
			UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, that is not passed.");

		//	extract all algebra indices for the element on level
			UG_ASSERT(vLevelDD[level].valid(), "DoF Distribution missing");
			vLevelDD[level]->inner_algebra_indices(elem, levelInd);

		//	check that index sets have same cardinality
			UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

		//	copy all elements of the vector
			for(size_t i = 0; i < surfaceInd.size(); ++i)
			{
			//	copy entries into level vector
				surfaceVector[surfaceInd[i]] = levelVector[levelInd[i]];
			}
		}
	}
}

/// projects surface function to level functions
template <typename TVector>
void ProjectLevelToSurface(TVector& surfVector,
                           ConstSmartPtr<DoFDistribution> surfDD,
                           const SurfaceView& surfView,
                           const std::vector<const TVector*>& vLevelVector,
                           std::vector<ConstSmartPtr<DoFDistribution> > vLevelDD,
                           const int baseLevel = 0)
{
	PROFILE_FUNC_GROUP("gmg");
//	check, that levelFuntions and level DoFDistributions are the same number
	if(vLevelVector.size() > vLevelDD.size())
		UG_THROW("ProjectLevelToSurface: Number of level Vectors ("
				<< vLevelVector.size() <<") greater than DoF Distributions ("
				<< vLevelDD.size() << "). Aborting.\n");

//	forward for all BaseObject types
	if(surfDD->max_dofs(VERTEX))
		ProjectLevelToSurface<VertexBase, TVector>
					(surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel);
	if(surfDD->max_dofs(EDGE))
		ProjectLevelToSurface<EdgeBase, TVector>
					(surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel);
	if(surfDD->max_dofs(FACE))
		ProjectLevelToSurface<Face, TVector>
					(surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel);
	if(surfDD->max_dofs(VOLUME))
		ProjectLevelToSurface<Volume, TVector>
					(surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel);

#ifdef UG_PARALLEL
//	copy storage type into surf vector
	if(vLevelVector.size() > 0)
	{
	//	flag if at least on vector given
		bool bFound = false;

		uint type = PST_UNDEFINED;
	//	find first valid vector and get its type
		for(size_t lev = baseLevel; lev < vLevelVector.size(); ++lev)
			if(vLevelVector[lev] != NULL && vLevelVector[lev]->size() > 0)
			{
				type = vLevelVector[lev]->get_storage_mask();
				bFound = true;
			}

	//	get intersection of types
		for(size_t lev = baseLevel; lev < vLevelVector.size(); ++lev)
			if(vLevelVector[lev] != NULL && vLevelVector[lev]->size() > 0)
				type = type & vLevelVector[lev]->get_storage_mask();

	//	check if union is defined
		if(type == PST_UNDEFINED && bFound == true)
		{
			UG_LOG("ERROR in 'ProjectLevelToSurface': storage type of level"
					" vectors is not ok. Must have at least on common type."
					" (e.g. additive/unique for all, or consistent for all)\n"
					" Types in levels are:\n");
 			for(size_t lev = baseLevel; lev < vLevelVector.size(); ++lev)
				if(vLevelVector[lev] != NULL)
					UG_LOG("  lev "<<lev<<": "<<vLevelVector[lev]->get_storage_mask()<<"\n");
 			UG_THROW("Cannot Project Level to Surface.")
		}

	//	set type of surface vector to common base
		surfVector.set_storage_type(type);
	}
#endif
}


////////////////////////////////////////////////////////////////////////////////
// Operation on Shadows/Shadowing
////////////////////////////////////////////////////////////////////////////////

/**
 * This functions adds the shadow values from a coarser grid to the shadowing
 * DoFs on the finer grid.
 *
 * \param[out]	fineVec			fine grid vector
 * \param[out]	coarseVec		coarse grid vector
 * \param[in]	scale			scaling, when adding
 * \param[in] 	ddFine			dof distribution on fine space
 * \param[in] 	ddCoarse		dof distribution on coarse space
 * \param[in]	surfView		surface view
 */
template <typename TBaseElem, typename TVector>
void AddProjectionOfShadows(const std::vector<TVector*>& vFineVector,
                            std::vector<ConstSmartPtr<DoFDistribution> > vDDFine,
                            const TVector& coarseVec,
                            ConstSmartPtr<DoFDistribution> ddCoarse,
                            const int level,
                            const number scale,
                            const SurfaceView& surfView)
{
	PROFILE_FUNC_GROUP("gmg");
	std::vector<size_t> fineInd, coarseInd;

// 	iterators
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator const_iterator;
	const_iterator iter, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < ddCoarse->num_subsets(); ++si)
	{
		iter = ddCoarse->begin<TBaseElem>(si);
		iterEnd = ddCoarse->end<TBaseElem>(si);

	// 	loop elements of coarse subset
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			TBaseElem* pElem = *iter;

		//	skip non-shadows
			if(!surfView.is_shadowed(pElem)) continue;

		// 	get child (i.e. shadow)
			TBaseElem* pShadowing = surfView.child_if_copy(pElem);

		//	offset to count which child currently handling
			int offset = 0;

		// 	get global indices
			ddCoarse->inner_algebra_indices(pElem, coarseInd);

		//	skip if not a copy
			while(pShadowing){
			//	increase offset
				++offset;

			// 	get global indices
				vDDFine[level+offset]->inner_algebra_indices(pShadowing, fineInd);

			//	add coarse vector entries to fine vector entries
				for(size_t i = 0; i < coarseInd.size(); ++i)
				{
					VecScaleAdd((*vFineVector[level+offset])[fineInd[i]],
								1.0, (*vFineVector[level+offset])[fineInd[i]],
								scale, coarseVec[coarseInd[i]]);
				}

			//	next child
				pElem = pShadowing;
				pShadowing = surfView.child_if_copy(pElem);
			}
		}
	}
}

template <typename TVector>
void AddProjectionOfShadows(const std::vector<TVector*>& vFineVector,
                            std::vector<ConstSmartPtr<DoFDistribution> > vDDFine,
                            const TVector& coarseVec,
                            ConstSmartPtr<DoFDistribution> ddCoarse,
                            const int level,
                            const number scale,
                            const SurfaceView& surfView)
{
	PROFILE_FUNC_GROUP("gmg");
//	forward for all BaseObject types
	if(ddCoarse->max_dofs(VERTEX))
		AddProjectionOfShadows<VertexBase, TVector>
					(vFineVector, vDDFine, coarseVec, ddCoarse, level, scale, surfView);
	if(ddCoarse->max_dofs(EDGE))
		AddProjectionOfShadows<EdgeBase, TVector>
					(vFineVector, vDDFine, coarseVec, ddCoarse, level, scale, surfView);
	if(ddCoarse->max_dofs(FACE))
		AddProjectionOfShadows<Face, TVector>
					(vFineVector, vDDFine, coarseVec, ddCoarse, level, scale, surfView);
	if(ddCoarse->max_dofs(VOLUME))
		AddProjectionOfShadows<Volume, TVector>
					(vFineVector, vDDFine, coarseVec, ddCoarse, level, scale, surfView);
}


/**
 * This functions sets the values of a vector to zero, where the index
 * corresponds to a refine-patch boundary (i.e. the geomeric object is a
 * shadowing object) for an element type
 *
 * \param[out]	vec					grid vector
 * \param[in] 	dd					DoFDistribution
 * \param[in]	surfView			SurfaceView
 * \param[in]	pmapGlobalToPatch	(optional) mapping of global indices to patch indices
 */
template <typename TBaseElem, typename TVector>
void SetZeroOnShadowing(TVector& vec,
                        ConstSmartPtr<DoFDistribution> dd,
                        const SurfaceView& surfView,
                        const std::vector<int>* pmapGlobalToPatch = NULL)
{
	PROFILE_FUNC_GROUP("gmg");
//	indices
	std::vector<size_t> ind;

// 	Vertex iterators
	typedef typename DoFDistribution::traits<TBaseElem>::const_iterator const_iterator;
	const_iterator iter, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < dd->num_subsets(); ++si)
	{
		iter = dd->begin<TBaseElem>(si);
		iterEnd = dd->end<TBaseElem>(si);

	// 	loop nodes of fine subset
		if(pmapGlobalToPatch){
			const std::vector<int>& mapGlobalToPatch = *pmapGlobalToPatch;

			for(; iter != iterEnd; ++iter)
			{
			//	get vertex
				TBaseElem* vrt = *iter;

				if(!surfView.is_shadowing(vrt))
					continue;

			// 	get global indices
				dd->inner_algebra_indices(vrt, ind);

			//	set vector entries to zero
				for(size_t i = 0; i < ind.size(); ++i)
				{
					UG_ASSERT(ind[i] < mapGlobalToPatch.size(),
							 "mapGlobalToPatch is too small on level " << dd->grid_level() << "."
							 << "size: " << mapGlobalToPatch.size() << ", "
							 << "index: " << ind[i]
							 << ", at " << GetGeometricObjectCenter(
									 *const_cast<MultiGrid*>(surfView.subset_handler()->multi_grid()), vrt)
							 << ", on level " << surfView.subset_handler()->multi_grid()->get_level(vrt));
					UG_ASSERT((mapGlobalToPatch[ind[i]] >= 0) && (mapGlobalToPatch[ind[i]] < (int)vec.size()),
							  "Some problem with mapGlobalToPatch... probably trying to set a ghost to zero? "
							  << "mapGlobalToPatch[ind[i]] = " << mapGlobalToPatch[ind[i]]
							  << ", num patch indices " << vec.size()
							  << ", at " << GetGeometricObjectCenter(
									  *const_cast<MultiGrid*>(surfView.subset_handler()->multi_grid()), vrt)
							  << ", on level " << surfView.subset_handler()->multi_grid()->get_level(vrt));
					vec[mapGlobalToPatch[ind[i]]] = 0.0;
				}
			}
		}
		else{
			for(; iter != iterEnd; ++iter)
			{
			//	get vertex
				TBaseElem* vrt = *iter;

				if(!surfView.is_shadowing(vrt))
					continue;

			// 	get global indices
				dd->inner_algebra_indices(vrt, ind);

			//	set vector entries to zero
				for(size_t i = 0; i < ind.size(); ++i)				{
					vec[ind[i]] = 0.0;
				}
			}
		}
	}
}

/**
 * This functions sets the values of a vector to zero, where the index
 * corresponds to a refine-patch boundary (i.e. the geomeric object is a
 * shadowing object)
 *
 * \param[out]	vec				grid vector
 * \param[in] 	dd				DoFDistribution
 * \param[in]	surfView		SurfaceView
 */
template <typename TVector>
void SetZeroOnShadowing(TVector& vec,
                        ConstSmartPtr<DoFDistribution> dd,
                        const SurfaceView& surfView,
                        const std::vector<int>* pmapGlobalToPatch = NULL)
{
//	forward for all BaseObject types
	if(dd->max_dofs(VERTEX))
		SetZeroOnShadowing<VertexBase, TVector>(vec, dd, surfView, pmapGlobalToPatch);
	if(dd->max_dofs(EDGE))
		SetZeroOnShadowing<EdgeBase, TVector>(vec, dd, surfView, pmapGlobalToPatch);
	if(dd->max_dofs(FACE))
		SetZeroOnShadowing<Face, TVector>(vec, dd, surfView, pmapGlobalToPatch);
	if(dd->max_dofs(VOLUME))
		SetZeroOnShadowing<Volume, TVector>(vec, dd, surfView, pmapGlobalToPatch);
}

////////////////////////////////////////////////////////////////////////////////
// Selections
////////////////////////////////////////////////////////////////////////////////

/// selects all non-shadows, that are adjacent to a shadow in the multigrid
void SelectNonShadowsAdjacentToShadows(BoolMarker& sel, const SurfaceView& surfView);

/// selects all non-shadows, that are adjacent to a shadow on a grid levels
void SelectNonShadowsAdjacentToShadowsOnLevel(BoolMarker& sel,
                                              const SurfaceView& surfView,
                                              int level);

#ifdef UG_PARALLEL
template <typename TElemBase, typename TIter>
void SelectNonGhosts(BoolMarker& sel,
                     DistributedGridManager& dstGrMgr,
                     TIter iter,
                     TIter iterEnd)
{
//	loop all base elems
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TElemBase* elem = *iter;

	//	select ghosts
		if(!dstGrMgr.is_ghost(elem)) sel.mark(elem);
	}
}
#endif

////////////////////////////////////////////////////////////////////////////////
// Matrix Copy operations
////////////////////////////////////////////////////////////////////////////////

/// copies a matrix from a larger one into a smaller one
/**
 * This function copies a matrix of a larger index set into a matrix with a
 * smaller (or equally sized) index set. The copying is performed using a
 * mapping between the index set, that returns smallIndex = vMap[largeIndex],
 * and a -1 if the largeIndex is dropped.
 */
template <typename TMatrix>
void CopyMatrixByMapping(TMatrix& smallMat,
                         const std::vector<int>& vMap,
                         const TMatrix& origMat)
{
	PROFILE_FUNC_GROUP("gmg");
//	check size
	UG_ASSERT(vMap.size() == origMat.num_rows(), "Size must match.");
	UG_ASSERT(vMap.size() == origMat.num_cols(), "Size must match.");

//	type of matrix row iterator
	typedef typename TMatrix::const_row_iterator const_row_iterator;

//	loop all mapped indices
	for(size_t origInd = 0; origInd < vMap.size(); ++origInd)
	{
	//	get mapped level index
		const int smallInd = vMap[origInd];

	//	skipped non-mapped indices (indicated by -1)
		if(smallInd < 0) continue;

	//	loop all connections of the surface dof to other surface dofs and copy
	//	the matrix coupling into the level matrix

		for(const_row_iterator conn = origMat.begin_row(origInd);
								conn != origMat.end_row(origInd); ++conn)
		{
		//	get corresponding level connection index
			const int smallConnIndex = vMap[conn.index()];

		//	check that index is in small matrix, too
			if(smallConnIndex < 0) continue;

		//	copy connection to smaller matrix
			smallMat(smallInd, smallConnIndex) = conn.value();
		}
	}

#ifdef UG_PARALLEL
	smallMat.set_storage_type(origMat.get_storage_mask());
#endif
}

/// copies a matrix into a equally sized second one using a mapping
/**
 * This function copies a matrix of a index set into a matrix with a
 * equally sized index set. The copying is performed using a
 * mapping between the index set, that returns newIndex = vMap[origIndex].
 */
template <typename TMatrix>
void CopyMatrixByMapping(TMatrix& newMat,
                         const std::vector<size_t>& vMap,
                         const TMatrix& origMat)
{
	PROFILE_FUNC_GROUP("gmg");
//	check size
	UG_ASSERT(vMap.size() <= newMat.num_rows(), "Size must match. Map:"<<vMap.size()<<", mat:"<<newMat.num_rows());
	UG_ASSERT(vMap.size() <= newMat.num_cols(), "Size must match. Map:"<<vMap.size()<<", mat:"<<newMat.num_cols());
	UG_ASSERT(vMap.size() <= origMat.num_rows(), "Size must match. Map:"<<vMap.size()<<", mat:"<<origMat.num_rows());
	UG_ASSERT(vMap.size() <= origMat.num_cols(), "Size must match. Map:"<<vMap.size()<<", mat:"<<origMat.num_cols());

	newMat.resize(0,0);
	newMat.resize(origMat.num_rows(), origMat.num_cols());
	newMat.set(0.0);

//	type of matrix row iterator
	typedef typename TMatrix::const_row_iterator const_row_iterator;

//	loop all mapped indices
	for(size_t origInd = 0; origInd < vMap.size(); ++origInd)
	{
	//	get mapped level index
		const size_t newInd = vMap[origInd];

	//	loop all connections of the surface dof to other surface dofs and copy
	//	the matrix coupling into the level matrix

		for(const_row_iterator conn = origMat.begin_row(origInd);
									conn != origMat.end_row(origInd); ++conn)
		{
			size_t newConnIndex = conn.index();

		//	get corresponding level connection index
			if(conn.index() < vMap.size())
				newConnIndex = vMap[conn.index()];

		//	copy connection to level matrix
			newMat(newInd, newConnIndex) = conn.value();
		}
	}

//	loop remaining indices as identity
	for(size_t origInd = vMap.size(); origInd < origMat.num_rows(); ++origInd)
	{
	//	get mapped level index
		const size_t newInd = origInd;

	//	loop all connections of the surface dof to other surface dofs and copy
	//	the matrix coupling into the level matrix

		for(const_row_iterator conn = origMat.begin_row(origInd);
									conn != origMat.end_row(origInd); ++conn)
		{
			size_t newConnIndex = conn.index();

		//	get corresponding level connection index
			if(conn.index() < vMap.size())
				newConnIndex = vMap[conn.index()];

		//	copy connection to level matrix
			newMat(newInd, newConnIndex) = conn.value();
		}
	}

#ifdef UG_PARALLEL
	newMat.set_storage_type(origMat.get_storage_mask());
#endif
}


} // end namespace ug
#endif
