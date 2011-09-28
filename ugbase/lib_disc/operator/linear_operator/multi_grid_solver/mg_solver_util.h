/*
 * mg_solver_util.h
 *
 *  Created on: 01.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_UTIL__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__MG_SOLVER_UTIL__

// extern header
#include <vector>

// library intern headers
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/distributed_grid.h"
#endif
#include "lib_algebra/operator/operator_iterator_interface.h"
#include "lib_algebra/operator/operator_inverse_interface.h"
#include "lib_algebra/operator/operator_interface.h"


namespace ug{

////////////////////////////////////////////////////////////////////////////////
// SurfaceToToplevelMap
////////////////////////////////////////////////////////////////////////////////

/// creates a mapping levIndex = vMap[surfIndex];
template <typename TElem, typename TDoFDistribution>
bool CreateSurfaceToToplevelMap(std::vector<size_t>& vMap,
                                const IDoFDistribution<TDoFDistribution>& surfDD,
                                const IDoFDistribution<TDoFDistribution>& topDD)
{
//	type of element iterator
	typedef typename geometry_traits<TElem>::const_iterator iter_type;

//	vector of indices
	typedef typename IDoFDistribution<TDoFDistribution>::algebra_index_vector_type ind_vec_type;
	ind_vec_type surfaceInd, levelInd;

	for(int si = 0; si < surfDD.num_subsets(); ++si)
	{
	//	iterators for subset
		iter_type iter = surfDD.template begin<TElem>(si);
		iter_type iterEnd = surfDD.template end<TElem>(si);

	//	loop all elements of type
		for( ; iter != iterEnd; ++iter)
		{
		//	get elem
			TElem* elem = *iter;

		//	extract all algebra indices for the element on surface
			surfDD.inner_algebra_indices(elem, surfaceInd);

		//	extract all algebra indices for the element on level
			topDD.inner_algebra_indices(elem, levelInd);

		//	check that index sets have same cardinality
			UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

		//	copy all elements of the vector
			for(size_t i = 0; i < surfaceInd.size(); ++i)
			{
			//	copy entries into level vector
				vMap[surfaceInd[i]] = levelInd[i];
			}
		}
	}

//	we're done
	return true;
}

template <typename TDoFDistribution>
bool CreateSurfaceToToplevelMap(std::vector<size_t>& vMap,
                                const IDoFDistribution<TDoFDistribution>& surfDD,
                                const IDoFDistribution<TDoFDistribution>& topDD)
{
//	check full refinement
	if(surfDD.num_indices() != topDD.num_indices())
	{
		UG_LOG("ERROR in 'CreateSurfaceToToplevelMap': This function can only"
				" be applied to a full refined grid, where the surface is the "
				" top level.\n");
		return false;
	}

//	success flag
	bool bRetVal = true;

//	resize mapping
	vMap.resize(surfDD.num_indices(), 10000555);

// 	add dofs on elements
	if(surfDD.has_indices_on(VERTEX))
		bRetVal &= CreateSurfaceToToplevelMap<VertexBase, TDoFDistribution>(vMap, surfDD, topDD);
	if(surfDD.has_indices_on(EDGE))
		bRetVal &= CreateSurfaceToToplevelMap<EdgeBase, TDoFDistribution>(vMap, surfDD, topDD);
	if(surfDD.has_indices_on(FACE))
		bRetVal &= CreateSurfaceToToplevelMap<Face, TDoFDistribution>(vMap, surfDD, topDD);
	if(surfDD.has_indices_on(VOLUME))
		bRetVal &= CreateSurfaceToToplevelMap<Volume, TDoFDistribution>(vMap, surfDD, topDD);

	return bRetVal;
}

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
template <typename TElem, typename TVector, typename TDoFDistributionImpl>
bool ProjectSurfaceToLevel(const std::vector<TVector*>& vLevelVector,
                           const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD,
                           const TVector& surfaceVector,
                           const IDoFDistribution<TDoFDistributionImpl>& surfaceDD,
                           const SurfaceView& surfaceView)
{
//	type of element iterator
	typedef typename geometry_traits<TElem>::const_iterator iter_type;

//	vector of indices
	typedef typename IDoFDistribution<TDoFDistributionImpl>::algebra_index_vector_type ind_vec_type;
	ind_vec_type surfaceInd, levelInd;

//	loop all elements of type
	for(int si = 0; si < surfaceDD.num_subsets(); ++si)
	{
	//	iterators for subset
		iter_type iter = surfaceDD.template begin<TElem>(si);
		iter_type iterEnd = surfaceDD.template end<TElem>(si);

		for( ; iter != iterEnd; ++iter)
		{
		//	get elem
			TElem* elem = *iter;

		//	extract all algebra indices for the element on surface
			surfaceDD.inner_algebra_indices(elem, surfaceInd);

		//	get level of element in hierarchy
			int level = surfaceView.get_level(elem);

		//	get corresponding level vector for element
			UG_ASSERT(vLevelVector[level] != NULL, "Vector missing");
			TVector& levelVector = *(vLevelVector[level]);

		//	check that level is correct
			UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, that is not passed.");

		//	extract all algebra indices for the element on level
			UG_ASSERT(vLevelDD[level] != NULL, "DoF Distribution missing");
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

//	we're done
	return true;
}

/// projects surface function to level functions
template <typename TVector, typename TDoFDistributionImpl>
bool ProjectSurfaceToLevel(const std::vector<TVector*>& vLevelVector,
                           const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD,
                           const TVector& surfVector,
                           const IDoFDistribution<TDoFDistributionImpl>& surfDD,
                           const SurfaceView& surfView)
{
//	check, that levelFuntions and level DoFDistributions are the same number
	if(vLevelVector.size() != vLevelDD.size())
	{
		UG_LOG("In ProjectSurfaceToLevel: Number of level Vectors ("
				<< vLevelVector.size() <<") and level DoF Distributions ("
				<< vLevelDD.size() << ") does not match. Aborting.\n");
		return false;
	}

//	return flag
	bool bRet = true;

//	forward for all BaseObject types
	if(surfDD.has_indices_on(VERTEX))
		bRet &= ProjectSurfaceToLevel<VertexBase, TVector, TDoFDistributionImpl>
					(vLevelVector, vLevelDD, surfVector, surfDD, surfView);
	if(surfDD.has_indices_on(EDGE))
		bRet &= ProjectSurfaceToLevel<EdgeBase, TVector, TDoFDistributionImpl>
					(vLevelVector, vLevelDD, surfVector, surfDD, surfView);
	if(surfDD.has_indices_on(FACE))
		bRet &= ProjectSurfaceToLevel<Face, TVector, TDoFDistributionImpl>
					(vLevelVector, vLevelDD, surfVector, surfDD, surfView);
	if(surfDD.has_indices_on(VOLUME))
		bRet &= ProjectSurfaceToLevel<Volume, TVector, TDoFDistributionImpl>
					(vLevelVector, vLevelDD, surfVector, surfDD, surfView);

#ifdef UG_PARALLEL
//	copy storage type into all vectors
	for(size_t lev = 0; lev < vLevelVector.size(); ++lev)
		if(vLevelVector[lev] != NULL)
			vLevelVector[lev]->copy_storage_type(surfVector);
#endif

//	we're done
	return bRet;
}

/// projects surface function to level functions
template <typename TElem, typename TVector, typename TDoFDistributionImpl>
bool ProjectLevelToSurface(TVector& surfaceVector,
                           const IDoFDistribution<TDoFDistributionImpl>& surfaceDD,
                           const SurfaceView& surfaceView,
						   const std::vector<const TVector*>& vLevelVector,
                           const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD,
                           const int baseLevel)
{
//	type of element iterator
	typedef typename geometry_traits<TElem>::const_iterator iter_type;

//	vector of indices
	typedef typename IDoFDistribution<TDoFDistributionImpl>::algebra_index_vector_type ind_vec_type;
	ind_vec_type surfaceInd, levelInd;

//	loop all elements of type
	for(int si = 0; si < surfaceDD.num_subsets(); ++si)
	{
	//	iterators for subset
		iter_type iter = surfaceDD.template begin<TElem>(si);
		iter_type iterEnd = surfaceDD.template end<TElem>(si);

		for( ; iter != iterEnd; ++iter)
		{
		//	get elem
			TElem* elem = *iter;

		//	skip shadows
			if(surfaceView.is_shadow(elem)) continue;

		//	extract all algebra indices for the element on surface
			surfaceDD.inner_algebra_indices(elem, surfaceInd);

		//	get level of element in hierarchy
			int level = surfaceView.get_level(elem);

		//	get corresponding level vector for element
			UG_ASSERT(vLevelVector[level] != NULL, "Vector missing");
			const TVector& levelVector = *(vLevelVector[level]);

		//	check that level is correct
			UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, that is not passed.");

		//	extract all algebra indices for the element on level
			UG_ASSERT(vLevelDD[level] != NULL, "DoF Distribution missing");
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

//	we're done
	return true;
}

/// projects surface function to level functions
template <typename TVector, typename TDoFDistributionImpl>
bool ProjectLevelToSurface(TVector& surfVector,
                           const IDoFDistribution<TDoFDistributionImpl>& surfDD,
                           const SurfaceView& surfView,
                           const std::vector<const TVector*>& vLevelVector,
                           const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD,
                           const int baseLevel)
{
//	check, that levelFuntions and level DoFDistributions are the same number
	if(vLevelVector.size() != vLevelDD.size())
	{
		UG_LOG("In ProjectLevelToSurface: Number of level Vectors ("
				<< vLevelVector.size() <<") and level DoF Distributions ("
				<< vLevelDD.size() << ") does not match. Aborting.\n");
		return false;
	}

//	return flag
	bool bRet = true;

//	forward for all BaseObject types
	if(surfDD.has_indices_on(VERTEX))
		bRet &= ProjectLevelToSurface<VertexBase, TVector, TDoFDistributionImpl>
					(surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel);
	if(surfDD.has_indices_on(EDGE))
		bRet &= ProjectLevelToSurface<EdgeBase, TVector, TDoFDistributionImpl>
					(surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel);
	if(surfDD.has_indices_on(FACE))
		bRet &= ProjectLevelToSurface<Face, TVector, TDoFDistributionImpl>
					(surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel);
	if(surfDD.has_indices_on(VOLUME))
		bRet &= ProjectLevelToSurface<Volume, TVector, TDoFDistributionImpl>
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
			return false;
		}

	//	set type of surface vector to common base
		surfVector.set_storage_type(type);
	}
#endif

//	we're done
	return bRet;
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
template <typename TBaseElem, typename TVector, typename TDoFDistributionImpl>
bool AddProjectionOfShadows(TVector& fineVec, const TVector& coarseVec,
                           const number scale,
                           const IDoFDistribution<TDoFDistributionImpl>& ddFine,
                           const IDoFDistribution<TDoFDistributionImpl>& ddCoarse,
                           const SurfaceView& surfView)
{

	typename TDoFDistributionImpl::algebra_index_vector_type fineInd;
	typename TDoFDistributionImpl::algebra_index_vector_type coarseInd;

// 	iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < ddCoarse.num_subsets(); ++si)
	{
		iter = ddCoarse.template begin<TBaseElem>(si);
		iterEnd = ddCoarse.template end<TBaseElem>(si);

	// 	loop elements of coarse subset
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			TBaseElem* pElem = *iter;

		//	skip non-shadows
			if(!surfView.is_shadow(pElem)) continue;

		// 	get child (i.e. shadow)
			TBaseElem* pShadowing = surfView.get_shadow_child(pElem);
			UG_ASSERT(pShadowing != NULL, "Shadow child does not exist. Error.");

		// 	get global indices
			ddCoarse.inner_algebra_indices(pElem, coarseInd);

		// 	get global indices
			ddFine.inner_algebra_indices(pShadowing, fineInd);

		//	add coarse vector entries to fine vector entries
			for(size_t i = 0; i < coarseInd.size(); ++i)
			{
				VecScaleAdd(fineVec[fineInd[i]],
				            1.0, fineVec[fineInd[i]],
				            scale, coarseVec[coarseInd[i]]);
			}
		}
	}

//	we're done
	return true;
}

template <typename TVector, typename TDoFDistributionImpl>
bool AddProjectionOfShadows(TVector& fineVec, const TVector& coarseVec,
                           const number scale,
                           const IDoFDistribution<TDoFDistributionImpl>& ddFine,
                           const IDoFDistribution<TDoFDistributionImpl>& ddCoarse,
                           const SurfaceView& surfView)
{
//	return flag
	bool bRet = true;

//	forward for all BaseObject types
	if(ddFine.has_indices_on(VERTEX))
		bRet &= AddProjectionOfShadows<VertexBase, TVector, TDoFDistributionImpl>
					(fineVec, coarseVec, scale, ddFine, ddCoarse, surfView);
	if(ddFine.has_indices_on(EDGE))
		bRet &= AddProjectionOfShadows<EdgeBase, TVector, TDoFDistributionImpl>
					(fineVec, coarseVec, scale, ddFine, ddCoarse, surfView);
	if(ddFine.has_indices_on(FACE))
		bRet &= AddProjectionOfShadows<Face, TVector, TDoFDistributionImpl>
					(fineVec, coarseVec, scale, ddFine, ddCoarse, surfView);
	if(ddFine.has_indices_on(VOLUME))
		bRet &= AddProjectionOfShadows<Volume, TVector, TDoFDistributionImpl>
					(fineVec, coarseVec, scale, ddFine, ddCoarse, surfView);

//	return success
	return bRet;
}


/**
 * This functions sets the values of a vector to zero, where the index
 * corresponds to a refine-patch boundary (i.e. the geomeric object is a
 * shadowing object) for an element type
 *
 * \param[out]	vec				grid vector
 * \param[in] 	dd				DoFDistribution
 * \param[in]	surfView		SurfaceView
 */
template <typename TBaseElem, typename TVector, typename TDoFDistributionImpl>
bool SetZeroOnShadowing(TVector& vec,
                        const IDoFDistribution<TDoFDistributionImpl>& dd,
                        const SurfaceView& surfView)
{
//	indices
	typename TDoFDistributionImpl::algebra_index_vector_type ind;

// 	Vertex iterators
	typename geometry_traits<TBaseElem>::const_iterator iter, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < dd.num_subsets(); ++si)
	{
		iter = dd.template begin<TBaseElem>(si);
		iterEnd = dd.template end<TBaseElem>(si);

	// 	loop nodes of fine subset
		for(; iter != iterEnd; ++iter)
		{
		//	get vertex
			TBaseElem* vrt = *iter;

		//	skip non-shadowing vertices
			if(!surfView.shadows(vrt)) continue;

		// 	get global indices
			dd.inner_algebra_indices(vrt, ind);

		//	set vector entries to zero
			for(size_t i = 0; i < ind.size(); ++i)
			{
				vec[ind[i]] = 0.0;
			}
		}
	}

//	we're done
	return true;
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
template <typename TVector, typename TDoFDistributionImpl>
bool SetZeroOnShadowing(TVector& vec,
                        const IDoFDistribution<TDoFDistributionImpl>& dd,
                        const SurfaceView& surfView)
{
//	return flag
	bool bRet = true;

//	forward for all BaseObject types
	if(dd.has_indices_on(VERTEX))
		bRet &= SetZeroOnShadowing<VertexBase, TVector, TDoFDistributionImpl>
					(vec, dd, surfView);
	if(dd.has_indices_on(EDGE))
		bRet &= SetZeroOnShadowing<EdgeBase, TVector, TDoFDistributionImpl>
					(vec, dd, surfView);
	if(dd.has_indices_on(FACE))
		bRet &= SetZeroOnShadowing<Face, TVector, TDoFDistributionImpl>
					(vec, dd, surfView);
	if(dd.has_indices_on(VOLUME))
		bRet &= SetZeroOnShadowing<Volume, TVector, TDoFDistributionImpl>
					(vec, dd, surfView);

//	return success
	return bRet;
}

////////////////////////////////////////////////////////////////////////////////
// Selections
////////////////////////////////////////////////////////////////////////////////

/// selects all non-shadows, that are adjacent to a shadow in the multigrid
bool SelectNonShadowsAdjacentToShadows(ISelector& sel, const SurfaceView& surfView);

/// selects all non-shadows, that are adjacent to a shadow on a grid levels
bool SelectNonShadowsAdjacentToShadowsOnLevel(ISelector& sel,
                                              const SurfaceView& surfView,
                                              int level);

#ifdef UG_PARALLEL
template <typename TElemBase>
void SelectNonGhosts(ISelector& sel,
                     DistributedGridManager& dstGrMgr,
                     typename geometry_traits<TElemBase>::iterator iter,
                     typename geometry_traits<TElemBase>::iterator iterEnd)
{
//	loop all base elems
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TElemBase* elem = *iter;

	//	select ghosts
		if(!dstGrMgr.is_ghost(elem)) sel.select(elem);
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
	smallMat.copy_storage_type(origMat);
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
//	check size
	UG_ASSERT(vMap.size() == newMat.num_rows(), "Size must match.");
	UG_ASSERT(vMap.size() == newMat.num_cols(), "Size must match.");
	UG_ASSERT(vMap.size() == origMat.num_rows(), "Size must match.");
	UG_ASSERT(vMap.size() == origMat.num_cols(), "Size must match.");

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
		//	get corresponding level connection index
			const size_t newConnIndex = vMap[conn.index()];

		//	copy connection to level matrix
			newMat(newInd, newConnIndex) = conn.value();
		}
	}

#ifdef UG_PARALLEL
	newMat.copy_storage_type(origMat);
#endif
}


} // end namespace ug
#endif
