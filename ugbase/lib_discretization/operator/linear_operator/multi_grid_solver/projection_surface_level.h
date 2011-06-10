/*
 * projection_surface_level.h
 *
 *  Created on: 01.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__PROJECTION_SURFACE_LEVEL__
#define __H__LIB_DISCRETIZATION__MULTI_GRID_SOLVER__PROJECTION_SURFACE_LEVEL__

// extern header
#include <vector>

// library intern headers
#include "lib_discretization/function_spaces/grid_function.h"

namespace ug{

/// creates a list of vectors from list of grid_functions
/**
 * This function extracts an std::vector of vectors from a std::vector of
 * Grid Functions. The out vector is cleared, if requested.
 *
 * \param[in,out]	vVectorInOut	vector of algebra Vectors
 * \param[in]		vGridFunc		vector of GridFunctions
 * \param[in]		clearCont		flag, if vector should be clear before extracting
 */
template <typename TGridFunction>
void ExtractVectorsFromGridFunction(std::vector<typename TGridFunction::vector_type*>& vVectorInOut,
                                    const std::vector<TGridFunction*>& vGridFunc,
                                    bool clearCont = true)
{
//	type of vector
	typedef typename TGridFunction::vector_type vector_type;

//	clear container iff required
	if(clearCont) vVectorInOut.clear();

//	loop grid functions
	for(size_t i = 0; i < vGridFunc.size(); ++i)
	{
	//	cast Grid Function and add to vector
		vVectorInOut.push_back(dynamic_cast<vector_type*>(vGridFunc[i]));
	}
}

/// creates a list of const vectors from list of const grid_functions
/**
 * This function extracts an std::vector of vectors from a std::vector of
 * Grid Functions. The out vector is cleared, if requested.
 *
 * \param[in,out]	vVectorInOut	vector of algebra Vectors
 * \param[in]		vGridFunc		vector of GridFunctions
 * \param[in]		clearCont		flag, if vector should be clear before extracting
 */
template <typename TGridFunction>
void ExtractVectorsFromGridFunction(std::vector<const typename TGridFunction::vector_type*>& vVectorInOut,
                                    const std::vector<TGridFunction*>& vGridFunc,
                                    bool clearCont = true)
{
//	type of vector
	typedef typename TGridFunction::vector_type vector_type;

//	clear container iff required
	if(clearCont) vVectorInOut.clear();

//	loop grid functions
	for(size_t i = 0; i < vGridFunc.size(); ++i)
	{
	//	cast Grid Function and add to vector
		vVectorInOut.push_back(dynamic_cast<const vector_type*>(vGridFunc[i]));
	}
}

/// projects surface function to level functions
template <typename TElem, typename TVector, typename TDoFDistributionImpl>
bool ProjectSurfaceToLevel(int si,
						   const std::vector<TVector*>& vLevelVector,
                           const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD,
                           const TVector& surfaceVector,
                           const IDoFDistribution<TDoFDistributionImpl>& surfaceDD,
                           const SurfaceView& surfaceView)
{
//	type of element iterator
	typedef typename geometry_traits<TElem>::const_iterator iter_type;

//	iterators for subset
	iter_type iter = surfaceDD.template begin<TElem>(si);
	iter_type iterEnd = surfaceDD.template end<TElem>(si);

//	vector of indices
	typedef typename IDoFDistribution<TDoFDistributionImpl>::algebra_index_vector_type ind_vec_type;
	ind_vec_type surfaceInd, levelInd;

//	loop all elements of type
	for( ; iter != iterEnd; ++iter)
	{
	//	get elem
		TElem* elem = *iter;

	//	extract all algebra indices for the element on surface
		surfaceDD.get_algebra_indices(elem, surfaceInd);

	//	get level of element in hierarchy
		int level = surfaceView.get_level(elem);

	//	get corresponding level vector for element
		UG_ASSERT(vLevelVector[level] != NULL, "Vector missing");
		TVector& levelVector = *(vLevelVector[level]);

	//	check that level is correct
		UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, that is not passed.");

	//	extract all algebra indices for the element on level
		UG_ASSERT(vLevelDD[level] != NULL, "DoF Distribution missing");
		vLevelDD[level]->get_algebra_indices(elem, levelInd);

	//	check that index sets have same cardinality
		UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

	//	copy all elements of the vector
		for(size_t i = 0; i < surfaceInd.size(); ++i)
		{
		//	copy entries into level vector
			levelVector[levelInd[i]] = surfaceVector[surfaceInd[i]];
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
		UG_LOG("In ProjectSurfaceToLevel: Number of level Vectors (" << vLevelVector.size() <<
		       ") and level DoF Distributions (" << vLevelDD.size() << ") does"
		       " not match. Aborting.\n");
		return false;
	}

//	return flag
	bool bRet = true;

//	loop all subsets
	for(int si = 0; si < (int)surfDD.num_subsets(); ++si)
	{
	//	skip empty subsets (they have no dimension)
		if(	surfDD.template num<VertexBase>(si) == 0 &&
			surfDD.template num<EdgeBase>(si) == 0 &&
			surfDD.template num<Face>(si) == 0 &&
			surfDD.template num<Volume>(si) == 0 ) continue;

	//	switch dimension of subset
		switch(surfDD.dim_subset(si))
		{
			case 0:
				bRet &= ProjectSurfaceToLevel<VertexBase, TVector, TDoFDistributionImpl>
							(si, vLevelVector, vLevelDD, surfVector, surfDD, surfView);
				break;
			case 1:
				bRet &= ProjectSurfaceToLevel<Edge, TVector, TDoFDistributionImpl>
							(si, vLevelVector, vLevelDD, surfVector, surfDD, surfView);
				break;
			case 2:
				bRet &= ProjectSurfaceToLevel<Triangle, TVector, TDoFDistributionImpl>
							(si, vLevelVector, vLevelDD, surfVector, surfDD, surfView);
				bRet &= ProjectSurfaceToLevel<Quadrilateral, TVector, TDoFDistributionImpl>
							(si, vLevelVector, vLevelDD, surfVector, surfDD, surfView);
				break;
			case 3:
				bRet &= ProjectSurfaceToLevel<Tetrahedron, TVector, TDoFDistributionImpl>
							(si, vLevelVector, vLevelDD, surfVector, surfDD, surfView);
				bRet &= ProjectSurfaceToLevel<Pyramid, TVector, TDoFDistributionImpl>
							(si, vLevelVector, vLevelDD, surfVector, surfDD, surfView);
				bRet &= ProjectSurfaceToLevel<Prism, TVector, TDoFDistributionImpl>
							(si, vLevelVector, vLevelDD, surfVector, surfDD, surfView);
				bRet &= ProjectSurfaceToLevel<Hexahedron, TVector, TDoFDistributionImpl>
							(si, vLevelVector, vLevelDD, surfVector, surfDD, surfView);
				break;
			default:
				UG_LOG("ERROR in 'ProjectSurfaceToLevel': Dimension of subset ("
						<< surfDD.dim_subset(si) << ") not supported.\n");
				return false;
		}
	}

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
bool ProjectLevelToSurface(int si,
                           TVector& surfaceVector,
                           const IDoFDistribution<TDoFDistributionImpl>& surfaceDD,
                           const SurfaceView& surfaceView,
						   const std::vector<const TVector*>& vLevelVector,
                           const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD,
			   const int baseLevel) // 'baseLevel' added (13042011ih)
{
//	type of element iterator
	typedef typename geometry_traits<TElem>::const_iterator iter_type;

//	iterators for subset
	iter_type iter = surfaceDD.template begin<TElem>(si);
	iter_type iterEnd = surfaceDD.template end<TElem>(si);

//	vector of indices
	typedef typename IDoFDistribution<TDoFDistributionImpl>::algebra_index_vector_type ind_vec_type;
	ind_vec_type surfaceInd, levelInd;

//	loop all elements of type
	for( ; iter != iterEnd; ++iter)
	{
	//	get elem
		TElem* elem = *iter;

	//	extract all algebra indices for the element on surface
		surfaceDD.get_algebra_indices(elem, surfaceInd);

	//	get level of element in hierarchy
		int level = surfaceView.get_level(elem);

	//	get corresponding level vector for element
		UG_ASSERT(vLevelVector[level] != NULL, "Vector missing");
		const TVector& levelVector = *(vLevelVector[level]);

	//	check that level is correct
		UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, that is not passed.");

	//	extract all algebra indices for the element on level
		UG_ASSERT(vLevelDD[level] != NULL, "DoF Distribution missing");
		vLevelDD[level]->get_algebra_indices(elem, levelInd);

	//	check that index sets have same cardinality
		UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

	//	copy all elements of the vector
		for(size_t i = 0; i < surfaceInd.size(); ++i)
		{
		//	copy entries into level vector
			surfaceVector[surfaceInd[i]] = levelVector[levelInd[i]];
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
			   const int baseLevel) // 'baseLevel' added (13042011ih)
{
//	check, that levelFuntions and level DoFDistributions are the same number
	if(vLevelVector.size() != vLevelDD.size())
	{
		UG_LOG("In ProjectLevelToSurface: Number of level Vectors (" << vLevelVector.size() <<
		       ") and level DoF Distributions (" << vLevelDD.size() << ") does"
		       " not match. Aborting.\n");
		return false;
	}

//	return flag
	bool bRet = true;

//	loop all subsets
	for(int si = 0; si < (int)surfDD.num_subsets(); ++si)
	{
	//	skip empty subsets (they have no dimension)
		if(	surfDD.template num<VertexBase>(si) == 0 &&
			surfDD.template num<EdgeBase>(si) == 0 &&
			surfDD.template num<Face>(si) == 0 &&
			surfDD.template num<Volume>(si) == 0 ) continue;

	//	switch dimension of subset
		switch(surfDD.dim_subset(si))
		{
			case 0:
				bRet &= ProjectLevelToSurface<VertexBase, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel); // 'baseLevel' added (13042011ih)
				break;
			case 1:
				bRet &= ProjectLevelToSurface<Edge, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel); // 'baseLevel' added (13042011ih)
				break;
			case 2:
				bRet &= ProjectLevelToSurface<Triangle, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel); // 'baseLevel' added (13042011ih)
				bRet &= ProjectLevelToSurface<Quadrilateral, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel); // 'baseLevel' added (13042011ih)
				break;
			case 3:
				bRet &= ProjectLevelToSurface<Tetrahedron, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel); // 'baseLevel' added (13042011ih)
				bRet &= ProjectLevelToSurface<Pyramid, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel); // 'baseLevel' added (13042011ih)
				bRet &= ProjectLevelToSurface<Prism, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel); // 'baseLevel' added (13042011ih)
				bRet &= ProjectLevelToSurface<Hexahedron, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD, baseLevel); // 'baseLevel' added (13042011ih)
				break;
			default:
				UG_LOG("ERROR in 'ProjectLevelToSurface': Dimension of subset ("
						<< surfDD.dim_subset(si) << ") not supported.\n");
				return false;
		}
	}

#ifdef UG_PARALLEL
//	copy storage type into surf vector
	if(vLevelVector.size() > 0)
	{
		ParallelStorageType type = vLevelVector[vLevelVector.size()-1]->get_storage_mask();

	//	get intersection of types
		for(size_t lev =  baseLevel /*0*/; lev < vLevelVector.size(); ++lev) // '0' changed to 'baseLevel' (13042011ih)
			if(vLevelVector[lev] != NULL && vLevelVector[lev]->size() > 0)
				type = type & vLevelVector[lev]->get_storage_mask();

	//	check if union is defined
		if(type == PST_UNDEFINED)
		{
			UG_LOG("ERROR in 'ProjectLevelToSurface': storage type of level"
					" vectors is not ok. Must have at least on common type."
					" (e.g. additive/unique for all, or consistent for all)\n"
					" Types in levels are:\n");
 			for(size_t lev =  baseLevel /*0*/; lev < vLevelVector.size(); ++lev) // '0' changed to 'baseLevel' (13042011ih)
				if(vLevelVector[lev] != NULL)
					UG_LOG("  lev " << lev << ": " << vLevelVector[lev]->get_storage_mask() << "\n");
			return false;
		}

	//	set type of surface vector to common base
		surfVector.set_storage_type(type);
	}
#endif

//	we're done
	return bRet;
}



/**
 * This functions adds the shadow values from a coarser grid to the shadowing
 * DoFs on the finer grid.
 *
 * \param[out]	fineVec			fine grid vector
 * \param[out]	coarseVec		coarse grid vector
 * \param[in] 	approxSpace		Approximation Space
 * \param[in]	coarseLevel		Coarse Level index
 * \param[in]	fineLevel		Fine Level index
 */
template <typename TApproximationSpace, typename TVector>
bool AddProjectionOfVertexShadows(TVector& fineVec, const TVector& coarseVec,
                                  TApproximationSpace& approxSpace,
                                  size_t coarseLevel, size_t fineLevel)
{
//	get DoFDistributions
	const typename TApproximationSpace::dof_distribution_type& coarseDoFDistr
		= approxSpace.get_level_dof_distribution(coarseLevel);
	const typename TApproximationSpace::dof_distribution_type& fineDoFDistr
		= approxSpace.get_level_dof_distribution(fineLevel);

//	get surface view
	const SurfaceView& surfView = *approxSpace.get_surface_view();

//  Allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDoFDistr.num_fct(); ++fct)
		if(fineDoFDistr.local_shape_function_set_id(fct)
				!= LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1))
		{
			UG_LOG("ERROR in 'AssembleVertexProjection': "
					"Interpolation only implemented for Lagrange P1 functions.\n");
			return false;
		}

// 	get MultiGrid
	MultiGrid& grid = approxSpace.get_domain().get_grid();

	typename TApproximationSpace::dof_distribution_type::algebra_index_vector_type fineInd;
	typename TApproximationSpace::dof_distribution_type::algebra_index_vector_type coarseInd;

// 	Vertex iterators
	geometry_traits<VertexBase>::const_iterator iter, iterBegin, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < fineDoFDistr.num_subsets(); ++si)
	{
		iterBegin = fineDoFDistr.template begin<VertexBase>(si);
		iterEnd = fineDoFDistr.template end<VertexBase>(si);

	// 	loop nodes of fine subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	skip non-shadowing vertices
			if(!surfView.shadows(*iter)) continue;

		// 	get father
			GeometricObject* geomObj = grid.get_parent(*iter);
			VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);

		//	Check if father is Vertex
			if(vert != NULL)
			{
				// get global indices
				coarseDoFDistr.get_inner_algebra_indices(vert, coarseInd);
			}
			else continue;

		// 	get global indices
			fineDoFDistr.get_inner_algebra_indices(*iter, fineInd);

		//	add coarse vector entries to fine vector entries
			for(size_t i = 0; i < coarseInd.size(); ++i)
			{
				fineVec[fineInd[i]] += coarseVec[coarseInd[i]];
			}
		}
	}

//	we're done
	return true;
}

/**
 * This functions sets the shadowing values from a finer grid to the shadow
 * DoFs on the coarser grid.
 *
 * \param[out]	coarseVec		coarse grid vector
 * \param[out]	fineVec			fine grid vector
 * \param[in] 	approxSpace		Approximation Space
 * \param[in]	coarseLevel		Coarse Level index
 * \param[in]	fineLevel		Fine Level index
 */
template <typename TApproximationSpace, typename TVector>
bool SetProjectionOfVertexShadowing(TVector& coarseVec, const TVector& fineVec,
                                    TApproximationSpace& approxSpace,
                                    size_t coarseLevel, size_t fineLevel)
{
//	get DoFDistributions
	const typename TApproximationSpace::dof_distribution_type& coarseDoFDistr
		= approxSpace.get_level_dof_distribution(coarseLevel);
	const typename TApproximationSpace::dof_distribution_type& fineDoFDistr
		= approxSpace.get_level_dof_distribution(fineLevel);

//	get surface view
	const SurfaceView& surfView = *approxSpace.get_surface_view();

//  Allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDoFDistr.num_fct(); ++fct)
		if(fineDoFDistr.local_shape_function_set_id(fct)
				!= LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1))
		{
			UG_LOG("ERROR in 'AssembleVertexProjection': "
					"Interpolation only implemented for Lagrange P1 functions.\n");
			return false;
		}

// 	get MultiGrid
	MultiGrid& grid = approxSpace.get_domain().get_grid();

	typename TApproximationSpace::dof_distribution_type::algebra_index_vector_type fineInd;
	typename TApproximationSpace::dof_distribution_type::algebra_index_vector_type coarseInd;

// 	Vertex iterators
	geometry_traits<VertexBase>::const_iterator iter, iterBegin, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < fineDoFDistr.num_subsets(); ++si)
	{
		iterBegin = fineDoFDistr.template begin<VertexBase>(si);
		iterEnd = fineDoFDistr.template end<VertexBase>(si);

	// 	loop nodes of fine subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	skip non-shadowing vertices
			if(!surfView.shadows(*iter)) continue;

		// 	get father
			GeometricObject* geomObj = grid.get_parent(*iter);
			VertexBase* vert = dynamic_cast<VertexBase*>(geomObj);

		//	Check if father is Vertex
			if(vert != NULL)
			{
				// get global indices
				coarseDoFDistr.get_inner_algebra_indices(vert, coarseInd);
			}
			else continue;

		// 	get global indices
			fineDoFDistr.get_inner_algebra_indices(*iter, fineInd);

		//	add coarse vector entries to fine vector entries
			for(size_t i = 0; i < coarseInd.size(); ++i)
			{
				 coarseVec[coarseInd[i]] = fineVec[fineInd[i]];
			}
		}
	}

//	we're done
	return true;
}

/**
 * This functions sets the values of a vector to zero, where the index
 * corresponds to a refine-patch boundary (i.e. the vertex is a shadowing
 * vertex)
 *
 * \param[out]	vec				grid vector
 * \param[in] 	approxSpace		Approximation Space
 * \param[in]	level			Level index
 */
template <typename TApproximationSpace, typename TVector>
bool SetZeroOnShadowingVertex(TVector& vec,
                            TApproximationSpace& approxSpace,
                            size_t level)
{
//	get DoFDistributions
	const typename TApproximationSpace::dof_distribution_type& dofDistr
		= approxSpace.get_level_dof_distribution(level);

//	get surface view
	const SurfaceView& surfView = *approxSpace.get_surface_view();

//  Allow only lagrange P1 functions
	for(size_t fct = 0; fct < dofDistr.num_fct(); ++fct)
		if(dofDistr.local_shape_function_set_id(fct)
				!= LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1))
		{
			UG_LOG("ERROR in 'AssembleVertexProjection': "
					"Interpolation only implemented for Lagrange P1 functions.\n");
			return false;
		}

//	indices
	typename TApproximationSpace::dof_distribution_type::algebra_index_vector_type ind;

// 	Vertex iterators
	geometry_traits<VertexBase>::const_iterator iter, iterBegin, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < dofDistr.num_subsets(); ++si)
	{
		iterBegin = dofDistr.template begin<VertexBase>(si);
		iterEnd = dofDistr.template end<VertexBase>(si);

	// 	loop nodes of fine subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;

		//	skip non-shadowing vertices
			if(!surfView.shadows(vrt)) continue;

		// 	get global indices
			dofDistr.get_inner_algebra_indices(vrt, ind);

		//	add coarse vector entries to fine vector entries
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
 * corresponds to a vertex shadowed by a refine-patch boundary
 *
 * \param[out]	vec				grid vector
 * \param[in] 	approxSpace		Approximation Space
 * \param[in]	level			Level index
 */
template <typename TApproximationSpace, typename TVector>
bool SetZeroOnVertexShadows(TVector& vec,
                            TApproximationSpace& approxSpace,
                            size_t level)
{
//	get DoFDistributions
	const typename TApproximationSpace::dof_distribution_type& dofDistr
		= approxSpace.get_level_dof_distribution(level);

//	get surface view
	const SurfaceView& surfView = *approxSpace.get_surface_view();

//  Allow only lagrange P1 functions
	for(size_t fct = 0; fct < dofDistr.num_fct(); ++fct)
		if(dofDistr.local_shape_function_set_id(fct)
				!= LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1))
		{
			UG_LOG("ERROR in 'AssembleVertexProjection': "
					"Interpolation only implemented for Lagrange P1 functions.\n");
			return false;
		}

//	indices
	typename TApproximationSpace::dof_distribution_type::algebra_index_vector_type ind;

// 	Vertex iterators
	geometry_traits<VertexBase>::const_iterator iter, iterBegin, iterEnd;

// 	loop subsets of fine level
	for(int si = 0; si < dofDistr.num_subsets(); ++si)
	{
		iterBegin = dofDistr.template begin<VertexBase>(si);
		iterEnd = dofDistr.template end<VertexBase>(si);

	// 	loop nodes of fine subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;

		//	skip non-shadowing vertices
			if(!surfView.is_shadow(vrt)) continue;

		// 	get global indices
			dofDistr.get_inner_algebra_indices(vrt, ind);

		//	add coarse vector entries to fine vector entries
			for(size_t i = 0; i < ind.size(); ++i)
			{
				vec[ind[i]] = 0.0;
			}
		}
	}

//	we're done
	return true;
}


template <typename TElemBase>
bool SelectNonShadowsAdjacentToShadows(ISelector& sel, const SurfaceView& surfView)
{
//	vectors for associated elements
	std::vector<VertexBase*> vAssVertex;
	std::vector<EdgeBase*> vAssEdge;
	std::vector<Face*> vAssFace;
	std::vector<Volume*> vAssVolume;

//	get grid
	Grid& grid = *sel.get_assigned_grid();

//	loop all subsets
	for(int si = 0; si < surfView.num_subsets(); ++si)
	{
	//	iterator type
		typename geometry_traits<TElemBase>::const_iterator iter, iterEnd;
		iterEnd = surfView.end<TElemBase>(si);

	//	loop all base elems
		for(iter = surfView.begin<TElemBase>(si); iter != iterEnd; ++iter)
		{
		//	get element
			TElemBase* elem = *iter;

		//	check if element is a shadow
			if(surfView.shadows(elem))
			{
			//	get shadow
				GeometricObject* shadow = surfView.get_parent(elem);

			//	get adjacent elemens
				CollectAssociated(vAssVertex, grid, shadow);
				CollectAssociated(vAssEdge, grid, shadow);
				CollectAssociated(vAssFace, grid, shadow);
				CollectAssociated(vAssVolume, grid, shadow);

			//	select associated elements
				for(size_t i = 0; i < vAssVertex.size(); ++i)
					if(surfView.is_contained(vAssVertex[i]))
						sel.select(vAssVertex[i]);
				for(size_t i = 0; i < vAssEdge.size(); ++i)
					if(surfView.is_contained(vAssEdge[i]))
						sel.select(vAssEdge[i]);
				for(size_t i = 0; i < vAssFace.size(); ++i)
					if(surfView.is_contained(vAssFace[i]))
						sel.select(vAssFace[i]);
				for(size_t i = 0; i < vAssVolume.size(); ++i)
					if(surfView.is_contained(vAssVolume[i]))
						sel.select(vAssVolume[i]);
			}
		}

	}

//	we're done
	return true;
}


inline bool SelectNonShadowsAdjacentToShadows(ISelector& sel, const SurfaceView& surfView)
{
//	clear all marks
	sel.clear();

//	get grid
	Grid& grid = *sel.get_assigned_grid();

//	select elements
	bool bRes = true;

//	note: the highest dimension of elements must not be loop, since there are
//		  no slaves of the highest dimension
	bRes &= SelectNonShadowsAdjacentToShadows<VertexBase>(sel, surfView);

	if(grid.num<Face>() > 0 || grid.num<Volume>() > 0)
		bRes &= SelectNonShadowsAdjacentToShadows<EdgeBase>(sel, surfView);

	if(grid.num<Volume>() > 0)
		bRes &= SelectNonShadowsAdjacentToShadows<Face>(sel, surfView);

//	we're done
	return bRes;
}

template <typename TElemBase>
bool SelectNonShadowsAdjacentToShadowsOnLevel(ISelector& sel,
                                              const SurfaceView& surfView,
                                              int level)
{
//	vectors for associated elements
	std::vector<VertexBase*> vAssVertex;
	std::vector<EdgeBase*> vAssEdge;
	std::vector<Face*> vAssFace;
	std::vector<Volume*> vAssVolume;

//	get grid
	Grid& grid = *sel.get_assigned_grid();

//	get multigrid
	MultiGrid& mg = *dynamic_cast<MultiGrid*>(&grid);

//	check multigrid
	if(&mg == NULL)
	{
		UG_LOG("ERROR in SelectNonShadowsAdjacentToShadowsOnLevel: No "
				"Multigrid given, selection ob level not possible.\n");
		return false;
	}

//	check level
	if(level >= (int) mg.num_levels() || level < 0)
	{
		UG_LOG("ERROR in SelectNonShadowsAdjacentToShadowsOnLevel: Requested "
				"level "<<level<<" does not exist in Multigrid.\n");
		return false;
	}

//	loop all subsets
	for(int si = 0; si < surfView.num_subsets(); ++si)
	{
	//	iterator type
		typename geometry_traits<TElemBase>::const_iterator iter, iterEnd;
		iterEnd = surfView.end<TElemBase>(si);

	//	loop all base elems
		for(iter = surfView.begin<TElemBase>(si); iter != iterEnd; ++iter)
		{
		//	get element
			TElemBase* elem = *iter;

		//	check if element is a shadow
			if(surfView.shadows(elem))
			{
			//	get shadow
				GeometricObject* shadow = surfView.get_parent(elem);

			//	check if this is the correct level
				if(mg.get_level(shadow) != level) continue;

			//	get adjacent elements
				CollectAssociated(vAssVertex, grid, shadow);
				CollectAssociated(vAssEdge, grid, shadow);
				CollectAssociated(vAssFace, grid, shadow);
				CollectAssociated(vAssVolume, grid, shadow);

			//	select associated elements
				for(size_t i = 0; i < vAssVertex.size(); ++i)
					if(surfView.is_contained(vAssVertex[i]))
						sel.select(vAssVertex[i]);
				for(size_t i = 0; i < vAssEdge.size(); ++i)
					if(surfView.is_contained(vAssEdge[i]))
						sel.select(vAssEdge[i]);
				for(size_t i = 0; i < vAssFace.size(); ++i)
					if(surfView.is_contained(vAssFace[i]))
						sel.select(vAssFace[i]);
				for(size_t i = 0; i < vAssVolume.size(); ++i)
					if(surfView.is_contained(vAssVolume[i]))
						sel.select(vAssVolume[i]);
			}
		}

	}

//	we're done
	return true;
}


inline bool SelectNonShadowsAdjacentToShadowsOnLevel(ISelector& sel,
                                              const SurfaceView& surfView,
                                              int level)
{
//	clear all marks
	sel.clear();

//	get grid
	Grid& grid = *sel.get_assigned_grid();

//	get multigrid
	MultiGrid& mg = *dynamic_cast<MultiGrid*>(&grid);

//	check multigrid
	if(&mg == NULL)
	{
		UG_LOG("ERROR in SelectNonShadowsAdjacentToShadowsOnLevel: No "
				"Multigrid given, selection ob level not possible.\n");
		return false;
	}

//	check level
	if(level >= (int) mg.num_levels() || level < 0)
	{
		UG_LOG("ERROR in SelectNonShadowsAdjacentToShadowsOnLevel: Requested "
				"level "<<level<<" does not exist in Multigrid.\n");
		return false;
	}

//	select elements
	bool bRes = true;

//	note: the highest dimension of elements must not be loop, since there are
//		  no slaves of the highest dimension
	bRes &= SelectNonShadowsAdjacentToShadowsOnLevel<VertexBase>(sel, surfView, level);

	if(grid.num<Face>() > 0 || grid.num<Volume>() > 0)
		bRes &= SelectNonShadowsAdjacentToShadowsOnLevel<EdgeBase>(sel, surfView, level);

	if(grid.num<Volume>() > 0)
		bRes &= SelectNonShadowsAdjacentToShadowsOnLevel<Face>(sel, surfView, level);

//	we're done
	return bRes;
}

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

/// projects surface function to level functions
template <typename TMatrix>
bool CopySmoothingMatrix(TMatrix& smoothMat,
                         const std::vector<int>& vMap,
                         const TMatrix& origMat)
{
//	type of matrix row iterator
	typedef typename TMatrix::const_row_iterator const_row_iterator;

//	loop all mapped indices
	for(size_t origInd = 0; origInd < vMap.size(); ++origInd)
	{
	//	get mapped level index
		const int smoothInd = vMap[origInd];

	//	skipped non-mapped indices (indicated by -1)
		if(smoothInd < 0) continue;

	//	loop all connections of the surface dof to other surface dofs and copy
	//	the matrix coupling into the level matrix

		for(const_row_iterator conn = origMat.begin_row(origInd);
						conn != origMat.end_row(origInd); ++conn)
		{
		//	get corresponding level connection index
			const int origConnIndex = vMap[conn.index()];

		//	check that index is from same level, too
			if(origConnIndex < 0) continue;

		//	copy connection to level matrix
			smoothMat(smoothInd, origConnIndex) = conn.value();
		}
	}

#ifdef UG_PARALLEL
	smoothMat.copy_storage_type(origMat);
#endif

	return true;
}


} // end namespace ug
#endif
