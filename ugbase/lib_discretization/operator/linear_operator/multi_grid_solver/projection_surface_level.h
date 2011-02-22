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
		TVector& levelVector = *(vLevelVector[level]);

	//	check that level is correct
		UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, that is not passed.");

	//	extract all algebra indices for the element on level
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
	{
		vLevelVector[lev]->copy_storage_type(surfVector);
	}
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
                           const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD)
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
		const TVector& levelVector = *(vLevelVector[level]);

	//	check that level is correct
		UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, that is not passed.");

	//	extract all algebra indices for the element on level
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
                           const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD)
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
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD);
				break;
			case 1:
				bRet &= ProjectLevelToSurface<Edge, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD);
				break;
			case 2:
				bRet &= ProjectLevelToSurface<Triangle, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD);
				bRet &= ProjectLevelToSurface<Quadrilateral, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD);
				break;
			case 3:
				bRet &= ProjectLevelToSurface<Tetrahedron, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD);
				bRet &= ProjectLevelToSurface<Pyramid, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD);
				bRet &= ProjectLevelToSurface<Prism, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD);
				bRet &= ProjectLevelToSurface<Hexahedron, TVector, TDoFDistributionImpl>
							(si, surfVector, surfDD, surfView, vLevelVector, vLevelDD);
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
		ParallelStorageType type = vLevelVector[0]->get_storage_mask();

	//	get intersection of types
		for(size_t lev = 1; lev < vLevelVector.size(); ++lev)
			type = type & vLevelVector[lev]->get_storage_mask();

	//	check if union is defined
		if(type == PST_UNDEFINED)
		{
			UG_LOG("ERROR in 'ProjectLevelToSurface': storage type of level"
					" vectors is not ok. Must have at least on common type."
					" (e.g. additive/unique for all, or consistent for all)\n"
					" Types in levels are:\n");
			for(size_t lev = 0; lev < vLevelVector.size(); ++lev)
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

} // end namespace ug
#endif
