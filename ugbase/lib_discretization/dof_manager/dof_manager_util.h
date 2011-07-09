/*
 * dof_manager_util.h
 *
 *  Created on: 19.04.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER_UTIL__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER_UTIL__

namespace ug{


/// projects surface function to level functions
template <typename TElem, typename TDoFDistributionImpl>
bool CreateSurfaceToLevelMapping(int si,
                                 std::vector<std::vector<int> >& vSurfLevelMapping,
                                 const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD,
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
		surfaceDD.algebra_indices(elem, surfaceInd);

	//	get level of element in hierarchy
		int level = surfaceView.get_level(elem);

	//	get corresponding level matrix for element
		UG_ASSERT(level < (int)vSurfLevelMapping.size(), "Level missing");
		std::vector<int>& levelMapping = vSurfLevelMapping[level];

	//	check that level is correct
		UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, "
				"									that is not passed.");

	//	extract all algebra indices for the element on level
		UG_ASSERT(vLevelDD[level] != NULL, "DoF Distribution missing");
		vLevelDD[level]->algebra_indices(elem, levelInd);

	//	check that index sets have same cardinality
		UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

	//	copy all elements of the matrix
		for(size_t i = 0; i < surfaceInd.size(); ++i)
		{
			UG_ASSERT(surfaceInd[i] < levelMapping.size(), "Index to large.");
			levelMapping[surfaceInd[i]] = levelInd[i];
		}
	}

//	we're done
	return true;
}

/// creates a mapping of indices from the surface dof distribution to the level dof distribution
template <typename TDoFDistributionImpl>
bool CreateSurfaceToLevelMapping(std::vector<std::vector<int> >& vSurfLevelMapping,
                                 const std::vector<const IDoFDistribution<TDoFDistributionImpl>*>& vLevelDD,
                                 const IDoFDistribution<TDoFDistributionImpl>& surfDD,
                                 const SurfaceView& surfView)
{
//	resize the mapping
	vSurfLevelMapping.clear();
	vSurfLevelMapping.resize(vLevelDD.size());
	for(size_t lev = 0; lev < vSurfLevelMapping.size(); ++lev)
		vSurfLevelMapping[lev].resize(surfDD.num_dofs(), -1);

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
				bRet &= CreateSurfaceToLevelMapping<VertexBase, TDoFDistributionImpl>
							(si, vSurfLevelMapping, vLevelDD, surfDD, surfView);
				break;
			case 1:
				bRet &= CreateSurfaceToLevelMapping<Edge, TDoFDistributionImpl>
							(si, vSurfLevelMapping, vLevelDD, surfDD, surfView);
				break;
			case 2:
				bRet &= CreateSurfaceToLevelMapping<Triangle, TDoFDistributionImpl>
							(si, vSurfLevelMapping, vLevelDD, surfDD, surfView);
				bRet &= CreateSurfaceToLevelMapping<Quadrilateral, TDoFDistributionImpl>
							(si, vSurfLevelMapping, vLevelDD, surfDD, surfView);
				break;
			case 3:
				bRet &= CreateSurfaceToLevelMapping<Tetrahedron, TDoFDistributionImpl>
							(si, vSurfLevelMapping, vLevelDD, surfDD, surfView);
				bRet &= CreateSurfaceToLevelMapping<Pyramid, TDoFDistributionImpl>
							(si, vSurfLevelMapping, vLevelDD, surfDD, surfView);
				bRet &= CreateSurfaceToLevelMapping<Prism, TDoFDistributionImpl>
							(si, vSurfLevelMapping, vLevelDD, surfDD, surfView);
				bRet &= CreateSurfaceToLevelMapping<Hexahedron, TDoFDistributionImpl>
							(si, vSurfLevelMapping, vLevelDD, surfDD, surfView);
				break;
			default:
				UG_LOG("ERROR in 'CreateSurfaceToLevelMapping': Dimension of subset ("
						<< surfDD.dim_subset(si) << ") not supported.\n");
				return false;
		}
	}

//	we're done
	return bRet;
}

/// projects surface function to level functions
template <typename TMatrix>
bool CopyMatrixSurfaceToLevel(TMatrix& levelMatrix,
                              const std::vector<int>& surfLevelMapping,
                              const TMatrix& surfMatrix)
{
//	type of matrix row iterator
	typedef typename TMatrix::const_row_iterator const_row_iterator;

//	loop all mapped indices
	for(size_t surfInd = 0; surfInd < surfLevelMapping.size(); ++surfInd)
	{
	//	get mapped level index
		const int levelInd = surfLevelMapping[surfInd];

	//	skipped non-mapped indices (indicated by -1)
		if(levelInd < 0) continue;

	//	loop all connections of the surface dof to other surface dofs and copy
	//	the matrix coupling into the level matrix

		for(const_row_iterator conn = surfMatrix.begin_row(surfInd);
						conn != surfMatrix.end_row(surfInd); ++conn)
		{
		//	get corresponding level connection index
			const int levelConnIndex = surfLevelMapping[conn.index()];

		//	check that index is from same level, too
			if(levelConnIndex < 0) continue;

		//	copy connection to level matrix
			levelMatrix(levelInd, levelConnIndex) = conn.value();
		}
	}

#ifdef UG_PARALLEL
	levelMatrix.copy_storage_type(surfMatrix);
#endif

	return true;
}

/// projects surface function to level functions
template <typename TMatrix>
bool CopyMatrixSurfaceToLevel(std::vector<TMatrix*>& vLevelMatrix,
                              const std::vector<std::vector<int> >& vSurfLevelMapping,
                              const TMatrix& surfMatrix)
{
//	check sizes
	if(vSurfLevelMapping.size() != vLevelMatrix.size())
	{
		UG_LOG("ERROR in CopyMatrixSurfaceToLevel: Number of matrices and "
				"mappings does not match. cannot proceed.\n");
		return false;
	}

//	loop all level
	bool bRes = true;
	for(size_t lev = 0; lev < vSurfLevelMapping.size(); ++lev)
		bRes &= CopyMatrixSurfaceToLevel(vSurfLevelMapping[lev], *vLevelMatrix[lev], surfMatrix);

//	we're done
	return bRes;
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__DOF_MANAGER_UTIL__ */
