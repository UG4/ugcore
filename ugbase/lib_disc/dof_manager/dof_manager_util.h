/*
 * dof_manager_util.h
 *
 *  Created on: 19.04.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__DOF_MANAGER_UTIL__
#define __H__UG__LIB_DISC__DOF_MANAGER__DOF_MANAGER_UTIL__

namespace ug{


/// projects surface function to level functions
template <typename TElem>
void CreateSurfaceToLevelMapping(std::vector<std::vector<int> >& vSurfLevelMapping,
                                 const std::vector<ConstSmartPtr<LevelDoFDistribution> >& vLevelDD,
                                 ConstSmartPtr<SurfaceDoFDistribution> surfaceDD,
                                 const SurfaceView& surfaceView)
{
//	type of element iterator
	typedef typename SurfaceDoFDistribution::traits<TElem>::const_iterator iter_type;

//	iterators for subset
	iter_type iter = surfaceDD->begin<TElem>();
	iter_type iterEnd = surfaceDD->end<TElem>();

//	vector of indices
	std::vector<size_t> surfaceInd, levelInd;

//	loop all elements of type
	for( ; iter != iterEnd; ++iter)
	{
	//	get elem
		TElem* elem = *iter;

	//	extract all algebra indices for the element on surface
		surfaceDD->inner_algebra_indices(elem, surfaceInd);

	//	get level of element in hierarchy
		int level = surfaceView.get_level(elem);

	//	get corresponding level matrix for element
		UG_ASSERT(level < (int)vSurfLevelMapping.size(), "Level missing");
		std::vector<int>& levelMapping = vSurfLevelMapping[level];

	//	check that level is correct
		UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, "
				"									that is not passed.");

	//	extract all algebra indices for the element on level
		UG_ASSERT(vLevelDD[level].valid(), "DoF Distribution missing");
		vLevelDD[level]->inner_algebra_indices(elem, levelInd);

	//	check that index sets have same cardinality
		UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

	//	copy all elements of the matrix
		for(size_t i = 0; i < surfaceInd.size(); ++i)
		{
			UG_ASSERT(surfaceInd[i] < levelMapping.size(), "Index to large.");
			levelMapping[surfaceInd[i]] = levelInd[i];
		}

	//	same for shadows
		TElem* parent = surfaceView.parent_if_copy(elem);
		while(parent && surfaceView.is_shadowed(parent))
		{
		//	get level of element in hierarchy
			int level = surfaceView.get_level(parent);

		//	get corresponding level matrix for element
			UG_ASSERT(level < (int)vSurfLevelMapping.size(), "Level missing");
			std::vector<int>& levelMapping = vSurfLevelMapping[level];

		//	check that level is correct
			UG_ASSERT(level < (int)vLevelDD.size(), "Element of level detected, "
					"									that is not passed.");

		//	extract all algebra indices for the element on level
			UG_ASSERT(vLevelDD[level].valid(), "DoF Distribution missing");
			vLevelDD[level]->inner_algebra_indices(parent, levelInd);

		//	check that index sets have same cardinality
			UG_ASSERT(surfaceInd.size() == levelInd.size(), "Number of indices does not match.");

		//	copy all elements of the matrix
			for(size_t i = 0; i < surfaceInd.size(); ++i)
			{
				UG_ASSERT(surfaceInd[i] < levelMapping.size(), "Index to large.");
				levelMapping[surfaceInd[i]] = levelInd[i];
			}

		//	next shadow
			elem = parent;
			parent = surfaceView.parent_if_copy(elem);
		}

	}
}

/// creates a mapping of indices from the surface dof distribution to the level dof distribution
inline void CreateSurfaceToLevelMapping(std::vector<std::vector<int> >& vSurfLevelMapping,
                                        const std::vector<ConstSmartPtr<LevelDoFDistribution> >& vLevelDD,
                                        ConstSmartPtr<SurfaceDoFDistribution> surfDD,
                                        const SurfaceView& surfView)
{
//	resize the mapping
	vSurfLevelMapping.clear();
	vSurfLevelMapping.resize(vLevelDD.size());
	for(size_t lev = 0; lev < vSurfLevelMapping.size(); ++lev)
		vSurfLevelMapping[lev].resize(surfDD->num_indices(), -1);

	if(surfDD->max_dofs(VERTEX))
		CreateSurfaceToLevelMapping<VertexBase>(vSurfLevelMapping, vLevelDD, surfDD, surfView);
	if(surfDD->max_dofs(EDGE))
		CreateSurfaceToLevelMapping<EdgeBase>(vSurfLevelMapping, vLevelDD, surfDD, surfView);
	if(surfDD->max_dofs(FACE))
		CreateSurfaceToLevelMapping<Face>(vSurfLevelMapping, vLevelDD, surfDD, surfView);
	if(surfDD->max_dofs(VOLUME))
		CreateSurfaceToLevelMapping<Volume>(vSurfLevelMapping, vLevelDD, surfDD, surfView);
}

/// projects surface function to level functions
template <typename TMatrix>
void CopyMatrixSurfaceToLevel(TMatrix& levelMatrix,
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
}

/// projects surface function to level functions
template <typename TMatrix>
void CopyMatrixSurfaceToLevel(std::vector<TMatrix*>& vLevelMatrix,
                              const std::vector<std::vector<int> >& vSurfLevelMapping,
                              const TMatrix& surfMatrix)
{
//	check sizes
	if(vSurfLevelMapping.size() != vLevelMatrix.size())
		UG_THROW("CopyMatrixSurfaceToLevel: Number of matrices and "
				"mappings does not match. cannot proceed.");

//	loop all level
	for(size_t lev = 0; lev < vSurfLevelMapping.size(); ++lev)
		CopyMatrixSurfaceToLevel(vSurfLevelMapping[lev], *vLevelMatrix[lev], surfMatrix);
}


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__DOF_MANAGER_UTIL__ */
