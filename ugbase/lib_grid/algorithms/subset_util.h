//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__SUBSET_UTIL__
#define __H__LIB_GRID__SUBSET_UTIL__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

/** \defgroup subsetUtil Subset Util
 * Several methods that ease subset-handling are grouped here.
 * @{
 */

////////////////////////////////////////////////////////////////////////
//	FindFirstFreeSubset
///	returns the index of the last subset, that contains elements of type TElem.
/**
 * returns -1 if no subset contains elements of type TElem.
 * All geometric objects and the 4 geometric base objects are valid types.
 *
 * call this method like this: GetMaxSubsetIndex<Face>(yourSubsetHandler);
 */
template <class TElem>
int GetMaxSubsetIndex(SubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
//	MakeSubsetsConsecutive
///	moves subsets so that no empty subset of type TElem is between filled ones.
/**
 * This algorithm does not change the order of filled subsets.
 * Since this algorithm only checks for elements of type TElem, it
 * is possible that after termination subsets are empty regarding
 * other element-types.
 */
template <class TElem>
void MakeSubsetsConsecutive(SubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
//	AssignFaceInterfaceEdgesToSubsets
///	assigns edges which belong to no subset and are adjacent to faces of different subsets to new subsets.
/**
 * the edges of each interface between two subset are put into
 * a seperate subset, starting at GetMaxSubsetIndex<EdgeBase>(sh) + 1.
 */
void AssignFaceInterfaceEdgesToSubsets(Grid& grid, SubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
//	AssignVolumeInterfaceFacesToSubsets
///	assigns faces which belong to no subset and are adjacent to volumes of different subsets to new subsets.
/**
 * the faces of each interface between two subset are put into
 * a seperate subset, starting at GetMaxSubsetIndex<Face>(sh) + 1.
 */
void AssignVolumeInterfaceFacesToSubsets(Grid& grid, SubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
//	AssignAssociatedVerticesToSubset
///	assigns vertices of the given elements to the subset at subsetIndex
/**
 * TIterator should be an stl-compatible iterator. Its value_type should be
 * a pointer to either EdgeBase, Face, Volume or an derived type of the three.
 */
template <class TIterator>
void AssignAssociatedVerticesToSubset(ISubsetHandler& sh, TIterator elemsBegin,
										TIterator elemsEnd, int subsetIndex);
										
////////////////////////////////////////////////////////////////////////
//	AdjustSubsetsForLgmNg
///	reorders subsets in a way that allows for easy export to lgm-ng.
/**
 * .lgm and .ng are the geometry-file-types of ug3. They require a very
 * special subset-structure. This method will try to adjust the
 * subsets in a way that makes export using \sa ExportGridToUG
 * easy.
 * Please note, that this method does not yet produce fully compatible
 * output. In order to ensure maximal compatibility you should avoid
 * to have empty subsets between filled ones.
 */
void AdjustSubsetsForLgmNg(Grid& grid, SubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
//	SeparateFaceSubsetsByNormal
///	separates faces by orthogonal axis-aligned normals.
/**
 * The faces of each existing subset are assigned to new subsets
 * based on the closest axis-aligned normal.
 * Faces with similar normals which are contained in different
 * subsets are assigned to different subsets.
 *
 * \param paNorm pointer to the normal attachment. NULL indicates
 * 				that normals shall be calculated on the fly (default).
 */
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								APosition aPos = aPosition,
								ANormal* paNorm = NULL);

////////////////////////////////////////////////////////////////////////
//	SeparateFaceSubsetsByNormal
///	separates subset by the given normals.
/**
 * The faces of each existing subset are assigned to new subsets
 * based on the closest matching normal.
 * Faces with similar normals which are contained in different
 * subsets are assigned to different subsets.
 *
 * \param paNorm pointer to the normal attachment. NULL indicates
 * 				that normals shall be calculated on the fly (default).
 */
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								std::vector<vector3> vNormals,
								APosition aPos = aPosition,
								ANormal* paNorm = NULL);

////////////////////////////////////////////////////////////////////////
//	SeparateVolumesByFaceSubsets
///	groups volumes by separating face-subsets.
/**
 * all volumes that are surrounded by the same face-subsets are
 * assigned to a common subset.
 * If material points are specified, the group of volumes is
 * assigned to the subset depending on the index of the
 * contained material-point.
 * \todo implement material-point support.
 */
void SeparateVolumesByFaceSubsets(Grid& grid, SubsetHandler& sh,
									vector3* pMaterialPoints = NULL,
									int numMaterialPoints = 0);

////////////////////////////////////////////////////////////////////////
//	AssignSubsetColors
///	assigns a different color to each subset
void AssignSubsetColors(ISubsetHandler& sh);

/**@}*/ // end of doxygen defgroup command
}//	end of namespace

////////////////////////////////////////////////
// include implementations of template methods.
#include "subset_util_impl.hpp"

#endif
