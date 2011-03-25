//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__SUBSET_UTIL__
#define __H__LIB_GRID__SUBSET_UTIL__

#include <vector>
#include "lib_grid/lg_base.h"
#include "callbacks/callbacks.h"

namespace ug
{

/**
 * Several methods that ease subset-handling are grouped here. 
 * \defgroup lib_grid_algorithms_subset_util subset util
 * \ingroup lib_grid_algorithms
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
//	AssignAssociatedVerticesToSubsets
///	Assigns associated vertices of elements of type TElem in sh to sh.
/**
 * Make sure that TElem is not of type VertexBase or any derivate.
 *
 * The method iterates over all elements of type TElem in sh and
 * assigns associated vertices to sh. The target subset-index is taken
 * from srcIndHandler.
 *
 * Valid types for TSubsetHandler are SubsetHandler and MGSubsetHandler
 * compatible types.
 *
 * This method is e.g. used for SurfaceView creation.
 */
template <class TElem, class TSubsetHandler>
void AssignAssociatedVerticesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler);

////////////////////////////////////////////////////////////////////////
//	AssignAssociatedEdgesToSubsets
///	Assigns associated edges of elements of type TElem in sh to sh.
/**
 * Make sure that TElem is not of type EdgeBase or any derivate.
 *
 * The method iterates over all elements of type TElem in sh and
 * assigns associated edges to sh. The target subset-index is taken
 * from srcIndHandler.
 *
 * Valid types for TSubsetHandler are SubsetHandler and MGSubsetHandler
 * compatible types.
 *
 * This method is e.g. used for SurfaceView creation.
 */
template <class TElem, class TSubsetHandler>
void AssignAssociatedEdgesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler);

////////////////////////////////////////////////////////////////////////
//	AssignAssociatedFacesToSubsets
///	Assigns associated faces of elements of type TElem in sh to sh.
/**
 * Make sure that TElem is not of type Face or any derivate.
 *
 * The method iterates over all elements of type TElem in sh and
 * assigns associated faces to sh. The target subset-index is taken
 * from srcIndHandler.
 *
 * Valid types for TSubsetHandler are SubsetHandler and MGSubsetHandler
 * compatible types.
 *
 * This method is e.g. used for SurfaceView creation.
 */
template <class TElem, class TSubsetHandler>
void AssignAssociatedFacesToSubsets(TSubsetHandler& sh,
									const ISubsetHandler& srcIndHandler);

////////////////////////////////////////////////////////////////////////
//	AssignAssociatedFacesToSubsets
///	Assigns associated sides of elements of type TElem in sh to sh.
/**
 * The method iterates over all elements of type TElem in sh and
 * assigns associated sides to sh. The target subset-index is taken
 * from srcIndHandler.
 *
 * Valid types for TSubsetHandler are SubsetHandler and MGSubsetHandler
 * compatible types.
 *
 * This method is e.g. used for SurfaceView creation.
 */
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedSidesToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler);

////////////////////////////////////////////////////////////////////////
//	AssignAssociatedLowerDimElems
///	Assigns associated elements of elements of type TElem in sh to sh.
/**
 * The method iterates over all elements of type TElem in sh and assigns
 * associated elements of lower dimension to sh. The subset-index to
 * which those elements are assigned are taken from srcIndHandler.
 *
 * Associated elements that are assigned to sh and have a subset-index
 * of -1 in srcIndHandler are assigned to the subset at alternativeSubsetIndex.
 *
 * Valid types for TSubsetHandler are SubsetHandler and MGSubsetHandler
 * compatible types.
 *
 * This method is e.g. used for SurfaceView creation.
 */
template <class TElem, class TSubsetHandlerDest, class TSubsetHandlerSrc>
void AssignAssociatedLowerDimElemsToSubsets(TSubsetHandlerDest& sh,
									const TSubsetHandlerSrc& srcIndHandler);

////////////////////////////////////////////////////////////////////////
//	CreateSurfaceView
///	Collects all elements between iterBegin and iterEnd that don't have any children.
/**
 * DEPRECIATED
 *
 * Elements which are on the surface of the multi-grid-hierarchy
 * (elements that don't have children) are assigned to a subset of the
 * shSurfaceViewOut. The subset-index is taken from sh.
 *
 * TIterator has to be an STL compatible iterator, whose value-type is a
 * pointer to a VertexBase, EdgeBase, Face, Volume or derived class.
 *
 * make sure that all elements between iterBegin and iterEnd are members
 * of the given MultiGrid.
 *
 * This method will extend the surface-view. The caller is responsible for
 * clearing it before calling this method.
 */
/*
template <class TIterator>
void CreateSurfaceView(SubsetHandler& shSurfaceViewOut, MultiGrid& mg,
						ISubsetHandler& sh, TIterator iterBegin,
						TIterator iterEnd);
*/
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
 *
 * Through keepExistingInterfaceSubsets you can decide whether interface
 * faces or edges which already are assigned to a subset shall be kept.
 * Note that the algorithm may not produce a fully compatible geometry
 * in this case (this depends on your initial subsets).
 */
void AdjustSubsetsForLgmNg(Grid& grid, SubsetHandler& sh,
							bool keepExistingInterfaceSubsets = false);

////////////////////////////////////////////////////////////////////////
//	SplitIrregularManifoldSubset
///	Keeps a regular part in the subset and assigns all other faces to another one.
/** THIS ALGORITHM USES Grid::mark
 *
 * If the faces in the given subset build an irregular manifold, then
 * this algorithm finds a regular part and assigns all other faces to
 * the given targetIndex.
 *
 * \return	true if the subset was splitted, false if not.
 */
bool SplitIrregularManifoldSubset(SubsetHandler& sh, int srcIndex,
								  int targetIndex);

////////////////////////////////////////////////////////////////////////
//	SeparateFaceSubsetsByNormal
///	separates faces by orthogonal axis-aligned normals.
/**
 * The faces of each existing subset are assigned to new subsets
 * based on the closest axis-aligned normal.
 * Faces with similar normals which are contained in different
 * subsets are assigned to different subsets.
 *
 * \param grid 	Grid
 * \param sh	Subset Handler
 * \param aPos	Position Attachment
 * \param paNorm pointer to the normal attachment. NULL indicates
 * 				that normals shall be calculated on the fly (default).
 * \param applyToSubset		Allows to specify which subset shall be separated.
 * 							-2: All subsets,
 * 							0, ..., numSubsets: The specified subset only.
 */
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								APosition aPos = aPosition,
								ANormal* paNorm = NULL,
								int applyToSubset = -2);

////////////////////////////////////////////////////////////////////////
//	SeparateFaceSubsetsByNormal
///	separates subset by the given normals.
/**
 * The faces of each existing subset are assigned to new subsets
 * based on the closest matching normal.
 * Faces with similar normals which are contained in different
 * subsets are assigned to different subsets.
 *
 * \param grid 		Grid
 * \param sh		Subset Handler
 * \param vNormals 	normals
 * \param aPos		Position Attachment
 * \param paNorm pointer to the normal attachment. NULL indicates
 * 				that normals shall be calculated on the fly (default).
 * \param applyToSubset		Allows to specify which subset shall be separated.
 * 							-2: All subsets,
 * 							0, ..., numSubsets: The specified subset only.
 */
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								std::vector<ug::vector3> vNormals,
								APosition aPos = aPosition,
								ANormal* paNorm = NULL,
								int applyToSubset = -2);

////////////////////////////////////////////////////////////////////////
///	assigns a region of volumes to a subset.
/**
 * Indirectly uses Grid::mark.
 *
 * All volumes that are connected directly or indirectly to proxyVol
 * (that means on can reach the volume without traversing a face that
 * is assigned to a subset) are considered to be in the same region
 * as proxyVol. All those volumes will be assigned to to the subset
 * in shVolsOut specified by newSubsetIndex.
 *
 * shFaces will be used to determine whether a side of the volume
 * is contained in a subset (!= -1) and is thus considered to be a
 * boundary of the region.
 *
 * shFaces and shVolsOut may refer to the same subset handler.
 */
void AssignRegionToSubset(Grid& grid, ISubsetHandler& shVolsOut,
						  const ISubsetHandler& shFaces,
						  Volume* proxyVol, int newSubsetIndex);

////////////////////////////////////////////////////////////////////////
///	finds regions by marker-points
/**
 * Indirectly uses Grid::mark.
 * If a region contains the i-th marker-point, it is assigned to the
 * subset (firstSubsetIndex + i).
 *
 * A region is a group of volumes that is surrounded by faces that lie
 * in subsets != -1.
 *
 * NOTE that this method currently only works for tetrahedral grids
 * (at least the markers are only associated with tetrahedrons).
 */
bool SeparateRegions(Grid& grid, ISubsetHandler& shVolsOut,
					 const ISubsetHandler& shFaces,
					 const MarkerPointManager& mpm,
					 int firstSubsetIndex);

////////////////////////////////////////////////////////////////////////
//	SeparateSubsetsByLowerDimSubsets
///	Assigns all elements of the given type to subsets.
/**	Different subsets are created for different regions. A region
 * is a set of elements of the given type, which are surrounded by
 * a closed set of lower dimensional elements, which are all assigned to
 * a subset.
 */
template <class TElem>
void SeparateSubsetsByLowerDimSubsets(Grid& grid, SubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
//	SeparateSubsetsByLowerDimSelection
///	Assigns all elements of the given type to subsets.
/**	Different subsets are created for different regions. A region
 * is a set of elements of the given type, which are surrounded by
 * a closed set of lower dimensional elements, which are all assigned to
 * a subset.
 */
template <class TElem>
void SeparateSubsetsByLowerDimSelection(Grid& grid, SubsetHandler& sh,
										Selector& sel);

////////////////////////////////////////////////////////////////////////
// 	SeparateSubsetsByLowerDimSeparators
///	Assigns all elements of the given type to subsets.
/**	Different subsets are created for different regions. A region
 * is a set of elements of the given type, which are surrounded by
 * a closed set of lower dimensional elements, which are all separators.
 *
 * Through a callback one can specify the elements which separate subsets.
 * If the callback returns true for an element, the element is regarded as
 * a separator.
 *
 * Note that the callback operates on elements which have one dimension less
 * than TElem. The callback is compatible with the CB_ConsiderVertex,
 * CB_ConsiderEdge and CB_ConsiderFace callbacks.
 */
template <class TElem>
void SeparateSubsetsByLowerDimSeparators(Grid& grid, SubsetHandler& sh,
					boost::function<bool (typename TElem::lower_dim_base_object*)>
						cbIsSeparator);

////////////////////////////////////////////////////////////////////////
//	AssignInnerAndBoundarySubsets
///	assigns objects to subsets depending on whether they are inner or boundary objects.
void AssignInnerAndBoundarySubsets(Grid& grid, ISubsetHandler& shOut,
									int inSubset, int bndSubset);
									
////////////////////////////////////////////////////////////////////////
//	AssignSubsetColors
///	assigns a different color to each subset
void AssignSubsetColors(ISubsetHandler& sh);



////////////////////////////////////////////////////////////////////////
///	Assigns all sides of elements of the given type to a separate subset
/**	All elements of type TElem::lower_dim_base_object are assigned to a
 * subset, depending on the subsets of neighbors of type TElem.
 * For each neighborhood constellation a separate subset index is chosen.
 * You may forbid assignment of sides, which already are in a subset through
 * the different preserve... parameters. Note that all have the
 * default parameter true.
 *
 * Valid parameters for TElem are: EdgeBase, Face, Volume
 *
 * \todo	Some performance improvements could be implemented, by checking
 * 			for the different preserve... states as early as possible.
 */
template <class TElem>
void AssignSidesToSubsets(ISubsetHandler& sh,
						bool preserveManifolds = true,
						bool preserveInterfaces = true,
						bool preserveInner = true);


////////////////////////////////////////////////////////////////////////
///	Assigns all elements of type TElem with subset index -1 to subset at index si
template <class TElem, class TSubsetHandler>
void AssignUnassignedElemsToSubset(TSubsetHandler& sh, int si);


////////////////////////////////////////////////////////////////////////
///	Adjust the grid so that it is ready for simulation with ug4
/**	For simulation it is crucial that all geometric objects of a
 * grid are assigned to a subset. Furthermore all lower dimensional
 * elements, which are only connected to higher dimensional elements
 * of one subset should be in the same subset as those elements.
 * Elements which form an interface between different subsets should
 * be in a separate subset.
 *
 * This method will erase all empty subsets before executing the
 * actual algorithm.
 *
 * You may forbid assignment of sides which already are in a subset through
 * the different preserve... parameters. Not that they all have true as
 * default parameter
 *
 * \todo	An additional parameter would be nice, with which one could
 * 			decide whether elements should only be added to existing
 * 			subsets (assignToExistingSubsets).
 */
template <class TSubsetHandler>
void AdjustSubsetsForSimulation(TSubsetHandler& sh,
								bool preserveManifolds = true,
								bool preserveInterfaces = true,
								bool preserveInner = true);

/**@}*/ // end of doxygen defgroup command
}//	end of namespace

////////////////////////////////////////////////
// include implementations of template methods.
#include "subset_util_impl.hpp"

#endif
