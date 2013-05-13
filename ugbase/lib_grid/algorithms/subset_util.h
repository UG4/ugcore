//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__SUBSET_UTIL__
#define __H__LIB_GRID__SUBSET_UTIL__

#include <vector>
#include "lib_grid/lg_base.h"
#include "common/ug_config.h"

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
///	returns the first subset, which does not contain any elements at all
UG_API
int GetFirstFreeSubset(const ISubsetHandler& sh);

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
///	Assigns all elements of the given grid to the given subset
/**	Make sure, that the given subset handler operates on the given grid.
 */
UG_API
void AssignGridToSubset(Grid& g, ISubsetHandler& sh, int subsetInd);

////////////////////////////////////////////////////////////////////////
///	Assigns all selected elements to the specified subset
/**	Make sure that the specified subset handler and the specified selector
 * operate on the same grid.
 */
UG_API
void AssignSelectionToSubset(ISelector& sel, ISubsetHandler& sh, int subsetInd);

////////////////////////////////////////////////////////////////////////
//	AssignFaceInterfaceEdgesToSubsets
///	assigns edges which belong to no subset and are adjacent to faces of different subsets to new subsets.
/**
 * the edges of each interface between two subset are put into
 * a seperate subset, starting at GetMaxSubsetIndex<EdgeBase>(sh) + 1.
 */
UG_API
void AssignFaceInterfaceEdgesToSubsets(Grid& grid, SubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
//	AssignVolumeInterfaceFacesToSubsets
///	assigns faces which belong to no subset and are adjacent to volumes of different subsets to new subsets.
/**
 * the faces of each interface between two subset are put into
 * a seperate subset, starting at GetMaxSubsetIndex<Face>(sh) + 1.
 */
UG_API
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
UG_API
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
UG_API
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
UG_API
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
UG_API
void SeparateFaceSubsetsByNormal(Grid& grid, SubsetHandler& sh,
								std::vector<vector3> vNormals,
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
UG_API
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
UG_API
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
void SeparateSubsetsByLowerDimSubsets(Grid& grid, SubsetHandler& sh,
									  bool appendAtEnd = false);

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
										Selector& sel, bool appendAtEnd = false);

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
					bool appendAtEnd,
					boost::function<bool (typename TElem::lower_dim_base_object*)>
						cbIsSeparator);

////////////////////////////////////////////////////////////////////////
//	AssignInnerAndBoundarySubsets
///	assigns objects to subsets depending on whether they are inner or boundary objects.
UG_API
void AssignInnerAndBoundarySubsets(Grid& grid, ISubsetHandler& shOut,
									int inSubset, int bndSubset);


////////////////////////////////////////////////////////////////////////
///	Returns an rgb vector (values ranging from 0 to 1), with the i-th default color.
//todo: Move this method to common/util or something like that.
UG_API
vector3 GetColorFromStandardPalette(int index);

////////////////////////////////////////////////////////////////////////
//	AssignSubsetColors
///	assigns a different color to each subset
UG_API
void AssignSubsetColors(ISubsetHandler& sh);



////////////////////////////////////////////////////////////////////////
///	Assigns all sides of elements of the given type to a separate subset
/**	All elements of type TElem::lower_dim_base_object, which are located in
 * subset -1, are assigned to a subset, depending on the subsets of neighbors
 * of type TElem. For each neighborhood constellation a separate subset index
 * is chosen.
 *
 * Valid parameters for TElem are: EdgeBase, Face, Volume
 *
 * You may specify a selector, which indicates which elements to check, when
 * looking for neighbors in different subsets. Please note, that the selector
 * is only used, if a selected side is encountered.
 */
template <class TElem>
void AssignSidesToSubsets(ISubsetHandler& sh, ISelector* psel = NULL);

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
 * You may specify, that existing subsets shall not be touched. This is
 * a good idea, if you already specified some boundaries etc.
 * However, you should in this case make sure, that all elements which are
 * not assigned to a subset intentionally, are assigned to subset -1.
 */
template <class TSubsetHandler>
void AdjustSubsetsForSimulation(TSubsetHandler& sh,
								bool preserveExistingSubsets);

////////////////////////////////////////////////////////////////////////
///	Assigns subset depending on the element type
UG_API
void AssignSubsetsByElementType(ISubsetHandler& sh);

////////////////////////////////////////////////////////////////////////
//  FaceArea
///	Returns the area sum of convex faces given by subset index and level.
/**
 * \param sh subset handler
 * \param si subset index
 * \param lvl grid level
 * \param aaPos position attachment
 *
 * \return \c number area sum of convex faces
 */
template <class TAAPosVRT>
UG_API
number FaceArea(ISubsetHandler& sh, int si, size_t lvl, TAAPosVRT& aaPos);

/**@}*/ // end of doxygen defgroup command
}//	end of namespace

////////////////////////////////////////////////
// include implementations of template methods.
#include "subset_util_impl.hpp"

#endif
