/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__LIB_GRID__SELECTION_UTIL__
#define __H__LIB_GRID__SELECTION_UTIL__

#include <stack>
#include "lib_grid/lg_base.h"
#include "common/ug_config.h"
#include "lib_grid/tools/selector_multi_grid.h"
#include "lib_grid/callbacks/basic_callbacks.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	selection util methods

/**
 * Several methods that ease selection-handling are grouped here.
 * \defgroup lib_grid_algorithms_selection_util selection util
 * \ingroup lib_grid_algorithms
 * @{
 */

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
///	calculates the center of selected objects
/**	This algorithm uses Grid::mark
 * The center is calculated by averaging the positions of all vertices, which
 * touch the selection.
 */
template <class TAAPosVRT>
bool CalculateCenter(typename TAAPosVRT::ValueType& centerOut,
					 Selector& sel, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
///	moves all vertices touching the selection by the specified offset.
/**	This algorithm uses Grid::mark
 */
template <class TAAPosVRT>
void TranslateSelection(Selector& sel, const typename TAAPosVRT::ValueType& offset,
						TAAPosVRT& aaPos);


////////////////////////////////////////////////////////////////////////
//	CollectVerticesTouchingSelection
///	Collects all vertices which are selected or which touch a selected element.
/**	This method uses Grid::mark
 * returns the number of collected vertices.
 */
UG_API
size_t CollectVerticesTouchingSelection(std::vector<Vertex*>& vrtsOut,
										ISelector& sel);

////////////////////////////////////////////////////////////////////////
//	EraseSelectedObjects
///	Erases selected objects from the associated grid
/**
 * TSelector has to either be of type Selector or MGSelector.
 */
template <class TSelector>
void EraseSelectedObjects(TSelector& sel);


////////////////////////////////////////////////////////////////////////
//	InvertSelection
///	Inverts the selection of the elements between begin and end.
/**
 * TSelector has to either be of type Selector or MGSelector.
 */
template <class TSelector, class TIterator>
void InvertSelection(TSelector& sel, TIterator begin, TIterator end);

////////////////////////////////////////////////////////////////////////
//	InvertSelection
///	Inverts the selection.
/**
 * TSelector has to either be of type Selector or MGSelector.
 */
template <class TSelector>
void InvertSelection(TSelector& sel);


////////////////////////////////////////////////////////////////////////
///	Selects all elements of type TElem, which touch an element between begin and end.
/**	To select all associated volumes of selected faces, call this method like this:
 * \code
 * SelectAssociated<Volume>(sel, sel.begin<Face>, sel.end<Face>());
 * \endcode
 */
template <class TElem, class TIterator>
void
SelectAssociated(ISelector& sel, TIterator begin, TIterator end,
				 ISelector::status_t status = ISelector::SELECTED);

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedVertices
///	selects all associated vertices of the elements between elemsBegin and elemsEnd
/**
 * TSelector has to feature a method select(TElemIterator::value_type&);
 *
 * TElemIterator has to be a stl-compatible iterator.
 * The underlying element-type has to be a pointer to a class that
 * features the following methods:
 * 
 * Vertex* vertex(int i);//returns the i-th vertex of the element.
 * uint num_vertices();//returns the number of vertices that the element holds.
 *
 * Valid classes are for example Edge, Face and Volume.
 *
 * Make sure that the elements only reference vertices that belong to the grid
 * at which the selector is registered.
 */
template <class TSelector, class TElemIterator>
void SelectAssociatedVertices(TSelector& sel, TElemIterator elemsBegin,
							  TElemIterator elemsEnd,
							  ISelector::status_t status = ISelector::SELECTED);
								
////////////////////////////////////////////////////////////////////////
//	SelectAssociatedEdges
///	selects all associated edges of the elements between elemsBegin and elemsEnd
/**
 * TSelector has to feature a method select(TElemIterator::value_type&);
 *
 * TElemIterator has to be a stl-compatible iterator.
 * The underlying element-type has to be a pointer to a class that
 * is supported by libGrid::CollectEdges(...)
 *
 * Valid classes are for example Face and Volume.
 *
 * Make sure that the elements only reference edges that belong to the grid
 * at which the selector is registered.
 */
template <class TSelector, class TElemIterator>
void SelectAssociatedEdges(TSelector& sel, TElemIterator elemsBegin,
						   TElemIterator elemsEnd,
						   ISelector::status_t status = ISelector::SELECTED);

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedFaces
///	selects all associated faces of the elements between elemsBegin and elemsEnd
/**
 * TSelector has to feature a method select(TElemIterator::value_type&);
 *
 * TElemIterator has to be a stl-compatible iterator.
 * The underlying element-type has to be a pointer to a class that
 * is supported by libGrid::CollectFaces(...)
 *
 * A valid classe is for example Volume.
 *
 * Make sure that the elements only reference faces that belong to the grid
 * at which the selector is registered.
 */
template <class TSelector, class TElemIterator>
void SelectAssociatedFaces(TSelector& sel, TElemIterator elemsBegin,
						   TElemIterator elemsEnd,
						   ISelector::status_t status = ISelector::SELECTED);

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedVolumes
///	selects all associated faces of the elements between elemsBegin and elemsEnd
/**
 * TSelector has to feature a method select(TElemIterator::value_type&);
 *
 * TElemIterator has to be a stl-compatible iterator.
 * The underlying element-type has to be a pointer to a class that
 * is supported by libGrid::CollectFaces(...)
 *
 * A valid classe is for example Volume.
 *
 * Make sure that the elements only reference faces that belong to the grid
 * at which the selector is registered.
 */
template <class TSelector, class TElemIterator>
void SelectAssociatedVolumes(TSelector& sel, TElemIterator elemsBegin,
						   TElemIterator elemsEnd,
						   ISelector::status_t status = ISelector::SELECTED);


//////////////////////////////////////////////////////////////////////////
/////	selects associated geometric objects of selected ones.
//UG_API
//void SelectAssociatedGridObjects(Selector& sel,
//							  ISelector::status_t status = ISelector::SELECTED);

////////////////////////////////////////////////////////////////////////
///	selects associated geometric objects of selected ones on each level.
template <class TSelector>
UG_API
void SelectAssociatedGridObjects(TSelector& sel,
							  ISelector::status_t status = ISelector::SELECTED);


////////////////////////////////////////////////////////////////////////
/// Selects all associated elements of lower dimensions
template <class TSelector>
UG_API
void CloseSelection (TSelector& sel);


////////////////////////////////////////////////////////////////////////
///	Assigns the selection state of selected elements to associated sides.
/**	If recursive is set to true, the method will recursively call itself, to
 * copy the state to sides of sides and so on.
 *
 * The new status of the side will be an or combination of the initial state
 * of the side and the states of adjacent elements.
 *
 * Valid types for TSelector are Selector and MGSelector.
 */
template <class TElem, class TSelector>
void AssignSelectionStateToSides(TSelector& sel, bool recursive);

////////////////////////////////////////////////////////////////////////
///	selects elements that lie on the associated grid's boundary
template <class TElemIterator>
void SelectBoundaryElements(ISelector& sel, TElemIterator elemsBegin,
						 TElemIterator elemsEnd);

////////////////////////////////////////////////////////////////////////
///	selects elements that do not lie on the associated grid's boundary
template <class TElemIterator>
void SelectInnerElements(ISelector& sel, TElemIterator elemsBegin,
						 TElemIterator elemsEnd);

////////////////////////////////////////////////////////////////////////
/// Selects edges which at which triangles meet in a large angle
template <class TEdgeIterator>
void SelectCreaseEdges(ISelector& sel, TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
						number minAngle, APosition aVrtPos,
						bool ignoreBoundaryEdges = true,
						ISelector::status_t state = ISelector::SELECTED);

////////////////////////////////////////////////////////////////////////
///	selects sides that are only adjacent to one of the given inner elements
/**
 * This algorithm uses Grid::mark.
 * selects the sides of the elements between begin and end
 * that are only adjacent to one of those elements.
 *
 * Edges that already are selected will stay selected, even if they are
 * inner edges.
 *
 * Please note that only existing sides are checked.
 */
template <class TIter>
void SelectAreaBoundary(ISelector& sel, const TIter begin, const TIter end);

////////////////////////////////////////////////////////////////////////
///	Selects elements which are adjacent to higher dimensional elements of different subsets
/**	Please note, that this method does not select boundary segments.
 * If regardSelectedNbrsOnly is set to true (default false), then only selected
 * neighbors are checked for different interfaces.*/
template <class TIter>
void SelectInterfaceElements(ISelector& sel, ISubsetHandler& sh,
							 const TIter begin, const TIter end,
							 bool regardSelectedNbrsOnly = false);

////////////////////////////////////////////////////////////////////////
/// selects all elements of the given type in the given subset
/**	If you want to deselect elements, pass ISelector::DESELECT to the status
 * argument.*/
template <class TElem>
void SelectSubsetElements(ISelector& sel, ISubsetHandler& sh, int subsetIndex,
						  ISelector::status_t status = ISelector::SELECTED);

////////////////////////////////////////////////////////////////////////
///	extends the selection to neighbours of selected elements.
/**	
 * This algorithm uses Grid::mark.
 *
 * Extension is performed extSize times.
 *
 * \todo: Performance can be improved. See implementation.
 */
template <class TSelector>
UG_API
void ExtendSelection(TSelector& sel, size_t extSize,
					 ISelector::status_t status = ISelector::SELECTED);


////////////////////////////////////////////////////////////////////////
///	extends the selection to neighbours of selected elements in the given direction.
/**	
 * This algorithm uses Grid::mark.
 *
 * Extension is performed extSize times.
 *
 * Only elements which can be reached in the given direction (centers are compared)
 * are selected.
 *
 * by setting minAngle to 0 and maxAngle to 10 you may e.g. extend only in
 * roughly the given direction.
 *
 * By setting minAngle to 89 and maxAngle to 91 you may e.g. extend only orthogonal
 * to the given direction
 *
 * \todo: Performance can be improved. See implementation.
 */
template <class TSelector, class TAAPos>
UG_API
void ExtendSelectionInDirection(
        TSelector& sel,
        size_t extSize,
        const typename TAAPos::ValueType& dir,
        number minAngle,
        number maxAngle,
      	const TAAPos& aaPos,
		ISelector::status_t status = ISelector::SELECTED);


////////////////////////////////////////////////////////////////////////
///	Extends the selection around selected objects until selected sides are reached.
/**	Selects all elements of the given base-type which are reachable
 * without traversing a bounding side.
 * Whether a side is bounding can be specified through the cbRegionBoundary callback
 *
 * - 	Use e.g. IsSelected(sel) to indicate that all selected edges are region boundaries.
 * -	Use e.g. IsNotInSubset(sh, -1) to indicate, that all faces which are in a subset
 * 		should be considered as region boundaryies.
 *
 * Those callbacks are declared in lib_grid/algorithms/callback_util.h"
 *
 * Valid types for TGeomBaseObj are Edge, Face and Volume.
 */
template <class TGeomObj>
void SelectionFill(Selector& sel,
			   	   typename Grid::traits<typename TGeomObj::side>::callback cbRegionBoundary);

////////////////////////////////////////////////////////////////////////
/// Selects the region which contains the given point
/**	Selects all elements of the given base-type which are reachable
 * without traversing a bounding side. The method starts at the element which
 * contains the given point.
 *
 * Whether a side is bounding can be specified through the cbRegionBoundary callback
 *
 * - 	Use e.g. IsSelected(sel) to indicate that all selected sides are region boundaries.
 * -	Use e.g. IsNotInSubset(sh, -1) to indicate, that all sides which are in a subset
 * 		should be considered as region boundaryies.
 *
 * Those callbacks are declared in lib_grid/algorithms/callback_util.h"
 *
 * The method tries to find the start element using a brute force approach with
 * runtime O(n). If no element which contains the given point was found, the method
 * returns false.
 *
 * Valid types for TGeomBaseObj are Edge, Face and Volume.
 */
template <class TGeomObj, class TAAPos>
bool SelectRegion(Selector& sel, const typename TAAPos::ValueType& p, TAAPos& aaPos,
			   	  typename Grid::traits<typename TGeomObj::side>::callback cbRegionBoundary);

////////////////////////////////////////////////////////////////////////
//	SelectAssociatedGenealogy
///	Selects the complete genealogy of all selected elements.
/**
 * After the method returns the selection in msel is complete
 * regarding the property that the parent of each selected
 * element is selected, too.
 *
 * If selectAssociatedElements is set to true, the selection will be
 * complete regarding the property that for each selected element all
 * elements of lower dimension are selected, too.
 *
 * It is assumed that the given selector references a valid multi-grid.
 * That means a grid whose elements only refence elements on the same
 * level.
 */
UG_API
void SelectAssociatedGenealogy(MGSelector& msel, bool selectAssociatedElements);

////////////////////////////////////////////////////////////////////////
///	Selects all edges that face a given direction
template <class TAAPos>
void SelectEdgesByDirection(
				Selector& sel,
				TAAPos& aaPos,
				const vector3& dir,
				number minDeviationAngle,
				number maxDeviationAngle,
				bool selectFlipped);


////////////////////////////////////////////////////////////////////////
//	SelectSmoothEdgePath
///	selects for each selected edge all edges that can be reached by a smooth path.
/**
 * \param sel: Selector
 * \param thresholdDegree defines the maximal degree at which the angle
 *			between two edges is regarded as smooth. Between 0 and 180.
 * \param normalWeight	defines how much the surrounding of an edge shall be
 *			taken into account when determining the path.
 *			1: The surrounding is the most important feature
 *			2: The direction of the edge is the most important feature
 *			]0,1[: Gradually shifting between the both extreme settings
 * \param stopAtSelVrts: If set to true, the edge-path will stop at selected
 *						vertices.
 * \param aPos: Position attachment
 * \todo: replace aPos by an template AttachmentAccessor TAAPosVrt.*/
UG_API
void SelectSmoothEdgePath(Selector& sel, number thresholdDegree,
							number normalWeight,
							bool stopAtSelVrts = true,
							APosition& aPos = aPosition);


////////////////////////////////////////////////////////////////////////
//	SelectShortPolychains
///	Selects regular polygonal chains which are shorter than a given threshold
template <class TAAPos>
void SelectShortPolychains(ISelector& sel, number maxLength, bool closedChainsOnly,
						   TAAPos aaPos);

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionVertices
///	selects all vertices that are only connected to selected elements.
/**	
 * This algorithm uses Grid::mark.
 *
 * TSelector either has to be of type Selector or MGSelector.
 * Only elements of higher dimension are regarded.
 */
template <class TSelector>
void SelectInnerSelectionVertices(TSelector& sel);

////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionEdges
///	selects all edges that are only connected to selected elements.
/**
 * This algorithm uses Grid::mark.
 *
 * TSelector either has to be of type Selector or MGSelector.
 * Only elements of higher dimension are regarded.
 */
template <class TSelector>
void SelectInnerSelectionEdges(TSelector& sel);

////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionFaces
///	selects all faces that are only connected to selected elements.
/**
 * This algorithm uses Grid::mark.
 *
 * TSelector either has to be of type Selector or MGSelector.
 * Only elements of higher dimension are regarded.
 */
template <class TSelector>
void SelectInnerSelectionFaces(TSelector& sel);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	DeselectBoundarySelectionVertices
///	deselects all vertices that are connected to unselected elements.
/**	
 * TSelector either has to be of type Selector or MGSelector.
 * Only elements of higher dimension are regarded.
 */
template <class TSelector>
void DeselectBoundarySelectionVertices(TSelector& sel);

////////////////////////////////////////////////////////////////////////
//	DeselectBoundarySelectionEdges
///	deselects all edges that are connected to unselected elements.
/**
 * TSelector either has to be of type Selector or MGSelector.
 * Only elements of higher dimension are regarded.
 */
template <class TSelector>
void DeselectBoundarySelectionEdges(TSelector& sel);

////////////////////////////////////////////////////////////////////////
//	DeselectBoundarySelectionFaces
///	deselects all faces that are connected to unselected elements.
/**
 * TSelector either has to be of type Selector or MGSelector.
 * Only elements of higher dimension are regarded.
 */
template <class TSelector>
void DeselectBoundarySelectionFaces(TSelector& sel);


////////////////////////////////////////////////////////////////////////
//	SelectLinkedElements
///	Repeatedly traverses sides of selected elements and selects associated elements
/**	The method extends the selection as long as it finds new candidates.
 * Given a selected element, the method first checks the elements sides whether
 * they may be traversed by calling cbIsTraversable on each. If a side may be traversed,
 * cbIsSelectable is called on the adjacent unselected elements to the given side.
 * If the callback returns true, the corresponding elements are selected and considered
 * as new starting points for the search. By default all elements and all sides are
 * considered to be selectable / traversable.
 */
template <class TElem>
void SelectLinkedElements(ISelector& sel,
		  typename Grid::traits<TElem>::callback
		  	  cbIsSelectable = ConsiderAll(),
		  typename Grid::traits<typename TElem::side>::callback
		  	  cbIsTraversable = ConsiderAll());

							   
////////////////////////////////////////////////////////////////////////
//	SelectLinkedFlatFaces
///	Extends the selection of faces to all neighbouring faces that have a similar normal.
/**
 * \param sel: Selector
 * \param maxDeviationAngle: in degree. Maximal angle between normals of faces considered as flat.
 * \param ignoreOrientation (default false): If true, neighboured faces which have
 * 							inverted orientation are traversed anyways.
 * \param aPos: Position attachment
 */
UG_API
void SelectLinkedFlatFaces(Selector& sel, number maxDeviationAngle,
						   bool ignoreOrientation = false,
						   bool stopAtSelectedEdges = false,
						   APosition& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
//	SelectLinkedFlatFaces
///	Extends the selection of faces to all neighbouring faces that have a similar normal.
/**
 * In contrast to SelectLinkedFlatFaces this method also traverses
 * degenerated faces over their non-degenerated sides.
 *
 * \param sel: Selector
 * \param maxDeviationAngle: in degree. Maximal angle between normals of faces considered as flat.
 * \param ignoreOrientation (default false): If true, neighboured faces which have
 * 							inverted orientation are traversed anyways.
 * \param aPos: Position attachment
 */
UG_API
void SelectLinkedFlatAndDegeneratedFaces(Selector& sel,
										 number maxDeviationAngle,
										 bool ignoreOrientation = false,
										 bool stopAtSelectedEdges = false,
										 number degThreshold = SMALL,
						   	   	   	     APosition& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
//	FaceArea
/**
 * Returns the area sum of convex faces selected by ISelector sel
 *
 * \param sel: Selector
 * \param aaPos: Position attachment
 *
 * \return \c area sum of convex faces
 */
template <class TAAPosVRT>
UG_API
number FaceArea(ISelector& sel, TAAPosVRT& aaPos);


///	Returns indices of selected elements in ascending order
/** \warning 	this method has linear complexity and should be avoided unless
 *				required to exchange selection states with external programs
 *				or scripts.
 * \{ */
template <class elem_t>
void GetSelectedElementIndices (const ISelector& sel, std::vector<int>& indsOut);

void GetSelectedElementIndices (const ISelector& sel,
                                std::vector<size_t>& vrtIndsOut,
                                std::vector<size_t>& edgeIndsOut,
                                std::vector<size_t>& faceIndsOut,
                                std::vector<size_t>& volIndsOut);
/** \} */


///	Selects elements with the specified indices
/** \warning 	this method has linear complexity and should be avoided
 *				unless required to exchange selection states with external
 *				programs or scripts.
 * \{ */
template <class elem_t>
void SelectElementsByIndex (ISelector& sel, const std::vector<size_t>& inds);

void SelectElementsByIndex (ISelector& sel,
                            const std::vector<size_t>& vrtInds,
                            const std::vector<size_t>& edgeInds,
                            const std::vector<size_t>& faceInds,
                            const std::vector<size_t>& volInds);

/** \} */

/**@}*/ // end of doxygen defgroup command
}// end of namespace


////////////////////////////////
//	include implementation
#include "selection_util_impl.hpp"

#endif
