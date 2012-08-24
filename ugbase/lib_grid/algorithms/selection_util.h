// created by Sebastian Reiter
// y09 m11 d06
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SELECTION_UTIL__
#define __H__LIB_GRID__SELECTION_UTIL__

#include <stack>
#include "lib_grid/lg_base.h"

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
 */
template <class TAAPosVRT>
bool CalculateCenter(typename TAAPosVRT::ValueType& centerOut,
					 Selector& sel, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
//	CollectVerticesTouchingSelection
///	Collects all vertices which are selected or which touch a selected element.
/**	This method uses Grid::mark
 * returns the number of collected vertices.
 */
size_t CollectVerticesTouchingSelection(std::vector<VertexBase*>& vrtsOut,
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
 * VertexBase* vertex(int i);//returns the i-th vertex of the element.
 * uint num_vertices();//returns the number of vertices that the element holds.
 *
 * Valid classes are for example EdgeBase, Face and Volume.
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


////////////////////////////////////////////////////////////////////////
///	selects associated geometric objects of selected ones.
void SelectAssociatedGeometricObjects(Selector& sel,
							  ISelector::status_t status = ISelector::SELECTED);

////////////////////////////////////////////////////////////////////////
///	selects associated geometric objects of selected ones on each level.
void SelectAssociatedGeometricObjects(MGSelector& msel,
							  ISelector::status_t status = ISelector::SELECTED);


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
void ExtendSelection(Selector& sel, size_t extSize);

////////////////////////////////////////////////////////////////////////
///	Extends the selection around selected objects until selected sides are reached.
/**	Selects all elements of the given base-type which are reachable
 * without traversing a selected side.
 *
 * Valid types for TGeomBaseObj are EdgeBase, Face and Volume.
 */
template <class TGeomObj>
void SelectionFill(Selector& sel);

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
void SelectAssociatedGenealogy(MGSelector& msel, bool selectAssociatedElements);

////////////////////////////////////////////////////////////////////////
//	SelectSmoothEdgePath
///	selects for each selected edge all edges that can be reached by a smooth path.
/**
 * \param sel: Selector
 * \param thresholdDegree defines the maximal degree at which the angle
 *			between two edges is regarded as smooth. Between 0 and 180.
 * \param stopAtSelVrts: If set to true, the edge-path will stop at selected
 *						vertices.
 * \param aPos: Position attachment
 * \todo: replace aPos by an template AttachmentAccessor TAAPosVrt.*/ 
void SelectSmoothEdgePath(Selector& sel, number thresholdDegree,
							bool stopAtSelVrts = true,
							APosition& aPos = aPosition);


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
//	SelectLinkedBoundaryElements
///	Selects linked boundary elements, starting from currently selected ones.
/** If stopAtSelectedSides is enabled, selected sides will not be
 * traversed during the search.
 *
 * Valid parameters for TElem are EdgeBase and Face.
 */
template <class TElem, class TAPos>
void SelectLinkedBoundaryElements(ISelector& sel, TAPos& aPos,
								  bool stopAtSelectedSides = true);
							   
////////////////////////////////////////////////////////////////////////
//	SelectLinkedFlatFaces
///	Extends the selection of faces to all neighbouring faces that have a similar normal.
/**
 * \param sel: Selector
 * \param maxDeviationAngle: in degree. Maximal angle between normals of faces considered as flat.
 * \param traverseFlipped (default false): If true, neighboured faces which have
 * 							inverted orientation are traversed anyways.
 * \param aPos: Position attachment
 */
void SelectLinkedFlatFaces(Selector& sel, number maxDeviationAngle,
						   bool traverseFlipped = false,
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
 * \param traverseFlipped (default false): If true, neighboured faces which have
 * 							inverted orientation are traversed anyways.
 * \param aPos: Position attachment
 */
void SelectLinkedFlatAndDegeneratedFaces(Selector& sel,
										 number maxDeviationAngle,
										 bool traverseFlipped = false,
										 bool stopAtSelectedEdges = false,
										 number degThreshold = SMALL,
						   	   	   	     APosition& aPos = aPosition);

/**@}*/ // end of doxygen defgroup command

}// end of namespace


////////////////////////////////
//	include implementation
#include "selection_util_impl.hpp"

#endif
