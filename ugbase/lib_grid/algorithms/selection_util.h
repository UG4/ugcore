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
								TElemIterator elemsEnd);
								
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
void SelectAssociatedEdges(TSelector& sel,
								TElemIterator elemsBegin,
								TElemIterator elemsEnd);

////////////////////////////////////////////////////////////////////////
///	selects edges that are only adjacent to one of the given faces
/**
 * This algorithm uses Grid::mark.
 * selects the edges of the faces between facesBegin and facesEnd
 * that are only adjacent to one of those faces.
 *
 * Edges that already are selected will stay selected, even if they are
 * inner edges.
 *
 * Please note that only existing edges are checked.
 */
void SelectAreaBoundaryEdges(ISelector& sel, FaceIterator facesBegin,
								  FaceIterator facesEnd);

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
void SelectAssociatedFaces(TSelector& sel,
								TElemIterator elemsBegin,
								TElemIterator elemsEnd);

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
 * \param thresholdDegree defines the maximal degree at which the angle
 *			between two edges is regarded as smooth. Between 0 and 180.
 * \param stopAtSelVrts: If set to true, the edge-path will stop at selected
 *						vertices.
 * \todo: replace aPos by an template AttachmentAccessor TAAPosVrt.*/ 
void SelectSmoothEdgePath(Selector& sel, number thresholdDegree,
							bool stopAtSelVrts = true,
							APosition& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionVertices
///	selects all vertices that are only connected to selected elements.
/**	
 * This algorithm uses Grid::mark.
 */
void SelectInnerSelectionVertices(Selector& sel);

////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionEdges
///	selects all edges that are only connected to selected elements.
/**
 * This algorithm uses Grid::mark.
 */
void SelectInnerSelectionEdges(Selector& sel);

////////////////////////////////////////////////////////////////////////
//	SelectInnerSelectionFaces
///	selects all faces that are only connected to selected elements.
/**
 * This algorithm uses Grid::mark.
 */
void SelectInnerSelectionFaces(Selector& sel);

}// end of namespace


////////////////////////////////
//	include implementation
#include "selection_util_impl.hpp"

#endif
