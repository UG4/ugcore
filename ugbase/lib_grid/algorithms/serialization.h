// created by Sebastian Reiter
// y09 m11 d04
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SERIALIZATION__
#define __H__LIB_GRID__SERIALIZATION__

#include <iostream>
#include "lib_grid/lg_base.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GRID

////////////////////////////////////////////////////////////////////////
///	Writes a part of the grids elements to a binary-stream.
/**
 * The passed GeometricObjectCollection goc may only reference
 * elements of the given grid. It is important, that the goc
 * is complete - that means that all referenced vertices are
 * contained in the goc.
 *
 * If you pack several different parts of your grid, you should use
 * this method, since it is faster than calling SerializeGridElements
 * without the attachment.
 *
 * the integer attachment aInt is used during this method to store
 * an index in each vertex of the goc. The initial content of
 * the referenced attachment is ignored.
 *
 * The caller is responsible to attach the aIntVRT attachment to the 
 * to the vertices of the grid before calling this method.
 * The caller is also responsible to detach aIntVRT from the grids
 * vertices when it is no longer required.
 *
 * After termination the attachment holds the indices at which
 * the respcetive vertices are stored in the pack.
 */
bool SerializeGridElements(Grid& grid, GeometricObjectCollection goc,
						   AInt& aIntVRT, std::ostream& out);

////////////////////////////////////////////////////////////////////////
///	Writes all grid elements into a binary-stream.
bool SerializeGridElements(Grid& grid, std::ostream& out);

////////////////////////////////////////////////////////////////////////
///	Writes a part of the grids elements to a binary-stream.
/**
 * The passed GeometricObjectCollection goc may only reference
 * elements of the given grid. It is important, that the goc
 * is complete - that means that all referenced vertices are
 * contained in the goc.
 *
 * If you're planning to serialize multiple parts of one grid, you
 * should consider to use the full-featured serialization method.
 */
bool SerializeGridElements(Grid& grid, GeometricObjectCollection goc,
						   std::ostream& out);

////////////////////////////////////////////////////////////////////////
///	Creates grid elements from a binary stream
bool DeserializeGridElements(Grid& grid, std::istream& in);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	MULTI-GRID

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
///	writes a part of the elements of a MultiGrid to a binary stream.
/**
 * The passed MultiLevelGeometricObjectCollection goc may only
 * reference elements of the given grid. It is important, that the goc
 * is complete - that means that all referenced vertices are
 * contained in the goc.
 * The goc has also to be complete in regard to the multi-grid hierarchy.
 * This means that for each element of the goc, the parent has to be
 * also part of the goc.
 *
 * If you pack several different parts of your grid, you should use
 * this method, since it is faster than calling SerializeGridElements
 * without the attachment.
 *
 * the integer attachment aInt is used during this method to store
 * an index in each element of the goc. The initial content of
 * the referenced attachment is ignored.
 *
 * The caller is responsible to attach the aIntVRT attachment to the 
 * to the vertices of the grid before calling this method.
 * The caller is also responsible to detach aIntVRT from the grids
 * vertices when it is no longer required.
 *
 * After termination the attachments hold the indices that were
 * assigned to the respective elements - starting from 0 for each
 * element type.
 */
bool SerializeMultiGridElements(MultiGrid& mg,
								MultiLevelGeometricObjectCollection mgoc,
								AInt& aIntVRT, AInt& aIntEDGE,
								AInt& aIntFACE, AInt& aIntVOL,
								std::ostream& out);

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
///	writes a part of the elements of a MultiGrid to a binary stream.
/**
 * The passed GeometricObjectCollection goc may only reference
 * elements of the given grid. It is important, that the goc
 * is complete - that means that all referenced vertices are
 * contained in the goc.
 * The goc has also to be complete in regard to the multi-grid hierarchy.
 * This means that for each element of the goc, the parent has to be
 * also part of the goc.
 *
 * If you're planning to serialize multiple parts of one grid, you
 * should consider to use the full-featured serialization method.
 */
bool SerializeMultiGridElements(MultiGrid& mg,
								MultiLevelGeometricObjectCollection mgoc,
								std::ostream& out);

////////////////////////////////////////////////////////////////////////
//	SerializeMultiGridElements
///	writes the elements of a MultiGrid to a binary stream.
bool SerializeMultiGridElements(MultiGrid& mg,
								std::ostream& out);
								
////////////////////////////////////////////////////////////////////////
///	Creates multi-grid elements from a binary stream
/**
 * If you pass a pointer to a std::vector using pvVrts, pvEdges,
 * pvFaces or pvVolumes, those vectors will contain the elements
 * of the grid in the order they were read.
 * Specifying those vectors does not lead to a performance loss.
 */
bool DeserializeMultiGridElements(MultiGrid& mg, std::istream& in,
									std::vector<VertexBase*>* pvVrts = NULL,
									std::vector<EdgeBase*>* pvEdges = NULL,
									std::vector<Face*>* pvFaces = NULL,
									std::vector<Volume*>* pvVols = NULL);



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	ATTACHMENTS

////////////////////////////////////////////////////////////////////////
//	copies attached values to a binary stream.
/**
 * copies attached values of the grids elements of the given type
 * to the binary stream.
 *
 * Make sure that attachment is attached to the specified elements.
 */
template <class TElem, class TAttachment>
bool SerializeAttachment(Grid& grid, TAttachment& attachment,
						 std::ostream& out);

////////////////////////////////////////////////////////////////////////
///	copies attached values to a binary stream.
/**
 * copies attached values of the elements between iterBegin and iterEnd
 * to the binary stream.
 *
 * Make sure that attachment is attached to the specified elements.
 */
template <class TElem, class TAttachment>
bool SerializeAttachment(Grid& grid, TAttachment& attachment,
						 typename geometry_traits<TElem>::iterator iterBegin,
						 typename geometry_traits<TElem>::iterator iterEnd,
						 std::ostream& out);

////////////////////////////////////////////////////////////////////////
///	copies attached values from a binary stream
/**
 * copies values from the given binary stream to the given attachment of
 * elements between iterBegin and iterEnd.
 * If attachment was not attached to the grid, then it will be attached
 * automatically.
 */
template <class TElem, class TAttachment>
bool DeserializeAttachment(Grid& grid, TAttachment& attachment,
						 std::istream& in);

////////////////////////////////////////////////////////////////////////
///	copies attached values from a binary stream
/**
 * copies values from the given binary stream to the given attachment of
 * elements between iterBegin and iterEnd.
 * If attachment was not attached to the grid, then it will be attached
 * automatically.
 */
template <class TElem, class TAttachment>
bool DeserializeAttachment(Grid& grid, TAttachment& attachment,
						 typename geometry_traits<TElem>::iterator iterBegin,
						 typename geometry_traits<TElem>::iterator iterEnd,
						 std::istream& in);


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	SUBSET-HANDLER

////////////////////////////////////////////////////////////////////////
///	writes the subset-indices of all elements in the goc to a stream.
bool SerializeSubsetHandler(Grid& grid, SubsetHandler& sh,
							GeometricObjectCollection goc,
							std::ostream& out);

////////////////////////////////////////////////////////////////////////
///	writes the subset-indices of all elements in the grid to a stream.
bool SerializeSubsetHandler(Grid& grid, SubsetHandler& sh,
							std::ostream& out);

////////////////////////////////////////////////////////////////////////
///	assigns subset-indices to all elements in the goc from a stream.
/**	One has to be very careful that the given goc only contains
 * the elements that were passed to the serialization routine.
 * Problems could be caused by automatic element creation.
 * consider to set grid.set_option(GRIDOPT_NONE) before loading
 * the grid.
 */
bool DeserializeSubsetHandler(Grid& grid, SubsetHandler& sh,
							GeometricObjectCollection goc,
							std::istream& in);

////////////////////////////////////////////////////////////////////////
///	assigns subset-indices to all elements in the grid from a stream.
/**	One has to be very careful that the given grid only contains
 * the elements that were passed to the serialization routine.
 * Problems could be caused by automatic element creation.
 * consider to set grid.set_option(GRIDOPT_NONE) before loading
 * the grid.
 */
bool DeserializeSubsetHandler(Grid& grid, SubsetHandler& sh,
							std::istream& in);

/*
bool SerializeSelector(Grid& grid, Selector& sel, std::ostream& out);

bool DeserializeSelector(Grid& grid, Selector& sel, std::istream& in);
*/
}

////////////////////////////////
//	include implementation
#include "serialization_impl.hpp"

#endif
