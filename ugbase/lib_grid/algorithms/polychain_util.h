// created by Sebastian Reiter
// y10 m12 d13
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__POLYCHAIN_UTIL__
#define __H__LIB_GRID__POLYCHAIN_UTIL__

#include <utility>
#include "lib_grid/lg_base.h"

namespace ug
{

/**
 * Methods that allow to interprete parts of a grid as a polygonal chain.
 * \defgroup lib_grid_algorithms_polychain_util polygonal-chain util
 * \ingroup lib_grid_algorithms
 * @{
 */

enum PolyChainTypes
{
	PCT_UNKNOWN = 0,
	PCT_CLOSED = 1,
	PCT_OPEN = 1 << 1,
	PCT_SEPARATED = 1 << 2,
	PCT_IRREGULAR = 1 << 3,
	PCT_EMPTY = 1 << 4
};

////////////////////////////////////////////////////////////////////////
///	returns an or combination of constants enumerated in PolyChainTypes.
template <class TEdgeIterator>
size_t
GetPolyChainType(Grid& grid, TEdgeIterator edgesBegin,
				  TEdgeIterator edgesEnd,
				  Grid::edge_traits::callback cbEdgeIsInPolyChain);
				  
////////////////////////////////////////////////////////////////////////
///	Returns the start-vertex and start-edge of a polygonal chain.
/**	The edges between iteratorBegin and iteratorEnd are considered to be
 *	edges of a polygonal chain. The algorithm then searches for a vertex that
 *	lies on the boundary of that chain. If a boundary vertex is found which
 *	is the first vertex in its associated chain-edge, it is preferred to
 *	a vertex, which is the second vertex of its associated chain-edge.
 *	
 *	Returned is an std::pair. The first component points to the vertex
 *	and the second to its associated edge.
 *
 *	If the polygonal chain is closed, the the edge at edgesBegin and its
 *	first vertex is returned.
 *
 *	The given callback cbConsiderEdge has to return true for all edges between
 *	edgesBegin and edgesEnd and false for all other edges.
 *
 *	If there are no edges between the given iterators, std::pair(NULL, NULL) is returned.
 */
template <class TEdgeIterator>
std::pair<VertexBase*, EdgeBase*>
GetFirstSectionOfPolyChain(Grid& grid, TEdgeIterator edgesBegin,
						  TEdgeIterator edgesEnd,
						  Grid::edge_traits::callback cbEdgeIsInPolyChain);
						  
////////////////////////////////////////////////////////////////////////
///	returns the next section in a polygonal chain.
/**	lastSection has to contain a vertex and an associated edge, which is
 *	part of the polygonal chain. You can retrieve such a section by a call
 *	to ug::GetFirstVertexOfPolyChain or to ug::GetNextSectionInPolyChain.
 *
 *	If no more sections can be found, std::pair(NULL, NULL) is returned.
 *
 *	Be careful with closed polychains! The algorithm will always find a
 *	next section in this case.
 *
 *	If the polychain is irregular (more than 2 chain-edges meet at one vertex)
 *	the algorithm still terminates. However assumtions on the outcome should not
 *	be made.
 */
std::pair<VertexBase*, EdgeBase*>
GetNextSectionOfPolyChain(Grid& grid, std::pair<VertexBase*, EdgeBase*> lastSection,
						  Grid::edge_traits::callback cbEdgeIsInPolyChain);

////////////////////////////////////////////////////////////////////////
///	Makes sure that the polychain at srcIndex is regular and not separated.
/**	This algorithm uses Grid::mark
 *
 * Makes sure that the polychain at srcIndex is regular and not separated.
 * This is achieved by assigning all problematic edges to targetIndex.
 *
 * Returns true if the chain was splitted and false if not.
 */
bool SplitIrregularPolyChain(SubsetHandler& sh, int srcIndex, int targetIndex);


////////////////////////////////////////////////////////////////////////
///	given a list of edges, this method collects associated vertices in a polychain
/**	This method uses Grid::mark.
 * edges between edgesBegin and edgesEnd should build a closed regular polygon.
 * \todo	add support for open chains.*/
template <class TEdgeIter>
bool CreatePolyChain(std::vector<VertexBase*>& polyChainOut, Grid& grid,
					TEdgeIter edgesBegin, TEdgeIter edgesEnd);

/**@}*/ // end of doxygen defgroup command

}//	end of namespace


////////////////////////////////
//	include implementation
#include "polychain_util_impl.hpp"

#endif
