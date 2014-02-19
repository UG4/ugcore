//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d18

#ifndef __H__LIB_GRID__NEIGHBORHOOD__
#define __H__LIB_GRID__NEIGHBORHOOD__

#include <vector>
#include "grid.h"

namespace ug
{

/**
 * Methods to access the neighborhood of geometric objects
 * \defgroup lib_grid_algorithms_neighborhood_util neighborhood util
 * \ingroup lib_grid_algorithms
 * @{
 */

///	Constants to specify a neighborhood
/**	Use arbitrary or combinations. */
enum NeighborhoodType
{
	NHT_DEFAULT = 0,
	NHT_VERTEX_NEIGHBORS = 1,
	NHT_EDGE_NEIGHBORS = 1<<1,
	NHT_FACE_NEIGHBORS = 1<<2,
	NHT_VOLUME_NEIGHBORS = 1<<3,
	NHT_ALL = NHT_VERTEX_NEIGHBORS
			| NHT_EDGE_NEIGHBORS
			| NHT_FACE_NEIGHBORS
			| NHT_VOLUME_NEIGHBORS
};

////////////////////////////////////////////////////////////////////////
///	Collects all vertices that are connected by elements of the specified type
/**
 * This method uses Grid::mark.
 *
 * You may specify the types of objects, which are regarded as connecting objects
 * through or-combinations of constants enumerated in NeighborhoodType.
 *
 * Through the consider... callbacks, you may specify whether a specific
 * edge/face/volume shall be considered, when collecting neighbors of a vertex.
 * By default, all geometric objects are considered.
 *
 * \param vNeighboursOut 	neighbor vertices
 * \param grid 				Grid
 * \param sh				Subset Handler
 * \param vrt				vertex
 * \param subsetIndex		subset index
 * \param nbhType: Accepts or-combinations of any NeighborhoodType
 *					enumerated constants.
 */
void CollectNeighbors(std::vector<Vertex*>& vNeighborsOut,
						Grid& grid, Vertex* vrt, uint nbhType = NHT_EDGE_NEIGHBORS,
						Grid::edge_traits::callback considerEdge	= Grid::edge_traits::cb_consider_all,
						Grid::face_traits::callback considerFace	= Grid::face_traits::cb_consider_all,
						Grid::volume_traits::callback considerVol	= Grid::volume_traits::cb_consider_all);

////////////////////////////////////////////////////////////////////////
//	CollectNeighbors
///	collects all edges that are connected to the given one.
/**
 * This method uses Grid::mark.
 *
 * if nbhType != NHT_VERTEX_NEIGHBORS, then vNeighborsOut will be empty.
 * This parameter is featured for compatibility reasons with the other
 * CollectNeighbors methods.
 */
void CollectNeighbors(std::vector<EdgeBase*>& vNeighborsOut, EdgeBase* e,
					   Grid& grid, NeighborhoodType nbhType = NHT_VERTEX_NEIGHBORS);


////////////////////////////////////////////////////////////////////////
//	CollectNeighbors
///	collects all faces that are connected to the given one.
/**
 * This method uses Grid::mark.
 *
 * Using nbhType, you can choose which neighborhood you want to receive:
 *	- NHT_EDGE_NEIGHBORS (default): All faces that share an edge with the given one are regarded as neighbours.
 *	- NHT_VERTEX_NEIGHBORS: All faces that share a vertex with the given one are regarded as neighbours.
 */
void CollectNeighbors(std::vector<Face*>& vNeighborsOut, Face* f,
					   Grid& grid, NeighborhoodType nbhType = NHT_EDGE_NEIGHBORS);

////////////////////////////////////////////////////////////////////////
//	CollectNeighbors
///	collects all volumes that are connected to the given one.
/**
 * This method uses Grid::mark.
 *
 * Using nbhType, you can choose which neighborhood you want to receive:
 *	- NHT_FACE_NEIGHBORS (default): All volumes that share a face with the given one are regarded as neighbors.
 *	- NHT_EDGE_NEIGHBORS: All volumes that share an edge with the given one are regarded as neighbors.
 *	- NHT_VERTEX_NEIGHBORS: All volumes that share a vertex with the given one are regarded as neighbors.
 */
void CollectNeighbors(std::vector<Volume*>& vNeighboursOut, Volume* v,
					   Grid& grid, NeighborhoodType nbhType = NHT_FACE_NEIGHBORS);

////////////////////////////////////////////////////////////////////////
///	Collects all neighbors in a given neighborhood of a vertex
/**	This algorithm uses Grid::mark
 */				   
void CollectNeighborhood(std::vector<Face*>& facesOut, Grid& grid,
						  Vertex* vrt, size_t range,
						  bool clearContainer = true);

/**@}*/ // end of doxygen defgroup command

}//	end of namespace

#endif
