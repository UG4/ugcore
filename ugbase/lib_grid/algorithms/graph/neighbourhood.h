//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d18

#ifndef __H__LIB_GRID__NEIGHBOURHOOD__
#define __H__LIB_GRID__NEIGHBOURHOOD__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

enum NeighbourhoodType
{
	NHT_DEFAULT = 0,
	NHT_VERTEX_NEIGHBOURS = 1,
	NHT_EDGE_NEIGHBOURS = 1<<1,
	NHT_FACE_NEIGHBOURS = 1<<2,
	NHT_VOLUME_NEIGHBOURS = 1<<3,
	NHT_ALL = NHT_VERTEX_NEIGHBOURS
			| NHT_EDGE_NEIGHBOURS
			| NHT_FACE_NEIGHBOURS
			| NHT_VOLUME_NEIGHBOURS
};

////////////////////////////////////////////////////////////////////////
///	Collects all vertices that are connected by edges in the specified subset.
/**
 * This method uses Grid::mark.
 *
 * \param vNeighboursOut 	neighbour vertices
 * \param grid 				Grid
 * \param sh				Subset Handler
 * \param vrt				vertex
 * \param subsetIndex		subset index
 * \param nbhType: Accepts or-combinations of any NeighbourhoodType
 *					enumerated constants.
 */
void CollectSubsetNeighbours(std::vector<VertexBase*>& vNeighboursOut,
							Grid& grid, SubsetHandler& sh,
							VertexBase* vrt, int subsetIndex,
							uint nbhType = NHT_ALL);

////////////////////////////////////////////////////////////////////////
//	CollectNeighbours
///	collects all edges that are connected to the given one.
/**
 * This method uses Grid::mark.
 *
 * if nbhType != NHT_VERTEX_NEIGHBOURS, then vNeighboursOut will be empty.
 * This parameter is featured for compatibility reasons with the other
 * CollectNeighbours methods.
 */
void CollectNeighbours(std::vector<EdgeBase*>& vNeighboursOut, EdgeBase* e,
					   Grid& grid, NeighbourhoodType nbhType = NHT_VERTEX_NEIGHBOURS);


////////////////////////////////////////////////////////////////////////
//	CollectNeighbours
///	collects all faces that are connected to the given one.
/**
 * This method uses Grid::mark.
 *
 * Using nbhType, you can choose which neighbourhood you want to receive:
 *	- NHT_EDGE_NEIGHBOURS (default): All faces that share an edge with the given one are regarded as neighbours.
 *	- NHT_VERTEX_NEIGHBOURS: All faces that share a vertex with the given one are regarded as neighbours.
 */
void CollectNeighbours(std::vector<Face*>& vNeighboursOut, Face* f,
					   Grid& grid, NeighbourhoodType nbhType = NHT_EDGE_NEIGHBOURS);

////////////////////////////////////////////////////////////////////////
//	CollectNeighbours
///	collects all volumes that are connected to the given one.
/**
 * This method uses Grid::mark.
 *
 * Using nbhType, you can choose which neighbourhood you want to receive:
 *	- NHT_FACE_NEIGHBOURS (default): All volumes that share a face with the given one are regarded as neighbours.
 *	- NHT_EDGE_NEIGHBOURS: All volumes that share an edge with the given one are regarded as neighbours.
 *	- NHT_VERTEX_NEIGHBOURS: All volumes that share a vertex with the given one are regarded as neighbours.
 */
void CollectNeighbours(std::vector<Volume*>& vNeighboursOut, Volume* v,
					   Grid& grid, NeighbourhoodType nbhType = NHT_FACE_NEIGHBOURS);

////////////////////////////////////////////////////////////////////////
///	Collects all neighbours in a given neighbourhood of a vertex
/**	This algorithm uses Grid::mark
 */				   
void CollectNeighbourhood(std::vector<Face*>& facesOut, Grid& grid,
						  VertexBase* vrt, size_t range,
						  bool clearContainer = true);

}//	end of namespace

#endif
