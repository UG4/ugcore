/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__NEIGHBORHOOD__
#define __H__LIB_GRID__NEIGHBORHOOD__

#include <vector>
#include "grid.h"
#include "lib_grid/callbacks/basic_callbacks.h"

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
enum class NeighborhoodType : uint8_t
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
using NeighborhoodType_t = uint8_t;

constexpr NeighborhoodType operator & (NeighborhoodType lhs, NeighborhoodType rhs) noexcept {
	return static_cast<NeighborhoodType>(
		static_cast<NeighborhoodType_t>(lhs) &
		static_cast<NeighborhoodType_t>(rhs)
	);
}

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
 * \param vNeighborsOut 	neighbor vertices
 * \param grid 				Grid
 * \param vrt				vertex
 * \param nbhType: Accepts or-combinations of any NeighborhoodType
 *					enumerated constants.
 * \param considerEdge todo
 * \param considerFace todo
 * \param considerVol todo
 */
void CollectNeighbors(std::vector<Vertex*>& vNeighborsOut,
						Grid& grid, Vertex* vrt, NeighborhoodType nbhType = NeighborhoodType::NHT_EDGE_NEIGHBORS,
						Grid::edge_traits::callback considerEdge = ConsiderAll(),
						Grid::face_traits::callback considerFace = ConsiderAll(),
						Grid::volume_traits::callback considerVol = ConsiderAll());

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
void CollectNeighbors(std::vector<Edge*>& vNeighborsOut, Edge* e,
					   Grid& grid, NeighborhoodType nbhType = NeighborhoodType::NHT_VERTEX_NEIGHBORS);


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
					   Grid& grid, NeighborhoodType nbhType = NeighborhoodType::NHT_EDGE_NEIGHBORS);

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
					   Grid& grid, NeighborhoodType nbhType = NeighborhoodType::NHT_FACE_NEIGHBORS);

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
