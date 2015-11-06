/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__prism_rules__
#define __H__UG__prism_rules__

//	only required for dummy-parameter ug::vector3*
#include "common/math/ugmath_types.h"

namespace ug{
namespace prism_rules{

////////////////////////////////////////////////////////////////////////////////
//	LOOKUP TABLES

const int NUM_VERTICES	= 6;
const int NUM_EDGES		= 9;
const int NUM_FACES		= 5;
const int NUM_TRIS		= 2;
const int NUM_QUADS		= 3;
const int MAX_NUM_INDS_OUT = 128;//todo: this is just an estimate!
const int MAX_NUM_CONVERT_TO_TETS_INDS_OUT = 15;
const int MAX_NUM_COLLAPSE_INDS_OUT = 6;

///	the local vertex indices of the given edge
const int EDGE_VRT_INDS[][2] = {	{0, 1}, {1, 2}, {2, 0},
									{0, 3}, {1, 4}, {2, 5},
									{3, 4}, {4, 5}, {5, 3}};

///	the local vertex indices of the given face
const int FACE_VRT_INDS[][4] = {	{0, 1, 2, -1},	{0, 3, 4, 1},
//									{1, 4, 5, 2},	{0, 2, 5, 3},
									{1, 4, 5, 2},	{2, 5, 3, 0},
									{3, 5, 4, -1}};

///	the index of the face opposed to the given one. -1 if no face is opposed
const int OPPOSED_FACE[NUM_FACES] = {4, -1, -1, -1, 0};


/** for each vertex, a pair containing the object type (0: vrt, 1: edge, 2: face)
 * and an index into the associated array, which describe the object which lies
 * on the opposite side of the prism, to a given vertex.*/
const int OPPOSED_OBJECT[][NUM_VERTICES] = {{1, 7}, {1, 8}, {1, 6}, {1, 1}, {1, 2}, {1, 0}};


////////////////////////////////////////////////////////////////////////////////
//	SOME HELPER TABLES
const int TOP_FACE =	4;
const int BOTTOM_FACE =	0;

const int IS_BOTTOM_EDGE[9] =	{1, 1, 1, 0, 0, 0, 0, 0, 0};
const int IS_SIDE_EDGE[9] =		{0, 0, 0, 1, 1, 1, 0, 0, 0};
const int IS_TOP_EDGE[9] =		{0, 0, 0, 0, 0, 0, 1, 1, 1};

const int TRIS[2] =		{0, 4};
const int QUADS[3] =	{1, 2, 3};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	NOTE: The lists below are all generated automatically

///	returns the j-th edge of the i-th face
const int FACE_EDGE_INDS[5][4] =	{{0, 1, 2, -1}, {3, 6, 4, 0}, {4, 7, 5, 1},
									 {5, 8, 3, 2}, {8, 7, 6, -1}};
//const int FACE_EDGE_INDS[5][4] =	{{0, 1, 2, -1}, {3, 6, 4, 0}, {4, 7, 5, 1},
//									 {2, 5, 8, 3}, {8, 7, 6, -1}};

///	tells whether the i-th face contains the j-th edge
const int FACE_CONTAINS_EDGE[][9] =
				{{1, 1, 1, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 1, 1, 0, 1, 0, 0},
				 {0, 1, 0, 0, 1, 1, 0, 1, 0}, {0, 0, 1, 1, 0, 1, 0, 0, 1},
				 {0, 0, 0, 0, 0, 0, 1, 1, 1}};

///	Associates the index of the connecting edge with each tuple of vertices.
/**	Use two vertex indices to index into this table to retrieve the index
 * of their connecting edge.
 */
const int EDGE_FROM_VRTS[6][6] =
				{{-1, 0, 2, 3, -1, -1}, {0, -1, 1, -1, 4, -1},
				 {2, 1, -1, -1, -1, 5}, {3, -1, -1, -1, 6, 8},
				 {-1, 4, -1, 6, -1, 7}, {-1, -1, 5, 8, 7, -1}};

///	Associates the index of the connecting face with each triple of vertices.
/**	Use three vertex indices to index into this table to retrieve the index
 * of their connecting face.
 */
const int FACE_FROM_VRTS[6][6][6] =
							{{{-1, -1, -1, -1, -1, -1}, {-1, -1, 0, 1, 1, -1},
							  {-1, 0, -1, 3, -1, 3}, {-1, 1, 3, -1, 1, 3},
							  {-1, 1, -1, 1, -1, -1}, {-1, -1, 3, 3, -1, -1}},
							 {{-1, -1, 0, 1, 1, -1}, {-1, -1, -1, -1, -1, -1},
							  {0, -1, -1, -1, 2, 2}, {1, -1, -1, -1, 1, -1},
							  {1, -1, 2, 1, -1, 2}, {-1, -1, 2, -1, 2, -1}},
							 {{-1, 0, -1, 3, -1, 3}, {0, -1, -1, -1, 2, 2},
							  {-1, -1, -1, -1, -1, -1}, {3, -1, -1, -1, -1, 3},
							  {-1, 2, -1, -1, -1, 2}, {3, 2, -1, 3, 2, -1}},
							 {{-1, 1, 3, -1, 1, 3}, {1, -1, -1, -1, 1, -1},
							  {3, -1, -1, -1, -1, 3}, {-1, -1, -1, -1, -1, -1},
							  {1, 1, -1, -1, -1, 4}, {3, -1, 3, -1, 4, -1}},
							 {{-1, 1, -1, 1, -1, -1}, {1, -1, 2, 1, -1, 2},
							  {-1, 2, -1, -1, -1, 2}, {1, 1, -1, -1, -1, 4},
							  {-1, -1, -1, -1, -1, -1}, {-1, 2, 2, 4, -1, -1}},
							 {{-1, -1, 3, 3, -1, -1}, {-1, -1, 2, -1, 2, -1},
							  {3, 2, -1, 3, 2, -1}, {3, -1, 3, -1, 4, -1},
							  {-1, 2, 2, 4, -1, -1}, {-1, -1, -1, -1, -1, -1}}};

///	given two edges, the table returns the face, which contains both (or -1)
const int FACE_FROM_EDGES[][9] =
				{{0, 0, 0, 1, 1, -1, 1, -1, -1}, {0, 0, 0, -1, 2, 2, -1, 2, -1},
				 {0, 0, 0, 3, -1, 3, -1, -1, 3}, {1, -1, 3, 1, 1, 3, 1, -1, 3},
				 {1, 2, -1, 1, 1, 2, 1, 2, -1}, {-1, 2, 3, 3, 2, 2, -1, 2, 3},
				 {1, -1, -1, 1, 1, -1, 1, 4, 4}, {-1, 2, -1, -1, 2, 2, 4, 2, 4},
				 {-1, -1, 3, 3, -1, 3, 4, 4, 3}};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**	returns an array of integers, which contains the indices of the objects
 * resulting from the refinement of a pyramid.
 *
 *
 * \param newIndsOut	Array which has to be of size MAX_NUM_INDS_OUT.
 * 						When the algorithm is done, the array will contain
 * 						sequences of integers: {{gridObjectID, ind1, ind2, ...}, ...}.
 *						gridObjectID is a constant enumerated in GridObjectID and
 *						describes the type of the grid-object that is
 *						built from the following set of corner indices.
 * 						Old vertices are referenced by their local index. Vertices
 * 						created on an edge are indexed by the index of the edge +
 * 						NUM_VERTICES.
 * 						Vertices created on a face are referenced by
 * 						NUM_VERTICES + NUM_EDGES + index_of_face.
 * 						If an inner vertex has to be created, it is referenced
 * 						by NUM_VERTICES + NUM_EDGES + NUM_FACES (in this case
 * 						newCenterOut is set to true).
 *
 * \param newEdgeVrts	Array of size NUM_EDGES, which has to contain 1 for each
 * 						edge, which shall be refined and 0 for each edge, which
 * 						won't be refined.
 *
 * \param newCenterOut	If the refinement-rule requires a center vertex, then
 * 						this parameter will be set to true. If not, it is set to
 * 						false.
 *
 * \param corners		Ignored.
 *
 * \returns	the number of entries written to newIndsOut or 0, if the refinement
 * 			could not be performed.
 */
int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut,
		   vector3* corners = NULL);


/// fills an array of integers describing tetrahedra that shall replace the prism
/**	The method requires a compare-operator that defines a strict (global) ordering on the
 * vertices of the prism. Note that this operator should return consistent results
 * for all vertices in a given grid. The ordering is used to decide along which
 * diagonal each quadrilateral is split. Each new diagonal will start at the
 * smallest vertex of the corresponding quadrilateral, regarding the given ordering.
 *
 * The specified compare function (or compare operator) 'bool cmp (int i0, int i1)'
 * will be called with two local corner indices and has to return true, if vertex
 * at the first corner shall be considered smaller (globally) than the vertex at
 * the second corner.
 *
 * The idea and implementation follows:
 * Dompierre et al., "How to Subdivide Pyramids, Prisms and Hexahedra into Tetrahedra"
 *
 * \param newIndsOut	Array which has to be of size MAX_NUM_CONVERT_TO_TETS_INDS_OUT.
 * 						When the algorithm is done, the array will contain
 * 						sequences of integers: {{gridObjectID, ind1, ind2, ...}, ...}.
 *						gridObjectID is a constant enumerated in GridObjectID and
 *						describes the type of the grid-object that is
 *						built from the following set of corner indices.
 *
 * \param cmp			A function object that induces a strict ordering on the
 *						corners of the prism. The method shall return true if
 *						the vertex at the first corner-index shall be considered
 *						smaller than the vertex at the second corner index.
 *						If multiple prisms shall be converted to tetrahedra, it
 *						is important that the given ordering is global, i.e., if
 *						a pair of vertices is present in two connected prisms,
 *						the 'cmp' operator has to return the same value.
 *
 * \returns	the number of entries written to newIndsOut or 0, if the refinement
 * 			could not be performed.*/
template <class TCmp>
int ConvertToTetrahedra(int* newIndsOut, TCmp cmp);


///	Creates new volume elements that result from collapsing the edge between v0 and v1 into v0
/**	Note that the returned array may be empty.
 *
 * \param newIndsOut	Array which has to be of size MAX_NUM_COLLAPSE_INDS_OUT.
 * 						When the algorithm is done, the array will contain
 * 						sequences of integers: {{gridObjectID, ind1, ind2, ...}, ...}.
 *						gridObjectID is a constant enumerated in GridObjectID and
 *						describes the type of the grid-object that is
 *						built from the following set of corner indices.
 *
 * \returns	the size of the resulting element list 'newIndsOut'.*/
int CollapseEdge(int* newIndsOut, int v0, int v1);

}//	end of namespace
}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "prism_rules_impl.h"

#endif
