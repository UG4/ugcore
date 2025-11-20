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

#ifndef __H__UG__pyramid_rules__
#define __H__UG__pyramid_rules__

//	only required for dummy-parameter ug::vector3*
#include "common/math/ugmath_types.h"

namespace ug{
namespace pyra_rules
{

////////////////////////////////////////////////////////////////////////////////
//	LOOKUP TABLES

constexpr int NUM_VERTICES = 5;
constexpr int NUM_EDGES = 8;
constexpr int NUM_FACES = 5;
constexpr int NUM_TRIS = 4;
constexpr int NUM_QUADS = 1;
constexpr int MAX_NUM_INDS_OUT = 128;//todo: this is just an estimate!
constexpr int MAX_NUM_CONVERT_TO_TETS_INDS_OUT = 10;
constexpr int MAX_NUM_COLLAPSE_INDS_OUT = 5;

///	the local vertex indices of the given edge
const int EDGE_VRT_INDS[][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0},
									{4, 0}, {4, 1}, {4, 2}, {4, 3}};

///	the local vertex indices of the given face
const int FACE_VRT_INDS[][4] = { {0, 1, 2, 3},	{0, 4, 1, -1},
									{1, 4, 2, -1}, {2, 4, 3, -1},
									{0, 3, 4, -1}};

///	the pyramids top
constexpr int TOP_VERTEX = 4;


/** for each vertex, a pair containing the object type (0: vrt, 1: edge, 2: face)
 * and an index into the associated array, which describe the object which lies
 * on the opposite side of the pyramid, to a given vertex.*/
const int OPPOSED_OBJECT[][NUM_VERTICES] = {{1, 6}, {1, 7}, {1, 4}, {1, 5}, {2, 0}};

constexpr int NUM_BOTTOM_EDGES = 4;
constexpr int BOTTOM_EDGE_INDS[NUM_BOTTOM_EDGES] = {0, 1, 2, 3};

constexpr int NUM_TOP_EDGES = 4;
constexpr int TOP_EDGE_INDS[NUM_BOTTOM_EDGES] = {4, 5, 6, 7};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	NOTE: The lists below are all generated automatically

///	returns the j-th edge of the i-th face
const int FACE_EDGE_INDS[5][4] = 	{{0, 1, 2, 3}, {4, 5, 0, -1}, {5, 6, 1, -1},
									 {6, 7, 2, -1}, {3, 7, 4, -1}};

///	tells whether the i-th face contains the j-th edge
const int FACE_CONTAINS_EDGE[][8] =
						{{1, 1, 1, 1, 0, 0, 0, 0}, {1, 0, 0, 0, 1, 1, 0, 0},
						 {0, 1, 0, 0, 0, 1, 1, 0}, {0, 0, 1, 0, 0, 0, 1, 1},
						 {0, 0, 0, 1, 1, 0, 0, 1}};

///	Associates the index of the connecting edge with each tuple of vertices.
/**	Use two vertex indices to index into this table to retrieve the index
 * of their connecting edge.
 */
const int EDGE_FROM_VRTS[5][5] =	{{-1, 0, -1, 3, 4}, {0, -1, 1, -1, 5},
									 {-1, 1, -1, 2, 6}, {3, -1, 2, -1, 7},
									 {4, 5, 6, 7, -1}};

///	Associates the index of the connecting face with each triple of vertices.
/**	Use three vertex indices to index into this table to retrieve the index
 * of their connecting face.
 */
const int FACE_FROM_VRTS[5][5][5] =
				{{{-1, -1, -1, -1, -1}, {-1, -1, 0, 0, 1}, {-1, 0, -1, 0, -1},
				  {-1, 0, 0, -1, 4}, {-1, 1, -1, 4, -1}},
				 {{-1, -1, 0, 0, 1}, {-1, -1, -1, -1, -1}, {0, -1, -1, 0, 2},
				  {0, -1, 0, -1, -1}, {1, -1, 2, -1, -1}},
				 {{-1, 0, -1, 0, -1}, {0, -1, -1, 0, 2}, {-1, -1, -1, -1, -1},
				  {0, 0, -1, -1, 3}, {-1, 2, -1, 3, -1}},
				 {{-1, 0, 0, -1, 4}, {0, -1, 0, -1, -1}, {0, 0, -1, -1, 3},
				  {-1, -1, -1, -1, -1}, {4, -1, 3, -1, -1}},
				 {{-1, 1, -1, 4, -1}, {1, -1, 2, -1, -1}, {-1, 2, -1, 3, -1},
				  {4, -1, 3, -1, -1}, {-1, -1, -1, -1, -1}}};

///	given two edges, the table returns the face, which contains both (or -1)
const int FACE_FROM_EDGES[][8] =
					{{0, 0, 0, 0, 1, 1, -1, -1}, {0, 0, 0, 0, -1, 2, 2, -1},
					 {0, 0, 0, 0, -1, -1, 3, 3}, {0, 0, 0, 0, 4, -1, -1, 4},
					 {1, -1, -1, 4, 1, 1, -1, 4}, {1, 2, -1, -1, 1, 1, 2, -1},
					 {-1, 2, 3, -1, -1, 2, 2, 3}, {-1, -1, 3, 4, 4, -1, 3, 3}};


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
 * \param isSnapPoint	(optional) An array of size NUM_VERTICES. If all entries
 *						are set to 'false' the behaviour is the same as if the
 *						array wasn't specified.
 *						If a corner of a quadrilateral is a snap-point and if
 *						edges of that quadrilateral are refined, then only new
 *						edges connected to the snap-point are introduced.
 *						Note that only special snap-point constellations
 *						are supported.
 *
 * \returns	the number of entries written to newIndsOut or 0, if the refinement
 * 			could not be performed.
 */
int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut,
		   vector3* corners = nullptr, bool* isSnapPoint = nullptr);


///	returns true if the specified edgeMarks would lead to a regular refinement
/**	A regular refinement leads to new elements which are all similar to the
 * original element. I.e. which are of the same type and which have similar
 * angles.
 *
 * \note	this method does not perform refinement. Use 'Refine' with the
 *			specified edges instead.
 *
 * \param edgeMarks	If the i-th edge shall be refined, the expression
 *					'edgeMarks & (1<<i) != 0' has to be true. You can
 *					specify multiple refine-edges using or-combinations. E.g.,
 *					'edgeMarks = (1<<i) | (1<<j)' would indicate that the
 *					i-th and the j-th edge shall be refined.*/
bool IsRegularRefRule(const int edgeMarks);


/// fills an array of integers describing tetrahedra that shall replace the pyramid
/**	The method requires a compare-operator that defines a strict (global) ordering on the
 * vertices of pyramid. Note that this operator should return consistent results
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
 *						corners of the pyramid. The method shall return true if
 *						the vertex at the first corner-index shall be considered
 *						smaller than the vertex at the second corner index.
 *						If multiple prisms shall be converted to tetrahedra, it
 *						is important that the given ordering is global, i.e., if
 *						a pair of vertices is present in two connected prisms,
 *						the 'cmp' operator has to return the same value.
 *
 * \returns	the number of entries written to newIndsOut or 0, if the refinement
 * 			could not be performed.*/
template <typename TCmp>
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

}//	end of namespace pyra_rules
}//	end of namespace ug


////////////////////////////////////////
//	include implementation
#include "pyramid_rules_impl.h"

#endif
