// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 31.05.2011 (m,d,y)

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
 * 						sequences of integers: {{numInds, ind1, ind2, ...}, ...}.
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

}//	end of namespace
}//	end of namespace

#endif
