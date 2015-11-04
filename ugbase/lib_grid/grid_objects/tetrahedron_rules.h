// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 27.05.2011 (m,d,y)

#ifndef __H__UG__tetrahedron_rules__
#define __H__UG__tetrahedron_rules__

#include "common/math/ugmath.h"

namespace ug{
namespace tet_rules
{

////////////////////////////////////////////////////////////////////////////////
//	REFINEMENT RULE

/// identification of refinement rule to be used
enum GlobalRefinementRule{
	STANDARD,
	HYBRID_TET_OCT
};

////////////////////////////////////////////////////////////////////////////////
//	LOOKUP TABLES

const int NUM_VERTICES	= 4;
const int NUM_EDGES		= 6;
const int NUM_FACES		= 4;
const int NUM_TRIS		= 4;
const int NUM_QUADS		= 0;
const int MAX_NUM_INDS_OUT = 64;//todo: this is just an estimate!

///	the local vertex indices of the given edge
const int EDGE_VRT_INDS[][2] = {	{0, 1}, {1, 2}, {2, 0},
									{0, 3}, {1, 3}, {2,3}};

///	the local vertex indices of the given face
const int FACE_VRT_INDS[][4] = {	{0, 1, 2, -1}, {1, 3, 2, -1},
									{0, 2, 3, -1}, {0, 3, 1, -1}};

///	for each edge the local index of the opposed edge
const int OPPOSED_EDGE[] = {5, 3, 4, 1, 2, 0};

/** for each vertex, a pair containing the object type (0: vrt, 1: edge, 2: face)
 * and an index into the associated array, which describe the object which lies
 * on the opposite side of the tetrahedron, to a given vertex.*/
const int OPPOSED_OBJECT[][NUM_VERTICES] = {{2, 1}, {2, 2}, {2, 3}, {2, 0}};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	NOTE: The lists below are all generated automatically

///	returns the j-th edge of the i-th face
const int FACE_EDGE_INDS[4][4] = 	{{0, 1, 2, -1}, {4, 5, 1, -1},
									 {2, 5, 3, -1}, {3, 4, 0, -1}};

///	tells whether the i-th face contains the j-th edge
const int FACE_CONTAINS_EDGE[][6] = {{1, 1, 1, 0, 0, 0}, {0, 1, 0, 0, 1, 1},
									 {0, 0, 1, 1, 0, 1}, {1, 0, 0, 1, 1, 0}};

///	Associates the index of the connecting edge with each tuple of vertices.
/**	Use two vertex indices to index into this table to retrieve the index
 * of their connecting edge.
 */
const int EDGE_FROM_VRTS[4][4] =	{{-1, 0, 2, 3}, {0, -1, 1, 4},
									 {2, 1, -1, 5}, {3, 4, 5, -1}};

///	Associates the index of the connecting face with each triple of vertices.
/**	Use three vertex indices to index into this table to retrieve the index
 * of their connecting face.
 */
const int FACE_FROM_VRTS[4][4][4] =
			{{{-1, -1, -1, -1}, {-1, -1, 0, 3}, {-1, 0, -1, 2}, {-1, 3, 2, -1}},
			 {{-1, -1, 0, 3}, {-1, -1, -1, -1}, {0, -1, -1, 1}, {3, -1, 1, -1}},
			 {{-1, 0, -1, 2}, {0, -1, -1, 1}, {-1, -1, -1, -1}, {2, 1, -1, -1}},
			 {{-1, 3, 2, -1}, {3, -1, 1, -1}, {2, 1, -1, -1}, {-1, -1, -1, -1}}};

///	given two edges, the table returns the face, which contains both (or -1)
const int FACE_FROM_EDGES[][6] =	{{0, 0, 0, 3, 3, -1}, {0, 0, 0, -1, 1, 1},
									 {0, 0, 0, 2, -1, 2}, {3, -1, 2, 2, 3, 2},
									 {3, 1, -1, 3, 1, 1}, {-1, 1, 2, 2, 1, 1}};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/**	returns an array of integers, which contains the indices of the objects
 * resulting from the refinement of a tetrahedron.
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
 * 						NUM_VERTICES. If an inner vertex has to be created, it is
 * 						referenced by NUM_VERTICES + NUM_EDGES + NUM_FACES.
 *
 * \param newEdgeVrts	Array of size NUM_EDGES, which has to contain 1 for each
 * 						edge, which shall be refined and 0 for each edge, which
 * 						won't be refined.
 *
 * \param newCenterOut	If the refinement-rule requires a center vertex, then
 * 						this parameter will be set to true. If not, it is set to
 * 						false.
 *
 * \param corners		(optional) List of the four corner positions of the
 * 						tetrahedron. If it is specified, it is used during full
 * 						refinement (all edges marked), to determine the best
 * 						diagonal along which inner tetrahedrons are created.
 * 						Corners are only considered during full refinement and are
 * 						thus irrelevant during recursive refinement of other elements.
 *
 * \returns	the number of entries written to newIndsOut or 0, if the refinement
 * 			could not be performed.
 */
int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut,
		   vector3* corners = NULL);



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void SetRefinementRule(GlobalRefinementRule refRule);
GlobalRefinementRule GetRefinementRule();

}//	end of namespace tet_rules
}//	end of namespace ug

#endif
