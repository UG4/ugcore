/*
 * orientation.cpp
 *
 *  Created on: 25.07.2013
 *      Author: andreasvogel
 */

#include "orientation.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug{


///////////////////////////////////////////////////////////////////////////////
//	Lagrange Offsets
///////////////////////////////////////////////////////////////////////////////

/*
 * Lagrange DoF Orientation of an Edge:
 * If DoFs are assigned to a lower-dimensional edge and we have a degree higher
 * than 2 (i.e. more than one DoF on the edge) orientation is required to
 * ensure continuity of the shape functions. This means, that each element
 * that has the edge as a subelement, must number the dofs on the edge equally
 * in global numbering.
 *
 * The idea is as follows: We induce a global ordering of dofs on the edge by
 * using the vertices of the edge itself. We define, that dofs are always
 * assigned in a line from the vertex 0 to the vertex 1.
 * Now, in the local ordering of dofs on the reference element, the edge may
 * have been a different numbering for the corners.
 * Thus, we have to distinguish two case:
 * a) Orientation matches: we can simply use the usual offset numbering
 * b) Orientation mismatches: we have to use the reverse order as offset numbering
 */
bool OrientationMatches(const EdgeVertices& e1, const EdgeVertices& e2)
{
	return e1.vertex(0) == e2.vertex(0);
}

void ComputeOrientationOffsetLagrange(std::vector<size_t>& vOrientOffset,
                                      EdgeDescriptor& ed, Edge* edge, const size_t p)
{
	vOrientOffset.reserve(p-1);
	vOrientOffset.clear();

//	the standard orientation is from co0 -> co1.
	if(OrientationMatches(ed, *edge))
	{
		for(size_t i = 0; i < p-1; ++i)
			vOrientOffset.push_back(i);
	}
//	... and for reverse order
	else
	{
		for(int i = ((int)p) - 2; i >= 0; --i)
			vOrientOffset.push_back(i);
	}
}

void ComputeOrientationOffsetLagrange(std::vector<size_t>& vOrientOffset,
                                      Face* face, Edge* edge, size_t nrEdge,
                                      const size_t p)
{
	EdgeDescriptor ed;
	face->edge_desc(nrEdge, ed);
	ComputeOrientationOffsetLagrange(vOrientOffset, ed, edge, p);
}

void ComputeOrientationOffsetLagrange(std::vector<size_t>& vOrientOffset,
                                      Volume* vol, Edge* edge, size_t nrEdge,
                                      const size_t p)
{
	EdgeDescriptor ed;
	vol->edge_desc(nrEdge, ed);
	ComputeOrientationOffsetLagrange(vOrientOffset, ed, edge, p);
}

/*
 * Lagrange DoF Orientation of a Face:
 * If DoFs are assigned to a lower-dimensional face and we have a degree higher
 * than 2 (i.e. more than one DoF on the face) orientation is required to
 * ensure continuity of the shape functions. This means, that each element
 * that has the face as a subelement, must number the dofs on the face equally
 * in global numbering.
 *
 * The idea of the ordering is as follows:
 * DoFs are always assigned to the face in a natural order, i.e. in an order such
 * that the numbering of the face itself gives the ordering. Now, give a 3d
 * element with that face, the reference elements provides us with a numbering
 * of the vertices of the face. This numbering must not match the natural ordering
 * as present in the face itself.
 * Now find the vertex-id (id0) of the face that matches the natural vertex 0 of
 * the face itself, and the vertex-id (id1) of the natural vertex 1.
 * On the natural face we have a situation like this:
 *
 *  	*						*--------*			^ j
 *  	|  \					|		 |			|
 *  	|    \				    |	     |			|
 * 		|      \				|		 |			|-----> i
 *  	*------ *				*--------*
 *  	id0     id1 			id0		id1
 *
 * We define that the DoFs on the face are always numbered in x direction first,
 * continuing in the next row in y direction, numbering in x again and
 * continuing in y, etc. E.g. this gives (showing only inner dofs):
 *
 *  	*						*-------*			^ j
 *  	|5 \					| 6	7 8	|			|
 *  	|3 4 \					| 3 4 5	|			|
 * 		|0 1 2 \				| 0 1 2 |			|-----> i
 *  	*-------*				*-------*
 *  	id0     id1 			id0		id1
 *        p = 5						p = 4
 *
 * Now all rotations and mirroring can appear. This are resolved constructing
 * a mapping.
 */

static void MapLagrangeMultiIndexQuad(std::vector<size_t>& vOrientOffset,
                                      const int id0, bool sameOrientation,
                                      const size_t pOuter)
{
//	in the inner, the number of dofs is as if it would be an element of order p-2.
	const size_t p = pOuter-2;

//	resize array
	vOrientOffset.clear();
	vOrientOffset.reserve((p+1)*(p+1));

//	loop mapped indices as required by rotated face
	for(size_t mj = 0; mj <= p; ++mj){
		for(size_t mi = 0; mi <= p; ++mi){

		//	get corresponding multiindex in "natural" numbering
			size_t i,j;
			switch(id0)
			{
				case 0: i = mi;   j = mj;   break;
				case 1: i = mj;   j = p-mi; break;
				case 2: i = p-mi; j = p-mj; break;
				case 3: i = p-mj; j = mi;   break;
				default: UG_THROW("Orientation Quad: Corner "<<id0<<" invalid.");
			}
			if(!sameOrientation) std::swap(i, j);

		//	linearize index
			const size_t naturalIndex = (p+1) * j + i;

		//	set mapping
			vOrientOffset.push_back(naturalIndex);
		}
	}
};

static void MapLagrangeMultiIndexTriangle(std::vector<size_t>& vOrientOffset,
                                          const int id0, bool sameOrientation,
                                          const size_t pOuter)
{
//	in the inner, the number of dofs is as if it would be an element of order p-3.
	const size_t p = pOuter-3;

//	resize array
	vOrientOffset.clear();

//	loop mapped indices as required by rotated face
	for(size_t mj = 0; mj <= p; ++mj){
		for(size_t mi = 0; mi <= p-mj; ++mi){

		//	get corresponding multiindex in "natural" numbering
			size_t i,j;
			switch(id0)
			{
				case 0: i = mi;      j = mj;      break;
				case 1: i = mj;      j = p-mi-mj; break;
				case 2: i = p-mi-mj; j = mi;      break;
				default: UG_THROW("Orientation Triangle: Corner "<<id0<<" invalid.");
			}
			if(!sameOrientation) std::swap(i, j);

		//	linearize index and mapped index
			size_t naturalIndex = i;
			for(size_t c = 0; c < j; ++c)
				naturalIndex += (p+1-c);

		//	set mapping
			vOrientOffset.push_back(naturalIndex);
		}
	}
};

void ComputeOrientationOffsetLagrange(std::vector<size_t>& vOrientOffset,
                                      Volume* volume, Face* face, size_t nrFace,
                                      const size_t p)
{
//	get face descriptor
	FaceDescriptor fd;
	volume->face_desc(nrFace, fd);

//	find id0 and orientation
	const int numCo = face->num_vertices();
	const int id0 = GetVertexIndex(&fd, face->vertex(0));
	const int id1 = GetVertexIndex(&fd, face->vertex(1));
	const bool sameOrientation = (id1 == (id0+1)%numCo);

	switch(numCo){
		case 3:
			MapLagrangeMultiIndexTriangle(vOrientOffset, id0, sameOrientation, p);
			break;
		case 4:
			MapLagrangeMultiIndexQuad(vOrientOffset, id0, sameOrientation, p);
			break;
		default: UG_THROW("No corner number "<<numCo<<" implemented.");
	}
};

void ComputeOrientationOffsetLagrange(std::vector<size_t>& vOrientOffset,
                                      GridObject* volume, GridObject* face, size_t nrFace,
                                      const size_t p)
{
	UG_THROW("Should never be called.")
}


////////////////////////////////////////////////////////////////////////////////
// General Implementation
////////////////////////////////////////////////////////////////////////////////

template <typename TBaseElem, typename TSubBaseElem>
void ComputeOrientationOffsetGeneric(std::vector<size_t>& vOrientOffset,
                                     TBaseElem* elem, TSubBaseElem* sub, size_t nrSub,
                                     const LFEID& lfeid)
{
	vOrientOffset.clear();

	// if subelem higher dim than elem, no orientation
	if(TSubBaseElem::dim >= TBaseElem::dim) return;

	switch(lfeid.type()){
		// lagrange: orientate
		case LFEID::LAGRANGE:
			// only orientate if sub-dim < fct-space dim
			if(!(TSubBaseElem::dim < lfeid.dim())) return;

			// only orientate for p > 2, (else only max 1 DoF on sub)
			if(lfeid.order() <= 2) return;

			// orientate
			ComputeOrientationOffsetLagrange(vOrientOffset, elem, sub, nrSub, lfeid.order());
			break;

		// other cases: no orientation
		default: return;
	}
}


void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Volume* vol, Face* face, size_t nrFace,
                              const LFEID& lfeid)
{
	ComputeOrientationOffsetGeneric(vOrientOffset, vol, face, nrFace, lfeid);
}

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Volume* vol, Edge* edge, size_t nrEdge,
                              const LFEID& lfeid)
{
	ComputeOrientationOffsetGeneric(vOrientOffset, vol, edge, nrEdge, lfeid);
}

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              Face* face, Edge* edge, size_t nrEdge,
                              const LFEID& lfeid)
{
	ComputeOrientationOffsetGeneric(vOrientOffset, face, edge, nrEdge, lfeid);
}

void ComputeOrientationOffset(std::vector<size_t>& vOrientOffset,
                              GridObject* Elem, GridObject* SubElem, size_t nrSub,
                              const LFEID& lfeid)
{
//	general case: no offset needed
	vOrientOffset.clear();
}

} // end namespace ug
