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

#include <vector>
#include <algorithm>
#include "grid_objects_2d.h"
//#include "../algorithms/geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	TOOLS
///	helpful if a local vertex-order is required
/**
 * cornersOut and cornersIn both have to be of size numCorners.
 * After termination cornersOut will contain the vertices of
 * cornersIn, starting from firstCorner, taking vertices modulo numCorners.
 * If cornersOut == cornersIn, the method will fail! This is ok since
 * the method is used locally and has been created for a special case.
 */
static inline
bool ReorderCornersCCW(Vertex** cornersOut, Vertex** const cornersIn,
					   int numCorners, int firstCorner)
{
	cornersOut[0] = cornersIn[firstCorner];
	for(int i = 1; i < numCorners; ++i)
		cornersOut[i] = cornersIn[(firstCorner + i) % numCorners];
	return true;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	FACES

////////////////////////////////////////////////////////////////////////
//	TriangleDescriptor
TriangleDescriptor::TriangleDescriptor(const TriangleDescriptor& td)
{
	m_vertex[0] = td.vertex(0);
	m_vertex[1] = td.vertex(1);
	m_vertex[2] = td.vertex(2);
}

TriangleDescriptor::TriangleDescriptor(Vertex* v1, Vertex* v2, Vertex* v3)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
}

////////////////////////////////////////////////////////////////////////
//	CustomTriangle
template <class ConcreteTriangleType, class BaseClass, class RefTriType, class RefQuadType>
CustomTriangle<ConcreteTriangleType, BaseClass, RefTriType, RefQuadType>::
CustomTriangle(const TriangleDescriptor& td)
{
	m_vertices[0] = td.vertex(0);
	m_vertices[1] = td.vertex(1);
	m_vertices[2] = td.vertex(2);
}

template <class ConcreteTriangleType, class BaseClass, class RefTriType, class RefQuadType>
CustomTriangle<ConcreteTriangleType, BaseClass, RefTriType, RefQuadType>::
CustomTriangle(Vertex* v1, Vertex* v2, Vertex* v3)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
}

template <class ConcreteTriangleType, class BaseClass, class RefTriType, class RefQuadType>
std::pair<GridBaseObjectId, int>
CustomTriangle<ConcreteTriangleType, BaseClass, RefTriType, RefQuadType>::
get_opposing_object(Vertex* vrt) const
{
	for(int i = 0; i < 3; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(EDGE, (i + 1) % 3);
		}
	}

	UG_THROW("The given vertex is not contained in the given face.");
}

template <class ConcreteTriangleType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomTriangle<ConcreteTriangleType, BaseClass, RefTriType, RefQuadType>::
refine(std::vector<Face*>& vNewFacesOut,
		Vertex** newFaceVertexOut,
		Vertex** newEdgeVertices,
		Vertex* newFaceVertex,
		Vertex** pSubstituteVertices,
		int snapPointIndex)
{
//TODO: complete triangle refine

	*newFaceVertexOut = newFaceVertex;
	vNewFacesOut.clear();

//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	get the number of new vertices.
	uint numNewVrts = 0;
	for(uint i = 0; i < 3; ++i)
	{
		if(newEdgeVertices[i] != NULL)
			++numNewVrts;
	}

//	if newFaceVertex is specified, then create three sub-triangles and
//	refine each. If not then refine the triangle depending
	if(newFaceVertex)
	{
		if(numNewVrts > 0){
			assert(!"Problem in CustomTriangle::refine: refine with newFaceVertex and newEdgeVertices is not yet supported.");
			return false;
		}

	//	create three new triangles
		vNewFacesOut.push_back(new RefTriType(vrts[0], vrts[1], newFaceVertex));
		vNewFacesOut.push_back(new RefTriType(vrts[1], vrts[2], newFaceVertex));
		vNewFacesOut.push_back(new RefTriType(vrts[2], vrts[0], newFaceVertex));
		return true;
	}
	else
	{
		switch(numNewVrts)
		{
			case 0: // this may happen when the triangle belongs to a prism being anisotropically refined
					// and the volume on the other side is not being refined
				vNewFacesOut.push_back(new RefTriType(vrts[0], vrts[1], vrts[2]));
				return true;

			case 1:
			{
			//	get the index of the edge that will be refined
				int iNew = -1;;
				for(int i = 0; i < 3; ++i){
					if(newEdgeVertices[i]){
						iNew = i;
						break;
					}
				}
				
			//	the corners. The first corner is the corner on the opposite side of iNew.
			//	Other follow in ccw order
				int iCorner[3];
				iCorner[0] = (iNew + 2) % 3;
				iCorner[1] = (iCorner[0] + 1) % 3;
				iCorner[2] = (iCorner[1] + 1) % 3;
					
			//	create the new triangles.
				vNewFacesOut.push_back(new RefTriType(vrts[iCorner[0]], vrts[iCorner[1]],
																newEdgeVertices[iNew]));
				vNewFacesOut.push_back(new RefTriType(vrts[iCorner[0]], newEdgeVertices[iNew],
																vrts[iCorner[2]]));
																
				return true;
			}

			case 2:
			{
			//	get the index of the edge that won't be refined
				int iFree = -1;
				for(int i = 0; i < 3; ++i){
					if(!newEdgeVertices[i]){
						iFree = i;
						break;
					}
				}
				
			//	the refined edges
				int iNew[2];
				iNew[0] = (iFree + 1) % 3;
				iNew[1] = (iFree + 2) % 3;
				
			//	the corners
				int iCorner[3];
				iCorner[0] = iFree;
				iCorner[1] = (iFree + 1) % 3;
				iCorner[2] = (iFree + 2) % 3;
				
			//	create the faces
				vNewFacesOut.push_back(new RefTriType(newEdgeVertices[iNew[0]],
																vrts[iCorner[2]],
																newEdgeVertices[iNew[1]]));
				vNewFacesOut.push_back(new RefQuadType(vrts[iCorner[0]], vrts[iCorner[1]],
														newEdgeVertices[iNew[0]], newEdgeVertices[iNew[1]]));
				return true;
			}

			case 3:
			{
			//	perform regular refine.
				vNewFacesOut.push_back(new RefTriType(vrts[0], newEdgeVertices[0], newEdgeVertices[2]));
				vNewFacesOut.push_back(new RefTriType(vrts[1], newEdgeVertices[1], newEdgeVertices[0]));
				vNewFacesOut.push_back(new RefTriType(vrts[2], newEdgeVertices[2], newEdgeVertices[1]));
				vNewFacesOut.push_back(new RefTriType(newEdgeVertices[0], newEdgeVertices[1], newEdgeVertices[2]));
				return true;
			}

		}
	}

	return false;
}

template <class ConcreteTriangleType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomTriangle<ConcreteTriangleType, BaseClass, RefTriType, RefQuadType>::
is_regular_ref_rule(int edgeMarks) const
{
	return edgeMarks == 7;
}

template <class ConcreteTriangleType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomTriangle<ConcreteTriangleType, BaseClass, RefTriType, RefQuadType>::
collapse_edge(std::vector<Face*>& vNewFacesOut,
				int edgeIndex, Vertex* newVertex,
				Vertex** pSubstituteVertices)
{
//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();
	return true;
}

template <class ConcreteTriangleType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomTriangle<ConcreteTriangleType, BaseClass, RefTriType, RefQuadType>::
collapse_edges(std::vector<Face*>& vNewFacesOut,
				std::vector<Vertex*>& vNewEdgeVertices,
				Vertex** pSubstituteVertices)
{
	if(vNewEdgeVertices.size() > BaseClass::num_edges())
	{
		assert(!"WARNING in Triangle::collapse_edges(...): bad number of newEdgeVertices.");
		LOG("WARNING in Triangle::collapse_edges(...): bad number of newEdgeVertices.");
		return false;
	}

//	check if there is a valid entry in vNewEdgeVertices
	bool bGotOne = false;
	for(uint i = 0; i < vNewEdgeVertices.size(); ++i)
	{
		if(vNewEdgeVertices[i] != NULL)
		{
			bGotOne = true;
			break;
		}
	}

	if(!bGotOne)
	{
		assert(!"WARNING in Triangle::collapse:edges(...): no new vertex was specified.");
		LOG("WARNING in Triangle::collapse:edges(...): no new vertex was specified.");
		return false;
	}

//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();
	return true;
}

//	BEGIN Depreciated
template <class ConcreteTriangleType, class BaseClass, class RefTriType, class RefQuadType>
void
CustomTriangle<ConcreteTriangleType, BaseClass, RefTriType, RefQuadType>::
create_faces_by_edge_split(int splitEdgeIndex,
								Vertex* newVertex,
								std::vector<Face*>& vNewFacesOut,
								Vertex** pSubstituteVertices)
{
	assert(((splitEdgeIndex >= 0) && (splitEdgeIndex < 3)) && "ERROR in Triangle::create_faces_by_edge_split(...): bad edge index!");

//	if pvSubstituteVertices is supplied, we will use the vertices in
//	pvSubstituteVertices instead the ones of 'this'.
//	If not, we will redirect the pointer to a local local vector,
//	that holds the vertices of 'this'
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	we have to find the indices ind0, ind1, ind2, where
//	ind0 is the index of the vertex on e before newVertex,
//	ind1 is the index of the vertex on e after newVertex
//	and ind2 is the index of the vertex not located on e.

	int ind0 = splitEdgeIndex;
	int ind1 = (ind0 + 1) % 3;
	int ind2 = (ind1 + 1) % 3;

	vNewFacesOut.push_back(new Triangle(vrts[ind0], newVertex, vrts[ind2]));
	vNewFacesOut.push_back(new Triangle(newVertex, vrts[ind1], vrts[ind2]));
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	QUADRILATERALS

////////////////////////////////////////////////////////////////////////
//	QuadrilateralDescriptor
QuadrilateralDescriptor::QuadrilateralDescriptor(const QuadrilateralDescriptor& qd)
{
	m_vertex[0] = qd.vertex(0);
	m_vertex[1] = qd.vertex(1);
	m_vertex[2] = qd.vertex(2);
	m_vertex[3] = qd.vertex(3);
}

QuadrilateralDescriptor::QuadrilateralDescriptor(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
{
	m_vertex[0] = v1;
	m_vertex[1] = v2;
	m_vertex[2] = v3;
	m_vertex[3] = v4;
}

////////////////////////////////////////////////////////////////////////
//	Quad

template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
CustomQuadrilateral(const QuadrilateralDescriptor& qd)
{
	m_vertices[0] = qd.vertex(0);
	m_vertices[1] = qd.vertex(1);
	m_vertices[2] = qd.vertex(2);
	m_vertices[3] = qd.vertex(3);
}

template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
CustomQuadrilateral(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4)
{
	m_vertices[0] = v1;
	m_vertices[1] = v2;
	m_vertices[2] = v3;
	m_vertices[3] = v4;
}

template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
get_opposing_side(EdgeVertices* e, EdgeDescriptor& edOut) const
{
	int localInd = Face::get_local_side_index(e);
	if(localInd == -1){
		return false;
	}

	edge_desc((localInd + 2) % 4, edOut);
	return true;
}

template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
std::pair<GridBaseObjectId, int>
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
get_opposing_object(Vertex* vrt) const
{
	for(int i = 0; i < 4; ++i){
		if(vrt == m_vertices[i]){
			return make_pair(VERTEX, (i + 2) % 4);
		}
	}

	UG_THROW("The given vertex is not contained in the given face.");
}

template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
void
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
create_faces_by_edge_split(int splitEdgeIndex,
							Vertex* newVertex,
							std::vector<Face*>& vNewFacesOut,
							Vertex** pSubstituteVertices)
{
	assert(((splitEdgeIndex >= 0) && (splitEdgeIndex < 4)) && "ERROR in Quadrilateral::create_faces_by_edge_split(...): bad edge index!");

//	if pvSubstituteVertices is supplied, we will use the vertices in
//	pvSubstituteVertices instead the ones of 'this'.
//	If not, we will redirect the pointer to a local local vector,
//	that holds the vertices of 'this'
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	we have to find the indices ind0, ind1, ind2, where
//	ind0 is the index of the vertex on e before newVertex,
//	ind1 is the index of the vertex on e after newVertex
//	and ind2 and ind3 are the indices of the vertices not located on e.
	int ind0 = splitEdgeIndex; //edge-index equals the index of its first vertex.
	int ind1 = (ind0 + 1) % 4;
	int ind2 = (ind0 + 2) % 4;
	int ind3 = (ind0 + 3) % 4;

	TriangleDescriptor td;

//	edge-split generates 3 triangles
	vNewFacesOut.push_back(new RefTriType(vrts[ind0], newVertex, vrts[ind3]));
	vNewFacesOut.push_back(new RefTriType(newVertex, vrts[ind1], vrts[ind2]));
	vNewFacesOut.push_back(new RefTriType(vrts[ind3], newVertex, vrts[ind2]));

//	we're done.
}

////////////////////////////////////////////////////////////////////////
//	Quadrilateral::refine
template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
refine(std::vector<Face*>& vNewFacesOut,
		Vertex** newFaceVertexOut,
		Vertex** edgeVrts,
		Vertex* newFaceVertex,
		Vertex** pSubstituteVertices,
		int snapPointIndex)
{
//TODO: complete quad refine
	*newFaceVertexOut = newFaceVertex;
	vNewFacesOut.clear();
	
//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	check which edges have to be refined and perform the required operations.
//	get the number of new vertices.
	uint numNewVrts = 0;
	for(uint i = 0; i < 4; ++i)
	{
		if(edgeVrts[i] != NULL)
			++numNewVrts;
	}

	switch(numNewVrts)
	{
		case 0:
			// refine with mid point if it exists
			if (newFaceVertex)
			{
			//	create four new triangles
				vNewFacesOut.push_back(new RefTriType(vrts[0], vrts[1], newFaceVertex));
				vNewFacesOut.push_back(new RefTriType(vrts[1], vrts[2], newFaceVertex));
				vNewFacesOut.push_back(new RefTriType(vrts[2], vrts[3], newFaceVertex));
				vNewFacesOut.push_back(new RefTriType(vrts[3], vrts[0], newFaceVertex));
				return true;
			}

			// in case the mid point does not exists, we need a simple copy
			// This may happen when the quad belongs to a hexahedron being anisotropically refined
			// and the volume on the other side is not being refined.
			vNewFacesOut.push_back(new RefQuadType(vrts[0], vrts[1], vrts[2], vrts[3]));
			return true;
			
		case 1:
		{
			if(newFaceVertex){
				assert(!"Problem in CustomQuadrilateral::refine: refine with newFaceVertex and one newEdgeVertex is not yet supported.");
				return false;
			}
			
			int iNew = -1;
			for(int i = 0; i < 4; ++i){
				if(edgeVrts[i]){
					iNew = i;
					break;
				}
			}

			if(snapPointIndex != -1 && (snapPointIndex == iNew || snapPointIndex == (iNew + 1) % 4)){
				UG_LOG("WARNING: Invalid snap-point distribution detected. Ignoring snap-points for this element.\n");
				snapPointIndex = -1;
			}
			
		//	the corners in a local ordering relative to iNew. Keeping ccw order.
			Vertex* corner[4];
			const int rot = (iNew + 3) % 4;
			ReorderCornersCCW(corner, vrts, 4, rot);

		//	create the new elements
			if(snapPointIndex == -1){
				vNewFacesOut.push_back(new RefTriType(corner[0], corner[1], edgeVrts[iNew]));
				vNewFacesOut.push_back(new RefTriType(corner[0], edgeVrts[iNew], corner[3]));
				vNewFacesOut.push_back(new RefTriType(corner[3], edgeVrts[iNew], corner[2]));
			}
			else{
				snapPointIndex = (snapPointIndex + 4 - rot) % 4;
				if(snapPointIndex == 0){
					vNewFacesOut.push_back(new RefTriType(corner[0], corner[1], edgeVrts[iNew]));
					vNewFacesOut.push_back(new RefQuadType(corner[0], edgeVrts[iNew], corner[2], corner[3]));
				}
				else if(snapPointIndex == 3){
					vNewFacesOut.push_back(new RefQuadType(corner[0], corner[1], edgeVrts[iNew], corner[3]));
					vNewFacesOut.push_back(new RefTriType(corner[3], edgeVrts[iNew], corner[2]));
				}
				else{
					UG_THROW("Unexpected snap-point index: " << snapPointIndex << ". This is an implementation error!");
				}
			}

			return true;
		}

		case 2:
		{
			if(newFaceVertex){
				assert(!"Problem in CustomQuadrilateral::refine: refine with newFaceVertex and two newEdgeVertices is not yet supported.");
				return false;
			}

		//	there are two cases (refined edges are adjacent or opposite to each other)
			int iNew[2];
			int counter = 0;
			for(int i = 0; i < 4; ++i){
				if(edgeVrts[i]){
					iNew[counter] = i;
					++counter;
				}
			}

		//	corners will be filled later on
			Vertex* corner[4];

		//	check which case applies
			if(iNew[1] - iNew[0] == 2){
			//	edges are opposite to each other
			//	the corners in a local ordering relative to iNew[0]. Keeping ccw order.
				ReorderCornersCCW(corner, vrts, 4, (iNew[0] + 3) % 4);
				
			//	create new faces
				vNewFacesOut.push_back(new RefQuadType(corner[0], corner[1],
													edgeVrts[iNew[0]], edgeVrts[iNew[1]]));

				vNewFacesOut.push_back(new RefQuadType(edgeVrts[iNew[1]], edgeVrts[iNew[0]],
														corner[2], corner[3]));
			}
			else{
			//	edges are adjacent
			//	we want that (iNew[0] + 1) % 4 == iNew[1]
				if((iNew[0] + 1) % 4 != iNew[1])
					swap(iNew[0], iNew[1]);
			
			//	the corners in a local ordering relative to iNew[0]. Keeping ccw order.
				const int rot = (iNew[0] + 3) % 4;
				ReorderCornersCCW(corner, vrts, 4, rot);

				if(snapPointIndex != -1 && ((snapPointIndex + 4 - rot) % 4) != 0){
					snapPointIndex = -1;
					UG_LOG("WARNING: Invalid snap-point distribution detected. Ignoring snap-points for this element.\n");
				}


			//	create new faces
				if(snapPointIndex == -1){
					vNewFacesOut.push_back(new RefTriType(corner[0], corner[1], edgeVrts[iNew[0]]));
					vNewFacesOut.push_back(new RefTriType(edgeVrts[iNew[0]], corner[2], edgeVrts[iNew[1]]));
					vNewFacesOut.push_back(new RefTriType(corner[3], corner[0], edgeVrts[iNew[1]]));
					vNewFacesOut.push_back(new RefTriType(corner[0], edgeVrts[iNew[0]], edgeVrts[iNew[1]]));
				}
				else{
					vNewFacesOut.push_back(new RefTriType(corner[0], corner[1], edgeVrts[iNew[0]]));
					vNewFacesOut.push_back(new RefTriType(corner[3], corner[0], edgeVrts[iNew[1]]));
					vNewFacesOut.push_back(new RefQuadType(corner[0], edgeVrts[iNew[0]], corner[2], edgeVrts[iNew[1]]));
				}
			}

			return true;
		}

		case 3:
		{
			if(newFaceVertex){
				assert(!"Problem in CustomQuadrilateral::refine: refine with newFaceVertex and three newEdgeVertices is not yet supported.");
				return false;
			}

			int iFree = -1;
			for(int i = 0; i < 4; ++i){
				if(!edgeVrts[i]){
					iFree = i;
					break;
				}
			}
			
		//	the vertices on the edges:
			Vertex* nvrts[3];
			nvrts[0] = edgeVrts[(iFree + 1) % 4];
			nvrts[1] = edgeVrts[(iFree + 2) % 4];
			nvrts[2] = edgeVrts[(iFree + 3) % 4];

		//	the corners in a local ordering relative to iNew. Keeping ccw order.
			Vertex* corner[4];
			ReorderCornersCCW(corner, vrts, 4, (iFree + 1) % 4);

		//	create the faces
			vNewFacesOut.push_back(new RefQuadType(corner[0], nvrts[0], nvrts[2], corner[3]));
			vNewFacesOut.push_back(new RefTriType(corner[1], nvrts[1], nvrts[0]));
			vNewFacesOut.push_back(new RefTriType(corner[2], nvrts[2], nvrts[1]));
			vNewFacesOut.push_back(new RefTriType(nvrts[0], nvrts[1], nvrts[2]));

			return true;
		}

		case 4:
		{
		//	we'll create 4 new quads. create a new center if required.
			if(!newFaceVertex)
				newFaceVertex = new RegularVertex;

			*newFaceVertexOut = newFaceVertex;
		
			vNewFacesOut.push_back(new RefQuadType(vrts[0], edgeVrts[0], newFaceVertex, edgeVrts[3]));
			vNewFacesOut.push_back(new RefQuadType(vrts[1], edgeVrts[1], newFaceVertex, edgeVrts[0]));
			vNewFacesOut.push_back(new RefQuadType(vrts[2], edgeVrts[2], newFaceVertex, edgeVrts[1]));
			vNewFacesOut.push_back(new RefQuadType(vrts[3], edgeVrts[3], newFaceVertex, edgeVrts[2]));
			return true;
		}
	}

	return false;
}

template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
is_regular_ref_rule(int edgeMarks) const
{
	static const int allEdges = 15;	// 1111
	static const int hEdges = 5;	// 0101
	static const int vEdges = 10;	// 1010

	return 		(edgeMarks == allEdges)
			||	(edgeMarks == hEdges)
			||	(edgeMarks == vEdges);
}


template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
collapse_edge(std::vector<Face*>& vNewFacesOut,
				int edgeIndex, Vertex* newVertex,
				Vertex** pSubstituteVertices)
{
//	if an edge of the triangle is collapsed, nothing remains
	vNewFacesOut.clear();

//	handle substitute vertices.
	Vertex** vrts;
	if(pSubstituteVertices)
		vrts = pSubstituteVertices;
	else
		vrts = m_vertices;

//	the collapsed edge connects vertices at ind0 and ind1.
	int ind0 = edgeIndex; //edge-index equals the index of its first vertex.
	int ind1 = (ind0 + 1) % 4;
	int ind2 = (ind1 + 1) % 4;
	int ind3 = (ind2 + 1) % 4;

//	ind0 and ind1 will be replaced by newVertex.
	vNewFacesOut.push_back(new RefTriType(newVertex, vrts[ind2], vrts[ind3]));
	return true;
}

template <class ConcreteQuadrilateralType, class BaseClass, class RefTriType, class RefQuadType>
bool
CustomQuadrilateral<ConcreteQuadrilateralType, BaseClass, RefTriType, RefQuadType>::
collapse_edges(std::vector<Face*>& vNewFacesOut,
				std::vector<Vertex*>& vNewEdgeVertices,
				Vertex** pSubstituteVertices)
{
	if(vNewEdgeVertices.size() > BaseClass::num_edges())
	{
		assert(!"WARNING in Quadrilateral::collapse_edges(...): bad number of newEdgeVertices.");
		LOG("WARNING in Quadrilateral::collapse_edges(...): bad number of newEdgeVertices.");
		return false;
	}

	vNewFacesOut.clear();

//	check if there is a valid entry in vNewEdgeVertices
	int collapseIndex = -1;
	int numCollapses = 0;
	for(uint i = 0; i < 4; ++i)
	{
		if(i < vNewEdgeVertices.size())
		{
			if(vNewEdgeVertices[i] != NULL)
			{
				++numCollapses;
				collapseIndex = i;
			}
		}
	}

//	if more than 1 edge is collapsed, nothing remains.
	if(numCollapses == 0)
	{
		assert(!"WARNING in Quadrilateral::collapse:edges(...): no new vertex was specified.");
		LOG("WARNING in Quadrilateral::collapse:edges(...): no new vertex was specified.");
		return false;
	}
	else if(numCollapses == 1)
	{
	//	call collapse_edge with the edges index.
		collapse_edge(vNewFacesOut, collapseIndex, vNewEdgeVertices[collapseIndex], pSubstituteVertices);
	}

	return true;
}

//	explicit instantiation
template class CustomTriangle<Triangle, Face, Triangle, Quadrilateral>;
template class CustomTriangle<ConstrainedTriangle, ConstrainedFace,
							  ConstrainedTriangle, ConstrainedQuadrilateral>;
template class CustomTriangle<ConstrainingTriangle, ConstrainingFace,
							  ConstrainingTriangle, ConstrainingQuadrilateral>;

template class CustomQuadrilateral<Quadrilateral, Face, Triangle, Quadrilateral>;
template class CustomQuadrilateral<ConstrainedQuadrilateral, ConstrainedFace,
								   ConstrainedTriangle, ConstrainedQuadrilateral>;
template class CustomQuadrilateral<ConstrainingQuadrilateral, ConstrainingFace,
								   ConstrainingTriangle, ConstrainingQuadrilateral>;


////////////////////////////////////////////////////////////////////////////////
//	CONSTRAINING FACE
template <> size_t
ConstrainingFace::
num_constrained<Vertex>() const
{
	return num_constrained_vertices();
}

template <> size_t
ConstrainingFace::
num_constrained<Edge>() const
{
	return num_constrained_edges();
}

template <> size_t
ConstrainingFace::
num_constrained<Face>() const
{
	return num_constrained_faces();
}


template <> Vertex*
ConstrainingFace::
constrained<Vertex>(size_t ind) const
{
	return constrained_vertex(ind);
}

template <> Edge*
ConstrainingFace::
constrained<Edge>(size_t ind) const
{
	return constrained_edge(ind);
}

template <> Face*
ConstrainingFace::
constrained<Face>(size_t ind) const
{
	return constrained_face(ind);
}

}//	end of namespace
