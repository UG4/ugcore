/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#include "face_util.h"

#include <vector>
//#include <stack>
#include <cassert>
#include "vertex_util.h"

//ø #include "../attachment_util.h"

using namespace std;

namespace ug {

////////////////////////////////////////////////////////////////////////
int GetFaceIndex(Volume* vol, Face* f)
{
	size_t numFaces = vol->num_faces();
	FaceDescriptor fd;
	for(uint i = 0; i < numFaces; ++i)
	{
		vol->face_desc(i, fd);
		if(CompareVertices(f, &fd))
			return static_cast<int>(i);
	}
	return -1;
}


////////////////////////////////////////////////////////////////////////
//	CalculateNormal
void CalculateNormal(vector3& vNormOut, const FaceVertices* face,
					Grid::AttachmentAccessor<Vertex, APosition>& aaPos)
{
	if(face->num_vertices() == 3)
	{
		CalculateTriangleNormal(vNormOut, aaPos[face->vertex(0)],
								aaPos[face->vertex(1)], aaPos[face->vertex(2)]);
		return;
	}
	else if(face->num_vertices() == 4)
	{
		vector3 n1, n2;
		CalculateTriangleNormalNoNormalize(n1, aaPos[face->vertex(0)],
								aaPos[face->vertex(1)], aaPos[face->vertex(2)]);
		CalculateTriangleNormalNoNormalize(n2, aaPos[face->vertex(2)],
								aaPos[face->vertex(3)], aaPos[face->vertex(0)]);
		VecAdd(vNormOut, n1, n2);
		VecNormalize(vNormOut, vNormOut);

		return;
	}
	else if(face->num_vertices() > 4)
	{
		CalculateTriangleNormal(vNormOut, aaPos[face->vertex(0)],
								aaPos[face->vertex(1)], aaPos[face->vertex(2)]);
		return;
	}

	vNormOut = vector3(0, 0, 0);
}

void CalculateNormalNoNormalize(vector3& vNormOut, FaceVertices* face,
								Grid::AttachmentAccessor<Vertex, APosition>& aaPos)
{
	if(face->num_vertices() == 3 || face->num_vertices() > 4)
	{
		CalculateTriangleNormalNoNormalize(vNormOut, aaPos[face->vertex(0)],
								aaPos[face->vertex(1)], aaPos[face->vertex(2)]);
		return;
	}
	else if(face->num_vertices() == 4)
	{
		vector3 n1, n2;
		CalculateTriangleNormalNoNormalize(n1, aaPos[face->vertex(0)],
								aaPos[face->vertex(1)], aaPos[face->vertex(2)]);
		CalculateTriangleNormalNoNormalize(n2, aaPos[face->vertex(2)],
								aaPos[face->vertex(3)], aaPos[face->vertex(0)]);
		VecAdd(vNormOut, n1, n2);
		VecScale(vNormOut, vNormOut, 0.5);

		return;
	}

	vNormOut = vector3(0, 0, 0);
}

////////////////////////////////////////////////////////////////////////
//	CalculateFaceNormals
void CalculateFaceNormals(Grid& grid, const FaceIterator& facesBegin,
						const FaceIterator& facesEnd,
						AVector3& aPos, AVector3& aNorm)
{
	if(!grid.has_vertex_attachment(aPos))
		return;

	if(!grid.has_face_attachment(aNorm))
		grid.attach_to_faces(aNorm);

	Grid::VertexAttachmentAccessor aaPos(grid, aPos);
	Grid::FaceAttachmentAccessor aaNorm(grid, aNorm);

//	iterate through the specified faces and calculate normals.
	for(FaceIterator iter = facesBegin; iter != facesEnd; ++iter)
		CalculateNormal(aaNorm[*iter], *iter, aaPos);
}

////////////////////////////////////////////////////////////////////////
int NumAssociatedVolumes(Grid& grid, Face* f)
{
//	check if FACEOPT_STORE_ASSOCIATED_VOLUMES is enabled.
//	if so, use it to count the number of adjacent volumes.
	int counter = 0;
	if(grid.option_is_enabled(FaceOptions::FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		for(auto iter = grid.associated_volumes_begin(f);
		    iter != grid.associated_volumes_end(f); ++iter)
		{
			++counter;
		}
	}
	else
	{
	//	iterate over all volumes which are connected to the first vertex
	//	and check if they contain the face...
		auto iterEnd = grid.associated_volumes_end(f->vertex(0));
		for(auto iter = grid.associated_volumes_begin(f->vertex(0));
		    iter != iterEnd; ++iter)
		{
			if(VolumeContains(*iter, f))
				++counter;
		}
	}

	return counter;
}

////////////////////////////////////////////////////////////////////////
bool IsVolumeBoundaryFace(Grid& grid, Face* f)
{
	return NumAssociatedVolumes(grid, f) == 1;
}

bool IsVolumeBoundaryFace(Grid& grid, ConstrainedFace* f)
{
	static vector<Volume*> vols;//avoid repeated reallocation
	CollectVolumes(vols, grid, f);
	if(vols.empty())
		return true;
	return false;
}

bool IsVolumeBoundaryFace(Grid& grid, ConstrainingFace* f)
{
	static vector<Volume*> vols;//avoid repeated reallocation
	CollectVolumes(vols, grid, f);
	if(vols.empty())
		return true;
	return false;
}


////////////////////////////////////////////////////////////////////////
//	FaceQuality
number FaceQuality(FaceVertices* f,
				Grid::VertexAttachmentAccessor<APosition> aaPos)
{
//	take dot-products at the corners.
	number quality = 1;
	uint numVrts = f->num_vertices();

	vector3 v1, v2;
	number len1, len2;

//	calculate the direction from first to last vertex.
	VecSubtract(v1, aaPos[f->vertex(numVrts-1)], aaPos[f->vertex(0)]);
	len1 = VecLength(v1);

	//VecNormalize(v1, v1);

	for(uint i = 0; i < numVrts; ++i)
	{
		VecSubtract(v2, aaPos[f->vertex((i+1)%numVrts)], aaPos[f->vertex(i)]);
		len2 = VecLength(v2);
		//VecNormalize(v2, v2);

		number nQual = 1. - fabs(VecDot(v1, v2) / (len1 * len2));
		if(nQual < quality)
			quality = nQual;
	//	v1 of the next iteration equals -v2 of this iteration
		VecScale(v1, v2, -1);
		len1 = len2;
	}

	if(numVrts == 3){
	//	since at least one angle is <= 60, we have to normalize the return value
		return quality * 2.;
	}
	return quality;
}

////////////////////////////////////////////////////////////////////////
//	TriangleQuality
number TriangleQuality(vector3& v1, vector3& v2, vector3& v3)
{
//	take dot-products at the corners.
	number quality = 1;

//	required for the iteration.
	vector3* pv[3];
	pv[0] = &v1; pv[1] = &v2; pv[2] = &v3;

	vector3 d1, d2;
//	calculate the direction from first to last vertex.
	VecSubtract(d1, v3, v1);
	VecNormalize(d1, d1);

	for(uint i = 0; i < 3; ++i)
	{
		VecSubtract(d2, *pv[i], *pv[(i+1)%3]);
		VecNormalize(d2, d2);
		number nQual = 1.f - fabs(VecDot(d1, d2));
		if(nQual < quality)
			quality = nQual;
	//	v1 of the next iteration equals -v2 of this iteration
		VecScale(d1, d2, -1);
	}

//	since at least one angle is <= 60, we have to normalize the return value
	return quality * 2.;
}

////////////////////////////////////////////////////////////////////////
//	Triangulate
void Triangulate(Grid& grid, Quadrilateral* q,
				Grid::VertexAttachmentAccessor<APosition>* paaPos)
{
//	there are two ways to triangulate a quadrilateral.
//	if quality is set to false, we simply choose the first way.
	if(!paaPos)
	{
		grid.create<Triangle>(TriangleDescriptor(q->vertex(0), q->vertex(1), q->vertex(2)), q);
		grid.create<Triangle>(TriangleDescriptor(q->vertex(2), q->vertex(3), q->vertex(0)), q);
	}
	else
	{
		Grid::VertexAttachmentAccessor<APosition>& aaPos = *paaPos;

		number q1 = 0;
		number q2 = 0;

	//	first check whether normals in the corners should affect the direction of the edge
		vector3 n[4];
		for(int i = 0; i < 4; ++i)
			CalculateVertexNormal(n[i], grid, q->vertex(i), aaPos);
		
		// vector3 dir[2];
		// VecSubtract(dir[0], aaPos[q->vertex(2)], aaPos[q->vertex(0)]);
		// VecNormalize(dir[0], dir[0]);
		// VecSubtract(dir[1], aaPos[q->vertex(3)], aaPos[q->vertex(1)]);
		// VecNormalize(dir[1], dir[1]);

		// number maxDot[2] = {0, 0};
		// for(int i = 0; i < 4; ++i){
		// 	for(int j = 0; j < 2; ++j){
		// 		maxDot[j] = std::max(maxDot[j], fabs(VecDot(dir[j], n[i])));
		// 	}
		// }

		// int imin = maxDot[0] < maxDot[1] ? 0 : 1;
		// int imax = maxDot[0] < maxDot[1] ? 1 : 0;

		// if(maxDot[imax] > 0 && maxDot[imax] - maxDot[imin] > SMALL){
		// 	q1 = 1. - maxDot[0];
		// 	q2 = 1. - maxDot[1];
		// }

		q1 = VecDot(n[0], n[2]);
		q2 = VecDot(n[1], n[3]);
		
		if(fabs(q1-q2) < SMALL)
		{
		//	if normals are equal we check which way produces the better triangles.
			q1 = std::min(
					TriangleQuality(aaPos[q->vertex(0)], aaPos[q->vertex(1)], aaPos[q->vertex(2)]),
					TriangleQuality(aaPos[q->vertex(2)], aaPos[q->vertex(3)], aaPos[q->vertex(0)]));

			q2 = std::min(
					TriangleQuality(aaPos[q->vertex(1)], aaPos[q->vertex(2)], aaPos[q->vertex(3)]),
					TriangleQuality(aaPos[q->vertex(3)], aaPos[q->vertex(0)], aaPos[q->vertex(1)]));
		}

		if(q1 >= q2)
		{
			grid.create<Triangle>(TriangleDescriptor(q->vertex(0), q->vertex(1), q->vertex(2)), q);
			grid.create<Triangle>(TriangleDescriptor(q->vertex(2), q->vertex(3), q->vertex(0)), q);
		}
		else
		{
			grid.create<Triangle>(TriangleDescriptor(q->vertex(1), q->vertex(2), q->vertex(3)), q);
			grid.create<Triangle>(TriangleDescriptor(q->vertex(3), q->vertex(0), q->vertex(1)), q);
		}
	}

	grid.erase(q);
}

////////////////////////////////////////////////////////////////////////
//	Triangulate
void Triangulate(Grid& grid,
				QuadrilateralIterator iterBegin,
				QuadrilateralIterator iterEnd,
				Grid::VertexAttachmentAccessor<APosition>* paaPos)
{
	while(iterBegin != iterEnd)
	{
		Quadrilateral* q = *iterBegin;
		++iterBegin;
		Triangulate(grid, q, paaPos);
	}
}

////////////////////////////////////////////////////////////////////////
//	GetNeighbours
void GetNeighbours(std::vector<Face*>& vFacesOut, Grid& grid, Face* f,
					int side, bool clearContainer)
{
//	in the current implementation this method requires, that all edges
//	are created for all faces.
//TODO: improve this!
	if(!grid.option_is_enabled(FaceOptions::FACEOPT_AUTOGENERATE_EDGES))
	{
		UG_LOG("WARNING: autoenabling FACEOPT_AUTOGENERATE_EDGES in GetNeighbours(Face).\n");
		grid.enable_options(FaceOptions::FACEOPT_AUTOGENERATE_EDGES);
	}
	
	if(clearContainer)
		vFacesOut.clear();
		
//	check if we can find an edge for the specified side
	Edge* e = grid.get_edge(f, side);
	assert(e && "edge not found though it should be there!");

//	get the connected faces
	vector<Face*> vFaces;
	CollectFaces(vFaces, grid, e);
	
//	push them to vFacesOut - except f
	for(size_t i = 0; i < vFaces.size(); ++i)
	{
		if(vFaces[i] != f)
			vFacesOut.push_back(vFaces[i]);
	}
}

void InsertCenterVertex(Grid& g, Face* f, Vertex* vrt, bool eraseOldFace)
{
//	get the sides of the face and create new elements
	EdgeDescriptor ed;
	for(size_t i = 0; i < f->num_edges(); ++i){
		f->edge_desc(i, ed);
		g.create<Triangle>(TriangleDescriptor(ed.vertex(0), ed.vertex(1), vrt), f);
	}

	if(eraseOldFace)
		g.erase(f);
}

}//	end of namespace
