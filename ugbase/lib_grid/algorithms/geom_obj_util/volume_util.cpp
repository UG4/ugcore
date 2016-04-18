/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Martin Stepniewski
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

#include "volume_util.h"
#include "lib_grid/lg_base.h"
#include "edge_util.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	GetNeighbours - sreiter
void GetNeighbours(std::vector<Volume*>& vVolsOut, Grid& grid, Volume* v,
					int side, bool clearContainer)
{
	if(clearContainer)
		vVolsOut.clear();

//	if VOLOPT_AUTOGENERATE_FACES and FACEOPT_STORE_ASSOCIATED_VOLUMES are
//	activated, we may use them to find the connected volume quite fast.
	if(grid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES
							| FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		Face* f = grid.get_face(v, side);
		Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(f);
		for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(f);
			iter != iterEnd; ++iter)
		{
			if(*iter != v)
				vVolsOut.push_back(*iter);
		}

		return;
	}

//	we can't assume that associated faces exist.
//	we have to find the neighbour by hand.
//	mark all vertices of the side
	grid.begin_marking();

	FaceDescriptor fd;
	v->face_desc(side, fd);
	uint numFaceVrts = fd.num_vertices();
	for(uint i = 0; i < numFaceVrts; ++ i)
		grid.mark(fd.vertex(i));

//	iterate over associated volumes of the first vertex and count
//	the number of marked vertices it contains.
	Vertex* vrt = fd.vertex(0);
	Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(vrt);
	for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(vrt);
		iter != iterEnd; ++iter)
	{
		Volume* vol = *iter;
		if(vol != v){
			size_t count = 0;
			uint numVrts = vol->num_vertices();
			for(uint i = 0; i < numVrts; ++i){
				if(grid.is_marked(vol->vertex(i)))
					++count;
			}

		//	if the number of marked vertices in vol matches the
		//	number of vertices of the specified side, we consider
		//	the volume to be a neighbout of that side.
			if(count == numFaceVrts)
				vVolsOut.push_back(vol);
		}
	}

	grid.end_marking();
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateVolumeMinHeight - mstepnie
number CalculateMinVolumeHeight(Tetrahedron* tet,
								Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	vector3& a = aaPos[tet->vertex(0)];
	vector3& b = aaPos[tet->vertex(1)];
	vector3& c = aaPos[tet->vertex(2)];
	vector3& d = aaPos[tet->vertex(3)];

	return CalculateMinTetrahedronHeight(a, b, c, d);
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateMinTetrahedronHeight - mstepnie
number CalculateMinTetrahedronHeight(const vector3& a, const vector3& b, 
									 const vector3& c, const vector3& d)
{
	number minHeight, tmpMinHeight;

//	Assume a tetrahedron with vertices a, b, c, d. Calculate its direction vectors
	vector3 ab;
	vector3 ac;
	vector3 ad;
	vector3 bd;
	vector3 bc;

	VecSubtract(ab, b, a);
	VecSubtract(ac, c, a);
	VecSubtract(ad, d, a);
	VecSubtract(bd, d, b);
	VecSubtract(bc, c, b);

//	calculate the 4 face normals
	vector3 nabc;
	vector3 nabd;
	vector3 nacd;
	vector3 nbcd;

	VecCross(nabc, ab, ac);
	VecCross(nabd, ad, ab);
	VecCross(nacd, ac, ad);
	VecCross(nbcd, bd, bc);

//	LOTFUSSVERFAHREN
	vector3 CutComb;

	///////////
	// FACE ABC
	//	set up matrix for calculating the orthogonal projection of vertex d on face abc
	matrix33 A;
	for(uint i = 0; i<3; ++i)
	{
		A[i][0] = ab[i];
	}
	for(uint i = 0; i<3; ++i)
	{
		A[i][1] = ac[i];
	}
	for(uint i = 0; i<3; ++i)
	{
		A[i][2] = -nabc[i];
	}

	// calculate the height of d respecting abc
	matrix33 A_inv;
	Inverse(A_inv, A);

	vector3 rhs;
	VecSubtract(rhs, d, a);

	MatVecMult(CutComb, A_inv, rhs);
	VecScale(nabc, nabc, CutComb[2]);
	minHeight = VecLength(nabc);


	///////////
	// FACE abd
	//	set up matrix for calculating the orthogonal projection of vertex c on face abd
	for(uint i = 0; i<3; ++i)
	{
		A[i][0] = ad[i];
	}

	for(uint i = 0; i<3; ++i)
	{
		A[i][1] = ab[i];
	}

	for(uint i = 0; i<3; ++i)
	{
		A[i][2] = -nabd[i];
	}

	// calculate the height of d respecting abc
	Inverse(A_inv, A);

	VecSubtract(rhs, c, a);

	MatVecMult(CutComb, A_inv, rhs);
	VecScale(nabd, nabd, CutComb[2]);
	tmpMinHeight = VecLength(nabd);

	if(tmpMinHeight > minHeight)
	{
		minHeight = tmpMinHeight;
	}


	///////////
	// FACE acd
	//	set up matrix for calculating the orthogonal projection of vertex b on face acd
	for(uint i = 0; i<3; ++i)
	{
		A[i][0] = ac[i];
	}
	for(uint i = 0; i<3; ++i)
	{
		A[i][1] = ad[i];
	}
	for(uint i = 0; i<3; ++i)
	{
		A[i][2] = -nacd[i];
	}

	// calculate the height of b respecting acd
	Inverse(A_inv, A);

	VecSubtract(rhs, b, a);

	MatVecMult(CutComb, A_inv, rhs);
	VecScale(nacd, nacd, CutComb[2]);
	tmpMinHeight = VecLength(nacd);

	if(tmpMinHeight < minHeight)
	{
		minHeight = tmpMinHeight;
	}


	///////////
	// FACE bcd
	//	set up matrix for calculating the orthogonal projection of vertex c on face bcd
	for(uint i = 0; i<3; ++i)
	{
		A[i][0] = bd[i];
	}
	for(uint i = 0; i<3; ++i)
	{
		A[i][1] = bc[i];
	}
	for(uint i = 0; i<3; ++i)
	{
		A[i][2] = -nbcd[i];
	}

	// calculate the height of a respecting bcd
	Inverse(A_inv, A);

	VecSubtract(rhs, a, b);

	MatVecMult(CutComb, A_inv, rhs);
	VecScale(nbcd, nbcd, CutComb[2]);
	tmpMinHeight = VecLength(nbcd);

	if(tmpMinHeight < minHeight)
	{
		minHeight = tmpMinHeight;
	}


	return minHeight;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronAspectRatio - mstepnie
number CalculateTetrahedronAspectRatio(Grid& grid, Tetrahedron* tet,
							Grid::VertexAttachmentAccessor<AVector3>& aaPos)
{
	/*
	 * optimal Aspect Ratio of a regular tetrahedron
	 * Q = sqrt(2/3) * a / a = 0.81...
	 */

	number AspectRatio;
	number maxEdgelength;
	number minTetrahedronHeight;

//	Collect tetrahedron edges, find longest edge and calculate its length
	vector<Edge*> edges;
	CollectAssociated(edges, grid, tet);
	Edge* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
	maxEdgelength = EdgeLength(longestEdge, aaPos);

//	Calculate the minimal tetrahedron height
	minTetrahedronHeight = CalculateMinVolumeHeight(tet, aaPos);

//	Calculate the aspect ratio
	AspectRatio =  maxEdgelength / minTetrahedronHeight;

	return AspectRatio;
}

////////////////////////////////////////////////////////////////////////
// IntersectPlaneWithTetrahedron
size_t IntersectPlaneWithTetrahedron
(
	vector3 intsOut[4],
	const vector3& planePoint,
	const vector3& planeNormal,
	const vector3 t[4]
)
{
    using namespace tet_rules;
    size_t numIntersections = 0;
    int intersectingEdges[4];
    for(size_t ie = 0; ie < (size_t) NUM_EDGES; ++ie)
    {
        vector3 v0 = t[EDGE_VRT_INDS[ie][0]];
        vector3 v1 = t[EDGE_VRT_INDS[ie][1]];
        vector3 dir;
        VecSubtract(dir, v1, v0);

        vector3 p;
        number s;
        if(RayPlaneIntersection(p, s, v0, dir, planePoint, planeNormal))
        {
            if(s > 0 && s < 1)
            {
            	if (numIntersections >= 4) // to avoid numerical artefacts
            		UG_THROW ("IntersectPlaneWithTetrahedron:"
            			" Illegal number of intersections of a plane with a tetrahedron.");
                intsOut[numIntersections] = p;
                intersectingEdges[numIntersections] = ie;
                ++numIntersections;
            }
        }
    }
    
    if(numIntersections == 4)
    {
        if(intersectingEdges[1] == OPPOSED_EDGE[intersectingEdges[0]])
            swap(intsOut[1], intsOut[2]);
        else if(intersectingEdges[2] == OPPOSED_EDGE[intersectingEdges[1]])
            swap(intsOut[0], intsOut[1]);
    }

    return numIntersections;
}

void InsertCenterVertex(Grid& g, Volume* vol, Vertex* vrt, bool eraseOldVol)
{
//	get the sides of the volume and create new elements
	FaceDescriptor fd;
	for(size_t i = 0; i < vol->num_faces(); ++i){
		vol->face_desc(i, fd);
		if(fd.num_vertices() == 3){
		//	create a tetrahedron
			g.create<Tetrahedron>(TetrahedronDescriptor(fd.vertex(2), fd.vertex(1),
														fd.vertex(0), vrt), vol);
		}
		else if(fd.num_vertices() == 4){
		//	create a pyramid
			g.create<Pyramid>(PyramidDescriptor(fd.vertex(3), fd.vertex(2),
												fd.vertex(1), fd.vertex(0), vrt), vol);
		}
		else{
			UG_THROW("Unsupported face type in InsertCenterVertex (#Corners "
					<< fd.num_vertices() << ")");
		}
	}

	if(eraseOldVol)
		g.erase(vol);
}

/*
double CMesh::calculate_volume_gauss() {
	
	double volume = 0.0;
	
	for(uint i = 0; i < number_of_triangles; i++) {
		
	 �uint a = triangles[i].a;
		uint b = triangles[i].b;
		uint c = triangles[i].c;
		
		//printf("%f %f %f\n", vertices[b].x(), vertices[b].y(), vertices[b].z());
		
		double x = (vertices[b].y() - vertices[a].y()) * (vertices[c].z() - vertices[a].z()) - (vertices[b].z() - vertices[a].z()) * (vertices[c].y() - vertices[a].y());
		double y = (vertices[b].z() - vertices[a].z()) * (vertices[c].x() - vertices[a].x()) - (vertices[b].x() - vertices[a].x()) * (vertices[c].z() - vertices[a].z());
		double z = (vertices[b].x() - vertices[a].x()) * (vertices[c].y() - vertices[a].y()) - (vertices[b].y() - vertices[a].y()) * (vertices[c].x() - vertices[a].x());
		
		double length = sqrt(x * x + y * y + z * z);
		
		double surface = 0.5 * length;
		
		if(length > 0.0) {
		
		 �x /= length;
	 �	y /= length; 
 �		z /= length;
			
			double sx = (vertices[a].x() + vertices[b].x() + vertices[c].x()) / 3.0;
			double sy = (vertices[a].y() + vertices[b].y() + vertices[c].y()) / 3.0;
			double sz = (vertices[a].z() + vertices[b].z() + vertices[c].z()) / 3.0;
			
			volume += 1.0 / 3.0 * surface * (sx * x + sy * y + sz * z);
		
		}
		
		
	}
	
	return volume;
	
	
	
}
*/
}//	end of namespace


