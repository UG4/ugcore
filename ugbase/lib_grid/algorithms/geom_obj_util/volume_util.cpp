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
	vector3 vfaceNormal;
	number minHeight = std::numeric_limits<double>::max();
	number tmpMinHeight;


//	Iterate over all tetrahedron vertices and calculate height to opposing face
	for(size_t i = 0; i < 4; ++i)
	{
		Vertex* vrt = tet->vertex(i);
		FaceDescriptor fd = tet->face_desc(tet->get_opposing_object(vrt).second);

		CalculateNormal(vfaceNormal, &fd, aaPos);
		tmpMinHeight = DistancePointToPlane(aaPos[vrt], aaPos[fd.vertex(0)], vfaceNormal);

		if(tmpMinHeight < minHeight)
		{
			minHeight = tmpMinHeight;
		}
	}

	return minHeight;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronAspectRatio - mstepnie
number CalculateTetrahedronAspectRatio(Grid& grid, Tetrahedron* tet,
									   Grid::VertexAttachmentAccessor<AVector3>& aaPos)
{
	/*
	 * optimal Aspect Ratio of a regular tetrahedron with edge lengths a:
	 * Q = hmin/lmax = sqrt(2/3)*a / a = 0.8164...
	 *
	 * Info: return value is normalized by factor sqrt(3/2)
	 * (s. Shewchuk 2002)
	 */

	number aspectRatio;
	number maxEdgeLength;
	number minTetrahedronHeight;

//	Collect tetrahedron edges, find longest edge and calculate its length
	vector<Edge*> edges;
	CollectAssociated(edges, grid, tet);
	Edge* longestEdge = FindLongestEdge(edges.begin(), edges.end(), aaPos);
	maxEdgeLength = EdgeLength(longestEdge, aaPos);

//	Calculate the minimal tetrahedron height
	minTetrahedronHeight = CalculateMinVolumeHeight(tet, aaPos);

//	Calculate the aspect ratio
	aspectRatio =  std::sqrt(3/2.0) * minTetrahedronHeight / maxEdgeLength;

	return aspectRatio;
}

////////////////////////////////////////////////////////////////////////////////////////////
// CalculateHexahedronAspectRatio
////////////////////////////////////////////////////////////////////////////////////////////
number CalculateHexahedronAspectRatio(Grid& grid, Hexahedron* hex,
									   Grid::VertexAttachmentAccessor<AVector3>& aaPos)
{
	/*
	 * Another important indicator of the mesh quality is the aspect ratio.
	 * The aspect ratio is a measure of the stretching of a cell. It is computed
	 * as the ratio of the maximum value to the minimum value of any of the
	 * following distances: the normal distances between the cell centroid and
	 * face centroids (computed as a dot product of the distance vector and the
	 * face normal), and the distances between the cell centroid and nodes.
	 * For a unit cube (see Figure 5.23: Calculating the Aspect Ratio for a Unit
	 * Cube), the maximum distance is 0.866, and the minimum distance is 0.5,
	 * so the aspect ratio is 1.732. This type of definition can be applied on
	 * any type of mesh, including polyhedral.
	 *
	 * \see https://www.sharcnet.ca/Software/Ansys/17.0/en-us/help/flu_ug/flu_ug_mesh_quality.html
	 *
	 */
	vector<Face*> faces;
	CollectAssociated(faces, grid, hex);
	ug::vector3 center = CalculateCenter(hex, aaPos);

	/// Max distance to face center
	std::vector<number> As;
	for (size_t i = 0; i < faces.size(); i++)
	{
		vector3 dir;
		vector3 faceCenter = CalculateCenter(faces[i], aaPos);
		VecSubtract(dir, faceCenter, center);
		vector3 n;
		CalculateNormal(n, faces[i], aaPos);
		number dist = VecDot(dir, n);
		As.push_back(dist);
	}
	number max = *std::max_element(As.begin(), As.end());

	/// Min distance to node of hexaeder
	std::vector<number> Bs;
	for (size_t i = 0; i < hex->num_vertices(); i++)
	{
		number dist = VecDistance(aaPos[hex->vertex(i)], center);
		Bs.push_back(dist);
	}
	number min = *std::min_element(Bs.begin(), Bs.end());

	UG_COND_WARNING(std::abs(min) < SMALL, "Near 0-length minimum distance detected.");

	return max / min;
}

////////////////////////////////////////////////////////////////////////////////////////////
// CalculateHexahedronAspectRatio
// order of vertices should be the same as described in \sa PyramidDescriptor
// v1, v2, v3, v4: bottom-vertices in counterclockwise order (if viewed from the top).
// v5: top-vertex.
////////////////////////////////////////////////////////////////////////////////////////////
number CalculatePyramidAspectRatio
(
	Grid& grid,
	Pyramid* pyr,
	Grid::VertexAttachmentAccessor<AVector3>& aaPos
)
{
	// average edge length of base of pyramid
	number avg_edge_length = 0;
	for (size_t i = 0; i < 4; i++)
	{
		avg_edge_length += VecDistance(aaPos[pyr->vertex(i%4)], aaPos[pyr->vertex((i+1)%4)]);
	}
	avg_edge_length /= 4;

	// distance from base to top of pyramid
	const Face* const face = grid.get_element(FaceDescriptor(pyr->vertex(0), pyr->vertex(1), pyr->vertex(2), pyr->vertex(3)));
	vector3 normal;
	CalculateNormal(normal, face, aaPos);
	const number distance = VecDot(normal, aaPos[pyr->vertex(4)]);

	// AR of pyramid is ratio between distance and average edge length of base
	return distance / avg_edge_length;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronRootMeanSquareFaceArea - mstepnie
////////////////////////////////////////////////////////////////////////////////////////////
UG_API
number CalculateTetrahedronRootMeanSquareFaceArea(Grid& grid,
												  Tetrahedron* tet,
												  Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
//	Collect tetrahedron faces
	vector<Face*> faces;
	CollectAssociated(faces, grid, tet);

	number A;
	number A_rms = 0.0;

	for(size_t i = 0; i < faces.size(); ++i)
	{
		A = FaceArea(faces[i], aaPos);
		A_rms += A*A;
	}

	A_rms *= 0.25;
	A_rms = std::sqrt(A_rms);

	return A_rms;
}


////////////////////////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronVolToRMSFaceAreaRatio - mstepnie
////////////////////////////////////////////////////////////////////////////////////////////
UG_API
number CalculateTetrahedronVolToRMSFaceAreaRatio(Grid& grid,
										  	  	 Tetrahedron* tet,
												 Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	/*
	 * optimal volume to root-mean-square face area ratio of a
	 * regular tetrahedron with edge lengths a:
	 * Q = V/A_rms^(3/2)
	 *
	 * Info: return value is normalized by factor pow(3, 7/4.0) / 2.0 / sqrt(2);
	 * (s. Shewchuk 2002)
	 */

	number vol = CalculateTetrahedronVolume(aaPos[tet->vertex(0)],
											aaPos[tet->vertex(1)],
											aaPos[tet->vertex(2)],
											aaPos[tet->vertex(3)]);

	number A_rms = CalculateTetrahedronRootMeanSquareFaceArea(grid, tet, aaPos);
	number normalization = std::pow(3, 7/4.0) / 2.0 / std::sqrt(2);

	return normalization * vol / std::pow(A_rms, 3/2.0);
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
		
	 uint a = triangles[i].a;
		uint b = triangles[i].b;
		uint c = triangles[i].c;
		
		//printf("%f %f %f\n", vertices[b].x(), vertices[b].y(), vertices[b].z());
		
		double x = (vertices[b].y() - vertices[a].y()) * (vertices[c].z() - vertices[a].z()) - (vertices[b].z() - vertices[a].z()) * (vertices[c].y() - vertices[a].y());
		double y = (vertices[b].z() - vertices[a].z()) * (vertices[c].x() - vertices[a].x()) - (vertices[b].x() - vertices[a].x()) * (vertices[c].z() - vertices[a].z());
		double z = (vertices[b].x() - vertices[a].x()) * (vertices[c].y() - vertices[a].y()) - (vertices[b].y() - vertices[a].y()) * (vertices[c].x() - vertices[a].x());
		
		double length = sqrt(x * x + y * y + z * z);
		
		double surface = 0.5 * length;
		
		if(length > 0.0) {
		
	 	 x /= length;
	 	y /= length;
 		z /= length;
			
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


