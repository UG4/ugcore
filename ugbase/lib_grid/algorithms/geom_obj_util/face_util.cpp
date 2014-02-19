//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m01 d15

#include <vector>
#include <stack>
#include <cassert>
#include "face_util.h"
#include "../attachment_util.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
int GetFaceIndex(Volume* vol, Face* f)
{
	size_t numFaces = vol->num_faces();
	FaceDescriptor fd;
	for(uint i = 0; i < numFaces; ++i)
	{
		vol->face_desc(i, fd);
		if(CompareVertices(f, &fd))
			return (int)i;
	}
	return -1;
}


////////////////////////////////////////////////////////////////////////
//	CalculateNormal
void CalculateNormal(vector3& vNormOut, FaceVertices* face,
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
	return;
}

void CalculateNormalNoNormalize(vector3& vNormOut, FaceVertices* face,
								Grid::AttachmentAccessor<Vertex, APosition>& aaPos)
{
	if(face->num_vertices() == 3)
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
	else if(face->num_vertices() > 4)
	{
		CalculateTriangleNormalNoNormalize(vNormOut, aaPos[face->vertex(0)],
								aaPos[face->vertex(1)], aaPos[face->vertex(2)]);
		return;
	}

	vNormOut = vector3(0, 0, 0);
	return;
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

	Grid::VertexAttachmentAccessor<AVector3> aaPos(grid, aPos);
	Grid::FaceAttachmentAccessor<AVector3> aaNorm(grid, aNorm);

//	iterate through the specified faces and calculate normals.
	for(FaceIterator iter = facesBegin; iter != facesEnd; ++iter)
		CalculateNormal(aaNorm[*iter], *iter, aaPos);
}

////////////////////////////////////////////////////////////////////////
bool IsVolumeBoundaryFace(Grid& grid, Face* f)
{
//	check if FACEOPT_STORE_ASSOCIATED_VOLUMES is enabled.
//	if so, use it to count the number of adjacent volumes.
	int counter = 0;
	if(grid.option_is_enabled(FACEOPT_STORE_ASSOCIATED_VOLUMES))
	{
		for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(f);
			iter != grid.associated_volumes_end(f); iter++)
		{
			counter++;
		}
	}
	else
	{
	//	iterate over all volumes which are connected to the first vertex
	//	and check if they contain the face...
		Grid::AssociatedVolumeIterator iterEnd = grid.associated_volumes_end(f->vertex(0));
		for(Grid::AssociatedVolumeIterator iter = grid.associated_volumes_begin(f->vertex(0));
			iter != iterEnd; iter++)
		{
			if(VolumeContains(*iter, f))
				counter++;
		}
	}

//	if there is only one adjacent volume, the triangle is a boundary triangle
	if(counter == 1)
		return true;
	return false;
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
number FaceQuality(Face* f,
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
	//	since at least one angle is <= 60�, we have to normalize the return value
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

//	since at least one angle is <= 60�, we have to normalize the return value
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
	//	check which way produces the better triangles.
		number q1 = std::min(
				TriangleQuality(aaPos[q->vertex(0)], aaPos[q->vertex(1)], aaPos[q->vertex(2)]),
				TriangleQuality(aaPos[q->vertex(2)], aaPos[q->vertex(3)], aaPos[q->vertex(0)]));

		number q2 = std::min(
				TriangleQuality(aaPos[q->vertex(1)], aaPos[q->vertex(2)], aaPos[q->vertex(3)]),
				TriangleQuality(aaPos[q->vertex(3)], aaPos[q->vertex(0)], aaPos[q->vertex(1)]));

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
		iterBegin++;
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
	if(!grid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
	{
		LOG("WARNING: autoenabling FACEOPT_AUTOGENERATE_EDGES in GetNeighbours(Face).\n");
		grid.enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}
	
	if(clearContainer)
		vFacesOut.clear();
		
//	check if we can find an edge for the specified side
	EdgeBase* e = grid.get_edge(f, side);
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

bool EdgeOrientationMatches(EdgeVertices* ev, Face* f)
{
//	find the first vertex of ed in f
	size_t i;
	for(i = 0; i < f->num_vertices(); ++i)
	{
		if(f->vertex(i) == ev->vertex(0))
			break;
	}

	if(i < f->num_vertices())
	{
	//	the first one has been found.
	//	check whether the second vertex of ed is the
	//	same as the next vertex of f
		if(ev->vertex(1) == f->vertex((i+1)%f->num_vertices()))
			return true;//	the orientation is the same
	}

//	the orientation is not the same.
	return false;
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
