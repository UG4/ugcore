//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m01 d15

#include <vector>
#include <stack>
#include <cassert>
#include "face_util.h"
#include "misc_util.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CalculateNormal
void CalculateNormal(vector3& vNormOut, Face* face,
					Grid::VertexAttachmentAccessor<APosition>& aaPos)
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
//	FaceQuality
number FaceQuality(Face* f,
				Grid::VertexAttachmentAccessor<APosition> aaPos)
{
//	take dot-products at the corners.
	number quality = 1;
	uint numVrts = f->num_vertices();

	vector3 v1, v2;
//	calculate the direction from first to last vertex.
	VecSubtract(v1, aaPos[f->vertex(numVrts-1)], aaPos[f->vertex(0)]);
	VecNormalize(v1, v1);

	for(uint i = 0; i < numVrts; ++i)
	{
		VecSubtract(v2, aaPos[f->vertex((i+1)%numVrts)], aaPos[f->vertex(i)]);
		VecNormalize(v2, v2);
		number nQual = 1.f - fabs(VecDot(v1, v2));
		if(nQual < quality)
			quality = nQual;
	//	v1 of the next iteration equals -v2 of this iteration
		VecScale(v1, v2, -1);
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

	return quality;
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
	EdgeBase* e = FindEdge(grid, f, side);
	assert(e && "edge not found though it should be there!");

//	get the connected faces
	vector<Face*> vFaces;
	CollectFaces(vFaces, grid, e);
	
//	push them to vFacesOut - except f
	for(int i = 0; i < vFaces.size(); ++i)
	{
		if(vFaces[i] != f)
			vFacesOut.push_back(vFaces[i]);
	}
}

bool EdgeOrientationMatches(EdgeDescriptor& ed, Face* f)
{
//	find the first vertex of ed in f
	int i;
	for(i = 0; i < f->num_vertices(); ++i)
	{
		if(f->vertex(i) == ed.vertex(0))
			break;
	}
	
	if(i < f->num_vertices())
	{
	//	the first one has been found.
	//	check whether the second vertex of ed is the
	//	same as the next vertex of f
		if(ed.vertex(1) == f->vertex((i+1)%f->num_vertices()))
			return true;//	the orientation is the same
	}
	
//	the orientation is not the same.
	return false;
}

////////////////////////////////////////////////////////////////////////
//	FixOrientation
void FixOrientation(Grid& grid, FaceIterator facesBegin, FaceIterator facesEnd)
{
//	we need a mark on each face that tells us
//	whether it has already been processed.
//	faces that shall not be orientated, will be marked with -1.
	AInt aMark;
	grid.attach_to_faces(aMark);
	Grid::FaceAttachmentAccessor<AInt> aaMark(grid, aMark);
	SetAttachmentValues(aaMark, grid.faces_begin(), grid.faces_end(), -1);
	SetAttachmentValues(aaMark, facesBegin, facesEnd, 0);
	
	
//	this edge descriptor will be used multiple times
	EdgeDescriptor ed;
	
//	containers to store neighbour elements.
	vector<Face*> vNeighbours;
	
//	stack that stores candidates
	stack<Face*> stkFaces;
	
//	we'll iterate through all faces
	while(facesBegin != facesEnd)
	{
	//	candidates are empty at this point.
	//	if the face is unprocessed it is a new candidate
		if(aaMark[*facesBegin] == 0)
		{
		//	mark it as candidate
			aaMark[*facesBegin] = 1;
			stkFaces.push(*facesBegin);
			
		//	while the stack is not empty
			while(!stkFaces.empty())
			{
			//	get the candidate
				Face* f = stkFaces.top();
				stkFaces.pop();
				
			//	get the neighbours for each side
				for(uint i = 0; i < f->num_edges(); ++i)
				{
					GetNeighbours(vNeighbours, grid, f, i);

				//	fix orientation of unprocessed neighbours.
					for(size_t j = 0; j < vNeighbours.size(); ++j)
					{
						Face* fn = vNeighbours[j];
						if(aaMark[fn] == 0)
						{
						//	check whether the orientation of f and fn differs.
							f->edge(i, ed);
							if(EdgeOrientationMatches(ed, fn))
							{
							//	the orientation of ed is the same as the orientation
							//	of an edge in fn.
							//	the faces thus have different orientation.
								grid.flip_orientation(fn);
							}
							
						//	mark the face as processed and add it to the stack
							aaMark[fn] = 1;
							stkFaces.push(fn);
						}
					}
				}
			}
		}
		
	//	check the next face
		++facesBegin;
	}
	
	grid.detach_from_faces(aMark);
}

}//	end of namespace
