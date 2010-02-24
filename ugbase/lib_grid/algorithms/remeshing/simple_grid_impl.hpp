//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d27

#ifndef __H__REMESHING__SIMPLE_GRID_IMPL__
#define __H__REMESHING__SIMPLE_GRID_IMPL__

#include <queue>

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	ObtainSimpleGrid
template <class TPosAcc, class TIntAcc, class TNormAcc>
bool ObtainSimpleGrid(SimpleGrid& sgOut, Grid& grid,
						VertexBase* vrt1, VertexBase* vrt2, size_t size,
						TPosAcc& aaPos, TNormAcc& aaNorm,
						TIntAcc& aaInt)
{
//	vVrts will be reused in each call. To avoid unnecessary allocations,
//	we'll reuse this vector.
	static std::vector<VertexBase*> vVrts;
	vVrts.clear();
	
//	clear the simple-grid
	sgOut.clear();
	
//	collect vertices and triangles
	grid.begin_marking();
	
//	mark the first two vertices and add them to simple-grid
	grid.mark(vrt1);
	grid.mark(vrt2);
	vVrts.push_back(vrt1);
	vVrts.push_back(vrt2);
	aaInt[vrt1] = 0;
	aaInt[vrt2] = 1;

//	this counter holds the next vertex for which we have to search for
//	associated triangles
	size_t nextVrt = 0;
//	this number holds the first index that is not checked for neighbours.
	size_t vrtsEnd = 2;

//	find the triangles that are adjacent to the edge between vrt1 and vrt2
//	at this point we assume that all associated faces are triangles.
//	If they are not they are simply treated as if they were some.
	FaceIterator iterEnd = grid.associated_faces_end(vrt1);
	for(FaceIterator iter = grid.associated_faces_begin(vrt1);
		iter != iterEnd; ++iter)
	{
		VertexBase* vUnmarked = NULL;
		Face* f = *iter;
		int counter = 0;
		
		for(uint j = 0; j < 3; ++j){
			if(grid.is_marked(f->vertex(j)))
				++counter;
			else
				vUnmarked = f->vertex(j);
		}
		
		if(counter > 1){
		//	we found an adjacent triangle. vUnmarked contains the connected vertex
			if(!vUnmarked) goto bail_out;
		//	push the connected vertex to vVrts and assign the index
			aaInt[vUnmarked] = vVrts.size();
			vVrts.push_back(vUnmarked);
		//	add the triangle
			sgOut.triangles.push_back(aaInt[f->vertex(0)]);
			sgOut.triangles.push_back(aaInt[f->vertex(1)]);
			sgOut.triangles.push_back(aaInt[f->vertex(2)]);
		//	mark the face
			grid.mark(f);
		}
	}
	
//	mark the vertices in vVrts that are not yet marked
	for(size_t i = nextVrt; i < vVrts.size(); ++i)
		grid.mark(vVrts[i]);

//	collect all faces in the neighbourhood
	for(size_t i = 0; i < size; ++i)
	{		
		for(; nextVrt < vrtsEnd; ++nextVrt)
		{
			VertexBase* vrt = vVrts[nextVrt];
		//	colelct neighbour faces
			FaceIterator iterEnd = grid.associated_faces_end(vrt);
			for(FaceIterator iter = grid.associated_faces_begin(vrt);
				iter != iterEnd; ++iter)
			{
				Face* f = *iter;
			//	if f is unmarked
				if(!grid.is_marked(f)){
				//	add unmarked vertices to vVrts
					for(uint j = 0; j < 3; ++j){
						if(!grid.is_marked(f->vertex(j))){
							aaInt[f->vertex(j)] = vVrts.size();
							grid.mark(f->vertex(j));
							vVrts.push_back(f->vertex(j));
						}
					}

				//	add the triangle
					grid.mark(f);
					sgOut.triangles.push_back(aaInt[f->vertex(0)]);
					sgOut.triangles.push_back(aaInt[f->vertex(1)]);
					sgOut.triangles.push_back(aaInt[f->vertex(2)]);
				}
			}
		}
	//	in the next iteration we'll check all vertices up to this point
		vrtsEnd = vVrts.size();
	}
	
//	copy the vertex-positions and the normals to the grid
	for(size_t i = 0; i < vVrts.size(); ++i)
	{
		sgOut.vertices.push_back(aaPos[vVrts[i]]);
		sgOut.vertexNormals.push_back(aaNorm[vVrts[i]]);
	}
	
//	calculate triangle normals
	CalculateTriangleNormals(sgOut);
	
	grid.end_marking();
	return true;
	
bail_out:
	grid.end_marking();
	return false;
}

}//	end of namespace

#endif
