//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d17

#include <vector>
#include "tetrahedralization.h"
#include "lib_grid/externals/include/tetgen/tetgen.h"
#include "../geom_obj_util/geom_obj_util.h"

using namespace std;

namespace ug
{

static bool PerformTetrahedralization(Grid& grid,
										SubsetHandler* pSH,
										APosition& aPos)
{
	if(!grid.has_vertex_attachment(aPos))
		return false;

	if(grid.num_faces() == 0)
		return false;

//	access position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
//	make sure that no double-points exist
	RemoveDoubles(grid, grid.vertices_begin(), grid.vertices_end(), aPos, 1e-12);
//	convert quadrilaterals to triangles
	if(grid.num<Quadrilateral>() > 0)
		Triangulate(grid, grid.begin<Quadrilateral>(), grid.end<Quadrilateral>(), &aaPos);

//	attach an index to the vertices
	AInt aInd;
	grid.attach_to_vertices(aInd);
	Grid::VertexAttachmentAccessor<AInt> aaInd(grid, aInd);

//	datastructures to communicate with tetgenio
	tetgenio in, out;

//	setup points
	{
		in.numberofpoints = grid.num_vertices();
		in.pointlist = new REAL[in.numberofpoints*3];

	//	copy position data
		int counter = 0;

		for(VertexBaseIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaInd[*iter] = counter;
			vector3& v = aaPos[*iter];
			in.pointlist[counter * 3] = v.x;
			in.pointlist[counter * 3 + 1] = v.y;
			in.pointlist[counter * 3 + 2] = v.z;
		}
	}

//	set up facets
	{
		in.numberoffacets = grid.num_faces();
		in.facetlist = new tetgenio::facet[in.numberoffacets];
		in.facetmarkerlist = new int[in.numberoffacets];

		int counter = 0;
		for(FaceIterator iter = grid.faces_begin();
			iter != grid.faces_end(); ++iter, ++counter)
		{
			Face* f = *iter;

		//	copy the face to the facetlist
			tetgenio::facet* tf = &in.facetlist[counter];
			tf->numberofpolygons = 1;
			tf->polygonlist = new tetgenio::polygon[tf->numberofpolygons];
			tf->numberofholes = 0;
			tf->holelist = NULL;
			tetgenio::polygon* p = &tf->polygonlist[0];
			p->numberofvertices = f->num_vertices();
			p->vertexlist = new int[p->numberofvertices];
			for(int i = 0; i < p->numberofvertices; ++i)
				p->vertexlist[i] = aaInd[f->vertex(i)];

		//	set the face mark
			if(pSH)
				in.facetmarkerlist[counter] = pSH->get_subset_index(f);
			else
				in.facetmarkerlist[counter] = 0;
		}
	}

//	call tetrahedralization
	tetrahedralize("pqYYQ", &in, &out);

	if(out.numberofpoints < (int)grid.num_vertices()){
		LOG("  aborting tetrahedralization - bad number of points\n");
		return false;
	}

//	add new vertices to the grid. store all vertices in a vector.
	vector<VertexBase*> vVrts(out.numberofpoints);
	{
		int counter = 0;
	//	add the old ones to the vector
		for(VertexBaseIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			vVrts[counter] = *iter;
		}
	//	create new ones and add them to the vector
		for(; counter < out.numberofpoints; ++counter)
		{
			Vertex* v = *grid.create<Vertex>();
			aaPos[v].x = out.pointlist[counter*3];
			aaPos[v].y = out.pointlist[counter*3+1];
			aaPos[v].z = out.pointlist[counter*3+2];
			vVrts[counter] = v;
		}
	}

/*
	if(out.numberoftrifaces > (int)grid.num<Triangle>()){
		LOG("aborting tetrahedralization - bad nuber of triangle faces.\n");
		return false;
	}
*/

/*
	grid.erase(grid.faces_begin(), grid.faces_end());

	LOG(out.numberoftrifaces << endl);
	LOG(out.numberoffacets << endl);

//	add new faces

	for(int i = 0; i < out.numberoftrifaces; ++i)
	{
		Triangle* tri = *grid.create<Triangle>(TriangleDescriptor(vVrts[out.trifacelist[i*3]],
												vVrts[out.trifacelist[i*3 + 1]],
												vVrts[out.trifacelist[i*3 + 2]]));

		if(pSH)
			pSH->assign_subset(tri, 0);
	}
*/

	if(out.numberoftetrahedra < 1)
		return false;
		
//	add new volumes
	for(int i = 0; i < out.numberoftetrahedra; ++i)
	{
		Tetrahedron* tet = *grid.create<Tetrahedron>(
								TetrahedronDescriptor(vVrts[out.tetrahedronlist[i*4]],
														vVrts[out.tetrahedronlist[i*4 + 1]],
														vVrts[out.tetrahedronlist[i*4 + 2]],
														vVrts[out.tetrahedronlist[i*4 + 3]]));
		if(pSH)
			pSH->assign_subset(tet, 0);
	}

	return true;
}


bool Tetrahedralize(Grid& grid, APosition& aPos)
{
	return PerformTetrahedralization(grid, NULL, aPos);
}

bool Tetrahedralize(Grid& grid, SubsetHandler& sh,
					APosition& aPos)
{
	return PerformTetrahedralization(grid, &sh, aPos);
}

}//	end of namespace
