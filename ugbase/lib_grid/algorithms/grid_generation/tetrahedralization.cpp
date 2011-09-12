//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d17

#include <vector>
#include <sstream>
#include "tetrahedralization.h"
#include "../geom_obj_util/geom_obj_util.h"

#ifdef UG_TETGEN
	#include "tetgen.h"
#endif

using namespace std;

namespace ug
{

static bool PerformTetrahedralization(Grid& grid,
										SubsetHandler* pSH,
										number quality,
										bool preserveBnds,
										bool preserveAll,
										APosition& aPos)
{
#ifdef UG_TETGEN
	if(!grid.has_vertex_attachment(aPos))
		return false;

	if(grid.num_faces() == 0)
		return false;

//	access position data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
/*
//	make sure that no double-points exist
	RemoveDoubles(grid, grid.vertices_begin(), grid.vertices_end(), aPos, 1e-12);
*/
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
		//	position types are casted to float, since this circumvents an
		//	error that occurs on some geometries. Somehow tetgen constructs
		//	selfintersecting facets otherwise (sometimes). I didn't really understand
		//	this behaviour yet.
		//TODO: Think about how the following code could be improved.
			/*
			in.pointlist[counter * 3] = (float)v.x;
			in.pointlist[counter * 3 + 1] = (float)v.y;
			in.pointlist[counter * 3 + 2] = (float)v.z;*/
			in.pointlist[counter * 3] = v.x;
			in.pointlist[counter * 3 + 1] = v.y;
			in.pointlist[counter * 3 + 2] = v.z;
		}
	}

//	set up facets
	{/*
		if(pSH){
		//	collect all edges that are assigned to subsets
			vector<EdgeBase*> edges;
			for(int si = 0; si < pSH->num_subsets(); ++si){
				for(EdgeBaseIterator iter = pSH->begin<EdgeBase>(si);
					iter != pSH->end<EdgeBase>(si); ++iter)
				{
					edges.push_back(*iter);
				}
			}
			
			in.numberofedges = (int)edges.size();
			in.edgelist = new int[in.numberofedges*2];
			in.edgemarkerlist = new int[in.numberofedges];
			
			for(size_t i = 0; i < edges.size(); ++i){
				in.edgelist[i*2] = aaInd[edges[i]->vertex(0)];
				in.edgelist[i*2 + 1] = aaInd[edges[i]->vertex(1)];
				in.edgemarkerlist[i] = pSH->get_subset_index(edges[i]);
			}
			UG_LOG("number of edges in: " << in.numberofedges << endl);
		}*/
		
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

//	the aInd attachment is no longer required
	grid.detach_from_vertices(aInd);
	aaInd.invalidate();

//	call tetrahedralization
	try{
		stringstream ss;
		ss << "p";
		if(quality > SMALL)
			ss << "qq" << quality;
		if(preserveBnds || preserveAll)
			ss << "Y";
		if(preserveAll)
			ss << "Y";	// if inner bnds shall be preserved "YY" has to be passed to tetgen

		ss << "Q";
		tetrahedralize(const_cast<char*>(ss.str().c_str()), &in, &out);
	}
	catch(int errCode){
		UG_LOG("  aborting tetrahedralization. Received error: " << errCode << endl);
		return false;
	}

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
			aaPos[*iter].x = out.pointlist[counter*3];
			aaPos[*iter].y = out.pointlist[counter*3+1];
			aaPos[*iter].z = out.pointlist[counter*3+2];
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

//	erase edges if boundary segments were not preserved
	if(!preserveAll){
		grid.erase(grid.edges_begin(), grid.edges_end());
	}

//	add new faces
	grid.erase(grid.faces_begin(), grid.faces_end());
	for(int i = 0; i < out.numberoftrifaces; ++i)
	{
		Triangle* tri = *grid.create<Triangle>(TriangleDescriptor(vVrts[out.trifacelist[i*3]],
												vVrts[out.trifacelist[i*3 + 1]],
												vVrts[out.trifacelist[i*3 + 2]]));

		if(pSH)
			pSH->assign_subset(tri, out.trifacemarkerlist[i]);
	}

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

#else
	UG_LOG("WARNING in PerformTetrahedralization: Tetgen is not available in the "
			"current build. Please consider recompiling with Tetgen support enabled.\n");
	return false;
#endif

}


bool Tetrahedralize(Grid& grid, number quality, bool preserveBnds,
					bool preserveAll, APosition& aPos)
{
	return PerformTetrahedralization(grid, NULL, quality, preserveBnds,
									preserveAll, aPos);
}

bool Tetrahedralize(Grid& grid, SubsetHandler& sh, number quality,
					bool preserveBnds, bool preserveAll, APosition& aPos)
{
	return PerformTetrahedralization(grid, &sh, quality, preserveBnds,
									preserveAll, aPos);
}

}//	end of namespace
