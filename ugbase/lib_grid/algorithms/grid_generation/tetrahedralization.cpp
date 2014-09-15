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

static const char* VerbosityToTetgenParam(int verbosity)
{
	if(verbosity <= 0)
		return "";
	else if(verbosity == 1)
		return "V";
	else if(verbosity == 2)
		return "VV";
	else
		return "VVV";
}

static bool PerformTetrahedralization(Grid& grid,
										ISubsetHandler* pSH,
										number quality,
										bool preserveBnds,
										bool preserveAll,
										APosition& aPos,
										int verbosity)
{
#ifdef UG_TETGEN
	if(!grid.has_vertex_attachment(aPos))
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

	if(grid.num_edges() > 0 || grid.num_faces() > 0){
	//	we have to make sure that all vertices in the grid are actually connected
	//	to a triangle. If this wouldn't be the case, tetgen would likely have problems.
		size_t numVrtsRemoved = 0;
		for(VertexIterator iter = grid.begin<Vertex>();
			iter != grid.end<Vertex>();)
		{
			Vertex* v = *iter;
			++iter;
			if((NumAssociatedEdges(grid, v) == 0) || (NumAssociatedFaces(grid, v) == 0)){
				grid.erase(v);
				++numVrtsRemoved;
			}
		}

		if(numVrtsRemoved > 0){
			UG_LOG("WARNING in Tetrahedralize: Removed " << numVrtsRemoved
					<< " vertices which were not connected to any edge or face\n");
		}
	}

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

		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaInd[*iter] = counter;
			vector3& v = aaPos[*iter];
		//	position types are casted to float, since this circumvents an
		//	error that occurs on some geometries. Somehow tetgen constructs
		//	selfintersecting facets otherwise (sometimes). I didn't really understand
		//	this behaviour yet.
		//TODO: Think about how the following code could be improved.

			in.pointlist[counter * 3] = (float)v.x();
			in.pointlist[counter * 3 + 1] = (float)v.y();
			in.pointlist[counter * 3 + 2] = (float)v.z();
		}
	}

//	set up facets
	{/*
		if(pSH){
		//	collect all edges that are assigned to subsets
			vector<Edge*> edges;
			for(int si = 0; si < pSH->num_subsets(); ++si){
				for(EdgeIterator iter = pSH->begin<Edge>(si);
					iter != pSH->end<Edge>(si); ++iter)
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

		if(grid.num_faces() > 0){
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
	}

//	the aInd attachment is no longer required
	grid.detach_from_vertices(aInd);
	aaInd.invalidate();

//	call tetrahedralization
	try{
		stringstream ss;
		ss << VerbosityToTetgenParam(verbosity);
		if(grid.num_faces() > 0){
			ss << "p";
			if(quality > SMALL)
				ss << "qq" << quality;
			if(preserveBnds || preserveAll)
				ss << "Y";
			if(preserveAll)
				ss << "Y";	// if inner bnds shall be preserved "YY" has to be passed to tetgen
		}

		ss << "Q";//"Q";
		tetrahedralize(const_cast<char*>(ss.str().c_str()), &in, &out);
	}
	catch(int errCode){
		UG_LOG("  aborting tetrahedralization. Received error: " << errCode << endl);
		return false;
	}
/*
	if(out.numberofpoints < (int)grid.num_vertices()){
		LOG("  aborting tetrahedralization - bad number of points\n");
		return false;
	}
*/
//	add new vertices to the grid. store all vertices in a vector.
	vector<Vertex*> vVrts(out.numberofpoints);
	{
		int counter = 0;
	//	add the old ones to the vector
		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaPos[*iter].x() = out.pointlist[counter*3];
			aaPos[*iter].y() = out.pointlist[counter*3+1];
			aaPos[*iter].z() = out.pointlist[counter*3+2];
			vVrts[counter] = *iter;
			if(counter == out.numberofpoints -1){
				UG_LOG("WARNING: Unused points may remain!\n");
				break;
			}
		}
	//	create new ones and add them to the vector
		for(; counter < out.numberofpoints; ++counter)
		{
			RegularVertex* v = *grid.create<RegularVertex>();
			aaPos[v].x() = out.pointlist[counter*3];
			aaPos[v].y() = out.pointlist[counter*3+1];
			aaPos[v].z() = out.pointlist[counter*3+2];
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

		if(pSH && out.trifacemarkerlist)
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


static bool PerformRetetrahedralization(Grid& grid,
										SubsetHandler& sh,
										number quality,
										bool preserveBnds,
										bool preserveAll,
										APosition& aPos,
										ANumber& aVolCon,
										bool applyVolumeConstraint,
										int verbosity)
{
#ifdef UG_TETGEN
	if(!grid.has_vertex_attachment(aPos))
		return false;

	if(!grid.has_volume_attachment(aVolCon))
		return false;

	if(grid.num<Tetrahedron>() == 0)
		return false;

//	access data
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
	Grid::VolumeAttachmentAccessor<ANumber> aaVolCon(grid, aVolCon);

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

		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaInd[*iter] = counter;
			vector3& v = aaPos[*iter];
		//	position types are casted to float, since this circumvents an
		//	error that occurs on some geometries. Somehow tetgen constructs
		//	selfintersecting facets otherwise (sometimes). I didn't really understand
		//	this behaviour yet.
		//TODO: Think about how the following code could be improved.

			in.pointlist[counter * 3] = (float)v.x();
			in.pointlist[counter * 3 + 1] = (float)v.y();
			in.pointlist[counter * 3 + 2] = (float)v.z();
			//in.pointlist[counter * 3] = v.x();
			//in.pointlist[counter * 3 + 1] = v.y();
			//in.pointlist[counter * 3 + 2] = v.z();
		}
	}

//	set up triangles
	{
	//	count the total number of triangles in subsets
		int numTris = 0;
		for(int i = 0; i < sh.num_subsets(); ++i)
			numTris += sh.num<Triangle>(i);

		in.numberoftrifaces = numTris;
		in.trifacelist = new int[in.numberoftrifaces * 3];
		in.trifacemarkerlist = new int[in.numberoftrifaces];

		int curTri = 0;
		for(TriangleIterator iter = grid.begin<Triangle>();
			iter != grid.end<Triangle>(); ++iter, ++curTri)
		{
			Triangle* t = *iter;

		//	ignore faces in subset -1
			if(sh.get_subset_index(t) == -1)
				continue;

		//	add the triangle
			in.trifacelist[curTri * 3] = aaInd[t->vertex(0)];
			in.trifacelist[curTri * 3 + 1] = aaInd[t->vertex(1)];
			in.trifacelist[curTri * 3 + 2] = aaInd[t->vertex(2)];

		//	set the face mark
			in.trifacemarkerlist[curTri] = sh.get_subset_index(t);
		}
	}

//	now fill the tetrahedron lists
	{
		in.numberoftetrahedra = grid.num<Tetrahedron>();
		in.tetrahedronlist = new int[in.numberoftetrahedra * 4];
		in.tetrahedronvolumelist = new REAL[in.numberoftetrahedra];
		in.numberoftetrahedronattributes = 1;
		in.tetrahedronattributelist = new REAL[in.numberoftetrahedra];

	//	fill the lists
		int curTet = 0;
		for(Grid::traits<Tetrahedron>::iterator iter = grid.begin<Tetrahedron>();
			iter != grid.end<Tetrahedron>(); ++iter, ++curTet)
		{
			Tetrahedron* t = *iter;
			in.tetrahedronlist[curTet * 4] = aaInd[t->vertex(0)];
			in.tetrahedronlist[curTet * 4 + 1] = aaInd[t->vertex(1)];
			in.tetrahedronlist[curTet * 4 + 2] = aaInd[t->vertex(2)];
			in.tetrahedronlist[curTet * 4 + 3] = aaInd[t->vertex(3)];

			in.tetrahedronvolumelist[curTet] = aaVolCon[t];
			in.tetrahedronattributelist[curTet] = (REAL)sh.get_subset_index(t);
		}
	}

//	the aInd attachment is no longer required
	grid.detach_from_vertices(aInd);
	aaInd.invalidate();

//	call tetrahedralization
	try{
		stringstream ss;
		ss << VerbosityToTetgenParam(verbosity);
		if(applyVolumeConstraint)
			ss << "a";
		ss << "r";
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
		UG_LOG("  aborting retetrahedralization. Received error: " << errCode << endl);
		return false;
	}

//	first we'll erase all existing tetrahedrons
	grid.erase(grid.begin<Tetrahedron>(), grid.end<Tetrahedron>());

//	add new vertices to the grid. store all vertices in a vector.
	vector<Vertex*> vVrts(out.numberofpoints);
	{
		int counter = 0;
	//	add the old ones to the vector
		for(VertexIterator iter = grid.vertices_begin();
			iter != grid.vertices_end(); ++iter, ++counter)
		{
			aaPos[*iter].x() = out.pointlist[counter*3];
			aaPos[*iter].y() = out.pointlist[counter*3+1];
			aaPos[*iter].z() = out.pointlist[counter*3+2];
			vVrts[counter] = *iter;
			if(counter == out.numberofpoints -1){
				UG_LOG("WARNING: Unused points may remain!\n");
				break;
			}
		}
	//	create new ones and add them to the vector
		for(; counter < out.numberofpoints; ++counter)
		{
			RegularVertex* v = *grid.create<RegularVertex>();
			aaPos[v].x() = out.pointlist[counter*3];
			aaPos[v].y() = out.pointlist[counter*3+1];
			aaPos[v].z() = out.pointlist[counter*3+2];
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

		sh.assign_subset(tri, out.trifacemarkerlist[i]);
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

		sh.assign_subset(tet, (int)out.tetrahedronattributelist[i]);
	}

	return true;

#else
	UG_LOG("WARNING in PerformRetetrahedralization: Tetgen is not available in the "
			"current build. Please consider recompiling with Tetgen support enabled.\n");
	return false;
#endif

}


bool Tetrahedralize(Grid& grid, number quality, bool preserveBnds,
					bool preserveAll, APosition& aPos, int verbosity)
{
	return PerformTetrahedralization(grid, NULL, quality, preserveBnds,
									preserveAll, aPos, verbosity);
}

bool Tetrahedralize(Grid& grid, ISubsetHandler& sh, number quality,
					bool preserveBnds, bool preserveAll, APosition& aPos,
					int verbosity)
{
	return PerformTetrahedralization(grid, &sh, quality, preserveBnds,
									preserveAll, aPos, verbosity);
}

bool Retetrahedralize(Grid& grid, SubsetHandler& sh,
					ANumber& aVolumeConstraint,
					number quality,
					bool preserveBnds,
					bool preserveAll,
					APosition& aPos,
					bool applyVolumeConstraint,
					int verbosity)
{
	return PerformRetetrahedralization(grid, sh, quality, preserveBnds,
									preserveAll, aPos, aVolumeConstraint,
									applyVolumeConstraint,
									verbosity);
}

}//	end of namespace
