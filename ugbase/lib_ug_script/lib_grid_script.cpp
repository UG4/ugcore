// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d19

#include "luabind/luabind.hpp"
#include "lib_grid/lib_grid.h"
#include "tmp_lib_grid_methods.h"
#include "lib_grid/tools/marker_points.h"
#include "ug_script.h"

using namespace std;
using namespace ug;

#define CHECK_POINTER(grid, retVal) {if(!grid){LOG("invalid pointer\n"); return retVal;}}

namespace lgscript
{
////////////////////////////////////////////////////////////////////////
//	marker points
MarkerPointManager* new_markers()
{
	return new MarkerPointManager;
}

void delete_markers(MarkerPointManager* mpm)
{
	CHECK_POINTER(mpm, );
	delete mpm;
}

bool load_markers(MarkerPointManager* mpm, const char* filename)
{
	CHECK_POINTER(mpm, false);
	return LoadMarkerPointsFromFile(*mpm, filename);
}

bool snap_markers_to_vertices(MarkerPointManager* mpm, Grid* grid,
									number normalOffset)
{
	CHECK_POINTER(mpm, false);
	CHECK_POINTER(grid, false);
	if(!grid->has_vertex_attachment(aPosition))
		return false;

	Grid::VertexAttachmentAccessor<APosition> aaPos(*grid, aPosition);
	
	for(size_t i = 0; i < mpm->num_markers(); ++i)
		SnapMarkerPointToGridVertex(mpm->get_marker(i), *grid, normalOffset, aaPos);
	
	return true;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
Grid* new_grid()
{
	Grid* grid = new Grid(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION);
	LOG("  grid created\n");
	return grid;
}

Grid* clone_grid(Grid* grid)
{
	CHECK_POINTER(grid, NULL);
	Grid* clone = new Grid(*grid);
	LOG("  grid cloned\n");
	return clone;
}

void copy_grid(Grid* dest, Grid* src)
{
	CHECK_POINTER(dest, );
	CHECK_POINTER(src, );
	
	*dest = *src;
	
	LOG("  grid copied\n");
}

void delete_grid(Grid* grid)
{
	CHECK_POINTER(grid,);
	LOG("  grid erased\n");
	delete grid;
}

MultiGrid* new_multi_grid()
{
	MultiGrid* multigrid = new MultiGrid(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION);
	LOG("  multi-grid created\n");
	return multigrid;
}

void delete_multi_grid(MultiGrid* multigrid)
{
	CHECK_POINTER(multigrid,);
	LOG("  multi-grid erased\n");
	delete multigrid;
}

SubsetHandler* new_subset_handler(Grid* grid)
{
	CHECK_POINTER(grid, NULL);
	LOG("  subset-handler created\n");
	return new SubsetHandler(*grid);
}

void copy_subset_handler(SubsetHandler* dest, SubsetHandler* src)
{
	CHECK_POINTER(dest, );
	CHECK_POINTER(src, );
	
	*dest = *src;
	LOG("  subset-handler copied\n");
}

void delete_subset_handler(SubsetHandler* sh)
{
	CHECK_POINTER(sh, );
	LOG("  subset-handler erased\n");
	delete sh;
}

MultiGridSubsetHandler* new_multi_grid_subset_handler(MultiGrid* multigrid)
{
	CHECK_POINTER(multigrid, NULL);
	LOG("  multi-grid-subset-handler created\n");
	return new MultiGridSubsetHandler(*multigrid);
}

void delete_multi_grid_subset_handler(MultiGridSubsetHandler* mg_sh)
{
	CHECK_POINTER(mg_sh, );
	LOG("  multi-grid-subset-handler erased\n");
	delete mg_sh;
}

bool load_grid(Grid* grid, SubsetHandler* sh, const char* filename)
{
	CHECK_POINTER(grid, false);

	LOG("  loading grid from file " << filename << " ... ");
	bool bSuccess = false;

	if(sh)
		bSuccess = LoadGridFromFile(*grid, filename, *sh);
	else
		bSuccess = LoadGridFromFile(*grid, filename);

	if(bSuccess)
	{
		LOG("done\n");
	}
	else
	{
		LOG("failed! Check the filenames!\n");
	}

	return bSuccess;
}

bool load_multi_grid(MultiGrid* multigrid, MultiGridSubsetHandler* mg_sh, const char* filename)
{
	CHECK_POINTER(multigrid, false);

	LOG("  loading grid from file " << filename << " ... ");
	bool bSuccess = false;

	if(mg_sh)
		bSuccess = LoadGridFromFile(*multigrid, filename, *mg_sh);
	else
		bSuccess = LoadGridFromFile(*multigrid, filename);

	if(bSuccess)
	{
		LOG("done\n");
	}
	else
	{
		LOG("failed! Check the filenames!\n");
	}

	return bSuccess;
}

bool save_grid(Grid* grid, SubsetHandler* sh, const char* filename)
{
	CHECK_POINTER(grid, false);

	LOG("  saving grid to file " << filename << " ... ");
	bool bSuccess = false;

	if(sh)
	{
		bSuccess = SaveGridToFile(*grid, filename, *sh);
	}
	else
		bSuccess = SaveGridToFile(*grid, filename);

	if(bSuccess)
	{
		LOG("done\n");
	}
	else
	{
		LOG("failed! Check the filenames!\n");
	}

	return bSuccess;
}

bool save_marked_edges_to_obj(Grid* grid, SubsetHandler* sh,
								const char* filename)
{
	CHECK_POINTER(grid, false);
	CHECK_POINTER(sh, false);
	LOG("  saving marked edges to " << filename << "...");
	bool retVal = SaveMarkedEdgesToObj(*grid, filename, *sh);
	LOG(" done\n");
	return retVal;
}

bool adjust_subsets_for_ug3(Grid* pGrid, SubsetHandler* pSH)
{
	CHECK_POINTER(pGrid, false);
	CHECK_POINTER(pSH, false);
	
	AdjustSubsetsForLgmNg(*pGrid, *pSH);
	
	return true;
}

bool export_grid_to_ug3(Grid* pGrid, SubsetHandler* pSH,
				const char* fileNamePrefix, const char* lgmName,
				const char* problemName, int convex)
{
	CHECK_POINTER(pGrid, false);
	CHECK_POINTER(pSH, false);

	UG_LOG("  exporting to ug3. writing " << fileNamePrefix
			<< " (.lgm, .ng)\n");

	Grid& grid = *pGrid;
	SubsetHandler& sh = *pSH;
	
//	we have to separate volume- and face-subsets
	SubsetHandler shFaces(grid, SHE_FACE);
	SubsetHandler shVols(grid, SHE_VOLUME);
	
//	assign subsets
	for(size_t i = 0; i < sh.num_subsets(); ++i){
		
		for(FaceIterator iter = sh.begin<Face>(i);
			iter != sh.end<Face>(i); ++iter)
		{
			shFaces.assign_subset(*iter, i);
		}
		
		for(VolumeIterator iter = sh.begin<Volume>(i);
			iter != sh.end<Volume>(i); ++iter)
		{
			shVols.assign_subset(*iter, i);
		}
	}

	return ExportGridToUG(grid, shFaces, shVols, fileNamePrefix,
						lgmName, problemName, convex);
}

bool triangulate(Grid* grid)
{
	CHECK_POINTER(grid, false);
	LOG("  triangulating...");
	Triangulate(*grid, grid->begin<Quadrilateral>(),
				grid->end<Quadrilateral>());
	LOG(" done\n");
	return true;
}

int remove_doubles(Grid* grid, number threshold)
{
	CHECK_POINTER(grid, 0);
	LOG("  removing doubles... ");
	size_t numVrtsAtStart = grid->num_vertices();
	
	RemoveDoubles(*grid, grid->vertices_begin(), grid->vertices_end(),
					aPosition, threshold);
	
	int numRemoved = int(numVrtsAtStart - grid->num_vertices());
	LOG(numRemoved << endl);
	
	return numRemoved;
}

void mark_crease_elements(Grid* grid, SubsetHandler* shMarks,
						float creaseAngle)
{
	CHECK_POINTER(grid, );
	CHECK_POINTER(shMarks, );

	LOG("  marking crease elements...");

	if(!grid->option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
	{
		LOG("  INFO: auto-enabling FACEOPT_AUTOGENERATE_EDGES in mark_crease_elements.\n");
		grid->enable_options(FACEOPT_AUTOGENERATE_EDGES);
	}

	MarkCreaseEdges(*grid, *shMarks, grid->edges_begin(),
					grid->edges_end(), RM_CREASE, creaseAngle);

	MarkFixedCreaseVertices(*grid, *shMarks, RM_CREASE, RM_FIXED);

	LOG(" got " << shMarks->num<EdgeBase>(RM_CREASE) << " crease edges.");
	LOG(" done\n");
}

bool adjust_edge_length(Grid* gridOut, SubsetHandler* shOut, SubsetHandler* shMarksOut,
						Grid* gridIn, SubsetHandler* shIn, SubsetHandler* shMarksIn,
						float minEdgeLen, float maxEdgeLen, int numIterations, bool bProjectPoints)
{
	CHECK_POINTER(gridOut, false);
	CHECK_POINTER(gridIn, false);

	bool bSuccess = false;

	LOG("  adjusting edge length...\n");

//	those handlers are required if a null-pointer was passed to one of the subset-handlers
	SubsetHandler tmpMarksIn;
	SubsetHandler tmpMarksOut;
	SubsetHandler tmpSHIn;
	SubsetHandler tmpSHOut;

	if(!shMarksIn){
		LOG("  no marks in. generating temporary mark-handler.\n");
		tmpMarksIn.assign_grid(*gridIn);
		shMarksIn = &tmpMarksIn;
	}
	if(!shMarksOut){
		LOG("  no marks out. generating temporary mark-handler.\n");
		tmpMarksOut.assign_grid(*gridOut);
		shMarksOut = &tmpMarksOut;
	}
	if(!shIn){
		LOG("  no subsets in. generating temporary subset-handler.\n");
		tmpSHIn.assign_grid(*gridIn);
		shIn = &tmpSHIn;
	}
	if(!shOut){
		LOG("  no subsets out. generating temporary subset-handler.\n");
		tmpSHOut.assign_grid(*gridOut);
		shOut = &tmpSHOut;
	}

	bSuccess = AdjustEdgeLength(*gridOut, *shOut, *shMarksOut, *gridIn, *shIn, *shMarksIn,
								minEdgeLen, maxEdgeLen, numIterations, bProjectPoints);

	LOG("  done\n");
	return bSuccess;
}

bool extrude_cylinder(Grid* grid, SubsetHandler* sh, double x, double y, double z,
					 double height, double radius)
{
	CHECK_POINTER(grid, false);
	
	bool bSuccess = false;

	Grid& g = *grid;

	if(g.num<VertexBase>() == 0)
		return false;

	Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
	if(!aaPos.valid())
		return false;

//	the specified point
	vector3 p(x, y, z);

//	find the closest point in the grid
	VertexBase* vrt;
	{
		VertexBaseIterator iter = g.begin<VertexBase>();
		vrt = *iter;
		number distSq = VecDistanceSq(aaPos[vrt], p);
		iter++;
		for(; iter != g.end<VertexBase>(); ++iter){
			number d = VecDistanceSq(aaPos[*iter], p);
			if(d < distSq){
				distSq = d;
				vrt = *iter;
			}
		}
	}

//	calculate the normal of the vertex
	vector3 n;
	CalculateVertexNormal(n, g, vrt, aaPos);

	if(sh)
	{
		int numSubs = sh->num_subsets();
		bSuccess = ExtrudeCylinder(g, *sh, vrt, n, height, radius, aaPos,
								numSubs, numSubs + 1);
	}
	else
		bSuccess = ExtrudeCylinder(g, vrt, n, height, radius, aaPos);
						
	return bSuccess;
}

bool extrude_cylinders(Grid* grid, SubsetHandler* sh, MarkerPointManager* mgm,
					 double height, double radius)
{
	CHECK_POINTER(grid, false);
	CHECK_POINTER(mgm, false);

//	iterate through the vertices of pPointGrid and extrude
//	a cylinder for each
	bool bSuccess = true;
	for(size_t i = 0; i < mgm->num_markers(); ++i)
	{
		vector3& v = mgm->get_marker(i).pos;
		bSuccess &= extrude_cylinder(grid, sh, v.x, v.y, v.z, height, radius);
	}

	return bSuccess;
}

bool tetrahedralize(Grid* grid, SubsetHandler* sh)
{
	CHECK_POINTER(grid, false);
	
	LOG("  tetrahedralizing...");
	bool bSuccess = false;

	if(sh)
		bSuccess = Tetrahedralize(*grid, *sh);
	else
		bSuccess = Tetrahedralize(*grid);

	if(bSuccess){
		LOG(" done\n");
	}
	else{
		LOG(" failed\n");
	}

	return bSuccess;
}

bool separate_regions(Grid* grid, SubsetHandler* shVolsOut,
					  SubsetHandler* shFaces, MarkerPointManager* mpm,
					  int firstSubsetIndex)
{
	CHECK_POINTER(grid, false);
	CHECK_POINTER(shVolsOut, false);
	CHECK_POINTER(shFaces, false);
	CHECK_POINTER(mpm, false);

	LOG("  separating regions...");
	if(SeparateRegions(*grid, *shVolsOut, *shFaces, *mpm, firstSubsetIndex))
	{
		LOG(" done\n");
		return true;
	}
	else{
		LOG(" failed\n");
		return false;
	}
}

bool fracture_to_subset(Grid* pGrid, SubsetHandler* pSH,
						int subsetIndex, number threshold)
{
//	make sure that FACEOPT_AUTOGENERATE_EDGES is enabled
	CHECK_POINTER(pGrid, false);
	CHECK_POINTER(pSH, false);
	
	Grid& grid = *pGrid;
	SubsetHandler& sh = *pSH;

//	compare squares
	number maxDistSq = threshold * threshold;
	
//	access the position attachment
	if(!grid.has_vertex_attachment(aPosition))
		return false;
		
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPosition);
	
//	iterate over all faces and check whether they contain
//	degenerated edges
	EdgeDescriptor ed;
	for(FaceIterator iter = grid.begin<Face>();
		iter != grid.end<Face>(); ++iter)
	{
		Face* f = *iter;
		
		for(size_t i = 0; i < f->num_edges(); ++i){
			f->edge(i, ed);
			if(VecDistanceSq(aaPos[ed.vertex(0)],
							aaPos[ed.vertex(1)]) < maxDistSq)
			{
			//	the element belongs to a fracture
				sh.assign_subset(f, subsetIndex);
			}
		}
	}
	
//	done
	return true;
}

////////////////////////////////////////////////////////////////////////
//	frac_extrude
bool frac_extrude(Grid* pGrid, SubsetHandler* pSH, int numSecs,
					number dirX, number dirY, number dirZ)
{
	CHECK_POINTER(pGrid, false);
	CHECK_POINTER(pSH, false);
	
	Grid& grid = *pGrid;
	SubsetHandler& sh = *pSH;
	
	if(!grid.option_is_enabled(GRIDOPT_STANDARD_INTERCONNECTION))
	{
		UG_LOG("  INFO in frac_extrude: auto-enabling GRIDOPT_STANDARD_INTERCONNECTION.\n");
		grid.enable_options(GRIDOPT_STANDARD_INTERCONNECTION);
	}
	
//	first of all we'll enable auto-selection of the subset-handler
	bool subsetInheritanceWasEnabled = sh.subset_inheritance_enabled();
	sh.enable_subset_inheritance(true);

//	perform the extrusion steps
//	this vector hold the faces that we want to extrude.
	vector<Face*> vExtrudeFaces;
	vExtrudeFaces.assign(grid.faces_begin(), grid.faces_end());

//	we need a selector that records all newly created edges, faces and volumes.
	Selector sel(grid);
	sel.enable_autoselection(true);

//	perform the extrusion of the cell
//	volumes derive their subset-index from the faces that they were extruded from
	for(int i = 0; i < numSecs; ++i){
		Extrude(grid, NULL, NULL, &vExtrudeFaces, vector3(dirX, dirY, dirZ));
	}

//	remove all edges and faces from the subset-handler that have been created by extrusion.
	sh.assign_subset(sel.begin<EdgeBase>(), sel.end<EdgeBase>(), -1);
	sh.assign_subset(sel.begin<Face>(), sel.end<Face>(), -1);
	
//	clean up
	sh.enable_subset_inheritance(subsetInheritanceWasEnabled);
	
	return true;
}

}//	end of namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	register the script-module at the module manager.
//	do this in this empty namespace
//	please note that this trick will not work if you link this file
//	into any form of library. Only if it is compiled and directly linked
//	to the executable, this code will actually be executed.
//	This is due the automatic removal of unreferenced code by the linker.
namespace ug{
namespace script
{

//	initialises lib-grids script methods
bool InitLibGridScript(lua_State* L)
{
	using namespace lgscript;
	
//	bind classes
	luabind::module(L)[luabind::class_<Grid>("Grid")];
	luabind::module(L)[luabind::class_<MultiGrid>("MultiGrid")];

	luabind::module(L)[luabind::class_<SubsetHandler>("SubsetHandler")];
	luabind::module(L)[luabind::class_<MultiGridSubsetHandler>("MultiGridSubsetHandler")];

	luabind::module(L)[luabind::class_<MarkerPointManager>("MarkerPointManager")];

//	bind methods to lua
	luabind::module(L)[luabind::def("new_grid", new_grid)];
	luabind::module(L)[luabind::def("clone_grid", clone_grid)];
	luabind::module(L)[luabind::def("copy_grid", copy_grid)];
	luabind::module(L)[luabind::def("delete_grid", delete_grid)];
	luabind::module(L)[luabind::def("new_multi_grid", new_multi_grid)];
	luabind::module(L)[luabind::def("delete_multi_grid", delete_multi_grid)];

	luabind::module(L)[luabind::def("new_subset_handler", new_subset_handler)];
	luabind::module(L)[luabind::def("copy_subset_handler", copy_subset_handler)];
	luabind::module(L)[luabind::def("delete_subset_handler", delete_subset_handler)];
	luabind::module(L)[luabind::def("new_multi_grid_subset_handler", new_multi_grid_subset_handler)];
	luabind::module(L)[luabind::def("delete_multi_grid_subset_handler", delete_multi_grid_subset_handler)];

	luabind::module(L)[luabind::def("load_grid", load_grid)];
	luabind::module(L)[luabind::def("save_grid", save_grid)];
	luabind::module(L)[luabind::def("load_multi_grid", load_multi_grid)];
	luabind::module(L)[luabind::def("adjust_subsets_for_ug3", adjust_subsets_for_ug3)];
	luabind::module(L)[luabind::def("export_grid_to_ug3", export_grid_to_ug3)];

	luabind::module(L)[luabind::def("save_marked_edges_to_obj", save_marked_edges_to_obj)];
	
	luabind::module(L)[luabind::def("remove_doubles", remove_doubles)];
	luabind::module(L)[luabind::def("convert_to_triangle_grid", triangulate)];
	luabind::module(L)[luabind::def("triangulate", triangulate)];
	luabind::module(L)[luabind::def("mark_crease_elements", mark_crease_elements)];
	luabind::module(L)[luabind::def("adjust_edge_length", adjust_edge_length)];

	luabind::module(L)[luabind::def("extrude_cylinder", extrude_cylinder)];
	luabind::module(L)[luabind::def("extrude_cylinders", extrude_cylinders)];
	
	luabind::module(L)[luabind::def("fracture_to_subset", fracture_to_subset)];
	luabind::module(L)[luabind::def("frac_extrude", frac_extrude)];
	
	luabind::module(L)[luabind::def("tetrahedralize", tetrahedralize)];

	luabind::module(L)[luabind::def("separate_regions", separate_regions)];
	
	luabind::module(L)[luabind::def("new_markers", new_markers)];
	luabind::module(L)[luabind::def("delete_markers", delete_markers)];
	luabind::module(L)[luabind::def("load_markers", load_markers)];
	luabind::module(L)[luabind::def("snap_markers_to_vertices", snap_markers_to_vertices)];

	return true;
}

void FinalizeLibGridScript(lua_State* L)
{
}

}}//	end of namespace


