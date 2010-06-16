// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d19

#include "luabind/luabind.hpp"
#include "lib_grid/lib_grid.h"
#include "tmp_lib_grid_methods.h"
#include "lib_grid/tools/marker_points.h"
#include "ug_script.h"
#include "lib_grid/file_io/file_io_ugx.h"

using namespace std;
using namespace ug;

#define CHECK_POINTER(grid, retVal) {if(!grid){LOG("invalid pointer\n"); return retVal;}}

///	contains all lib_grid related script-functions.
namespace lgscript
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup objectCreation Object Creation
 * @{
 */
 
////////////////////////////////////////////////////////////////////////
//	marker points
///	Creates a MarkerPointManager.
MarkerPointManager* new_markers()
{
	return new MarkerPointManager;
}

///	Deletes a previously created MarkerPointManager.
void delete_markers(MarkerPointManager* mpm)
{
	CHECK_POINTER(mpm, );
	delete mpm;
}

////////////////////////////////////////////////////////////////////////
///	Creates an empty grid.
Grid* new_grid()
{
	Grid* grid = new Grid(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION);
	LOG("  grid created\n");
	return grid;
}

///	Clones a grid.
Grid* clone_grid(Grid* grid)
{
	CHECK_POINTER(grid, NULL);
	Grid* clone = new Grid(*grid);
	LOG("  grid cloned\n");
	return clone;
}

///	Copies a grid to another.
void copy_grid(Grid* dest, Grid* src)
{
	CHECK_POINTER(dest, );
	CHECK_POINTER(src, );
	
	*dest = *src;
	
	LOG("  grid copied\n");
}

///	Deletes a previously created grid.
void delete_grid(Grid* grid)
{
	CHECK_POINTER(grid,);
	LOG("  grid erased\n");
	delete grid;
}

///	Creates a MultiGrid.
MultiGrid* new_multi_grid()
{
	MultiGrid* multigrid = new MultiGrid(GRIDOPT_VERTEXCENTRIC_INTERCONNECTION);
	LOG("  multi-grid created\n");
	return multigrid;
}

///	Deletes a previously created MultiGrid.
void delete_multi_grid(MultiGrid* multigrid)
{
	CHECK_POINTER(multigrid,);
	LOG("  multi-grid erased\n");
	delete multigrid;
}

///	Creates a SubsetHandler.
SubsetHandler* new_subset_handler(Grid* grid)
{
	CHECK_POINTER(grid, NULL);
	LOG("  subset-handler created\n");
	return new SubsetHandler(*grid);
}

///	Copies a SubsetHandler to another SubsetHandler.
void copy_subset_handler(SubsetHandler* dest, SubsetHandler* src)
{
	CHECK_POINTER(dest, );
	CHECK_POINTER(src, );
	
	*dest = *src;
	LOG("  subset-handler copied\n");
}

///	Deletes a SubsetHandler.
void delete_subset_handler(SubsetHandler* sh)
{
	CHECK_POINTER(sh, );
	LOG("  subset-handler erased\n");
	delete sh;
}

///	Creates a MultiGridSubsetHandler.
MultiGridSubsetHandler* new_multi_grid_subset_handler(MultiGrid* multigrid)
{
	CHECK_POINTER(multigrid, NULL);
	LOG("  multi-grid-subset-handler created\n");
	return new MultiGridSubsetHandler(*multigrid);
}

///	Deletes a MultiGridSubsetHandler.
void delete_multi_grid_subset_handler(MultiGridSubsetHandler* mg_sh)
{
	CHECK_POINTER(mg_sh, );
	LOG("  multi-grid-subset-handler erased\n");
	delete mg_sh;
}

/**@}*/ // end of objectCreation


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup loadAndSave Load and Save
 * @{
 */
 
///	Loads marker points from a file.
bool load_markers(MarkerPointManager* mpm, const char* filename)
{
	CHECK_POINTER(mpm, false);
	return LoadMarkerPointsFromFile(*mpm, filename);
}

///	Loads a grid from a file.
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

///	Loads a MultiGrid from a file.
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

///	Saves a grid to a file.
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

///	Saves marked edges to a 'Wavefront .obj' file.
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

///	Adjusts subsets so that they are compatible with ug3.
bool adjust_subsets_for_ug3(Grid* pGrid, SubsetHandler* pSH)
{
	CHECK_POINTER(pGrid, false);
	CHECK_POINTER(pSH, false);
	
	AdjustSubsetsForLgmNg(*pGrid, *pSH);
	
	return true;
}

///	Exports a grid to ug3 (.lgm / .ng).
/**	Make sure that the subsets in pSH have been adjusted for export.
 *	\sa adjust_subsets_for_ug3*/
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

/**@}*/ // end of loadAndSave



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup geometry Geometry
 * @{
 */
///	Converts quadrilaterals to triangles.
bool triangulate(Grid* grid)
{
	CHECK_POINTER(grid, false);
	LOG("  triangulating...");
	Triangulate(*grid, grid->begin<Quadrilateral>(),
				grid->end<Quadrilateral>());
	LOG(" done\n");
	return true;
}

///	Merges vertices which are closer to each other than a given threshold.
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

///	fills a closed triangular surface with tetrahedrons.
/**	Make sure that the given surface-grid does not contain any
 *	holes or degenerated faces and make sure that the surface
 *	does not intersect itself.*/
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

/**@}*/ // end of geometry

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup marks Marks
 * @{
 */
///	Moves each marker-point to its closest vertex.
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

///	Marks edges with a high dihedral angle as crease-edges.
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

///	marks selected edges as crease elements.
void mark_selected_edges(Selector* sel, SubsetHandler* shMarks)
{
	CHECK_POINTER(sel, );
	CHECK_POINTER(shMarks, );
	LOG("  marking selected edges... ");
	shMarks->assign_subset(sel->begin<EdgeBase>(), sel->end<EdgeBase>(),
							RM_CREASE);
							
	MarkFixedCreaseVertices(*sel->get_assigned_grid(), *shMarks,
							RM_CREASE, RM_FIXED);
	
	LOG("done\n");
}

/**@}*/ // end of marks

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup optimization Optimization
 * @{
 */
///	Retriangulates the grid so that all edges approximatly have a certain legth.
/**	Please note that this algorithm only works for triangular grids.*/
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

/**@}*/ // end of optimization

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup neuro Neuro
 * @{
 */
 
///	Extrudes a cylinder around a given point.
/**	The algorithm first moves the given point to the closest
 *	vertex, then adapts the grid around that vertex to the
 *	given cylinder and then extrudes the triangles which lie
 *	inside the cylinder.*/
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

///	extrudes multiple cylinders.
/**	\sa extrude_cylinder.*/
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
/**@}*/ // end of neuro


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup fractures Fractures
 * @{
 */
///	Collects degenerated faces and moves them to the specified subset.
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
///	extrudes a 2d geometry. Suited for fractured geometries.
/**	Extrudes a geometry along the given direction (dirX, dirY, dirZ).
 *	\param numSecs: defines how many extrusion steps are performed.
 *	\sa fracture to subset.*/
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

/**@}*/ // end of fractures

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup subsets Subsets
 * @{
 */
///	Separates volume-regions.
/**	A volume region is a set of volumes, which is surrounded by
 *	a set of faces.*/
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

/**@}*/ // end of subsets

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/** \defgroup selection Selection
 * @{
 */
 

/**@}*/ // end of selection


// only for tests

void test()
{
//	create a new grid and load a file
	Grid grid;
	SubsetHandler sh(grid);
	if(!LoadGridFromFile(grid, "/Users/sreiter/Projects/ug4/trunk/data/grids/unit_square_quads_1x1.obj", sh)){
		UG_LOG("  file-load failed. aborting test\n");
		return;
	}
	
	UG_LOG("  grid loaded... attempting xml-write\n");
	FileAccessorUGX ugx;
	ugx.add_grid(grid, "grid", aPosition);
	ugx.write_to_stream(cout);
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
	using namespace luabind;
	
	module(L)[
	//	bind classes
		class_<Grid>("Grid"),
		class_<MultiGrid>("MultiGrid"),
		
		class_ <Selector>("Selector"),
		
		class_<SubsetHandler>("SubsetHandler"),
		class_<MultiGridSubsetHandler>("MultiGridSubsetHandler"),

		class_<MarkerPointManager>("MarkerPointManager"),

	//	bind methods to lua
		def("new_grid", new_grid),
		def("clone_grid", clone_grid),
		def("copy_grid", copy_grid),
		def("delete_grid", delete_grid),
		def("new_multi_grid", new_multi_grid),
		def("delete_multi_grid", delete_multi_grid),

		def("new_subset_handler", new_subset_handler),
		def("copy_subset_handler", copy_subset_handler),
		def("delete_subset_handler", delete_subset_handler),
		def("new_multi_grid_subset_handler", new_multi_grid_subset_handler),
		def("delete_multi_grid_subset_handler", delete_multi_grid_subset_handler),

		def("load_grid", load_grid),
		def("save_grid", save_grid),
		def("load_multi_grid", load_multi_grid),
		def("adjust_subsets_for_ug3", adjust_subsets_for_ug3),
		def("export_grid_to_ug3", export_grid_to_ug3),

		def("save_marked_edges_to_obj", save_marked_edges_to_obj),
		
		def("remove_doubles", remove_doubles),
		def("convert_to_triangle_grid", triangulate),
		def("triangulate", triangulate),
		def("mark_crease_elements", mark_crease_elements),
		def("adjust_edge_length", adjust_edge_length),

		def("extrude_cylinder", extrude_cylinder),
		def("extrude_cylinders", extrude_cylinders),
		
		def("fracture_to_subset", fracture_to_subset),
		def("frac_extrude", frac_extrude),
		
		def("tetrahedralize", tetrahedralize),

		def("separate_regions", separate_regions),
		
		def("new_markers", new_markers),
		def("delete_markers", delete_markers),
		def("load_markers", load_markers),
		def("snap_markers_to_vertices", snap_markers_to_vertices),
		def("mark_selected_edges", mark_selected_edges),
		
		def("test", test)
	];

	return true;
}

void FinalizeLibGridScript(lua_State* L)
{
}

}}//	end of namespace


