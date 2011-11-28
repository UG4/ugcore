//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include <iostream>
#include <sstream>
#include "registry/registry.h"
#include "bridge.h"
#include "lib_grid/lib_grid.h"
#include "common/profiler/profiler.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/tools/partition_map.h"

using namespace std;

namespace ug
{
namespace bridge
{
///	Wrapper object that simplifies script creation
class GridObject : public Grid
{
	public:
		GridObject() : Grid(GRIDOPT_STANDARD_INTERCONNECTION), m_sh(*this)	{}
		inline Grid& grid()	{return *this;}
		inline SubsetHandler& subset_handler()	{return m_sh;}
		
	protected:
		SubsetHandler m_sh;
};

bool LoadGridObject(GridObject& go, const char* filename)
{
	return LoadGridFromFile(go.grid(), go.subset_handler(), filename);
}

bool SaveGridObject(GridObject& go, const char* filename)
{
	return SaveGridToFile(go.grid(), go.subset_handler(), filename);
}

GridObject* CreateGridObject(const char* filename)
{
	GridObject* go = new GridObject;
	if(!LoadGridObject(*go, filename)){
		delete go;
		return NULL;
	}
	return go;
}

bool CreateFractal(Grid& grid, HangingNodeRefiner_Grid& href,
					number scaleFac, size_t numIterations)
{
	PROFILE_FUNC();
//	HangingNodeRefiner_IR1 href(grid);
	return CreateFractal_NormalScale(grid, href, scaleFac, numIterations);
//	return true;
}

///	Saves a grid hierarchy by offsetting levels along the z-axis.
/**	Note that this method might better be implemented for domains.*/
bool SaveGridHierarchyTransformed(MultiGrid& mg, const SubsetHandler& csh,
								  const char* filename, number offset)
{
//	cast away constness
	SubsetHandler& sh = *const_cast<SubsetHandler*>(&csh);

	APosition aPos;
//	uses auto-attach
	Grid::AttachmentAccessor<VertexBase, APosition> aaPos(mg, aPos, true);

//	copy the existing position to aPos. We take care of dimension differences.
//	Note:	if the method was implemented for domains, this could be implemented
//			in a nicer way.
	if(mg.has_vertex_attachment(aPosition))
		ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition, aPos);
	else if(mg.has_vertex_attachment(aPosition2))
		ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition2, aPos);
	else if(mg.has_vertex_attachment(aPosition1))
		ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition1, aPos);

//	iterate through all vertices and apply an offset depending on their level.
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(VertexBaseIterator iter = mg.begin<VertexBase>(lvl);
			iter != mg.end<VertexBase>(lvl); ++iter)
		{
			aaPos[*iter].z += (number)lvl * offset;
		}
	}

//	finally save the grid
	//SaveGridToFile(mg, sh, filename, aPos);
	bool writeSuccess = SaveGridToUGX(mg, sh, filename, aPos);

//	clean up
	mg.detach_from_vertices(aPos);

	return writeSuccess;
}

///	Saves a grid hierarchy by offsetting levels along the z-axis.
/**	Note that this method might better be implemented for domains.*/
bool SaveGridHierarchyTransformed(MultiGrid& mg, const char* filename, number offset)
{
//	cast away constness
	SubsetHandler& sh = mg.get_hierarchy_handler();

	APosition aPos;
//	uses auto-attach
	Grid::AttachmentAccessor<VertexBase, APosition> aaPos(mg, aPos, true);

//	copy the existing position to aPos. We take care of dimension differences.
//	Note:	if the method was implemented for domains, this could be implemented
//			in a nicer way.
	if(mg.has_vertex_attachment(aPosition))
		ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition, aPos);
	else if(mg.has_vertex_attachment(aPosition2))
		ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition2, aPos);
	else if(mg.has_vertex_attachment(aPosition1))
		ConvertMathVectorAttachmentValues<VertexBase>(mg, aPosition1, aPos);

//	iterate through all vertices and apply an offset depending on their level.
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(VertexBaseIterator iter = mg.begin<VertexBase>(lvl);
			iter != mg.end<VertexBase>(lvl); ++iter)
		{
			aaPos[*iter].z += (number)lvl * offset;
		}
	}

//	finally save the grid
	//SaveGridToFile(mg, sh, filename, aPos);
	bool writeSuccess = SaveGridToUGX(mg, sh, filename, aPos);

//	clean up
	mg.detach_from_vertices(aPos);

	return writeSuccess;
}

bool LoadGrid(Grid& grid, ISubsetHandler& sh, const char* filename)
{
	return LoadGridFromFile(grid, sh, filename);
}

bool LoadGrid(Grid& grid, const char* filename)
{
	return LoadGridFromFile(grid, filename);
}

bool SaveGrid(Grid& grid, const char* filename)
{
	return SaveGridToFile(grid, filename);
}

bool SaveGrid(Grid& grid, SubsetHandler& sh, const char* filename)
{
	return SaveGridToFile(grid, sh, filename);
}

bool SaveGrid(Grid& grid, const SubsetHandler& sh, const char* filename)
{
	return SaveGridToFile(grid, *const_cast<SubsetHandler*>(&sh), filename);
}

bool SaveGridHierarchy(MultiGrid& mg, const char* filename)
{
	return SaveGridToFile(mg, mg.get_hierarchy_handler(), filename);
}


void TestSubdivision(const char* fileIn, const char* fileOut, int numRefs)
{
//todo: Callbacks have to make sure that their attachment is accessible in the grid.
//		even if they were initialized before the attachment was attached to the grid.
	MultiGrid mg;
	SubsetHandler sh(mg);
	RefinementCallbackSubdivisionLoop<APosition> refCallback(mg, aPosition, aPosition);
	GlobalMultiGridRefiner ref(mg, &refCallback);
	
	if(LoadGridFromFile(mg, sh, fileIn)){
		for(int lvl = 0; lvl < numRefs; ++lvl){
			ref.refine();
		}

		ProjectToLimitPLoop(mg, aPosition, aPosition);
		SaveGridToFile(mg, mg.get_hierarchy_handler(), fileOut);

	}
	else{
		UG_LOG("Load failed. aborting...\n");
	}
}

bool CreateSmoothHierarchy(MultiGrid& mg, size_t numRefs)
{
	IRefinementCallback* refCallback = NULL;
//	we're only checking for the main attachments here.
//todo: improve this - add a domain-based hierarchy creator.
	if(mg.has_vertex_attachment(aPosition1))
		refCallback = new RefinementCallbackSubdivisionLoop<APosition1>(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		refCallback = new RefinementCallbackSubdivisionLoop<APosition2>(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		refCallback = new RefinementCallbackSubdivisionLoop<APosition>(mg, aPosition, aPosition);
		
	if(!refCallback){
		UG_LOG("No standard position attachment found. Aborting.\n");
		return false;
	}
	
	GlobalMultiGridRefiner ref(mg, refCallback);

	for(size_t lvl = 0; lvl < numRefs; ++lvl){
		ref.refine();
	}

	if(mg.has_vertex_attachment(aPosition1))
		ProjectToLimitPLoop(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		ProjectToLimitPLoop(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		ProjectToLimitPLoop(mg, aPosition, aPosition);

	delete refCallback;
	return true;
}

bool CreateSemiSmoothHierarchy(MultiGrid& mg, size_t numRefs)
{
	IRefinementCallback* refCallback = NULL;
//	we're only checking for the main attachments here.
//todo: improve this - add a domain-based hierarchy creator.
	if(mg.has_vertex_attachment(aPosition1))
		refCallback = new RefinementCallbackSubdivBoundary<APosition1>(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		refCallback = new RefinementCallbackSubdivBoundary<APosition2>(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		refCallback = new RefinementCallbackSubdivBoundary<APosition>(mg, aPosition, aPosition);
		
	if(!refCallback){
		UG_LOG("No standard position attachment found. Aborting.\n");
		return false;
	}
	
	GlobalMultiGridRefiner ref(mg, refCallback);

	for(size_t lvl = 0; lvl < numRefs; ++lvl){
		ref.refine();
	}

	if(mg.has_vertex_attachment(aPosition1))
		ProjectToLimitSubdivBoundary(mg, aPosition1, aPosition1);
	else if(mg.has_vertex_attachment(aPosition2))
		ProjectToLimitSubdivBoundary(mg, aPosition2, aPosition2);
	else if(mg.has_vertex_attachment(aPosition))
		ProjectToLimitSubdivBoundary(mg, aPosition, aPosition);

	delete refCallback;
	return true;
}

template <class TElem>
void MarkForRefinement(MultiGrid& mg,
					  IRefiner& refiner,
					  float percentage)
{

	typedef typename geometry_traits<TElem>::iterator iterator;
	for(iterator iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter)
	{
		if(urand<float>(0, 99) < percentage){
			refiner.mark(*iter);
		}
	}

/*
	Grid::VertexAttachmentAccessor<APosition> aaPos(mg, aPosition);
	TElem* elem = FindByCoordinate<TElem>(vector3(-0.00001, -0.00001, -0.00001),
									mg.begin<TElem>(mg.num_levels()-1),
									mg.end<TElem>(mg.num_levels()-1),
									aaPos);

	if(elem)
		refiner.mark(elem);
	else{
		UG_LOG("No element found for refinement.\n");
	}*/


}

bool TestHangingNodeRefiner_MultiGrid(const char* filename,
									  const char* outFilename,
									  int numIterations,
									  float percentage)
{
	MultiGrid mg;
	MGSubsetHandler sh(mg);
	HangingNodeRefiner_MultiGrid refiner(mg);

	PROFILE_BEGIN(PROFTEST_loading);
	if(!LoadGridFromFile(mg, sh, filename)){
		UG_LOG("  could not load " << filename << endl);
		return false;
	}
	PROFILE_END();

	{
		PROFILE_BEGIN(PROFTEST_refining);
		for(int i = 0; i < numIterations; ++i){
			UG_LOG("refinement step " << i+1 << endl);

			if(mg.num<Volume>() > 0)
				MarkForRefinement<Volume>(mg, refiner, percentage);
			else if(mg.num<Face>() > 0)
				MarkForRefinement<Face>(mg, refiner, percentage);
			else
				MarkForRefinement<EdgeBase>(mg, refiner, percentage);

			refiner.refine();
		}
	}

	UG_LOG("saving to " << outFilename << endl;)
	SaveGridHierarchy(mg, outFilename);

	UG_LOG("grid element numbers:\n");
	PrintGridElementNumbers(mg);

//	create a surface view
	SurfaceView surfView(mg);

	CreateSurfaceView(surfView, mg, sh);

	SaveGridToFile(mg, surfView, "surface_view.ugx");

	UG_LOG("surface view element numbers:\n");
	PrintGridElementNumbers(surfView);

	return true;
}

////////////////////////////////////////////////////////////////////////
///	This method only makes sense during the development of the grid redistribution
/**	Note that the source grid is completely distributed (no vertical interfaces).
 *
 * The method loads a grid on proc 0 and distributes half of it to proc 1.
 * It then creates a new partition map on each process and redistributes
 * half of the grid on process 0 to process 2 and half of the grid on
 * process 1 to process 3.
 */
void TestGridRedistribution(const char* filename)
{
#ifndef UG_PARALLEL
	UG_LOG("WARNING in TestGridRedistribution: ");
	UG_LOG("This method only works in a parallel environment.\n");
#else
	UG_LOG("Executing TestGridRedistribution...\n");
	MultiGrid mg;
	SubsetHandler sh(mg);
	DistributedGridManager distGridMgr(mg);
	GridLayoutMap& glm = distGridMgr.grid_layout_map();

	int numProcs = pcl::GetNumProcesses();
	if(numProcs != 4){
		UG_LOG(" This test-method only runs with exactly 4 processes.\n");
		return;
	}

UG_LOG(" performing initial distribution\n");
	if(pcl::GetProcRank() == 0){
	//	use this for loading and distribution.
		//MultiGrid tmg;
		//SubsetHandler tsh(tmg);
		SubsetHandler shPart(mg);

		bool success = true;
		string strError;

		if(!LoadGridFromFile(mg, sh, filename)){
			strError = "File not found.";
			success = false;
		}
		else{
		//	partition the grid once (only 2 partitions)
			if(mg.num_volumes() > 0){
				PartitionElementsByRepeatedIntersection<Volume, 3>(
												shPart, mg,
												mg.num_levels() - 1,
												2, aPosition);
			}
			else if(mg.num_faces() > 0){
				PartitionElementsByRepeatedIntersection<Face, 2>(
												shPart, mg,
												mg.num_levels() - 1,
												2, aPosition);
			}
			else{
				strError = "This test-method only runs on geometries which contain "
							"faces or volumes.\n";
				success = false;
			}
		}

	//	make sure that all processes exit if something went wrong.
		if(!pcl::AllProcsTrue(success)){
			UG_LOG(strError << endl);
			return;
		}

	//	now distribute the grid
		//if(!DistributeGrid(tmg, tsh, shPart, 0, &mg, &sh, &glm)){
		if(!DistributeGrid_KeepSrcGrid(mg, sh, glm, shPart, 0, false)){
			UG_LOG("Distribution failed\n");
			return;
		}
	}
	else if(pcl::GetProcRank() == 1){
		if(!pcl::AllProcsTrue(true)){
			UG_LOG("Problems occured on process 0. Aborting.\n");
			return;
		}

		if(!ReceiveGrid(mg, sh, glm, 0, true)){
			UG_LOG("Receive failed.\n");
			return;
		}
	}
	else{
		if(!pcl::AllProcsTrue(true)){
			UG_LOG("Problems occured on process 0. Aborting.\n");
			return;
		}
	}

	distGridMgr.grid_layouts_changed(true);

UG_LOG(" done\n");
UG_LOG(" preparing for redistribution\n");
////////////////////////////////
//	The grid is now distributed to two processes (proc 1 and proc 2).
//	note that it was distributed completely (no source grid kept).

////////////////////////////////
//	Now lets redistribute the grid.
	SubsetHandler shPart(mg);

//	at this point it is clear that we either have volumes or faces.
//	note that we now cut the geometry along the second coordinate (start counting at 0).
	int rank = pcl::GetProcRank();
	if(rank < 2){
		if(mg.num_volumes() > 0){
			PartitionElementsByRepeatedIntersection<Volume, 3>(
											shPart, mg,
											mg.num_levels() - 1,
											2, aPosition, 1);
		}
		else if(mg.num_faces() > 0){
			PartitionElementsByRepeatedIntersection<Face, 2>(
											shPart, mg,
											mg.num_levels() - 1,
											2, aPosition, 1);
		}
		else{
			throw UGError("TestGridRedistribution: Method should have been aborted earlier.");
		}


	//	since currently no proc-map is supported during redistribution, we will now
	//	adjust the subsets so that they match the process-ids.
	//	We do this by hand here. If a real process map would exist this could be automated.
		switch(rank){
			case 0:	shPart.move_subset(1, 2);
					break;

			case 1:	shPart.move_subset(1, 3);
					shPart.move_subset(0, 1);
					break;

			default:
					break;
		}

	//	save the partition maps
		{
			stringstream ss;
			ss << "partition_map_" << pcl::GetProcRank() << ".ugx";
			SaveGridToFile(mg, shPart, ss.str().c_str());
		}
	}

UG_LOG(" done\n");
UG_LOG(" redistributing grid\n");

//	data serialization
	GeomObjAttachmentSerializer<VertexBase, APosition>
		posSerializer(mg, aPosition);
	SubsetHandlerSerializer shSerializer(sh);

	GridDataSerializationHandler serializer;
	serializer.add(&posSerializer);
	serializer.add(&shSerializer);

//	now call redistribution
	RedistributeGrid(distGridMgr, shPart, serializer, serializer, false);
UG_LOG("done\n");

//	save the hierarchy on each process
	{
		stringstream ss;
		ss << "hierarchy_" << pcl::GetProcRank() << ".ugx";
		SaveGridHierarchy(mg, ss.str().c_str());
	}

//	save the domain on each process
	{
		stringstream ss;
		ss << "domain_" << pcl::GetProcRank() << ".ugx";
		SaveGridToFile(mg, sh, ss.str().c_str());
	}

#endif // UG_PARALLEL
}


////////////////////////////////////////////////////////////////////////
///	A helper class for ExpandLayers.
/**	This class should never be publicly available, especially since
 * deriving from std::vector is a bad idea (compare 'Effective C++').
 * However, it is very useful in this situation.
 *
 * The class simply extends std::vector<FractureInfo> by an add_layer method.
 */
class ExpandLayersDesc : public std::vector<FractureInfo>
{
	public:
		ExpandLayersDesc() {}

		void add_layer(int subsetInd, int newSubsetInd, number width)
		{
			push_back(FractureInfo(subsetInd, newSubsetInd, width));
		}
};



////////////////////////////////////////////////////////////////////////
///	test tetrakaidekahedron generator
void TestTKDGenerator(const char* outfile, number height, number baseEdgeLength, number diameter)
{
	Grid g;
	SubsetHandler sh(g);
	sh.set_default_subset_index(0);

	g.attach_to_vertices(aPosition);

	tkdGenerator::GenerateTetrakaidecahedron(g, height, baseEdgeLength, diameter);
	SaveGridToFile(g, sh, outfile);
}

////////////////////////////////////////////////////////////////////////
bool RegisterLibGridInterface(Registry& reg, string parentGroup)
{
	try
	{
	//	get group string
		stringstream groupString; groupString << parentGroup << "/Grid";
		string grp = groupString.str();

	//	Grid
		reg.add_class_<Grid>("Grid", grp)
			.add_constructor()
			.add_method("clear", static_cast<void (Grid::*)()>(&Grid::clear))
			.add_method("num_vertices", &Grid::num_vertices)
			.add_method("num_edges", &Grid::num_edges)
			.add_method("num_faces", &Grid::num_faces)
			.add_method("num_triangles", &Grid::num<Triangle>)
			.add_method("num_quadrilaterals", &Grid::num<Quadrilateral>)
			.add_method("num_volumes", &Grid::num_volumes)
			.add_method("num_tetrahedrons", &Grid::num<Tetrahedron>)
			.add_method("num_pyramids", &Grid::num<Pyramid>)
			.add_method("num_prisms", &Grid::num<Prism>)
			.add_method("num_hexahedrons", &Grid::num<Hexahedron>)
			.add_method("reserve_vertices", &Grid::reserve<VertexBase>)
			.add_method("reserve_edges", &Grid::reserve<EdgeBase>)
			.add_method("reserve_faces", &Grid::reserve<Face>)
			.add_method("reserve_volumes", &Grid::reserve<Volume>);

	//	MultiGrid
		reg.add_class_<MultiGrid, Grid>("MultiGrid", grp)
			.add_constructor()
			.add_method("num_levels", &MultiGrid::num_levels)

			.add_method("num_vertices", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<VertexBase>)
			.add_method("num_edges", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<EdgeBase>)
			.add_method("num_faces", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Face>)
			.add_method("num_triangles", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Triangle>)
			.add_method("num_quadrilaterals", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Quadrilateral>)
			.add_method("num_volumes", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Volume>)
			.add_method("num_tetrahedrons", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Tetrahedron>)
			.add_method("num_pyramids", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Pyramid>)
			.add_method("num_prisms", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Prism>)
			.add_method("num_hexahedrons", (size_t (MultiGrid::*)(int) const) &MultiGrid::num<Hexahedron>);


	////////////////////////
	//	SUBSET HANDLERS

	//  ISubsetHandler
		reg.add_class_<ISubsetHandler>("ISubsetHandler", grp)
			.add_method("num_subsets", &ISubsetHandler::num_subsets)
			.add_method("get_subset_name", &ISubsetHandler::get_subset_name)
			.add_method("set_subset_name", &ISubsetHandler::set_subset_name)
			.add_method("get_subset_index", static_cast<int (ISubsetHandler::*)(const char*) const>(
											&ISubsetHandler::get_subset_index));
		
	//	SubsetHandler
		reg.add_class_<SubsetHandler, ISubsetHandler>("SubsetHandler", grp)
			.add_constructor()
			.add_method("assign_grid", &SubsetHandler::assign_grid);

	//	MGSubsetHandler
		reg.add_class_<MGSubsetHandler, ISubsetHandler>("MGSubsetHandler", grp)
			.add_constructor()
			.add_method("assign_grid", &MGSubsetHandler::assign_grid);

	//	SurfaceView
		reg.add_class_<SurfaceView, SubsetHandler>("SurfaceView", grp)
			.add_constructor()
			.add_method("assign_grid", static_cast<void (SurfaceView::*)(MultiGrid&)>(&SurfaceView::assign_grid));

		reg.add_function("CheckSurfaceView", &CheckSurfaceView, grp);


	////////////////////////
	//	REFINEMENT

	//	IRefiner
		reg.add_class_<IRefiner>("IRefiner", grp)
			.add_method("refine", &IRefiner::refine)
			.add_method("clear_marks", &IRefiner::clear_marks);

	//	HangingNodeRefiner
		reg.add_class_<HangingNodeRefiner_Grid, IRefiner>("HangingNodeRefiner_Grid", grp)
			.add_constructor()
			.add_method("assign_grid", &HangingNodeRefiner_Grid::assign_grid);

		reg.add_class_<HangingNodeRefiner_MultiGrid, IRefiner>("HangingNodeRefiner_MultiGrid", grp)
			.add_constructor()
			.add_method("assign_grid", &HangingNodeRefiner_MultiGrid::assign_grid);

	//	GlobalMultiGridRefiner
		reg.add_class_<GlobalMultiGridRefiner, IRefiner>("GlobalMultiGridRefiner", grp)
			.add_constructor()
			.add_method("assign_grid", static_cast<void (GlobalMultiGridRefiner::*)(MultiGrid&)>(&GlobalMultiGridRefiner::assign_grid));
	
	//	parallel refinement
	#ifdef UG_PARALLEL
		reg.add_class_<ParallelHangingNodeRefiner_MultiGrid, HangingNodeRefiner_MultiGrid>
			("ParallelHangingNodeRefiner_MultiGrid", grp)
			.add_constructor();
	#endif

	//	GridObject
		reg.add_class_<GridObject, Grid>("GridObject", grp)
			.add_constructor()
			.add_method("grid", &GridObject::grid)
			.add_method("subset_handler", &GridObject::subset_handler);

	//	Grid functions
		reg.add_function("CreateFractal", &CreateFractal, grp)
			.add_function("PrintAttachmentInfo", &PrintAttachmentInfo, grp);

	//  GridObject functions
		reg.add_function("LoadGrid", static_cast<bool (*)(Grid&, ISubsetHandler&, const char*)>(&LoadGrid), grp)
			.add_function("LoadGrid", static_cast<bool (*)(Grid&, const char*)>(&LoadGrid), grp)
			.add_function("SaveGrid", static_cast<bool (*)(Grid&, const SubsetHandler&, const char*)>(&SaveGrid), grp)
			.add_function("SaveGrid", static_cast<bool (*)(Grid&, SubsetHandler&, const char*)>(&SaveGrid), grp)
			.add_function("SaveGrid", static_cast<bool (*)(Grid&, const char*)>(&SaveGrid), grp)
			.add_function("LoadGridObject", &LoadGridObject, grp)
			.add_function("SaveGridObject", &SaveGridObject, grp)
			.add_function("SaveGridHierarchyTransformed",
						  static_cast<bool (*)(MultiGrid&, const SubsetHandler&, const char*, number)>(
								  &SaveGridHierarchyTransformed),
						  grp)
			.add_function("SaveGridHierarchyTransformed",
						  static_cast<bool (*)(MultiGrid&, const char*, number)>(
								  &SaveGridHierarchyTransformed),
						  grp)
			.add_function("CreateGridObject", &CreateGridObject, grp)
			.add_function("PrintGridElementNumbers", static_cast<void (*)(MultiGrid&)>(&PrintGridElementNumbers), grp)
			.add_function("PrintGridElementNumbers", static_cast<void (*)(Grid&)>(&PrintGridElementNumbers), grp);

	//	refinement
		reg.add_function("TestSubdivision", &TestSubdivision, grp)
			.add_function("TestHangingNodeRefiner_MultiGrid", &TestHangingNodeRefiner_MultiGrid, grp)
			.add_function("CreateSmoothHierarchy", &CreateSmoothHierarchy, grp)
			.add_function("CreateSemiSmoothHierarchy", &CreateSemiSmoothHierarchy, grp)
			.add_function("SaveGridHierarchy", &SaveGridHierarchy, grp)
			.add_function("TestGridRedistribution", &TestGridRedistribution, grp);
		
	//	subset util
		reg.add_function("AdjustSubsetsForSimulation",
						static_cast<void (*)(SubsetHandler&, bool, bool, bool)>(
						&AdjustSubsetsForSimulation<SubsetHandler>), grp)
			.add_function("AdjustSubsetsForSimulation",
						static_cast<void (*)(MGSubsetHandler&, bool, bool, bool)>(
						&AdjustSubsetsForSimulation<MGSubsetHandler>), grp);

	//	PartitionMap
		reg.add_class_<PartitionMap>("PartitionMap", grp)
			.add_constructor()
			.add_method("clear", &PartitionMap::clear)
			.add_method("get_partition_handler", &PartitionMap::get_partition_handler)
			.add_method("add_target_proc", &PartitionMap::add_target_proc)
			.add_method("add_target_procs", &PartitionMap::add_target_procs)
			.add_method("num_target_procs", &PartitionMap::num_target_procs)
			.add_method("get_target_proc", &PartitionMap::get_target_proc)
			.add_method("shift_target_procs", &PartitionMap::shift_target_procs);

	//	ExpandLayers
		typedef std::vector<FractureInfo> FracInfoVec;
		reg.add_class_<FracInfoVec>("FractureInfoVec", grp);

		reg.add_class_<ExpandLayersDesc, FracInfoVec>("ExpandLayersDesc", grp)
			.add_constructor()
			.add_method("add_layer", &ExpandLayersDesc::add_layer);

		reg.add_function("ExpandLayers2d", &ExpandFractures2d, grp)
			.add_function("ExpandLayers3d", &ExpandFractures3d, grp);

	//	add TKD-Generator method
		reg.add_function("TestTKDGenerator", &TestTKDGenerator, grp);
	}
	catch(UG_REGISTRY_ERROR_RegistrationFailed ex)
	{
		UG_LOG("### ERROR in RegisterLibGridInterface: "
				"Registration failed (using name " << ex.name << ").\n");
		return false;
	}

	return true;
}

}//	end of namespace 
}//	end of namespace 
