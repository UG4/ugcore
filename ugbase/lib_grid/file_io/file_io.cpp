//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d13

#include <string>
#include "file_io.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "common/util/path_provider.h"
#include "common/util/file_util.h"
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/tools/surface_view.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
	#include "pcl/pcl_util.h"
#endif

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
bool FindFileInStandardGridPaths(std::string& filenameOut, const char* filename)
{
	filenameOut = filename;
	if(FileExists(filenameOut.c_str()))
		return true;


//	Now check whether the file was specified relative to the current
//	working directory
	filenameOut = PathProvider::get_current_path();
	filenameOut.append("/").append(filename);

	if(FileExists(filenameOut.c_str()))
		return true;

//	now check the grid path
	filenameOut = PathProvider::get_path(GRID_PATH);
	filenameOut.append("/").append(filename);

	if(FileExists(filenameOut.c_str()))
		return true;

//	filename couldn't be located
	filenameOut = "";
	return false;
}


////////////////////////////////////////////////////////////////////////////////
//	this method performs the actual loading.
static bool LoadGrid3d_IMPL(Grid& grid, ISubsetHandler* pSH,
					   const char* filename, AVector3& aPos)
{
	string strName = filename;
	bool bAutoassignFaces = false;
	bool bSuccess = false;
	if(strName.find(".txt") != string::npos)
	{
		bAutoassignFaces = true;
		bSuccess = LoadGridFromTXT(grid, filename, aPos);
	}
	else if(strName.find(".obj") != string::npos)
		bSuccess = LoadGridFromOBJ(grid, filename, aPos, NULL, pSH);
	else if(strName.find(".lgb") != string::npos)
	{
		int numSHs = 0;
		if(pSH)
			numSHs = 1;

		bSuccess = LoadGridFromLGB(grid, filename, &pSH, numSHs, aPos);
	}
	else if(strName.find(".net") != string::npos)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strName.find(".art") != string::npos)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strName.find(".dat") != string::npos)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strName.find(".lgm") != string::npos)
		bSuccess = ImportGridFromLGM(grid, filename, aPos, pSH);
	else if(strName.find(".ng") != string::npos)
		bSuccess = ImportGridFromNG(grid, filename, aPos, pSH);
	else if(strName.find(".dump") != string::npos)
	{
		bSuccess = LoadGridFromDUMP(grid, filename, pSH, aPos);
	}
	else if(strName.find(".ele") != string::npos)
		return LoadGridFromELE(grid, filename, pSH, aPos);
	else if(strName.find(".msh") != string::npos)
		bSuccess = LoadGridFromMSH(grid, filename, pSH, aPos);

	if(bAutoassignFaces && pSH)
		pSH->assign_subset(grid.faces_begin(), grid.faces_end(), 0);

	return bSuccess;
}


static bool LoadGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition1& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
	if(LoadGrid3d_IMPL(grid, psh, filename, aPosTMP)){
	//	convert the position data from 3d to the required dimension.
		ConvertMathVectorAttachmentValues<VertexBase>(grid, aPosTMP, aPos);
		grid.detach_from_vertices(aPosTMP);
		return true;
	}
	grid.detach_from_vertices(aPosTMP);
	return false;
}

static bool LoadGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition2& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
	if(LoadGrid3d_IMPL(grid, psh, filename, aPosTMP)){
	//	convert the position data from 3d to the required dimension.
		ConvertMathVectorAttachmentValues<VertexBase>(grid, aPosTMP, aPos);
		grid.detach_from_vertices(aPosTMP);
		return true;
	}
	grid.detach_from_vertices(aPosTMP);
	return false;
}

static bool LoadGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition3& aPos)
{
	return LoadGrid3d_IMPL(grid, psh, filename, aPos);
}

////////////////////////////////////////////////////////////////////////////////
///	This method calls specific load routines or delegates loading to LoadGrid3d
template <class TAPos>
static bool LoadGrid(Grid& grid, ISubsetHandler* psh,
					 const char* filename, TAPos& aPos,
					 int procId)
{
//	For convenience, we support multiple different standard paths, from which
//	grids may be loaded. We thus first check, where the specified file is
//	located and load it from that location afterwards.
	bool loadingGrid = true;
	#ifdef UG_PARALLEL
		if((procId != -1) && (pcl::GetProcRank() != procId))
			loadingGrid = false;
	#endif

	grid.message_hub()->post_message(GridMessage_Creation(GMCT_CREATION_STARTS, procId));
	bool retVal = false;
	if(loadingGrid){
	//	Now perform the actual loading.
	//	first all load methods, which do accept template position types are
	//	handled. Then all those which only work with 3d position types are processed.
		string tfile;
		if(FindFileInStandardGridPaths(tfile, filename)){
			if(tfile.find(".ugx") != string::npos){
				if(psh)
					retVal = LoadGridFromUGX(grid, *psh, tfile.c_str(), aPos);
				else{
				//	we have to create a temporary subset handler, since
					SubsetHandler shTmp(grid);
					retVal = LoadGridFromUGX(grid, shTmp, tfile.c_str(), aPos);
				}
			}
			else{
			//	now we'll handle those methods, which only support 3d position types.
				retVal = LoadGrid3d(grid, psh, tfile.c_str(), aPos);
			}
		}
	}

	grid.message_hub()->post_message(GridMessage_Creation(GMCT_CREATION_STOPS, procId));

	#ifdef UG_PARALLEL
		pcl::ProcessCommunicator procCom;
		if(procId == -1)
			retVal = pcl::AllProcsTrue(retVal, procCom);
		else
			retVal = pcl::OneProcTrue(retVal, procCom);
	#endif

	return retVal;
}

////////////////////////////////////////////////////////////////////////////////
template <class TAPos>
bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh,
					  const char* filename, TAPos& aPos, int procId)
{
	return LoadGrid(grid, &sh, filename, aPos, procId);
}

template <class TAPos>
bool LoadGridFromFile(Grid& grid, const char* filename, TAPos& aPos, int procId)
{
	return LoadGrid(grid, NULL, filename, aPos, procId);
}

bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh, const char* filename, int procId)
{
	return LoadGrid(grid, &sh, filename, aPosition, procId);
}

bool LoadGridFromFile(Grid& grid, const char* filename, int procId)
{
	return LoadGrid(grid, NULL, filename, aPosition, procId);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	this method performs the actual save.
static bool SaveGrid3d_IMPL(Grid& grid, ISubsetHandler* pSH,
							const char* filename, AVector3& aPos)
{
	string strName = filename;
	if(strName.find(".txt") != string::npos)
		return SaveGridToTXT(grid, filename, aPos);
	else if(strName.find(".obj") != string::npos)
		return SaveGridToOBJ(grid, filename, aPos, NULL, pSH);
	else if(strName.find(".lgb") != string::npos)
	{
		int numSHs = 0;
		if(pSH)
			numSHs = 1;

		return SaveGridToLGB(grid, filename, &pSH, numSHs, aPos);
	}
	else if(strName.find(".ele") != string::npos)
		return SaveGridToELE(grid, filename, pSH, aPos);
	else if(strName.find(".net") != string::npos)
		return SaveGridToART(grid, filename, pSH, aPos);
	else if(strName.find(".art") != string::npos)
		return SaveGridToART(grid, filename, pSH, aPos);
	else if(strName.find(".ncdf") != string::npos)
		return SaveGridToNCDF(grid, filename, pSH, aPos);
	else if(strName.find(".stl") != string::npos)
		return SaveGridToSTL(grid, filename, pSH, aPos);

	return false;
}

////////////////////////////////////////////////////////////////////////////////
static bool SaveGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition1& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
//	convert the position data from the given dimension to 3d
	ConvertMathVectorAttachmentValues<VertexBase>(grid, aPos, aPosTMP);
	if(SaveGrid3d_IMPL(grid, psh, filename, aPosTMP)){
		grid.detach_from_vertices(aPosTMP);
		return true;
	}
	grid.detach_from_vertices(aPosTMP);
	return false;
}

////////////////////////////////////////////////////////////////////////////////
static bool SaveGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition2& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
//	convert the position data from the given dimension to 3d
	ConvertMathVectorAttachmentValues<VertexBase>(grid, aPos, aPosTMP);
	if(SaveGrid3d_IMPL(grid, psh, filename, aPosTMP)){
		grid.detach_from_vertices(aPosTMP);
		return true;
	}
	grid.detach_from_vertices(aPosTMP);
	return false;
}

////////////////////////////////////////////////////////////////////////////////
static bool SaveGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition3& aPos)
{
	return SaveGrid3d_IMPL(grid, psh, filename, aPos);
}

////////////////////////////////////////////////////////////////////////////////
template <class TAPos>
static bool SaveGrid(Grid& grid, ISubsetHandler* psh,
					 const char* filename, TAPos& aPos)
{
	string strName = filename;
	if(strName.find(".ugx") != string::npos){
		if(psh)
			return SaveGridToUGX(grid, *psh, filename, aPos);
		else {
			SubsetHandler shTmp(grid);
			return SaveGridToUGX(grid, shTmp, filename, aPos);
		}
	}
	else
		return SaveGrid3d(grid, psh, filename, aPos);
}


////////////////////////////////////////////////////////////////////////////////
template <class TAPos>
bool SaveGridToFile(Grid& grid, ISubsetHandler& sh,
					const char* filename, TAPos& aPos)
{
	return SaveGrid(grid, &sh, filename, aPos);
}

template <class TAPos>
bool SaveGridToFile(Grid& grid, const char* filename, TAPos& aPos)
{
	return SaveGrid(grid, NULL, filename, aPos);
}

bool SaveGridToFile(Grid& grid, ISubsetHandler& sh, const char* filename)
{
//	check whether one of the standard attachments is attached and call
//	SaveGrid with that attachment
	if(grid.has_vertex_attachment(aPosition))
		return SaveGrid(grid, &sh, filename, aPosition);
	if(grid.has_vertex_attachment(aPosition2))
		return SaveGrid(grid, &sh, filename, aPosition2);
	if(grid.has_vertex_attachment(aPosition1))
		return SaveGrid(grid, &sh, filename, aPosition1);

	return false;
}

bool SaveGridToFile(Grid& grid, const char* filename)
{
//	check whether one of the standard attachments is attached and call
//	SaveGrid with that attachment
	if(grid.has_vertex_attachment(aPosition))
		return SaveGrid(grid, NULL, filename, aPosition);
	if(grid.has_vertex_attachment(aPosition2))
		return SaveGrid(grid, NULL, filename, aPosition2);
	if(grid.has_vertex_attachment(aPosition1))
		return SaveGrid(grid, NULL, filename, aPosition1);
	return false;
}

bool SaveGridHierarchyTransformed(MultiGrid& mg, ISubsetHandler& sh,
								  const char* filename, number offset)
{
	PROFILE_FUNC_GROUP("grid");
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
			aaPos[*iter].z() += (number)lvl * offset;
		}
	}

//	finally save the grid
	bool writeSuccess = SaveGridToFile(mg, sh, filename, aPos);

//	clean up
	mg.detach_from_vertices(aPos);

	return writeSuccess;
}

bool SaveGridHierarchyTransformed(MultiGrid& mg, const char* filename,
									  number offset)
{
	PROFILE_FUNC_GROUP("grid");
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
			aaPos[*iter].z() += (number)lvl * offset;
		}
	}

//	finally save the grid
	bool writeSuccess = SaveGridToFile(mg, sh, filename, aPos);

//	clean up
	mg.detach_from_vertices(aPos);

	return writeSuccess;
}

template <class TElem>
static void AssignSubsetsByInterfaceType(SubsetHandler& sh, MultiGrid& mg)
{
	const int siNormal = 0;
	const int siHMaster = 1;
	const int siHSlave = 1 << 1;
	const int siVMaster = 1 << 2;
	const int siVSlave = 1 << 3;

	const char* subsetNames[] = {"normal", "hmaster", "hslave", "hslave+hmaster",
						  "vmaster", "vmaster+hmaster", "vmaster+hslave",
						  "vmaster+hslave+hmaster", "vslave", "vslave+hmaster",
						  "vslave+hslave", "vslave+hslave+hmaster",
						  "vslave+vmaster", "vslave+vmaster+hmaster",
						  "vslave+vmaster+hslave", "vslave+vmaster+hmaster+hslave"};

	for(int i = 0; i < 16; ++i)
		sh.subset_info(i).name = subsetNames[i];

	typedef typename Grid::traits<TElem>::iterator TIter;
	for(TIter iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter){
		int status = ES_NONE;

		#ifdef UG_PARALLEL
			DistributedGridManager* distGridMgr = mg.distributed_grid_manager();
			if(distGridMgr)
				status = distGridMgr->get_status(*iter);
		#endif

		int index = siNormal;
		if(status & ES_H_MASTER)
			index |= siHMaster;
		if(status & ES_H_SLAVE)
			index |= siHSlave;
		if(status & ES_V_MASTER)
			index |= siVMaster;
		if(status & ES_V_SLAVE)
			index |= siVSlave;

		sh.assign_subset(*iter, index);
	}
}

bool SaveParallelGridLayout(MultiGrid& mg, const char* filename, number offset)
{
	PROFILE_FUNC_GROUP("grid");

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
			aaPos[*iter].z() += (number)lvl * offset;
		}
	}

//	create a subset handler which holds different subsets for the different interface types
	SubsetHandler sh(mg);

	AssignSubsetsByInterfaceType<VertexBase>(sh, mg);
	AssignSubsetsByInterfaceType<EdgeBase>(sh, mg);
	AssignSubsetsByInterfaceType<Face>(sh, mg);
	AssignSubsetsByInterfaceType<Volume>(sh, mg);

	AssignSubsetColors(sh);
	EraseEmptySubsets(sh);

//	finally save the grid
	bool writeSuccess = SaveGridToFile(mg, sh, filename, aPos);

//	clean up
	mg.detach_from_vertices(aPos);

	return writeSuccess;
}

//template <class TElem>
//static void AssignSubsetsBySurfaceViewState(SubsetHandler& sh, const SurfaceView& sv,
//											MultiGrid& mg)
//{
//	const int siSurface = 0;
//	const int siShadow = 1;
//	const int siShadowing = 2;
//	const int siHidden = 3;
//
//	const char* subsetNames[] = {"surface", "shadow", "shadowing", "not_in_surface"};
//
//	for(int i = 0; i < 4; ++i)
//		sh.subset_info(i).name = subsetNames[i];
//
//	typedef typename Grid::traits<TElem>::iterator TIter;
//	for(TIter iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter){
//		TElem* e = *iter;
//		int si = siHidden;
//
//		if(sv.is_shadowing(e))
//			si = siShadowing;
//		else if(sv.is_shadowed(e))
//			si = siShadow;
//		else if(sv.is_surface_element(e))
//			si = siSurface;
//
//		sh.assign_subset(e, si);
//	}
//}

template <class TElem>
static void AssignSubsetsBySurfaceViewState(SubsetHandler& sh, const SurfaceView& sv,
											MultiGrid& mg)
{
	typedef typename Grid::traits<TElem>::iterator TIter;
	for(TIter iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter){
		TElem* e = *iter;

		sh.assign_subset(e, sv.get_surface_state(e).get());
	}
	for(int i = 0; i < sh.num_subsets(); ++i)
		sh.subset_info(i).name = "unknown";
	sh.subset_info(SurfaceView::UNASSIGNED).name = "unassigned";
	sh.subset_info(SurfaceView::PURE_SURFACE).name = "pure-surface";
	sh.subset_info(SurfaceView::SHADOWING).name = "shadowing";
	sh.subset_info(SurfaceView::SHADOW_COPY).name = "shadow-copy";
	sh.subset_info(SurfaceView::SHADOW_NONCOPY).name = "shadow-noncopy";
	//EraseEmptySubsets(sh);
}

bool SaveSurfaceViewTransformed(MultiGrid& mg, const SurfaceView& sv,
								const char* filename, number offset)
{
	PROFILE_FUNC_GROUP("grid");

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
			aaPos[*iter].z() += (number)lvl * offset;
		}
	}

//	create a subset handler which holds different subsets for the different interface types
	SubsetHandler sh(mg);

	AssignSubsetsBySurfaceViewState<VertexBase>(sh, sv, mg);
	AssignSubsetsBySurfaceViewState<EdgeBase>(sh, sv, mg);
	AssignSubsetsBySurfaceViewState<Face>(sh, sv, mg);
	AssignSubsetsBySurfaceViewState<Volume>(sh, sv, mg);

	AssignSubsetColors(sh);
	EraseEmptySubsets(sh);

//	finally save the grid
	bool writeSuccess = SaveGridToFile(mg, sh, filename, aPos);

//	clean up
	mg.detach_from_vertices(aPos);

	return writeSuccess;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiation
template bool LoadGridFromFile(Grid&, ISubsetHandler&, const char*, AVector1&, int);
template bool LoadGridFromFile(Grid&, ISubsetHandler&, const char*, AVector2&, int);
template bool LoadGridFromFile(Grid&, ISubsetHandler&, const char*, AVector3&, int);

template bool LoadGridFromFile(Grid&, const char*, AVector1&, int);
template bool LoadGridFromFile(Grid&, const char*, AVector2&, int);
template bool LoadGridFromFile(Grid&, const char*, AVector3&, int);

template bool SaveGridToFile(Grid&, ISubsetHandler&, const char*, AVector1&);
template bool SaveGridToFile(Grid&, ISubsetHandler&, const char*, AVector2&);
template bool SaveGridToFile(Grid&, ISubsetHandler&, const char*, AVector3&);

template bool SaveGridToFile(Grid&, const char*, AVector1&);
template bool SaveGridToFile(Grid&, const char*, AVector2&);
template bool SaveGridToFile(Grid&, const char*, AVector3&);

}//	end of namespace
