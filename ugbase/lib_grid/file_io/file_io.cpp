/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <string>
#include "file_io.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "common/util/path_provider.h"
#include "common/util/file_util.h"
#include "common/util/string_util.h"
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/tools/surface_view.h"

#include "file_io_2df.h"
#include "file_io_art.h"
#include "file_io_asc.h"
#include "file_io_txt.h"
#include "file_io_tetgen.h"
#include "file_io_obj.h"
#include "file_io_lgm.h"
#include "file_io_lgb.h"
#include "file_io_ng.h"
#include "file_io_ug.h"
#include "file_io_dump.h"
#include "file_io_ncdf.h"
#include "file_io_ugx.h"
#include "file_io_msh.h"
#include "file_io_stl.h"
#include "file_io_tikz.h"
#include "file_io_vtu.h"
#include "file_io_swc.h"
#include "file_io_grdecl.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl_process_communicator.h"
	#include "pcl/pcl_util.h"
#endif

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
//	this method performs the actual loading.
static bool LoadGrid3d_IMPL(Grid& grid, ISubsetHandler* pSH,
					   const char* filename, AVector3& aPos)
{
	string strExt = GetFilenameExtension(string(filename));
	strExt = ToLower(strExt);

	bool bAutoassignFaces = false;
	bool bSuccess = false;
	if(strExt.compare("txt") == 0)
	{
		bAutoassignFaces = true;
		bSuccess = LoadGridFromTXT(grid, filename, aPos);
	}
	else if(strExt.compare("grdecl") == 0)
	{
		bAutoassignFaces = true;
		bSuccess = LoadGridFromGRDECL(grid, filename, aPos);
	}
	else if(strExt.compare("obj") == 0)
		bSuccess = LoadGridFromOBJ(grid, filename, aPos, NULL, pSH);
	else if(strExt.compare("lgb") == 0)
	{
		int numSHs = 0;
		if(pSH)
			numSHs = 1;

		bSuccess = LoadGridFromLGB(grid, filename, &pSH, numSHs, NULL, aPos);
	}
	else if(strExt.compare("2df") == 0)
		bSuccess = LoadGridFrom2DF(grid, filename, pSH, aPos);
	else if(strExt.compare("stl") == 0)
		bSuccess = LoadGridFromSTL(grid, filename, pSH, aPos);
	else if(strExt.compare("net") == 0)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strExt.compare("art") == 0)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strExt.compare("dat") == 0)
		bSuccess = LoadGridFromART(grid, filename, pSH, aPos);
	else if(strExt.compare("lgm") == 0)
		bSuccess = ImportGridFromLGM(grid, filename, aPos, pSH);
	else if(strExt.compare("ng") == 0)
		bSuccess = ImportGridFromNG(grid, filename, aPos, pSH);
	else if(strExt.compare("dump") == 0)
	{
		bSuccess = LoadGridFromDUMP(grid, filename, pSH, aPos);
	}
	else if(strExt.compare("ele") == 0)
		return LoadGridFromELE(grid, filename, pSH, aPos);
	else if(strExt.compare("msh") == 0)
		bSuccess = LoadGridFromMSH(grid, filename, pSH, aPos);
	else if(strExt.compare("smesh") == 0)
		bSuccess = LoadGridFromSMESH(grid, filename, aPos, pSH);
	else if(strExt.compare("asc") == 0){
		bSuccess = LoadGridFromASC(grid, filename, aPos);
		bAutoassignFaces = true;
	}
	else if(strExt.compare("swc") == 0){
		bSuccess = LoadGridFromSWC(grid, pSH, filename, aPos);
	}

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
		ConvertMathVectorAttachmentValues<Vertex>(grid, aPosTMP, aPos);
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
		ConvertMathVectorAttachmentValues<Vertex>(grid, aPosTMP, aPos);
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
		if((procId != -1) && (pcl::ProcRank() != procId))
			loadingGrid = false;
	#endif

	grid.message_hub()->post_message(GridMessage_Creation(GMCT_CREATION_STARTS, procId));
	bool retVal = false;
	if(loadingGrid){
	//	Now perform the actual loading.
	//	first all load methods, which do accept template position types are
	//	handled. Then all those which only work with 3d position types are processed.
		string tfile = FindFileInStandardPaths(filename);
		if(!tfile.empty()){
			if(tfile.find(".ugx") != string::npos){
				if(psh)
					retVal = LoadGridFromUGX(grid, *psh, tfile.c_str(), aPos);
				else{
				//	we have to create a temporary subset handler
					SubsetHandler shTmp(grid);
					retVal = LoadGridFromUGX(grid, shTmp, tfile.c_str(), aPos);
				}
			}
			else if(tfile.find(".vtu") != string::npos){
				if(psh)
					retVal = LoadGridFromVTU(grid, *psh, tfile.c_str(), aPos);
				else{
				//	we have to create a temporary subset handler
					SubsetHandler shTmp(grid);
					retVal = LoadGridFromVTU(grid, shTmp, tfile.c_str(), aPos);
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

template <class TAPos>
static bool LoadGrid(Grid& grid, SPProjectionHandler* ph, size_t& num_ph, ISubsetHandler* psh, std::vector<std::string> additionalSHNames,
						std::vector<SmartPtr<ISubsetHandler>> ash, const char* filename, TAPos& aPos, int procId)
{
//	For convenience, we support multiple different standard paths, from which
//	grids may be loaded. We thus first check, where the specified file is
//	located and load it from that location afterwards.
	bool loadingGrid = true;
	#ifdef UG_PARALLEL
		if((procId != -1) && (pcl::ProcRank() != procId))
			loadingGrid = false;
	#endif

	grid.message_hub()->post_message(GridMessage_Creation(GMCT_CREATION_STARTS, procId));

	bool retVal = false;
	if(loadingGrid){
	//	Now perform the actual loading.
	//	first all load methods, which do accept template position types are
	//	handled. Then all those which only work with 3d position types are processed.
		string tfile = FindFileInStandardPaths(filename);
		if(!tfile.empty()){
			if(tfile.find(".ugx") != string::npos){
				if(psh)
					retVal = LoadGridFromUGX(grid, *ph, num_ph, *psh, additionalSHNames, ash, tfile.c_str(), aPos);
				else{
				//	we have to create a temporary subset handler
					SubsetHandler shTmp(grid);
					retVal = LoadGridFromUGX(grid, *ph, num_ph, shTmp, additionalSHNames, ash, tfile.c_str(), aPos);
				}
			}

			else if(tfile.find(".vtu") != string::npos){
				if(psh)
					retVal = LoadGridFromVTU(grid, *psh, tfile.c_str(), aPos);
				else{
				//	we have to create a temporary subset handler
					SubsetHandler shTmp(grid);
					retVal = LoadGridFromVTU(grid, shTmp, tfile.c_str(), aPos);
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
bool LoadGridFromFile(Grid& grid, SPProjectionHandler& ph, size_t& num_ph, ISubsetHandler& sh, vector<string> additionalSHNames,
						vector<SmartPtr<ISubsetHandler>> ash, const char* filename, TAPos& aPos, int procId)
{
	return LoadGrid(grid, &ph, num_ph, &sh, additionalSHNames, ash, filename, aPos, procId);
}

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
	if(strName.find(".2df") != string::npos)
		return SaveGridTo2DF(grid, filename, pSH, aPos);
	else if(strName.find(".obj") != string::npos)
		return SaveGridToOBJ(grid, filename, aPos, NULL, pSH);
	else if(strName.find(".lgb") != string::npos)
	{
		int numSHs = 0;
		if(pSH)
			numSHs = 1;

		return SaveGridToLGB(grid, filename, &pSH, numSHs, NULL, aPos);
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
	else if(strName.find(".smesh") != string::npos)
		return ExportGridToSMESH(grid, filename, aPos, pSH);
	else if((strName.find(".tikz") != string::npos)
			|| (strName.find(".tex") != string::npos))
	{
		return ExportGridToTIKZ(grid, filename, pSH, aPos);
	}
	else if (strName.find(".swc") != string::npos)
		return ExportGridToSWC(grid, pSH, filename, aPos);

	return false;
}

////////////////////////////////////////////////////////////////////////////////
static bool SaveGrid3d(Grid& grid, ISubsetHandler* psh,
					   const char* filename, APosition1& aPos)
{
	APosition aPosTMP;
	grid.attach_to_vertices(aPosTMP);
//	convert the position data from the given dimension to 3d
	ConvertMathVectorAttachmentValues<Vertex>(grid, aPos, aPosTMP);
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
	ConvertMathVectorAttachmentValues<Vertex>(grid, aPos, aPosTMP);
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
		 #if (defined UG_PARALLEL && defined UG_DEBUG)
		 std::size_t found=strName.find(".ugx");
		 strName=strName.replace(found, 4, "");
		 strName=strName.append(std::to_string(pcl::ProcRank()));
		 strName.append(".ugx");   
		#endif
		
		if(psh)
			return SaveGridToUGX(grid, *psh, strName.c_str(), aPos);
		else {
			SubsetHandler shTmp(grid);
			return SaveGridToUGX(grid, shTmp, strName.c_str(), aPos);
		}
	}
	else if(strName.find(".vtu") != string::npos){
		#if (defined UG_PARALLEL && defined UG_DEBUG)
                 std::size_t found=strName.find(".vtu");
	         	 strName=strName.replace(found, 4, "");
                 strName=strName.append(std::to_string(pcl::ProcRank()));
                 strName.append(".vtu");
		#endif
		return SaveGridToVTU(grid, psh, strName.c_str(), aPos);
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
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(mg, aPos, true);

//	copy the existing position to aPos. We take care of dimension differences.
//	Note:	if the method was implemented for domains, this could be implemented
//			in a nicer way.
	if(mg.has_vertex_attachment(aPosition))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition, aPos);
	else if(mg.has_vertex_attachment(aPosition2))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition2, aPos);
	else if(mg.has_vertex_attachment(aPosition1))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition1, aPos);

//	iterate through all vertices and apply an offset depending on their level.
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(VertexIterator iter = mg.begin<Vertex>(lvl);
			iter != mg.end<Vertex>(lvl); ++iter)
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
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(mg, aPos, true);

//	copy the existing position to aPos. We take care of dimension differences.
//	Note:	if the method was implemented for domains, this could be implemented
//			in a nicer way.
	if(mg.has_vertex_attachment(aPosition))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition, aPos);
	else if(mg.has_vertex_attachment(aPosition2))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition2, aPos);
	else if(mg.has_vertex_attachment(aPosition1))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition1, aPos);

//	iterate through all vertices and apply an offset depending on their level.
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(VertexIterator iter = mg.begin<Vertex>(lvl);
			iter != mg.end<Vertex>(lvl); ++iter)
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
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(mg, aPos, true);

//	copy the existing position to aPos. We take care of dimension differences.
//	Note:	if the method was implemented for domains, this could be implemented
//			in a nicer way.
	if(mg.has_vertex_attachment(aPosition))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition, aPos);
	else if(mg.has_vertex_attachment(aPosition2))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition2, aPos);
	else if(mg.has_vertex_attachment(aPosition1))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition1, aPos);

//	iterate through all vertices and apply an offset depending on their level.
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(VertexIterator iter = mg.begin<Vertex>(lvl);
			iter != mg.end<Vertex>(lvl); ++iter)
		{
			aaPos[*iter].z() += (number)lvl * offset;
		}
	}

//	create a subset handler which holds different subsets for the different interface types
	SubsetHandler sh(mg);

	AssignSubsetsByInterfaceType<Vertex>(sh, mg);
	AssignSubsetsByInterfaceType<Edge>(sh, mg);
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

template <class TElem>
static void AssignSubsetsBySurfaceViewState(SubsetHandler& sh, const SurfaceView& sv,
											MultiGrid& mg)
{
	typedef typename Grid::traits<TElem>::iterator TIter;
	for(TIter iter = mg.begin<TElem>(); iter != mg.end<TElem>(); ++iter){
		TElem* e = *iter;

		sh.assign_subset(e, sv.surface_state(e).get());
	}
	for(int i = 0; i < sh.num_subsets(); ++i)
		sh.subset_info(i).name = "unknown";
	sh.subset_info(SurfaceView::MG_SHADOW_PURE).name = "shadow-pure";
	sh.subset_info(SurfaceView::MG_SURFACE_PURE).name = "surface-pure";
	sh.subset_info(SurfaceView::MG_SURFACE_RIM).name = "surface-rim";
	sh.subset_info(SurfaceView::MG_SHADOW_RIM_COPY).name = "shadow-rim-copy";
	sh.subset_info(SurfaceView::MG_SHADOW_RIM_NONCOPY).name = "shadow-rim-noncopy";
	//EraseEmptySubsets(sh);
}

bool SaveSurfaceViewTransformed(MultiGrid& mg, const SurfaceView& sv,
								const char* filename, number offset)
{
	PROFILE_FUNC_GROUP("grid");

	APosition aPos;
//	uses auto-attach
	Grid::AttachmentAccessor<Vertex, APosition> aaPos(mg, aPos, true);

//	copy the existing position to aPos. We take care of dimension differences.
//	Note:	if the method was implemented for domains, this could be implemented
//			in a nicer way.
	if(mg.has_vertex_attachment(aPosition))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition, aPos);
	else if(mg.has_vertex_attachment(aPosition2))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition2, aPos);
	else if(mg.has_vertex_attachment(aPosition1))
		ConvertMathVectorAttachmentValues<Vertex>(mg, aPosition1, aPos);

//	iterate through all vertices and apply an offset depending on their level.
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		for(VertexIterator iter = mg.begin<Vertex>(lvl);
			iter != mg.end<Vertex>(lvl); ++iter)
		{
			aaPos[*iter].z() += (number)lvl * offset;
		}
	}

//	create a subset handler which holds different subsets for the different interface types
	SubsetHandler sh(mg);

	AssignSubsetsBySurfaceViewState<Vertex>(sh, sv, mg);
	AssignSubsetsBySurfaceViewState<Edge>(sh, sv, mg);
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

template<class TElem>
void CopyGridLevelElements(MultiGrid& srcMG, Grid& destGrid,
				           ISubsetHandler& srcSH, ISubsetHandler& destSH,
						   int lvl, AVertex& aNewVrt)
{
	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcMG, aNewVrt);
	GridObjectCollection goc = srcMG.get_grid_objects();
	CustomVertexGroup vrts;

	typedef typename Grid::traits<TElem>::iterator iter_t;

	for(iter_t eIter = goc.begin<TElem>(lvl); eIter != goc.end<TElem>(lvl); ++eIter)
	{
		TElem* e = *eIter;
		vrts.resize(e->num_vertices());

		for(size_t iv = 0; iv < e->num_vertices(); ++iv)
		{
			vrts.set_vertex(iv, aaNewVrt[e->vertex(iv)]);
		}

		TElem* ne = *destGrid.create_by_cloning(e, vrts);
		destSH.assign_subset(ne, srcSH.get_subset_index(e));
	}
}

template <class TAPos>
void CopyGridLevel(MultiGrid& srcMG, Grid& destGrid,
				   ISubsetHandler& srcSH, ISubsetHandler& destSH,
				   int lvl, TAPos aPos)
{
	Grid::VertexAttachmentAccessor<TAPos> aaPos(destGrid, aPos);
	Grid::VertexAttachmentAccessor<TAPos> aaSrcPos(srcMG, aPos);
	GridObjectCollection goc = srcMG.get_grid_objects();

	AVertex aNewVrt;
	srcMG.attach_to_vertices(aNewVrt);
	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcMG, aNewVrt);

	for(int si = destSH.num_subsets(); si < srcSH.num_subsets(); ++si)
	{
		destSH.subset_info(si) = srcSH.subset_info(si);
	}

	for(VertexIterator vrtIter = goc.begin<Vertex>(lvl); vrtIter != goc.end<Vertex>(lvl); ++vrtIter)
	{
		Vertex* srcVrt  = *vrtIter;
		Vertex* destVrt = *destGrid.create_by_cloning(srcVrt);

		aaNewVrt[srcVrt] = destVrt;
		aaPos[destVrt] = aaSrcPos[srcVrt];
		destSH.assign_subset(destVrt, srcSH.get_subset_index(srcVrt));
	}

	CopyGridLevelElements<Edge>(srcMG, destGrid, srcSH, destSH, lvl, aNewVrt);
	CopyGridLevelElements<Face>(srcMG, destGrid, srcSH, destSH, lvl, aNewVrt);
	CopyGridLevelElements<Volume>(srcMG, destGrid, srcSH, destSH, lvl, aNewVrt);

	srcMG.detach_from_vertices(aNewVrt);
}

template <class TAPos>
void CopyGrid(Grid& srcGrid, Grid& destGrid,
			  ISubsetHandler& srcSH, ISubsetHandler& destSH,
			  TAPos aPos)
{
	Grid::VertexAttachmentAccessor<TAPos> aaPos(destGrid, aPos);
	Grid::VertexAttachmentAccessor<TAPos> aaSrcPos(srcGrid, aPos);
	GridObjectCollection goc = srcGrid.get_grid_objects();

	AVertex aNewVrt;
	srcGrid.attach_to_vertices(aNewVrt);
	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	for(int si = destSH.num_subsets(); si < srcSH.num_subsets(); ++si)
	{
		destSH.subset_info(si) = srcSH.subset_info(si);
	}

	for(VertexIterator vrtIter = goc.begin<Vertex>(); vrtIter != goc.end<Vertex>(); ++vrtIter)
	{
		Vertex* srcVrt  = *vrtIter;
		Vertex* destVrt = *destGrid.create_by_cloning(srcVrt);
		aaNewVrt[srcVrt] = destVrt;
		aaPos[destVrt] = aaSrcPos[srcVrt];
		destSH.assign_subset(destVrt, srcSH.get_subset_index(srcVrt));
	}

	CopyGridElements<Edge>(srcGrid, destGrid, srcSH, destSH, aNewVrt);
	CopyGridElements<Face>(srcGrid, destGrid, srcSH, destSH,  aNewVrt);
	CopyGridElements<Volume>(srcGrid, destGrid, srcSH, destSH, aNewVrt);

	srcGrid.detach_from_vertices(aNewVrt);
}

template <class TAPos>
bool SaveGridLevel(MultiGrid& srcMG, ISubsetHandler& srcSH,
				   int lvl, const char* filename, TAPos aPos)
{
	Grid destGrid;
	SubsetHandler destSH(destGrid);

	destGrid.attach_to_vertices(aPos);

	CopyGridLevel(srcMG, destGrid, srcSH, destSH, lvl, aPos);
	SaveGridToFile(destGrid, destSH, filename);

	return true;
}

template <typename TAPos>
void MergeGrids
(
	Grid& mrgGrid,
	Grid& grid,
	ISubsetHandler& mrgSH,
	ISubsetHandler& sh,
	TAPos aPos,
	bool joinSubsets
)
{
	// add offset to not join subsets with same index
	int subsetBaseInd = joinSubsets ? 0 : mrgSH.num_subsets();

	// attach data
	AVertex aVrt;
	grid.attach_to_vertices(aVrt);

	// attachments accessors for position and vertex index
	Grid::AttachmentAccessor<Vertex, TAPos> aaPosMRG(mrgGrid, aPos, true);
	Grid::AttachmentAccessor<Vertex, TAPos> aaPos(grid, aPos, true);
	Grid::AttachmentAccessor<Vertex, AVertex> aaVrt(grid, aVrt);

	// copy vertices
	for (VertexIterator iter = grid.begin<Vertex>();
			iter != grid.end<Vertex>(); ++iter)
	{
		Vertex* nvrt = *mrgGrid.create_by_cloning(*iter);
		aaPosMRG[nvrt] = aaPos[*iter];
		aaVrt[*iter] = nvrt;
		mrgSH.assign_subset(nvrt, subsetBaseInd + sh.get_subset_index(*iter));
	}

	//	copy edges
	EdgeDescriptor ed;
	for (EdgeIterator iter = grid.begin<Edge>();
			iter != grid.end<Edge>(); ++iter)
	{
		Edge* eSrc = *iter;
		ed.set_vertices(aaVrt[eSrc->vertex(0)], aaVrt[eSrc->vertex(1)]);
		Edge* e = *mrgGrid.create_by_cloning(eSrc, ed);
		mrgSH.assign_subset(e, subsetBaseInd + sh.get_subset_index(eSrc));
	}

	//	copy faces
	FaceDescriptor fd;
	for (FaceIterator iter = grid.begin<Face>();
		iter != grid.end<Face>(); ++iter)
	{
		Face* fSrc = *iter;
		fd.set_num_vertices((uint)fSrc->num_vertices());
		for (size_t i = 0; i < fd.num_vertices(); ++i) {
			fd.set_vertex((uint)i, aaVrt[fSrc->vertex(i)]);
		}
		Face* f = *mrgGrid.create_by_cloning(fSrc, fd);
		mrgSH.assign_subset(f, subsetBaseInd + sh.get_subset_index(fSrc));
	}

	// copy volumes
	VolumeDescriptor vd;
	for (VolumeIterator iter = grid.begin<Volume>();
			iter != grid.end<Volume>(); ++iter)
	{
		Volume* vSrc = *iter;
		vd.set_num_vertices((uint)vSrc->num_vertices());
		for (size_t i = 0; i < vd.num_vertices(); ++i) {
			vd.set_vertex((uint)i, aaVrt[vSrc->vertex(i)]);
		}

		Volume* v = *mrgGrid.create_by_cloning(vSrc, vd);
		mrgSH.assign_subset(v, subsetBaseInd + sh.get_subset_index(vSrc));
	}

	// remove the temporary attachment
	mrgGrid.detach_from_vertices(aVrt);

	// overwrite subset names
	for (int i_sub = 0; i_sub < sh.num_subsets(); ++i_sub){
		mrgSH.subset_info(subsetBaseInd + i_sub) = sh.subset_info(i_sub);
	}
}

bool SaveGridLevelToFile(MultiGrid& srcMG, ISubsetHandler& srcSH, int lvl, const char* filename)
{
//	check whether one of the standard attachments is attached and call
//	SaveGridLevel with that attachment
	/*#ifdef UG_PARALLEL
	    std::size_t found=filename.find(".ugx")
            if(found != string::npos){
		filename=filename.replace(found, 4, "");    
	    }
	    filename=filename.append(std::to_string(pcl::ProRank()));
	    if(found != string::npos){
	    		filename.append(".ugx");
	    }
	#endif*/
	if(srcMG.has_vertex_attachment(aPosition))
		return SaveGridLevel(srcMG, srcSH, lvl, filename, aPosition);
	if(srcMG.has_vertex_attachment(aPosition2))
		return SaveGridLevel(srcMG, srcSH, lvl, filename, aPosition2);
	if(srcMG.has_vertex_attachment(aPosition1))
		return SaveGridLevel(srcMG, srcSH, lvl, filename, aPosition1);

	return false;
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

template bool LoadGridFromFile(Grid&, SPProjectionHandler&, size_t&, ISubsetHandler&, std::vector<std::string>, std::vector<SmartPtr<ISubsetHandler>>, const char*, AVector1&, int);
template bool LoadGridFromFile(Grid&, SPProjectionHandler&, size_t&, ISubsetHandler&, std::vector<std::string>, std::vector<SmartPtr<ISubsetHandler>>, const char*, AVector2&, int);
template bool LoadGridFromFile(Grid&, SPProjectionHandler&, size_t&, ISubsetHandler&, std::vector<std::string>, std::vector<SmartPtr<ISubsetHandler>>, const char*, AVector3&, int);

template bool SaveGridToFile(Grid&, ISubsetHandler&, const char*, AVector1&);
template bool SaveGridToFile(Grid&, ISubsetHandler&, const char*, AVector2&);
template bool SaveGridToFile(Grid&, ISubsetHandler&, const char*, AVector3&);

template bool SaveGridToFile(Grid&, const char*, AVector1&);
template bool SaveGridToFile(Grid&, const char*, AVector2&);
template bool SaveGridToFile(Grid&, const char*, AVector3&);

template void CopyGrid(Grid&, Grid&, ISubsetHandler&, ISubsetHandler&, APosition1);
template void CopyGrid(Grid&, Grid&, ISubsetHandler&, ISubsetHandler&, APosition2);
template void CopyGrid(Grid&, Grid&, ISubsetHandler&, ISubsetHandler&, APosition3);

template void MergeGrids(Grid&, Grid&, ISubsetHandler&, ISubsetHandler&, APosition1, bool);
template void MergeGrids(Grid&, Grid&, ISubsetHandler&, ISubsetHandler&, APosition2, bool);
template void MergeGrids(Grid&, Grid&, ISubsetHandler&, ISubsetHandler&, APosition3, bool);

}//	end of namespace
