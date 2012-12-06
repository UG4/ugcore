//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d30

#include <sstream>
#include "parallelization_util.h"
#include "load_balancing.h"
#include "grid_distribution.h"
#include "lib_grid/algorithms/refinement/global_multi_grid_refiner.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/algorithms/subset_util.h"
#include "pcl/pcl_layout_tests.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
int GetAssociatedInterfaceType(int interfaceType)
{
	switch(interfaceType){
		case INT_H_MASTER:	return INT_H_SLAVE;
		case INT_H_SLAVE:	return INT_H_MASTER;
		case INT_V_MASTER:	return INT_V_SLAVE;
		case INT_V_SLAVE:	return INT_V_MASTER;
		default: return INT_NONE;
	}
}

////////////////////////////////////////////////////////////////////////
bool PartitionGrid_Bisection(SubsetHandler& partitionOut,
							  MultiGrid& mg, ISubsetHandler& sh,
							  size_t numProcs)
{
	if(mg.num<Volume>() > 0)
		PartitionElementsByRepeatedIntersection<Volume, 3>(
									partitionOut, mg,
									mg.num_levels() - 1,
									numProcs, aPosition);
	else if(mg.num<Face>() > 0)
		PartitionElementsByRepeatedIntersection<Face, 2>(
											partitionOut, mg,
											mg.num_levels() - 1,
											numProcs, aPosition);
	else if(mg.num<EdgeBase>() > 0)
		PartitionElementsByRepeatedIntersection<EdgeBase, 1>(
											partitionOut, mg,
											mg.num_levels() - 1,
											numProcs, aPosition);
	else{
		LOG("partitioning could not be performed - "
			<< "grid neither containes edges nor faces nor volumes. aborting...\n");
		return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
AdjustGrid_GlobalRefinement::
AdjustGrid_GlobalRefinement(size_t numRefinements) :
	m_numRefinements(numRefinements)
{
}

void AdjustGrid_GlobalRefinement::
operator()(MultiGrid& mg, ISubsetHandler& sh)
{
	if(m_numRefinements > 0){
		GlobalMultiGridRefiner refiner(mg);
		for(size_t i = 0; i < m_numRefinements; ++i)
			refiner.refine();
	}
}

////////////////////////////////////////////////////////////////////////
AdjustGrid_AutoAssignSubsets::
AdjustGrid_AutoAssignSubsets(int innerSubsetIndex,
							 int boundarySubsetIndex) :
	m_innerSubsetIndex(innerSubsetIndex),
	m_boundarySubsetIndex(boundarySubsetIndex)
{
}

void AdjustGrid_AutoAssignSubsets::
operator()(MultiGrid& mg, ISubsetHandler& sh)
{
	if(m_innerSubsetIndex > -2 && m_boundarySubsetIndex > -2)
		AssignInnerAndBoundarySubsets(mg, sh, m_innerSubsetIndex, m_boundarySubsetIndex);
}


////////////////////////////////////////////////////////////////////////
AdjustGrid_AutoAssignSubsetsAndRefine::
AdjustGrid_AutoAssignSubsetsAndRefine(int innerSubsetIndex,
									  int outerSubsetIndex,
									  size_t numRefinements) :
	m_assignSubsets(innerSubsetIndex, outerSubsetIndex),
	m_globalRefinement(numRefinements)
{
}

void AdjustGrid_AutoAssignSubsetsAndRefine::
operator()(MultiGrid& mg, ISubsetHandler& sh)
{
	m_assignSubsets(mg, sh);
	m_globalRefinement(mg, sh);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
bool LoadAndDistributeGrid(DistributedGridManager& distGridMgrOut,
						   ISubsetHandler& shOut, int numProcs,
						   const char* filename,
						   bool keepSrcGrid,
						   FuncAdjustGrid funcAdjustGrid,
						   FuncPartitionGrid funcPartitionGrid)
{

	UG_ASSERT(dynamic_cast<MultiGrid*>(distGridMgrOut.get_assigned_grid()),
				"Error during LoadAndDistributeGrid: DistributedGridManager has to operate "
				"on a MultiGrid in the current implementation.");

	MultiGrid& mgOut = *dynamic_cast<MultiGrid*>(distGridMgrOut.get_assigned_grid());
	GridLayoutMap& glmOut = distGridMgrOut.grid_layout_map();

//	make sure that each grid has a position attachment - even if no data
//	will be received.
	if(!mgOut.has_vertex_attachment(aPosition))
		mgOut.attach_to_vertices(aPosition);

//	process 0 loads and distributes the grid. The others receive it.
	if(pcl::GetProcRank() == 0)
	{
	//	load and distribute the grid
		{
			MultiGrid tmpMG(GRIDOPT_FULL_INTERCONNECTION);
			SubsetHandler tmpSH(tmpMG);

		//	depending on whether the source-grid should be kept, we either
		//	operate on mgOut directly or else load the grid into a temporary grid.
			MultiGrid* pmg = &tmpMG;
			ISubsetHandler* psh = &tmpSH;
			if(keepSrcGrid){
				pmg = &mgOut;
				psh = &shOut;
			}
			MultiGrid& mg = *pmg;
			ISubsetHandler& sh = *psh;

		//	load
			LOG("  loading... ");
			if(!LoadGridFromFile(mg, sh, filename))
			{
				LOG(" failed. Bad or missing file: " << filename << endl);
				return 0;
			}
			LOG("done\n");

		//	adjust the grid
			funcAdjustGrid(mg, sh);

			LOG("  performing load balancing\n");

		//	perform load-balancing
		//TODO: if grid partitioning fails, the whole process should be aborted.
		//		this has to be communicated to the other processes.
			SubsetHandler shPartition(mg);
			if(!funcPartitionGrid(shPartition, mg, sh, numProcs)){
				UG_LOG("  grid partitioning failed. proceeding anyway...\n");
			}

		//	get min and max num elements
			int maxElems = 0;
			int minElems = 0;
			if(mg.num<Volume>() > 0){
				minElems = mg.num<Volume>();
				for(int i = 0; i < shPartition.num_subsets(); ++i){
					minElems = min(minElems, (int)shPartition.num<Volume>(i));
					maxElems = max(maxElems, (int)shPartition.num<Volume>(i));
				}
			}
			else if(mg.num<Face>() > 0){
				minElems = mg.num<Face>();
				for(int i = 0; i < shPartition.num_subsets(); ++i){
					minElems = min(minElems, (int)shPartition.num<Face>(i));
					maxElems = max(maxElems, (int)shPartition.num<Face>(i));
				}
			}

            LOG("  Element Distribution - min: " << minElems << ", max: " << maxElems << endl);
/*
			const char* partitionMapFileName = "partitionMap.obj";
			LOG("saving partition map to " << partitionMapFileName << endl);
			SaveGridToFile(mg, partitionMapFileName, shPartition);
*/
			//	distribute the grid.
			LOG("  distributing grid... ");
			bool bSuccess = false;
			if(keepSrcGrid)
				bSuccess = DistributeGrid_KeepSrcGrid(mg, sh, glmOut, shPartition, 0);
			else
				bSuccess = DistributeGrid(mg, sh, shPartition, 0,
										   &mgOut, &shOut, &glmOut);
			if(!bSuccess)
			{
				UG_LOG("failed!\n");
			}
			else{
				UG_LOG("done\n");
			}
		}
	}
	else
	{
	//	a grid will only be received, if the process-rank is smaller than numProcs
		if(pcl::GetProcRank() < numProcs){
			if(!ReceiveGrid(mgOut, shOut, glmOut, 0, keepSrcGrid)){
				UG_LOG("  ReceiveGrid failed on process " << pcl::GetProcRank() <<
				". Aborting...\n");
				return false;
			}
		}
	}

//	tell the distGridMgr that the associated layout changed.
	distGridMgrOut.grid_layouts_changed(true);

	return true;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
bool AdjustAndDistributeGrid(DistributedGridManager& distGridMgrOut,
						    ISubsetHandler& shOut,
							MultiGrid& srcGrid,
							ISubsetHandler& srcSh,
							int numProcs,
							bool keepSrcGrid,
							FuncAdjustGrid funcAdjustGrid,
							FuncPartitionGrid funcPartitionGrid,
							int rootProc)
{

	UG_ASSERT(dynamic_cast<MultiGrid*>(distGridMgrOut.get_assigned_grid()),
				"Error during LoadAndDistributeGrid: DistributedGridManager has to operate "
				"on a MultiGrid in the current implementation.");

	MultiGrid& mgOut = *dynamic_cast<MultiGrid*>(distGridMgrOut.get_assigned_grid());
	GridLayoutMap& glmOut = distGridMgrOut.grid_layout_map();

	if(keepSrcGrid){
		if(&srcGrid != &mgOut){
			throw(UGError("if keepSrcGrid is enabled, srcGrid and the grid in distGridMgrOut have to be identical instances."));
		}
	}
	else{
		if(&srcGrid == &mgOut){
			throw(UGError("if keepSrcGrid is disabled, srcGrid and the grid in distGridMgrOut have to be different instances."));
		}
	}

	if(shOut.grid() != &mgOut){
		throw(UGError("shOut has to operate on the grid in distGridMgrOut."));
	}

//	make sure that each grid has a position attachment - even if no data
//	will be received.
	if(!mgOut.has_vertex_attachment(aPosition))
		mgOut.attach_to_vertices(aPosition);

//	process 0 adjusts and distributes the grid. The others receive it.
	if(pcl::GetProcRank() == rootProc)
	{
	//	adjust and distribute the grid
		{
		//	adjust the grid
			funcAdjustGrid(srcGrid, srcSh);

			LOG("  performing load balancing\n");

		//	perform load-balancing
		//TODO: if grid partitioning fails, the whole process should be aborted.
		//		this has to be communicated to the other processes.
			SubsetHandler shPartition(srcGrid);
			if(!funcPartitionGrid(shPartition, srcGrid, srcSh, numProcs)){
				UG_LOG("  grid partitioning failed. proceeding anyway...\n");
			}

		//	get min and max num elements
			int maxElems = 0;
			int minElems = 0;
			if(srcGrid.num<Volume>() > 0){
				minElems = srcGrid.num<Volume>();
				for(int i = 0; i < shPartition.num_subsets(); ++i){
					minElems = min(minElems, (int)shPartition.num<Volume>(i));
					maxElems = max(maxElems, (int)shPartition.num<Volume>(i));
				}
			}
			else if(srcGrid.num<Face>() > 0){
				minElems = srcGrid.num<Face>();
				for(int i = 0; i < shPartition.num_subsets(); ++i){
					minElems = min(minElems, (int)shPartition.num<Face>(i));
					maxElems = max(maxElems, (int)shPartition.num<Face>(i));
				}
			}

            LOG("  Element Distribution - min: " << minElems << ", max: " << maxElems << endl);
/*
			const char* partitionMapFileName = "partitionMap.ugx";
			LOG("saving partition map to " << partitionMapFileName << endl);
			SaveGridToFile(srcGrid, partitionMapFileName, shPartition);
*/
			//	distribute the grid.
			LOG("  distributing grid... ");
			bool bSuccess = false;
			if(keepSrcGrid)
			//todo	currently the genealogy is distributed. This should be
			//		controllable from the outside.
				bSuccess = DistributeGrid_KeepSrcGrid(mgOut, shOut, glmOut,
													shPartition, rootProc, false);
			else
				bSuccess = DistributeGrid(srcGrid, srcSh, shPartition, rootProc,
										   &mgOut, &shOut, &glmOut);
			if(!bSuccess)
			{
				UG_LOG("failed!\n");
			}
			else{
				UG_LOG("done\n");
			}
		}
	}
	else
	{
	//	a grid will only be received, if the process-rank is smaller than numProcs
		if(pcl::GetProcRank() < numProcs){
			if(!ReceiveGrid(mgOut, shOut, glmOut, rootProc, keepSrcGrid)){
				UG_LOG("  ReceiveGrid failed on process " << pcl::GetProcRank() <<
				". Aborting...\n");
				return false;
			}
		}
	}

//	tell the distGridMgr that the associated layout changed.
	distGridMgrOut.grid_layouts_changed(true);

	return true;
}



template <class TElem, class TAVrtPos>
class ToElementPosition
{
	public:
		typedef typename TAVrtPos::ValueType	TValue;

		ToElementPosition(Grid& g, TAVrtPos& aPos)
		{
			if(g.has_vertex_attachment(aPos))
				m_aaPos.access(g, aPos);
		}

		TValue operator() (VertexBase* e)	{return m_aaPos[e];}
		TValue operator() (EdgeBase* e)		{return CalculateCenter(e, m_aaPos);}
		TValue operator() (Face* e)			{return CalculateCenter(e, m_aaPos);}
		TValue operator() (Volume* e)		{return CalculateCenter(e, m_aaPos);}

	private:
		Grid::VertexAttachmentAccessor<TAVrtPos>	m_aaPos;
};


template <class TAPos>
bool TestGridLayoutMap(MultiGrid& mg, GridLayoutMap& glm, TAPos& aPos)
{
	typedef typename TAPos::ValueType	TValue;
	typedef VertexLayout::LevelLayout	VrtLevelLayout;

	bool bSuccess = true;

//	check the interfaces
	pcl::InterfaceCommunicator<VertexLayout::LevelLayout> com;
	pcl::ProcessCommunicator procCom;

	ToElementPosition<VertexBase, TAPos> toPos(mg, aPos);
	boost::function<TValue (VertexBase*)> cbToPos = toPos;

	UG_LOG("Testing horizontal layouts...\n");
	{
		VertexLayout& masterLayout = glm.get_layout<VertexBase>(INT_H_MASTER);
		VertexLayout& slaveLayout = glm.get_layout<VertexBase>(INT_H_SLAVE);

	//	we have to retrieve the max level of all layouts
		int locMaxLevel = max(slaveLayout.num_levels(), masterLayout.num_levels());
		int globMaxLevel = locMaxLevel;
		procCom.allreduce(&locMaxLevel, &globMaxLevel, 1, PCL_DT_INT, PCL_RO_MAX);

		for(int i = 0; i < globMaxLevel; ++i){
			UG_LOG("Testing VertexLayout on level " << i << ":" << endl);
			bSuccess &= pcl::TestLayout<VrtLevelLayout, TValue>(procCom, com, masterLayout.layout_on_level(i),
											slaveLayout.layout_on_level(i), true, toPos, true);
		}
	}

	UG_LOG("\nTesting vertical layouts...\n");
	{
		VertexLayout& masterLayout = glm.get_layout<VertexBase>(INT_V_MASTER);
		VertexLayout& slaveLayout = glm.get_layout<VertexBase>(INT_V_SLAVE);
		int locMaxLevel = max(slaveLayout.num_levels(), masterLayout.num_levels());
		int globMaxLevel = locMaxLevel;
		procCom.allreduce(&locMaxLevel, &globMaxLevel, 1, PCL_DT_INT, PCL_RO_MAX);

		for(int i = 0; i < globMaxLevel; ++i){
			UG_LOG("Testing VertexLayout on level " << i << ":" << endl);
			bSuccess &= pcl::TestLayout<VrtLevelLayout, TValue>(procCom, com, masterLayout.layout_on_level(i),
											slaveLayout.layout_on_level(i), true, toPos, true);
		}
	}

	return bSuccess;
}

bool TestGridLayoutMap(MultiGrid& mg, GridLayoutMap& glm)
{
	if(mg.has_vertex_attachment(aPosition))
		return TestGridLayoutMap(mg, glm, aPosition);
	else if(mg.has_vertex_attachment(aPosition2))
		return TestGridLayoutMap(mg, glm, aPosition2);
	else if(mg.has_vertex_attachment(aPosition1))
		return TestGridLayoutMap(mg, glm, aPosition1);
	else
		UG_LOG("ERROR in TestGridLayoutMap: A standard position attachment"
				" is required.\n");
	return false;
}

}//	end of namespace

