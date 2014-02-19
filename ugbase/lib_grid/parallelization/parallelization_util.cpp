//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d30

#include <sstream>
#include "parallelization_util.h"
#include "load_balancing.h"
#include "distribution.h"
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

		TValue operator() (Vertex* e)	{return m_aaPos[e];}
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
	typedef EdgeLayout::LevelLayout		EdgeLevelLayout;
	typedef FaceLayout::LevelLayout		FaceLevelLayout;

	bool bSuccess = true;
	const bool verbose = true;

//	check the interfaces
	pcl::InterfaceCommunicator<VertexLayout::LevelLayout> comVrt;
	pcl::ProcessCommunicator procCom;

	ToElementPosition<Vertex, TAPos> toPosVrt(mg, aPos);

	UG_LOG("Testing horizontal vertex layouts...\n");
	{
		VertexLayout& masterLayout = glm.get_layout<Vertex>(INT_H_MASTER);
		VertexLayout& slaveLayout = glm.get_layout<Vertex>(INT_H_SLAVE);

	//	we have to retrieve the max level of all layouts
		int locMaxLevel = max(slaveLayout.num_levels(), masterLayout.num_levels());
		int globMaxLevel = locMaxLevel;
		procCom.allreduce(&locMaxLevel, &globMaxLevel, 1, PCL_DT_INT, PCL_RO_MAX);

		for(int i = 0; i < globMaxLevel; ++i){
			if(verbose){
				UG_LOG("Testing VertexLayout on level " << i << ":" << endl);
			}
			bSuccess &= pcl::TestLayout<VrtLevelLayout, TValue>(procCom, comVrt, masterLayout.layout_on_level(i),
											slaveLayout.layout_on_level(i), verbose, toPosVrt, true);
		}
	}

	UG_LOG("Testing vertical vertex layouts...\n");
	{
		VertexLayout& masterLayout = glm.get_layout<Vertex>(INT_V_MASTER);
		VertexLayout& slaveLayout = glm.get_layout<Vertex>(INT_V_SLAVE);
		int locMaxLevel = max(slaveLayout.num_levels(), masterLayout.num_levels());
		int globMaxLevel = locMaxLevel;
		procCom.allreduce(&locMaxLevel, &globMaxLevel, 1, PCL_DT_INT, PCL_RO_MAX);

		for(int i = 0; i < globMaxLevel; ++i){
			if(verbose){
				UG_LOG("Testing VertexLayout on level " << i << ":" << endl);
			}
			bSuccess &= pcl::TestLayout<VrtLevelLayout, TValue>(procCom, comVrt, masterLayout.layout_on_level(i),
											slaveLayout.layout_on_level(i), verbose, toPosVrt, true);
		}
	}


	pcl::InterfaceCommunicator<EdgeLayout::LevelLayout> comEdge;
	ToElementPosition<EdgeBase, TAPos> toPosEdge(mg, aPos);

	UG_LOG("Testing horizontal edge layouts...\n");
	{
		EdgeLayout& masterLayout = glm.get_layout<EdgeBase>(INT_H_MASTER);
		EdgeLayout& slaveLayout = glm.get_layout<EdgeBase>(INT_H_SLAVE);

	//	we have to retrieve the max level of all layouts
		int locMaxLevel = max(slaveLayout.num_levels(), masterLayout.num_levels());
		int globMaxLevel = locMaxLevel;
		procCom.allreduce(&locMaxLevel, &globMaxLevel, 1, PCL_DT_INT, PCL_RO_MAX);

		for(int i = 0; i < globMaxLevel; ++i){
			if(verbose){
				UG_LOG("Testing EdgeLayout on level " << i << ":" << endl);
			}
			bSuccess &= pcl::TestLayout<EdgeLevelLayout, TValue>(procCom, comEdge, masterLayout.layout_on_level(i),
											slaveLayout.layout_on_level(i), verbose, toPosEdge, true);
		}
	}

	UG_LOG("Testing vertical edge layouts...\n");
	{
		EdgeLayout& masterLayout = glm.get_layout<EdgeBase>(INT_V_MASTER);
		EdgeLayout& slaveLayout = glm.get_layout<EdgeBase>(INT_V_SLAVE);
		int locMaxLevel = max(slaveLayout.num_levels(), masterLayout.num_levels());
		int globMaxLevel = locMaxLevel;
		procCom.allreduce(&locMaxLevel, &globMaxLevel, 1, PCL_DT_INT, PCL_RO_MAX);

		for(int i = 0; i < globMaxLevel; ++i){
			if(verbose){
				UG_LOG("Testing EdgeLayout on level " << i << ":" << endl);
			}
			bSuccess &= pcl::TestLayout<EdgeLevelLayout, TValue>(procCom, comEdge, masterLayout.layout_on_level(i),
											slaveLayout.layout_on_level(i), verbose, toPosEdge, true);
		}
	}

	pcl::InterfaceCommunicator<FaceLayout::LevelLayout> comFace;
	ToElementPosition<Face, TAPos> toPosFace(mg, aPos);

	UG_LOG("Testing horizontal face layouts...\n");
	{
		FaceLayout& masterLayout = glm.get_layout<Face>(INT_H_MASTER);
		FaceLayout& slaveLayout = glm.get_layout<Face>(INT_H_SLAVE);

	//	we have to retrieve the max level of all layouts
		int locMaxLevel = max(slaveLayout.num_levels(), masterLayout.num_levels());
		int globMaxLevel = locMaxLevel;
		procCom.allreduce(&locMaxLevel, &globMaxLevel, 1, PCL_DT_INT, PCL_RO_MAX);

		for(int i = 0; i < globMaxLevel; ++i){
			if(verbose){
				UG_LOG("Testing FaceLayout on level " << i << ":" << endl);
			}
			bSuccess &= pcl::TestLayout<FaceLevelLayout, TValue>(procCom, comFace, masterLayout.layout_on_level(i),
											slaveLayout.layout_on_level(i), verbose, toPosFace, true);
		}
	}

	UG_LOG("Testing vertical edge layouts...\n");
	{
		FaceLayout& masterLayout = glm.get_layout<Face>(INT_V_MASTER);
		FaceLayout& slaveLayout = glm.get_layout<Face>(INT_V_SLAVE);
		int locMaxLevel = max(slaveLayout.num_levels(), masterLayout.num_levels());
		int globMaxLevel = locMaxLevel;
		procCom.allreduce(&locMaxLevel, &globMaxLevel, 1, PCL_DT_INT, PCL_RO_MAX);

		for(int i = 0; i < globMaxLevel; ++i){
			if(verbose){
				UG_LOG("Testing FaceLayout on level " << i << ":" << endl);
			}
			bSuccess &= pcl::TestLayout<FaceLevelLayout, TValue>(procCom, comFace, masterLayout.layout_on_level(i),
											slaveLayout.layout_on_level(i), verbose, toPosFace, true);
		}
	}
	return pcl::AllProcsTrue(bSuccess);
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

