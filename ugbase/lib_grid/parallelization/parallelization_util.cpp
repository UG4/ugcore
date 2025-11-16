/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include <sstream>
#include "parallelization_util.h"
#include "distribution.h"
#include "lib_grid/refinement/global_multi_grid_refiner.h"
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

template <class TElem, class TAVrtPos>
class ToElementPosition
{
	public:
		using TValue = typename TAVrtPos::ValueType;

		ToElementPosition(Grid& g, TAVrtPos& aPos)
		{
			if(g.has_vertex_attachment(aPos))
				m_aaPos.access(g, aPos);
		}

		TValue operator() (Vertex* e)	{return m_aaPos[e];}
		TValue operator() (Edge* e)		{return CalculateCenter(e, m_aaPos);}
		TValue operator() (Face* e)			{return CalculateCenter(e, m_aaPos);}
		TValue operator() (Volume* e)		{return CalculateCenter(e, m_aaPos);}

	private:
		Grid::VertexAttachmentAccessor<TAVrtPos>	m_aaPos;
};


template <class TAPos>
bool TestGridLayoutMap(MultiGrid& mg, GridLayoutMap& glm, TAPos& aPos, bool verbose = true)
{
	using TValue = typename TAPos::ValueType;
	using VrtLevelLayout = VertexLayout::LevelLayout;
	using EdgeLevelLayout = EdgeLayout::LevelLayout;
	using FaceLevelLayout = FaceLayout::LevelLayout;

	bool bSuccess = true;

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
	ToElementPosition<Edge, TAPos> toPosEdge(mg, aPos);

	UG_LOG("Testing horizontal edge layouts...\n");
	{
		EdgeLayout& masterLayout = glm.get_layout<Edge>(INT_H_MASTER);
		EdgeLayout& slaveLayout = glm.get_layout<Edge>(INT_H_SLAVE);

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
		EdgeLayout& masterLayout = glm.get_layout<Edge>(INT_V_MASTER);
		EdgeLayout& slaveLayout = glm.get_layout<Edge>(INT_V_SLAVE);
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

	UG_LOG("Testing vertical face layouts...\n");
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

bool TestGridLayoutMap(MultiGrid& mg, GridLayoutMap& glm, bool verbose)
{
	if(mg.has_vertex_attachment(aPosition))
		return TestGridLayoutMap(mg, glm, aPosition, verbose);
	else if(mg.has_vertex_attachment(aPosition2))
		return TestGridLayoutMap(mg, glm, aPosition2, verbose);
	else if(mg.has_vertex_attachment(aPosition1))
		return TestGridLayoutMap(mg, glm, aPosition1, verbose);
	else
		UG_LOG("ERROR in TestGridLayoutMap: A standard position attachment"
				" is required.\n");
	return false;
}

}//	end of namespace

