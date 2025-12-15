/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

#ifndef __H__UG__REFERENCE_MAPPING_TEST__
#define __H__UG__REFERENCE_MAPPING_TEST__

#include <vector>
#include <cassert>

#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/local_finite_element/local_dof_set.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
//#include "lib_disc/reference_element/reference_element_util.h"

namespace ug {

//
//	Experimental reference mapping test procedures for octahedrons, tetrahedrons and edges
//

void OctReferenceMappingTest(std::vector<number> vCornerCoord0, std::vector<number> vCornerCoord1, std::vector<number> vCornerCoord2,
						  	 std::vector<number> vCornerCoord3, std::vector<number> vCornerCoord4, std::vector<number> vCornerCoord5,
							 std::vector<number> vGlobPos)
{
	UG_LOG(">> Starting OctReferenceMappingTest: " << std::endl);

	std::vector<MathVector<3, number> > vCornerCoords;
	std::vector<MathVector<3, number> > vLocPos(1, MathVector<3, number>(0.0));
	std::vector<MathVector<3, number> > vGlobPositions;

	MathVector<3> GlobPos(vGlobPos[0], vGlobPos[1], vGlobPos[2]);
	vGlobPositions.push_back(GlobPos);

	MathVector<3> CornerCoord0(vCornerCoord0[0], vCornerCoord0[1], vCornerCoord0[2]);
	MathVector<3> CornerCoord1(vCornerCoord1[0], vCornerCoord1[1], vCornerCoord1[2]);
	MathVector<3> CornerCoord2(vCornerCoord2[0], vCornerCoord2[1], vCornerCoord2[2]);
	MathVector<3> CornerCoord3(vCornerCoord3[0], vCornerCoord3[1], vCornerCoord3[2]);
	MathVector<3> CornerCoord4(vCornerCoord4[0], vCornerCoord4[1], vCornerCoord4[2]);
	MathVector<3> CornerCoord5(vCornerCoord5[0], vCornerCoord5[1], vCornerCoord5[2]);

	vCornerCoords.push_back(CornerCoord0);
	vCornerCoords.push_back(CornerCoord1);
	vCornerCoords.push_back(CornerCoord2);
	vCornerCoords.push_back(CornerCoord3);
	vCornerCoords.push_back(CornerCoord4);
	vCornerCoords.push_back(CornerCoord5);

	try
	{
		DimReferenceMapping<3, 3>& map = ReferenceMappingProvider::get<3, 3>(ROID_OCTAHEDRON, vCornerCoords);
		map.global_to_local(vLocPos, vGlobPositions);
	}
	UG_CATCH_THROW("OctReferenceMappingTest() could not map global to local.");

	UG_LOG("Calculated vLocPos: " << vLocPos[0] << std::endl);
}


void TetReferenceMappingTest(std::vector<number> vCornerCoord0, std::vector<number> vCornerCoord1, std::vector<number> vCornerCoord2,
						  	 std::vector<number> vCornerCoord3, std::vector<number> vGlobPos)
{
	UG_LOG(">> Starting TetReferenceMappingTest: " << std::endl);

	std::vector<MathVector<3, number> > vCornerCoords;
	std::vector<MathVector<3, number> > vLocPos(1, MathVector<3, number>(0.0));
	std::vector<MathVector<3, number> > vGlobPositions;

	MathVector<3> GlobPos(vGlobPos[0], vGlobPos[1], vGlobPos[2]);
	vGlobPositions.push_back(GlobPos);

	MathVector<3> CornerCoord0(vCornerCoord0[0], vCornerCoord0[1], vCornerCoord0[2]);
	MathVector<3> CornerCoord1(vCornerCoord1[0], vCornerCoord1[1], vCornerCoord1[2]);
	MathVector<3> CornerCoord2(vCornerCoord2[0], vCornerCoord2[1], vCornerCoord2[2]);
	MathVector<3> CornerCoord3(vCornerCoord3[0], vCornerCoord3[1], vCornerCoord3[2]);

	vCornerCoords.push_back(CornerCoord0);
	vCornerCoords.push_back(CornerCoord1);
	vCornerCoords.push_back(CornerCoord2);
	vCornerCoords.push_back(CornerCoord3);

	try
	{
		DimReferenceMapping<3, 3>& map = ReferenceMappingProvider::get<3, 3>(ROID_TETRAHEDRON, vCornerCoords);
		map.global_to_local(vLocPos, vGlobPositions);
	}
	UG_CATCH_THROW("TetReferenceMappingTest() could not map global to local.");

	UG_LOG("Calculated vLocPos: " << vLocPos[0] << std::endl);
	UG_LOG("PointIsInsideTetrahedron: " << PointIsInsideTetrahedron(GlobPos, CornerCoord0, CornerCoord1, CornerCoord2, CornerCoord3) << std::endl);
}


void EdgeReferenceMappingTest(std::vector<number> vCornerCoord0, std::vector<number> vCornerCoord1,
							  std::vector<number> vGlobPos)
{
	UG_LOG(">> Starting EdgeReferenceMappingTest: " << std::endl);

	std::vector<MathVector<3, number> > vCornerCoords;
	std::vector<MathVector<1, number> > vLocPos(1, MathVector<1, number>(0.0));
	std::vector<MathVector<3, number> > vGlobPositions;

	MathVector<3> GlobPos(vGlobPos[0], vGlobPos[1], vGlobPos[2]);
	vGlobPositions.push_back(GlobPos);

	MathVector<3> CornerCoord0(vCornerCoord0[0], vCornerCoord0[1], vCornerCoord0[2]);
	MathVector<3> CornerCoord1(vCornerCoord1[0], vCornerCoord1[1], vCornerCoord1[2]);

	vCornerCoords.push_back(CornerCoord0);
	vCornerCoords.push_back(CornerCoord1);

	try
	{
		DimReferenceMapping<1, 3>& map = ReferenceMappingProvider::get<1, 3>(ROID_EDGE, vCornerCoords);
		map.global_to_local(vLocPos, vGlobPositions);
	}
	UG_CATCH_THROW("EdgeReferenceMappingTest() could not map global to local.");

	UG_LOG("Calculated vLocPos: " << vLocPos[0] << std::endl);
}

}//	end of namespace

#endif
