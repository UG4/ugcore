/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#include "grid_statistics.h"
#include "geom_obj_util/volume_util.h"

using namespace std;

namespace ug
{


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	AssignTetrahedronAttributesByAspectRatio
/// assigns attributes to tetrahedral elements of a grid respecting their aspect ratio
bool AssignTetrahedronAttributesByAspectRatio(	Grid& grid,
												SubsetHandler& shVolume,
												AInt& aTetrahedronAspectRatioClass,
												vector<double>& offsets,
												Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	VolumeSelector VolumeSel(grid);
	Grid::VolumeAttachmentAccessor<AInt> aaTARC(grid, aTetrahedronAspectRatioClass);

//	order the volume subsets correctly, regarding user-specified aElementAttribute
	for(int i = 1; i < (int)offsets.size(); ++i)
	{
		int moveTo = i;
		for(int j = i-1; j >= 0; --j)
		{
			if(offsets[i] > offsets[j])
				moveTo = j;
			else
				break;
		}

		if(moveTo != i)
		{
			swap(offsets[i], offsets[moveTo]);
		}
	}

//	add zero as last offset
	offsets.push_back(0.0);

//	iterate through all tetrahedrons and assign attributes by their aspect ratio
	for(TetrahedronIterator vIter = grid.begin<Tetrahedron>(); vIter != grid.end<Tetrahedron>(); ++vIter)
	{
		Tetrahedron* tet = *vIter;
		number AspectRatio = CalculateTetrahedronAspectRatio(grid, tet, aaPos);

		for(size_t i = 0; i<offsets.size(); ++i)
		{
			if(AspectRatio >= offsets[i])
			{
				shVolume.assign_subset(tet, i);
				aaTARC[tet] = i;
				break;
			}
		}
	}

	return true;
}

}//	end of namespace
