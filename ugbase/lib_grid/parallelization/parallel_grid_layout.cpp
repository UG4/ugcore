/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "parallel_grid_layout.h"
#include "pcl/pcl_layout_util.h"

namespace ug{
AGeomObjID 	aGeomObjID("globalID", false);

template <>
size_t hash_key<GeomObjID>(const GeomObjID& key)
{
//	of course this hash does not completly avoid collisions.
//	One should check whether the chosen key is fine.
	return (unsigned long)(99971 * key.first + key.second * key.second);
}

std::ostream& operator<<(std::ostream& out, const GeomObjID& goId)
{
	out << "(" << goId.first << ", " << goId.second << ")";
	return out;
}

bool operator<(const GeomObjID& gid1, const GeomObjID& gid2)
{
	if(gid1.first < gid2.first)
		return true;
	if(gid1.first > gid2.first)
		return false;
	return gid1.second < gid2.second;
}

///	A helper method for GridLayoutMap::remove_empty_interfaces()
template <class TGeomObj>
static void RemoveEmptyInterfaces(
		typename GridLayoutMap::Types<TGeomObj>::Map& map)
{
	using TLayout = typename GridLayoutMap::Types<TGeomObj>::Layout;

	for(auto layoutIter = map.begin(); layoutIter != map.end(); ++layoutIter)
	{
		TLayout& layout = layoutIter->second;
		RemoveEmptyInterfaces(layout);
	}
}

void GridLayoutMap::remove_empty_interfaces()
{
	RemoveEmptyInterfaces<Vertex>(m_vertexLayoutMap);
	RemoveEmptyInterfaces<Edge>(m_edgeLayoutMap);
	RemoveEmptyInterfaces<Face>(m_faceLayoutMap);
	RemoveEmptyInterfaces<Volume>(m_volumeLayoutMap);
}

}// end of namespace
