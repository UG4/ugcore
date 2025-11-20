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

#include "grid_object_collection.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	GridObjectCollection
GridObjectCollection::
GridObjectCollection(size_t levelEstimate)
{
	m_levels.reserve(levelEstimate);
}

GridObjectCollection::
GridObjectCollection(ElementStorage<Vertex>::SectionContainer* vrtCon,
							ElementStorage<Edge>::SectionContainer* edgeCon,
							ElementStorage<Face>::SectionContainer* faceCon,
							ElementStorage<Volume>::SectionContainer* volCon)
{
	m_levels.reserve(1);
	add_level(vrtCon, edgeCon, faceCon, volCon);
}

GridObjectCollection::
GridObjectCollection(const GridObjectCollection& mgoc)
{
	assign(mgoc);
}

GridObjectCollection&
GridObjectCollection::
operator = (const GridObjectCollection& mgoc)
{
	assign(mgoc);
	return *this;
}

void
GridObjectCollection::
assign(const GridObjectCollection& mgoc)
{
	m_levels.resize(mgoc.num_levels());
	for(size_t i = 0; i < m_levels.size(); ++i)
		m_levels[i] = mgoc.m_levels[i];
}

void
GridObjectCollection::
add_level(ElementStorage<Vertex>::SectionContainer* vrtCon,
			ElementStorage<Edge>::SectionContainer* edgeCon,
			ElementStorage<Face>::SectionContainer* faceCon,
			ElementStorage<Volume>::SectionContainer* volCon)
{
	m_levels.emplace_back(	vrtCon,
											edgeCon,
											faceCon,
											volCon);
}


GridObjectCollection::ContainerCollection::
ContainerCollection(ElementStorage<Vertex>::SectionContainer* vrtCon,
					ElementStorage<Edge>::SectionContainer* edgeCon,
					ElementStorage<Face>::SectionContainer* faceCon,
					ElementStorage<Volume>::SectionContainer* volCon)
{
	vrtContainer = vrtCon;
	edgeContainer = edgeCon;
	faceContainer = faceCon;
	volContainer = volCon;
}

}//	end of namespace
