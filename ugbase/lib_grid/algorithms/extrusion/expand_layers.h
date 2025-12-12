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

#ifndef __H__UG__expand_layers__
#define __H__UG__expand_layers__

#include <vector>

#include "lib_grid/lg_base.h"

namespace ug {

/// Used to tell ExpandLayers_... which subsets should be regarded as layers.
/**
 * 	- subsetIndex defines the subset of the source (low dimensional) layer.
 *	- newSubsetIndex defines the subset into which the newly generated elements will go.
 *	- width describes the width to which a layer shall be expanded.
 */
struct FractureInfo{
	FractureInfo(int subsetInd, int newSubsetInd, double w) :
		subsetIndex(subsetInd), newSubsetIndex(newSubsetInd), width(w)	{}

	int subsetIndex;
	int newSubsetIndex;
	double width;
};

/**
 * This algorithm indirectly uses Grid::mark.
 *
 * 1 dimensional fractures specified in fracInfos are expanded to 2 dimensional subsets.
 * the resulting fractures will then consist of 2 layers of quadrilaterals. On the
 * boundaries triangles are inserted.
 *
 * Through expandFracBoundaries you can tell the algorithm whether inner fracture
 * boundaries shall be expanded. Note that this means that an additional node is
 * introduced at each inner fracture boundary vertex and that the associated
 * fracture elements are connected at two sides.
 * Note that fractures are always expanded at boundaries which lie on the geometries
 * boundary.
 *
 *	This algorithm requires the option FACEOPT_AUTOGENERATE_EDGES.
 *	The option is automatically enabled if required.
 */
bool ExpandFractures2d(Grid& grid, SubsetHandler& sh,
						const std::vector<FractureInfo>& fracInfos,
						bool expandInnerFracBnds, bool expandOuterFracBnds);


/**
 * This algorithm indirectly uses Grid::mark.
 *
 * 2 dimensional fractures specified in fracInfos are expanded to 3 dimensional subsets.
 * the resulting fractures will then consist of 2 layers of hexahedrons. On the
 * boundaries tetrahedrons, prisms and pyramids are inserted.
 *
 * Through expandFracBoundaries you can tell the algorithm whether inner fracture
 * boundaries shall be expanded. Note that this means that an additional node is
 * introduced at each inner fracture boundary vertex and that the associated
 * fracture elements are connected at two sides.
 * Note that fractures are always expanded at boundaries which lie on the geometries
 * boundary.
 *
 *	This algorithm requires the option FACEOPT_AUTOGENERATE_EDGES.
 *	The option is automatically enabled if required.
 *
 *	This algorithm requires the option VOLOPT_AUTOGENERATE_FACES.
 *	The option is automatically enabled if required.
 */
bool ExpandFractures3d(Grid& grid, SubsetHandler& sh,
						const std::vector<FractureInfo>& fracInfos,
						bool expandInnerFracBnds, bool expandOuterFracBnds);

}//	end of namespace

#endif
