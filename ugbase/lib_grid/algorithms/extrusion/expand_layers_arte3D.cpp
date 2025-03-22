/*
 *  * expand fractures using the Arte algorithm, 3D case
 *
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Knodel, inspired by Arte from Fuchs and Sebastian Reiters code for fracture expansion without Arte
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

#include <boost/function.hpp>
#include <stack>
#include <vector>
#include "lib_grid/lg_base.h"
#include "expand_layers.h"
#include "expand_layers_arte.h"
#include "expand_layers_arte3D.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/grid/grid_util.h"
//#include "lib_grid/util/simple_algebra/least_squares_solver.h"

#include <vector>

#include "lib_grid/algorithms/extrusion/ArteExpandFracs3D.h"

using namespace std;

namespace ug{


bool ExpandFractures3dArte( Grid& grid, SubsetHandler& sh,
						    std::vector<FractureInfo> const & fracInfos,
							bool useTrianglesInDiamonds, bool establishDiamonds )
{

	bool need2Restart = false;

	bool runResult = false;

	do
	{

		ArteExpandFracs3D ef3dA ( grid, sh, fracInfos,
							  useTrianglesInDiamonds, establishDiamonds );

		runResult = ef3dA.run( need2Restart );

//		return runResult;

		if( runResult )
		{
			return true;
		}

	} while( need2Restart );

	return runResult;

}



}// end of namespace


