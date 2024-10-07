/*
 * ArteExpandFracs3D.h
 *
 *  Created on: 06.10.2024
 *      Author: Markus M. Knodel
 *
 *  * expand fractures using the Arte algorithm, 3D case
 *
 * Author: Markus Knodel, inspired by Arte from Fuchs and Sebastian Reiters code for fracture expansion without Arte
 *
 * implementing a class that gives the basic tools for Arte in 3D
 * might be templated at a later stage to fulfill 2D and 3D purposes, if suitable
 * ( so far 2D case one entire function, not so perfect, but running....)
 *
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



#ifndef UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_ARTEEXPANDFRACS3D_H_
#define UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_ARTEEXPANDFRACS3D_H_

namespace ug {

class ArteExpandFracs3D
{

public:

	ArteExpandFracs3D( Grid & grid, SubsetHandler & sh,
		    std::vector<FractureInfo> const & fracInfos,
			bool useTrianglesInDiamonds, bool establishDiamonds );

	virtual ~ArteExpandFracs3D();


public:

	bool run();

private:

	Grid & m_grid;
	SubsetHandler & m_sh;
	std::vector<FractureInfo> m_fracInfos;
	bool m_useTrianglesInDiamonds, m_establishDiamonds;

	Grid::VertexAttachmentAccessor<APosition> m_aaPos;

	//	objects for temporary results
	FaceDescriptor m_facDescr;
	VolumeDescriptor m_volDescr;

	std::vector<Edge*> m_tmpEdges; // used for temporary results.
	std::vector<Face*> m_tmpFaces; // used for temporary results.
	std::vector<Volume*> m_tmpVols; // used for temporary results.



	bool initialize();

};

} /* namespace ug */

#endif /* UGCORE_UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_ARTEEXPANDFRACS3D_H_ */
