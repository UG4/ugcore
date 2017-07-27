/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_projectors
#define __H__UG_projectors

#include <boost/mpl/pair.hpp>
#include <boost/mpl/string.hpp>
#include <boost/mpl/vector.hpp>

#include "refinement_projector.h"
#include "cylinder_cut_projector.h"
#include "cylinder_projector.h"
#include "cylinder_volume_projector.h"
#include "plane_cut_projector.h"
#include "projection_handler.h"
#include "raster_layers_projector.h"
#include "smooth_projector.h"
#include "sphere_projector.h"
#include "subdivision_projector.h"
#include "neurite_projector.h"

namespace boost {
namespace mpl {
namespace tmp {

//	ATTENTION:	unfortunately, multi-chars may have a length of at most 4 to be
//				compatible with different compilers
//
//	NOTE:		ug::ProjectionHandler is not contained in this list on purpose,
//				since it isn't used for serialization.
//
//	NOTE:		ug::CylinderCutProjector and ug::PlaneCutProjector are not contained
//				in this list, since they are only usable in specialized algorithms.
typedef vector<
			pair <ug::RefinementProjector,	string<'defa','ult'> >,
            pair <ug::CylinderProjector,    string<'cyli','nder'> >,
            pair <ug::CylinderVolumeProjector,    string<'c','y','l','v','o','l'> >,
			pair <ug::SphereProjector,		string<'sphe','re'> >,
			pair <ug::SubdivisionProjector,	string<'subd','ivis', 'ion'> >,
			pair <ug::SmoothProjector,		string<'smoo','th'> >,
			pair <ug::RasterLayersProjector,string<'rast','er'> >,
			pair <ug::NeuriteProjector,     string<'neur','ite'> >
			>
	ProjectorTypes;	
}
}
}


namespace ug{

typedef boost::mpl::tmp::ProjectorTypes	ProjectorTypes;

}//	end of namespace

#endif	//__H__UG_projectors
