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

#ifndef __H__LIBGRID__EXTRUDE__
#define __H__LIBGRID__EXTRUDE__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{
/// \addtogroup lib_grid_algorithms_extrusion
///	@{

///	used to specify the behavior of ug::Extrude.
enum ExtrusionOptions
{
	EO_CREATE_FACES 	= 1,
	EO_CREATE_VOLUMES	= 1 << 1,

	EO_DEFAULT 			= 	EO_CREATE_FACES |
							EO_CREATE_VOLUMES
};

////////////////////////////////////////////////////////////////////////
//	Extrude
///	extrudes geometry and creates new edges, faces and volumes.
/**
 * For each vertex in pvVerticesInOut a new edge will be created.
 * NULL is a valid parameter for pvVerticesInOut.
 *
 * For each edge in pvEdgesInOut a new edge.
 * if EO_CREATE_FACES is enabled, a new face that connects the two
 * edges will be created too.
 * NULL is a valid parameter for pvEdgesInOut.
 *
 * For each face in pvFacesInOut a new face will be created.
 * if EO_CREATE_VOLUMES is enabled a new volume that connects the
 * two faces will be created too.
 * NULL is a valid parameter for pvFacesInOut.
 *
 * The element from which new elements are extruded is passed to
 * each new element as the parent element.
 *
 * Pass any or-combination of const in ExtrusionOptions as extrusionOptions.
 *
 * After the algorithm has finished, the in-out-vectors hold the elements
 * which have been created from the input-elements of those vectors.
 * The in-out-vectors thus can be directly used in a new call to Extrude.
 *
 * All newly created volume-elements will be pushed to pvVolsOut if that vector
 * was specified (NULL by default). Volume elements are created in the order in
 * which faces are specified in pvFacesInOut.
 *
 * If you need to have access to all newly created elements you could use
 * a ug::Selector with enabled autoselection.
 * \{ */
template <class vector_t>
void Extrude(Grid& grid,
			std::vector<Vertex*>* pvVerticesInOut,
			std::vector<Edge*>* pvEdgesInOut,
			std::vector<Face*>* pvFacesInOut,
			const vector_t& direction,
			uint extrusionOptions = EO_DEFAULT,
			Attachment<vector_t>& aPos = aPosition,
			std::vector<Volume*>* pvVolsOut = NULL);

template <class TAAPos>
void Extrude(Grid& grid,
			std::vector<Vertex*>* pvVerticesInOut,
			std::vector<Edge*>* pvEdgesInOut,
			std::vector<Face*>* pvFacesInOut,
			const typename TAAPos::ValueType& direction,
			TAAPos aaPos,
			uint extrusionOptions = EO_DEFAULT,
			std::vector<Volume*>* pvVolsOut = NULL);
/** \} */

/*
void Extrude(Grid& grid,
			Selector& sel,
			const vector3& direction,
			uint extrusionOptions = EO_DEFAULT,
			APosition& aPos = aPosition);
*/

/// @}

}//	end of namespace


#endif
