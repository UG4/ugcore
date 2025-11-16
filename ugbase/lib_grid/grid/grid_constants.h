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

#ifndef __H__LIB_GRID__GRID_CONSTANTS__
#define __H__LIB_GRID__GRID_CONSTANTS__

namespace ug
{

/// \addtogroup lib_grid
/// @{

////////////////////////////////////////////////////////////////////////
//	VertexOptions
///	Used to specify the way in which Grid manages vertex-specific data.
enum VertexOptions
{
	VRTOPT_NONE = 0x00000000,
	VRTOPT_STORE_ASSOCIATED_EDGES = 0x00000001,
	VRTOPT_STORE_ASSOCIATED_FACES = 0x00000002,
	VRTOPT_STORE_ASSOCIATED_VOLUMES = 0x00000004,
};

////////////////////////////////////////////////////////////////////////
//	EdgeOptions
///	Used to specify the way in which Grid manages edge-specific data.
enum EdgeOptions
{
	EDGEOPT_NONE = 0x00000000,
	EDGEOPT_STORE_ASSOCIATED_FACES = 0x00000100,
	EDGEOPT_STORE_ASSOCIATED_VOLUMES = 0x00000200,
};

////////////////////////////////////////////////////////////////////////
//	FaceOptions
///	Used to specify the way in which Grid manages face-specific data.
enum FaceOptions
{
	FACEOPT_NONE = 0x00000000,
	FACEOPT_STORE_ASSOCIATED_EDGES = 0x00010000,	///< minor speed-improvement for grid.get_edge(Face*, int)
	FACEOPT_STORE_ASSOCIATED_VOLUMES = 0x00020000,
	FACEOPT_AUTOGENERATE_EDGES = 0x00080000
};

////////////////////////////////////////////////////////////////////////
//	VolumeOptions
///	Used to specify the way in which Grid manages volume-specific data.
enum VolumeOptions
{
	VOLOPT_NONE = 0x00000000,
	VOLOPT_STORE_ASSOCIATED_EDGES = 0x01000000,		///< minor speed-improvement for grid.get_edge(Volume*, int)
	VOLOPT_STORE_ASSOCIATED_FACES = 0x02000000,		///< speed-improvement for grid.get_face(Face*, int) ~15%
	VOLOPT_AUTOGENERATE_EDGES = 0x08000000,
	VOLOPT_AUTOGENERATE_FACES = 0x10000000
};

////////////////////////////////////////////////////////////////////////
//	GridOptions
///	Specify how references between associated objects are stored in a grid.
enum GridOptions
{
	GRIDOPT_NONE = 0x00000000,
	GRIDOPT_NO_INTERCONNECTION = 0x00000000,
///	vertices store lists of associated geometric objects.
/**	Note that this is the minimal required interconnection for many dynamic
 * algorithms (i.e. deleting an object from a grid automatically enables this
 * option.*/
	GRIDOPT_VERTEXCENTRIC_INTERCONNECTION = VRTOPT_STORE_ASSOCIATED_EDGES
											| VRTOPT_STORE_ASSOCIATED_FACES
											| VRTOPT_STORE_ASSOCIATED_VOLUMES,

///	sides are automatically created
	GRIDOPT_AUTOGENERATE_SIDES = FACEOPT_AUTOGENERATE_EDGES
									| VOLOPT_AUTOGENERATE_FACES,

///	All elements store references to associated lower dimensional geometric objects
/**	Additionally GRIDOPT_VERTEXCENTRIC_INTERCONNECTION is used.*/
	GRIDOPT_STANDARD_INTERCONNECTION = GRIDOPT_VERTEXCENTRIC_INTERCONNECTION
										| GRIDOPT_AUTOGENERATE_SIDES
										| FACEOPT_STORE_ASSOCIATED_EDGES
										| VOLOPT_STORE_ASSOCIATED_EDGES
										| VOLOPT_STORE_ASSOCIATED_FACES,

///	All elements store references to all associated elements
/**	This includes GRIDOPT_VERTEXCENTRIC_INTERCONNECTION and
 * GRIDOPT_STANDARD_INTERCONNECTION.*/
	GRIDOPT_FULL_INTERCONNECTION = GRIDOPT_STANDARD_INTERCONNECTION
									| EDGEOPT_STORE_ASSOCIATED_FACES
									| EDGEOPT_STORE_ASSOCIATED_VOLUMES
									| FACEOPT_STORE_ASSOCIATED_VOLUMES,

	GRIDOPT_DEFAULT = GRIDOPT_VERTEXCENTRIC_INTERCONNECTION
					  | GRIDOPT_AUTOGENERATE_SIDES
};

/// @}
}

#endif
