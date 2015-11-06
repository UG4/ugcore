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

#include <string>
#include <vector>
#include "common/math/ugmath.h"
#include "../grid/grid.h"
#include "../common_attachments.h"
#include "../subset_handler.h"

#ifndef __H__LIB_GRID__FILE_IO_OBJ__
#define __H__LIB_GRID__FILE_IO_OBJ__

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	OBJMaterial
///	holds information about the materials in the obj-file.
struct OBJMaterial
{
	std::string 	m_strName;
	std::string 	m_strTextureDiffuse;
	vector4			m_vDiffuse;
	float			m_fAlpha;
};

////////////////////////////////////////////////////////////////////////
//	LoadGridFromOBJ
///	Loads a file from a wavefront '.obj' file. Fills optional subset-infos.
/**
 * paTex, pSubsetHandler and pvMaterials may be specified optionally.
 * aPos and paTex are used as vertex-attachments.
 *
 * \param grid: Grid
 * \param filename: Filename
 * \param aPos: position attachment
 * \param paTexCoord: vertex attachment. If specified texture coords will be stored in the
 * 		grids vertices.
 *
 * \param pSubsetHandler: if specified then faces are assigned to subsets depending
 * 		on their associated objects.
 *
 * \param pvMaterials: Holds material data.
 */
bool LoadGridFromOBJ(Grid& grid, const char* filename, AVector3& aPos = aPosition,
		AVector2* paTexCoord = NULL,
		ISubsetHandler* pSubsetHandler = NULL,
		std::vector<OBJMaterial>* pvMaterials = NULL);

////////////////////////////////////////////////////////////////////////
//	saves a file to obj
///	Saves a file from a wavefront '.obj' file. Writes optional subset-infos.
bool SaveGridToOBJ(Grid& grid, const char* filename, AVector3& aPos = aPosition,
		AVector2* paTexCoord = NULL,
		ISubsetHandler* pSubsetHandler = NULL,
		std::vector<OBJMaterial>* pvMaterials = NULL);
};

#endif
