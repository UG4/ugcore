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

#ifndef __H__LIBGRID__FILE_IO_TETGEN__
#define __H__LIBGRID__FILE_IO_TETGEN__

#include <vector>
#include "../grid/grid.h"
#include "../subset_handler.h"
#include "../common_attachments.h"

namespace ug
{

bool LoadGridFromELE(Grid& grid, const char* filename, ISubsetHandler* pSH = nullptr,
					APosition& aPos = aPosition);
					
bool SaveGridToELE(Grid& grid, const char* filename, ISubsetHandler* pSH = nullptr,
					APosition& aPos = aPosition,
					ANumber* paVolumeConstraint = nullptr);
					
					
bool ImportGridFromTETGEN(Grid& grid,
						const char* nodesFilename, const char* facesFilename,
						const char* elemsFilename, AVector3& aPos,
						std::vector<AFloat>* pvNodeAttributes = nullptr,
						AInt* paNodeBoundaryMarker = nullptr,
						AInt* paFaceBoundaryMarker = nullptr,
						AInt* paElementAttribute = nullptr);

bool ImportGridFromTETGEN(Grid& grid,
						const char* nodesFilename, const char* facesFilename,
						const char* elemsFilename, AVector3& aPos,
						ISubsetHandler* psh = nullptr,
						std::vector<AFloat>* pvNodeAttributes = nullptr);
					
////////////////////////////////////////////////////////////////////////
//	ExportGridToSMESH
///	writes an SMESH file.
/**
 * Be sure that grid consists of closed, 2-dimensional surfaces.
 * except of grid, filename and aPos, all parameters are optional.
 * Be sure that all attachments that you pass to the method have been properly attached.
 * If you don't want to specify a parameter simply pass nullptr.
 * meaning of optional parameters:
 * pvNodeAttributes: a vector that holds attachments of type AFloat,
					 which store vertex-attributes.
 * paNodeBoundaryMarker: specifies the boundary-index of a vertex.
 * paFaceBoundaryMarker: specifies the boundary-index of a face.
 * pvHoles: specifies where holes are located.
 * pvRegionPositions: specifies where regions are located.
 * vRegionAttributes: specifies the attribute of a region.
						Has to be of the same size as vRegionPositions.
						Will only be used if pvRegionAttributes has been specified.
 * vRegionVolumeConstraints: specifies the maximal volume of a tetrahedron in each region.
							 Has to be of the same size as vRegionPositions.
							 Will only be used if pvRegionAttributes has been specified.
 */
bool ExportGridToSMESH(Grid& grid, const char* filename, AVector3& aPos,
						std::vector<AFloat>* pvNodeAttributes = nullptr,
						AInt* paNodeBoundaryMarker = nullptr,
						AInt* paFaceBoundaryMarker = nullptr,
						std::vector<vector3>* pvHoles = nullptr,
						std::vector<vector3>* pvRegionPositions = nullptr,
						std::vector<int>* pvRegionAttributes = nullptr,
						std::vector<float>* pvRegionVolumeConstraints = nullptr);

////////////////////////////////////////////////////////////////////////
//	ExportGridToSMESH
///	writes an SMESH file.
/**
 * Be sure that grid consists of closed, 2-dimensional surfaces.
 * except of grid, filename and aPos, all parameters are optional.
 * Be sure that all attachments that you pass to the method have been properly attached.
 * If you don't want to specify a parameter simply pass nullptr.
 * meaning of optional parameters:
 * pvNodeAttributes: a vector that holds attachments of type AFloat,
					 which store vertex-attributes.
 * psh: subset-handler specifying boundary-indices of vertices and faces.
 * pvHoles: specifies where holes are located.
 * pvRegionPositions: specifies where regions are located.
 * vRegionAttributes: specifies the attribute of a region.
						Has to be of the same size as vRegionPositions.
						Will only be used if pvRegionAttributes has been specified.
 * vRegionVolumeConstraints: specifies the maximal volume of a tetrahedron in each region.
							 Has to be of the same size as vRegionPositions.
							 Will only be used if pvRegionAttributes has been specified.
 */
bool ExportGridToSMESH(Grid& grid, const char* filename, AVector3& aPos,
						ISubsetHandler* psh = nullptr,
						std::vector<AFloat>* pvNodeAttributes = nullptr,
						std::vector<vector3>* pvHoles = nullptr,
						std::vector<vector3>* pvRegionPositions = nullptr,
						std::vector<int>* pvRegionAttributes = nullptr,
						std::vector<float>* pvRegionVolumeConstraints = nullptr);

bool LoadGridFromSMESH(Grid& grid, const char* filename, AVector3& aPos,
						ISubsetHandler* psh = nullptr);

bool ExportGridToTETGEN(Grid& grid, const char* filename,
						AVector3& aPos, std::vector<AFloat>* pvNodeAttributes = nullptr,
						AInt* paNodeBoundaryMarker = nullptr,
						AInt* paFaceBoundaryMarker = nullptr,
						AInt* paElementAttribute = nullptr,
						ANumber* paVolumeConstraint = nullptr);

}//	end of namespace

#endif
