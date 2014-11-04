//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d18

#ifndef __H__LIBGRID__FILE_IO_TETGEN__
#define __H__LIBGRID__FILE_IO_TETGEN__

#include <vector>
#include "../grid/grid.h"
#include "../subset_handler.h"
#include "../common_attachments.h"

namespace ug
{

bool LoadGridFromELE(Grid& grid, const char* filename, ISubsetHandler* pSH = NULL,
					APosition& aPos = aPosition);
					
bool SaveGridToELE(Grid& grid, const char* filename, ISubsetHandler* pSH = NULL,
					APosition& aPos = aPosition);
					
					
bool ImportGridFromTETGEN(Grid& grid,
						const char* nodesFilename, const char* facesFilename,
						const char* elemsFilename, AVector3& aPos,
						std::vector<AFloat>* pvNodeAttributes = NULL,
						AInt* paNodeBoundaryMarker = NULL,
						AInt* paFaceBoundaryMarker = NULL,
						AInt* paElementAttribute = NULL);

bool ImportGridFromTETGEN(Grid& grid,
						const char* nodesFilename, const char* facesFilename,
						const char* elemsFilename, AVector3& aPos,
						ISubsetHandler* psh = NULL,
						std::vector<AFloat>* pvNodeAttributes = NULL);
					
////////////////////////////////////////////////////////////////////////
//	ExportGridToSMESH
///	writes an SMESH file.
/**
 * Be sure that grid consists of closed, 2-dimensional surfaces.
 * except of grid, filename and aPos, all parameters are optional.
 * Be sure that all attachments that you pass to the method have been properly attached.
 * If you don't want to specify a parameter simply pass NULL.
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
						std::vector<AFloat>* pvNodeAttributes = NULL,
						AInt* paNodeBoundaryMarker = NULL,
						AInt* paFaceBoundaryMarker = NULL,
						std::vector<vector3>* pvHoles = NULL,
						std::vector<vector3>* pvRegionPositions = NULL,
						std::vector<int>* pvRegionAttributes = NULL,
						std::vector<float>* pvRegionVolumeConstraints = NULL);

////////////////////////////////////////////////////////////////////////
//	ExportGridToSMESH
///	writes an SMESH file.
/**
 * Be sure that grid consists of closed, 2-dimensional surfaces.
 * except of grid, filename and aPos, all parameters are optional.
 * Be sure that all attachments that you pass to the method have been properly attached.
 * If you don't want to specify a parameter simply pass NULL.
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
						ISubsetHandler* psh = NULL,
						std::vector<AFloat>* pvNodeAttributes = NULL,
						std::vector<vector3>* pvHoles = NULL,
						std::vector<vector3>* pvRegionPositions = NULL,
						std::vector<int>* pvRegionAttributes = NULL,
						std::vector<float>* pvRegionVolumeConstraints = NULL);

bool ExportGridToTETGEN(Grid& grid, const char* nodesFilename,
						const char* facesFilename, const char* elemsFilename,
						AVector3& aPos, std::vector<AFloat>* pvNodeAttributes = NULL,
						AInt* paNodeBoundaryMarker = NULL,
						AInt* paFaceBoundaryMarker = NULL,
						AInt* paElementAttribute = NULL);

}//	end of namespace

#endif
