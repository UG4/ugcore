//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d02

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
