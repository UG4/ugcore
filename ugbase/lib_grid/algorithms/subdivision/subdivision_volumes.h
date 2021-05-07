/*
 * Copyright (c) 2014-2021:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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

#ifndef __H__UG__SUBDIVISION_VOLUMES__
#define __H__UG__SUBDIVISION_VOLUMES__

#include <vector>
#include <cassert>
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/grid_objects/tetrahedron_rules.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/algorithms/subdivision/subdivision_loop.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/util/compol_attachment_reduce.h"
	#include "lib_grid/parallelization/util/compol_copy_attachment.h"
	#include "lib_grid/parallelization/util/attachment_operations.hpp"
	#include "lib_grid/parallelization/distributed_grid.h"
	#include "pcl/pcl_interface_communicator.h"
#endif

#include  "common/profiler/profiler.h"

namespace ug
{


////////////////////////////////////////////////////////////////////////////////
//	Global boundary refinement rule setup
////////////////////////////////////////////////////////////////////////////////

/// enumeration for identification of global boundary refinement rule to be used
enum GlobalBoundaryRefinementRule
{
	LINEAR, // the default value
	SUBDIV_SURF_LOOP_SCHEME,
	SUBDIV_SURF_AVERAGING_SCHEME,
	SUBDIV_SURF_BUTTERFLY_SCHEME,
	SUBDIV_VOL
};

/// setting procedure for global boundary refinement rule variable
void SetBoundaryRefinementRule(GlobalBoundaryRefinementRule refRule);

///	get procedure for boundary refinement rule variable
GlobalBoundaryRefinementRule GetBoundaryRefinementRule();


/// Function for checking the number of associated volumes of all edges
/** This function calculates the number of associated volumes
 * 	for all edges.
 *
 * 	@param mg			reference to MultiGrid
 * 	@param sh			reference to SubsetHandler
 * 	@param sh			filename specification for output
**/
void CheckValences(MultiGrid& mg, MGSubsetHandler& markSH, const char* filename);


/// Function for printing the subdivision volumes tet-oct refinement mask
/** This function prints the subdivision volumes
 *  tet-oct refinement mask.
**/
void PrintSubdivisionVolumesRefinementMask();


/// Function for splitting an octahedron to 4 sub-tetrahedrons
/** Recall the refinement of a tetrahedron (s. tetrahdron_rules.cpp). A tetrahedron is
 *	refined into 4 outer tetrahedrons and 4 inner tetrahedrons. After the 4 outer
 *	tetrahedrons are created the remaining inner cavity corresponds to an octahedron.
 *	This octahedron can be split into 4 tetrahedrons in 3 different ways, depending
 *	on the length of the following diagonals:
 *	Based on the original tetrahedron we look at the three diagonals between the
 *	following edge-centers: 0-5, 1-3, 2-4
 *
 *	The diagonal between edge-centers 0-5 of the tetrahedron equals
 *	a segment between vertices 1 and 3 of the octahedron
 *
 *	The diagonal between edge-centers 1-3 of the tetrahedron equals
 *	a segment between vertices 0 and 5 of the octahedron
 *
 *	the diagonal between edge-centers 2-4 of the tetrahedron equals
 *	a segment between vertices 2 and 4 of the octahedron
 *
 *	HINT: preferably use bestDiag = 0, as it is the inherent diagonal
 *		  along which the octahedron was adaptively orientated according to
 *		  tetrahedron_rules.cpp
 *
 *
 *	@param grid			reference to grid
 * 	@param oct			pointer to octahedron
 * 	@param parentVol	pointer to parent volume
 * 	@param vTetsOut		reference to vector with pointers to new tetrahedrons

**/
void SplitOctahedronToTetrahedrons(	Grid& grid, Octahedron* oct, Volume* parentVol,
									std::vector<Tetrahedron*>& vTetsOut, int bestDiag);


/// Conversion function for hybrid tetra-/octahedral multigrids
/** This function converts each octahedron in all levels to
 * 	four tetrahedra and deletes the original octahedra
 * 	from the multigrid.
 *
 * 	WARNING: correct parent <-> childhood relationships won't persist
 *
 * 	@param mg			reference to MultiGrid
 * 	@param bestDiag		specify fixed or adaptive diagonal for octahedral split
**/
void TetrahedralizeHybridTetOctGrid(MultiGrid& mg, int bestDiag);


/// Projection function for smooth subdivision (volumes+surface) refinement
/** This function projects the vertices of all levels to their smooth limit
 * 	positions determined by the subdivision volumes refinement.
 *
 * 	@param mg			reference to MultiGrid
 * 	@param aPos			reference to position attachment
**/
template <class TAPosition>
void ProjectHierarchyToSubdivisionLimit(MultiGrid& mg, TAPosition& aPos);


/// Parent level vertex smoothing function for subdivision crease surfaces refinement (by C. Loop, 1987)
/** This function calculates the smoothed positions of all parent level vertices
 * 	determined by the subdivision surfaces refinement for surfaces that are not closed
 * 	handling the boundary polylines separately
 *
 * 	@param mg						reference to MultiGrid
 * 	@param aPos						reference to position attachment
 * 	@param markSH					reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param linearManifoldSH			reference to user-specified linearManifoldSubsets SubsetHandler
 * 	@param aSmoothBndPosEvenVrt		reference to aSmoothBndPosEvenVrt
 * 	@param aSmoothBndPosOddVrt		reference to aSmoothBndPosOddVrt
 * 	@param aNumManifoldEdges		reference to aNumManifoldEdges
**/
template <class TAPosition>
void CalculateSmoothCreaseManifoldPosInParentLevelLoopScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
											 	 	   	   	 MGSubsetHandler& linearManifoldSH,
															 TAPosition& aSmoothBndPosEvenVrt,
															 TAPosition& aSmoothBndPosOddVrt,
															 AInt& aNumManifoldEdges);


/// Parent level vertex smoothing function for subdivision surfaces refinement (by C. Loop, 1987)
/** This function calculates the smoothed positions of all parent level vertices
 * 	determined by the subdivision surfaces refinement.
 *
 * 	@param mg						reference to MultiGrid
 * 	@param aPos						reference to position attachment
 * 	@param markSH					reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param linearManifoldSH			reference to user-specified linearManifoldSubsets SubsetHandler
 * 	@param aSmoothBndPosEvenVrt		reference to aSmoothBndPosEvenVrt
 * 	@param aSmoothBndPosOddVrt		reference to aSmoothBndPosOddVrt
 * 	@param aNumManifoldEdges		reference to aNumManifoldEdges
**/
template <class TAPosition>
void CalculateSmoothManifoldPosInParentLevelLoopScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
											 	 	   MGSubsetHandler& linearManifoldSH,
													   TAPosition& aSmoothBndPosEvenVrt,
													   TAPosition& aSmoothBndPosOddVrt,
													   AInt& aNumManifoldEdges);


/// Toplevel vertex smoothing function for subdivision surfaces refinement (Averaging scheme)
/** This function calculates the smoothed positions of all toplevel vertices
 * 	determined by the subdivision surfaces refinement.
 *
 * 	@param mg					reference to MultiGrid
 * 	@param aPos					reference to position attachment
 * 	@param markSH				reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param linearManifoldSH		reference to user-specified linearManifoldSubsets SubsetHandler
 * 	@param aSmoothBndPos_tri	reference to aSmoothBndPos_tri
 *	@param aSmoothBndPos_quad	reference to aSmoothBndPos_quad
**/
template <class TAPosition>
void CalculateSmoothManifoldPosInTopLevelAveragingScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
														 MGSubsetHandler& linearManifoldSH,
														 TAPosition& aSmoothBndPos_tri,
														 TAPosition& aSmoothBndPos_quad);


/// Parent level vertex smoothing function for subdivision surfaces refinement (Butterfly scheme)
/** This function calculates the smoothed positions of all parent level vertices
 * 	determined by the subdivision surfaces refinement.
 *
 * 	@param mg						reference to MultiGrid
 * 	@param aPos						reference to position attachment
 * 	@param markSH					reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param linearManifoldSH			reference to user-specified linearManifoldSubsets SubsetHandler
 * 	@param aSmoothBndPosOddVrt		reference to aSmoothBndPosOddVrt
 * 	@param aNumManifoldEdges		reference to aNumManifoldEdges
**/
template <class TAPosition>
void CalculateSmoothManifoldPosInParentLevelButterflyScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
											 	 	   MGSubsetHandler& linearManifoldSH,
													   TAPosition& aSmoothBndPosOddVrt,
													   AInt& aNumManifoldEdges);


/// Toplevel vertex smoothing function for subdivision volumes refinement
/** This function calculates the smoothed positions of all toplevel vertices
 * 	determined by the subdivision volumes refinement.
 *
 * 	@param mg					reference to MultiGrid
 * 	@param markSH				reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param aSmoothVolPos_toc	reference to aSmoothVolPos_toc
 * 	@param aSmoothVolPos_prism	reference to aSmoothVolPos_prism
 * 	@param aSmoothVolPos_hex	reference to aSmoothVolPos_hex
**/
void CalculateSmoothVolumePosInTopLevel(MultiGrid& mg, MGSubsetHandler& markSH,
										APosition& aSmoothVolPos_toc,
										APosition& aSmoothVolPos_prism,
										APosition& aSmoothVolPos_hex);


/// Toplevel vertex smoothing function for subdivision volumes refinement
/** This function calculates the smoothed positions of all toplevel vertices
 * 	determined by the constrained subdivision volumes refinement.
 *
 * 	@param mg					reference to MultiGrid
 * 	@param markSH				reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param aSmoothVolPos_toc	reference to aSmoothVolPos_toc
**/
void CalculateConstrainedSmoothVolumePosInTopLevel(MultiGrid& mg, MGSubsetHandler& markSH,
												   APosition& aSmoothVolPos_toc);


/// Function for calculating the number of associated volumes of all toplevel vertices
/** This function calculates the number of associated volumes
 * 	for all toplevel vertices.
 *
 * 	@param mg				reference to MultiGrid
 * 	@param aNumElems_toc	reference to aNumElems_toc
 * 	@param aNumElems_prism	reference to aNumElems_prism
 * 	@param aNumElems_hex	reference to aNumElems_hex
**/
void CalculateNumElemsVertexAttachmentInTopLevel(MultiGrid& mg, AInt& aNumElems_toc, AInt& aNumElems_prism, AInt& aNumElems_hex);


/// Function for calculating the number of associated manifold edges of all parent level vertices
/** This function calculates the number of associated manifold edges
 * 	for all parent level vertices.
 *
 * 	@param mg					reference to MultiGrid
 * 	@param markSH				reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param aNumManifoldEdges	reference to aNumManifoldEdges
**/
void CalculateNumManifoldEdgesVertexAttachmentInParentLevel(MultiGrid& mg, MGSubsetHandler& markSH,
															AInt& aNumManifoldEdges);


/// Function for calculating the number of associated manifold faces (triangles) of all parent level manifold edges
/** This function calculates the number of associated manifold faces (triangles)
 * 	for all parent level edges.
 *
 * 	@param mg					reference to MultiGrid
 * 	@param markSH				reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param aNumManifoldFaces	reference to aNumManifoldFaces
**/
void CalculateNumManifoldFacesEdgeAttachmentInParentLevel(MultiGrid& mg, MGSubsetHandler& markSH, AInt& aNumManifoldFaces);


/// Function for calculating the number of associated manifold faces of all toplevel manifold vertices
/** This function calculates the number of of associated manifold faces (triangles and quadrilaterals)
 * 	for all toplevel vertices.
 *
 * 	@param mg						reference to MultiGrid
 * 	@param markSH					reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param aNumManifoldFaces_tri	reference to aNumManifoldFaces_tri
 * 	@param aNumManifoldFaces_quad	reference to aNumManifoldFaces_quad
**/
void CalculateNumManifoldFacesVertexAttachmentInTopLevel(MultiGrid& mg, MGSubsetHandler& markSH, AInt& aNumManifoldFaces_tri, AInt& aNumManifoldFaces_quad);


/// Procedure to initialize the crease boundary feature vertex attachment of all parent level vertices
/** This procedure initializes crease boundary feature vertex attachment
 * 	for all parent level vertices.
 *
 * 	@param mg					reference to MultiGrid
 * 	@param aNumManifoldFaces	reference to aNumManifoldFaces
 * 	@param aBoolIsBndVrt		reference to aBoolIsCreaseBndVrt
**/
void InitIsCreaseBoundaryVertexAttachmentInParentLevel(MultiGrid& mg, AInt& aNumManifoldFaces, ABool& aBoolIsCreaseBndVrt);


/// Procedure to initialize the linear boundary manifold subsets SubsetHandler with user-specified subsets
/** This procedure initializes the referenced linear boundary manifold subsets SubsetHandler
 * 	s.t. user-specified subsets
 *
 * 	@param dom						reference to Domain
 * 	@param linearManifoldSH			reference to user-specified linearManifoldSubsets SubsetHandler
**/
void InitLinearManifoldSubsetHandler(MultiGrid& mg, MGSubsetHandler& sh,
											   MGSubsetHandler& linearManifoldSH,
											   const char* linearManifoldSubsets);


/// Toplevel vertex repositioning function for subdivision surfaces refinement (by C. Loop, 1987)
/** This function repositions all toplevel manifold vertices to their smoothed positions
 * 	determined by the subdivision surfaces refinement.
 *
 * 	@param mg						reference to MultiGrid
 * 	@param aPos						reference to position attachment
 * 	@param markSH					reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param linearManifoldSH			reference to user-specified linearManifoldSubsets SubsetHandler
 * 	@param bClosedSurface 			bool switch for surfaces that are not closed in case of subdivision
 * 									refinement by C. Loop, 1987 handling the boundary polylines separately
**/
template <class TAPosition>
void ApplySmoothManifoldPosToTopLevelLoopScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
												MGSubsetHandler& linearManifoldSH, bool bClosedSurface);


/// Toplevel vertex repositioning function for subdivision surfaces refinement (Butterfly scheme)
/** This function repositions all toplevel manifold vertices to their smoothed positions
 * 	determined by the subdivision surfaces refinement.
 *
 * 	@param mg						reference to MultiGrid
 * 	@param aPos						reference to position attachment
 * 	@param markSH					reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param linearManifoldSH			reference to user-specified linearManifoldSubsets SubsetHandler
**/
template <class TAPosition>
void ApplySmoothManifoldPosToTopLevelButterflyScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
												     MGSubsetHandler& linearManifoldSH);


/// Toplevel vertex repositioning function for subdivision surfaces refinement (Averaging scheme)
/** This function repositions all toplevel manifold vertices to their smoothed positions
 * 	determined by the subdivision surfaces refinement.
 *
 * 	@param mg					reference to MultiGrid
 * 	@param aPos					reference to position attachment
 * 	@param markSH				reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param linearManifoldSH		reference to user-specified linearManifoldSubsets SubsetHandler
**/
template <class TAPosition>
void ApplySmoothManifoldPosToTopLevelAveragingScheme(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& markSH,
													 MGSubsetHandler& linearManifoldSH);


/// Toplevel vertex repositioning function for subdivision volumes refinement
/** This function repositions all toplevel inner vertices to their smoothed positions
 * 	determined by the subdivision volumes refinement.
 *
 * 	@param mg					reference to MultiGrid
 * 	@param markSH				reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param bConstrained			bool switch for constrained smooth subdivision volumes scheme
**/
void ApplySmoothVolumePosToTopLevel(MultiGrid& mg, MGSubsetHandler& markSH,
									bool bConstrained);


/// Function to create a smooth subdivision volumes hierarchy
/** This function transforms a linearly refined hybrid tetra-/octahedral volume
 * 	grid hierarchy into a hierarchy with smoothed boundary manifold
 * 	(s. Schaefer et al, "Smooth subdivision of tetrahedral meshes", 2004)
 *
 * 	@param mg						reference to MultiGrid
 * 	@param aPos						reference to position attachment
 * 	@param sh						reference to standard SubsetHandler
 * 	@param markSH					reference to SubsetHandler containing marked (inner) boundary manifold
 * 	@param linearManifoldSH 		reference to user-specified linearManifoldSubsets SubsetHandler
 * 	@param bClosedSurface 			bool switch for surfaces that are not closed in case of subdivision
 * 									refinement by C. Loop, 1987 handling the boundary polylines separately
**/
template <class TAPosition>
void ApplySmoothSubdivisionSurfacesToTopLevel(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& sh,
											  MGSubsetHandler& markSH, MGSubsetHandler& linearManifoldSH,
											  bool bClosedSurface);


/// Function to create a smooth subdivision volumes hierarchy
/** This function transforms a linearly refined hybrid tetra-/octahedral volume
 * 	grid hierarchy into a smoothed subdivision volumes hierarchy
 * 	(s. Schaefer et al, "Smooth subdivision of tetrahedral meshes", 2004)
 *
 * 	@param mg						reference to MultiGrid
 * 	@param sh						reference to standard SubsetHandler
 * 	@param markSH					reference to SubsetHandler containing marked (inner) boundary manifold
 * 	@param linearManifoldSH 		reference to user-specified linearManifoldSubsets SubsetHandler
 * 	@param bConstrained				bool switch for constrained smooth subdivision volumes scheme
**/
void ApplySmoothSubdivisionVolumesToTopLevel(MultiGrid& mg, MGSubsetHandler& sh, MGSubsetHandler& markSH,
											 MGSubsetHandler& linearManifoldSH, bool bConstrained);


/// Wrapper procedures for smooth subdivision surfaces/volumes hierarchy creation
/** These functions call the actual ApplySmoothSubdivisionVolumesToTopLevel procedure
 *
 * 	@param mg						reference to MultiGrid
 * 	@param aPos						reference to position attachment
 * 	@param markSH					reference to SubsetHandler markSH containing marked (inner) boundary manifold
 * 	@param linearManifoldSubsets 	user-specified linearManifoldSubsets
 * 	@param bConstrained				bool switch for constrained smooth subdivision volumes scheme
 * 	@param bClosedSurface 			bool switch for surfaces that are not closed in case of subdivision
 * 									refinement by C. Loop, 1987 handling the boundary polylines separately
**/
template <class TAPosition>
void ApplySmoothSubdivisionSurfacesToTopLevel(MultiGrid& mg, TAPosition& aPos, MGSubsetHandler& sh,
											  MGSubsetHandler& markSH, const char* linearManifoldSubsets,
											  bool bClosedSurface);

void ApplySmoothSubdivisionVolumesToTopLevel(MultiGrid& mg, MGSubsetHandler& sh, MGSubsetHandler& markSH,
											 const char* linearManifoldSubsets);


void ApplyConstrainedSmoothSubdivisionVolumesToTopLevel(MultiGrid& mg, MGSubsetHandler& sh, MGSubsetHandler& markSH,
														const char* linearManifoldSubsets);


}//	end of namespace

#endif
