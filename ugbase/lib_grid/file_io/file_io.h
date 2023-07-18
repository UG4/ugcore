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

#ifndef __H__LIB_GRID__FILE_IO__
#define __H__LIB_GRID__FILE_IO__

#include "common/types.h"
#include "common/ug_config.h"
#include "common/error.h"
#include <lib_grid/attachments/attachment_pipe.h>
#include <lib_grid/grid/grid_base_objects.h>
#include "lib_grid/refinement/projectors/projection_handler.h"

namespace ug
{

//	predeclaration of surface view (declared in lib_grid/tools/surface_view.h)
class Grid;
class MultiGrid;
class SurfaceView;
class ISubsetHandler;

////////////////////////////////////////////////////////////////////////////////
///	Loads a grid from a file. Position data is written to the specified attachment.
/**
 * Make sure that the given position attachment is either of type AVector1,
 * AVector2 or AVector3.
 *
 * If the given file can't be found, LoadGridFromFile will looks for it reative
 * to the following additional places:
 * 	- PathProvider::get_current_path()
 * 	- PathProvider::get_path(GRID_PATH)
 *
 * The method posts the following messages at the message hub of the given grid:
 * 	- GridMessage_Creation(GMCT_CREATION_STARTS)
 * 	- GridMessage_Creation(GMCT_CREATION_STOPS)
 *
 * procId can be used to only load a grid on one process in a parallel environment.
 * By default, procId is set to -1, which means that the domain is loaded on all
 * processes.
 * \{
 */
template <class TAPos>
UG_API
bool LoadGridFromFile(Grid& grid, SPProjectionHandler& ph, size_t& num_ph, ISubsetHandler& sh, std::vector<std::string> additionalSHNames,
						std::vector<SmartPtr<ISubsetHandler>> ash, const char* filename, TAPos& aPos, int procId = -1);

template <class TAPos>
UG_API
bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh, const char* filename, TAPos& aPos, int procId = -1);

template <class TAPos>
UG_API
bool LoadGridFromFile(Grid& grid, const char* filename, TAPos& aPos, int procId = -1);
/**	\} */

////////////////////////////////////////////////////////////////////////////////
///	Loads a grid from a file. Position data is written to aPosition.
/**
 * If the given file can't be found, LoadGridFromFile will looks for it reative
 * to the following additional places:
 * 	- PathProvider::get_current_path()
 * 	- PathProvider::get_path(GRID_PATH)
 *
 * The method posts the following messages at the message hub of the given grid:
 * 	- GridMessage_Creation(GMCT_CREATION_STARTS)
 * 	- GridMessage_Creation(GMCT_CREATION_STOPS)
 *
 * procId can be used to only load a grid on one process in a parallel environment.
 * By default, procId is set to -1, which means that the domain is loaded on all
 * processes.
 * \{ */
UG_API
bool LoadGridFromFile(Grid& grid, ISubsetHandler& sh, const char* filename, int procId = -1);

UG_API
bool LoadGridFromFile(Grid& grid, const char* filename, int procId = -1);
/** \} */


////////////////////////////////////////////////////////////////////////////////
///	Saves a grid to a file. Position data is read from the specified attachment.
/**
 * Make sure that the given position attachment is either of type AVector1,
 * AVector2 or AVector3.
 *
 * \{
 */
template <class TAPos>
UG_API
bool SaveGridToFile(Grid& grid, ISubsetHandler& sh,
					const char* filename, TAPos& aPos);

template <class TAPos>
UG_API
bool SaveGridToFile(Grid& grid, const char* filename, TAPos& aPos);
/**	\} */

////////////////////////////////////////////////////////////////////////////////
///	Saves a grid to a file. Position data is read from aPosition.
/** \{ */
UG_API
bool SaveGridToFile(Grid& grid, ISubsetHandler& sh, const char* filename);

UG_API
bool SaveGridToFile(Grid& grid, const char* filename);
/** \} */


///	Saves a grid hierarchy by offsetting levels along the z-axis.
/**	\todo:	Support any type of position attachment*/
bool SaveGridHierarchyTransformed(MultiGrid& mg, ISubsetHandler& sh,
								  const char* filename, number offset);

///	Saves a grid hierarchy by offsetting levels along the z-axis.
/**	\todo:	Support any type of position attachment*/
bool SaveGridHierarchyTransformed(MultiGrid& mg, const char* filename,
									  number offset);

///	Saves the grid-layout of parallel multi-grids
/**	\todo:	Support any type of position attachment*/
bool SaveParallelGridLayout(MultiGrid& mg, const char* filename,
							number offset = 0.1);


///	Saves a grid and assigns elements to subsets based on their surface-view-state.
/**	\todo:	Support any type of position attachment*/
bool SaveSurfaceViewTransformed(MultiGrid& mg, const SurfaceView& sv,
								const char* filename, number offset = 0.1);


////////////////////////////////////////////////////////////////////////////////
///	Copies a grid level to a separate grid.
/**
 * Extracts a specified level from a given MultiGrid with a given SubsetHandler
 * and copies it to a new grid. Position data is read from aPos.
 */
template <class TAPos>
void CopyGridLevel(MultiGrid& srcMG, Grid& destGrid,
				   ISubsetHandler& srcSH, ISubsetHandler& destSH,
				   int lvl, TAPos aPos);


////////////////////////////////////////////////////////////////////////////////
///	Copies a grid to a separate grid.
/**
 * Extracts the grid from a given Grid with a given SubsetHandler
 * and copies it to a new grid. Position data is read from aPos.
 * Note: To copy a grid w/o subset information can use Grid's copy constructor
 * \param[in] srcGrid the grid to be copied
 * \param[out] destGrid the copied grid
 * \param[in] srcSh the SubsetHandler for the srcGrid
 * \param[out] destSh the SubsetHandler for the destGrid
 * \tparam[in] aPos position attachment
 */
template <class TAPos>
void CopyGrid
(
	Grid& srcGrid,
	Grid& destGrid,
	ISubsetHandler& srcSH,
	ISubsetHandler& destSH,
	TAPos aPos
);

/// Merges to grids into one grid.
/**
 * Merges all grid elements and subset information from grid to mrgGrid
 * \param[in,out] mrgGrid merged grid
 * \param[in] grid input to be merged into merged grid
 * \param[in] mrgSH the SubsetHandler for mrgGrid
 * \param[in,out] sh the SubsetHandler for grid
 * \tparam[in] aPos position attachment
 * \param[in] joinSubsets if true then subsets will be joined otherwise not
 */
template <typename TAPos>
void MergeGrids
(
	Grid& mrgGrid,
	Grid& grid,
	ISubsetHandler& mrgSh,
	ISubsetHandler& sh,
	TAPos aPos,
	bool joinSubsets
);

///	Saves a grid level to a file.
/**
 * Extracts a specified level from a given MultiGrid with a given SubsetHandler
 * and saves it to a file. Position data is read from aPos.
 */
template <class TAPos>
bool SaveGridLevel(MultiGrid& srcMG, ISubsetHandler& srcSH, int lvl,
				   const char* filename, TAPos aPos);


///	Saves a grid level to a file.
/**
 * Extracts a specified level from a given MultiGrid with a given SubsetHandler
 * and saves it to a file. Position data is read from MultiGrid.
 */
bool SaveGridLevelToFile(MultiGrid& srcMG, ISubsetHandler& srcSH, int lvl,
		  	  	  	     const char* filename);

/// Copy grid elements of type TElem from srcGrid to destGrid
template<class TElem>
void CopyGridElements(Grid& srcGrid, Grid& destGrid,
				      ISubsetHandler& srcSH, ISubsetHandler& destSH,
					  Attachment<Vertex*>& aNewVrt);
/// end namespace ug
}

/// include implementation
#include "file_io_impl.h"

#endif
