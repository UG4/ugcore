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

#ifndef __H__PCL__PCL_BASE__
#define __H__PCL__PCL_BASE__

namespace pcl
{

/// \addtogroup pcl
/// \{

/// call this method before 'Init' to avoid a call to MPI_Init.
/** This method may be useful if you use this program together with another
 * program which calls MPI_Init.
 * \note This will also stop MPI_Finalize from being called.*/
void DisableMPIInit ();

///	call this method before any other pcl-operations.
void Init(int *argcp, char ***argvp);

///	call this method to abort all mpi processes
void Abort(int errorcode = 1);

///	call this method right before quitting your application
void Finalize();

///	returns the number of processes
int NumProcs();

///	returns the rank of the process
int ProcRank();

/// sets error handler
void SetErrHandler();

// end group pcl
/// \}

}//	end of namespace
#endif
