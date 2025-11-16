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




#include <mpi.h>

#include "common/log.h"

#include "pcl_base.h"
#include "pcl_profiling.h"
#include "pcl_comm_world.h"


namespace pcl
{

static bool PERFORM_MPI_INITIALIZATION = true;

////////////////////////////////////////////////////////////////////////
void DisableMPIInit ()
{
	PERFORM_MPI_INITIALIZATION = false;
}

////////////////////////////////////////////////////////////////////////
void Init(int *argcp, char ***argvp)
{
	PCL_PROFILE(MPI_Init);

	if(PERFORM_MPI_INITIALIZATION){
		//	init mpi /* was: MPI_Init(argcp, argvp);*/
		int flag;
		MPI_Initialized(&flag);
		if (!flag) MPI_Init(argcp, argvp);
		SetErrHandler();
	}

	// default to MPI_COMM_WORLD
	if  (PCL_COMM_WORLD == MPI_COMM_NULL) {PCL_COMM_WORLD= MPI_COMM_WORLD;}
}

////////////////////////////////////////////////////////////////////////
void Abort(int errorcode)
{
	PCL_PROFILE(pclAbort);
	MPI_Abort(PCL_COMM_WORLD, 1);
}

////////////////////////////////////////////////////////////////////////
void Finalize()
{
	PCL_PROFILE(pclFinalize);
	if(PERFORM_MPI_INITIALIZATION)
		MPI_Finalize();
}

////////////////////////////////////////////////////////////////////////
int ProcRank()
{
	int rank;
	MPI_Comm_rank(PCL_COMM_WORLD, &rank);
	return rank;
}

////////////////////////////////////////////////////////////////////////
int NumProcs()
{
	int numProcesses;
	MPI_Comm_size(PCL_COMM_WORLD, &numProcesses);
	return numProcesses;
}

////////////////////////////////////////////////////////////////////////

/*
 * this function is a MPI Error Handler which displays the error
 * using UG_LOG and UG_ASSERT
 */
void MPIErrorHandler( MPI_Comm *comm, int *err, ... )
{
	char error_string[256];
	int length_of_error_string=256, error_class;

	MPI_Error_class(*err, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);
	UG_LOG_ALL_PROCS("MPI ERROR: " << error_string << std::endl);
	ug::GetLogAssistant().flush();
	Abort();
}

void SetErrHandler()
{
	MPI_Errhandler newerr;
	MPI_Comm_create_errhandler( MPIErrorHandler, &newerr );
	MPI_Comm_set_errhandler( PCL_COMM_WORLD, newerr );
}

}//	end of namespace
