/*
 * pcl_errhandler.cpp
 *
 *  Created on: Jun 14, 2012
 *      Author: Martin Rupp
 */




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "mpi.h"
#include "common/log.h"
#include "common/assert.h"

using namespace ug;

namespace pcl{

/**
 * this function is a MPI Error Handler which displays the error
 * using UG_LOG and UG_ASSERT
 */
void MPIErrorHandler( MPI_Comm *comm, int *err, ... )
{
	UG_LOG("MPI ERROR!\n");
	char error_string[256];
	int length_of_error_string=256, error_class;

	MPI_Error_class(*err, &error_class);
	MPI_Error_string(error_class, error_string, &length_of_error_string);
	UG_LOG("MPI Error: " << error_string << "\n");
	GetLogAssistant().flush();
	UG_ASSERT(0, "MPI Error: " << error_string << "\n" );
	return;
}

/// setting the mpi error handler for ug
void SetMPIErrHandler()
{
	MPI_Errhandler newerr;
#if MPI_VERSION > 1
	MPI_Comm_create_errhandler( MPIErrorHandler, &newerr );
	MPI_Comm_set_errhandler( MPI_COMM_WORLD, newerr );
#else
	MPI_Errhandler_create( MPIErrorHandler, &newerr );
	MPI_Errhandler_set( MPI_COMM_WORLD, newerr );
#endif

}

} //namespace ug
