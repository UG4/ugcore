//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d4

#include "mpi.h"
#include "pcl_base.h"
#include "pcl_profiling.h"
#include "common/log.h"

namespace pcl
{

////////////////////////////////////////////////////////////////////////
void Init(int *argcp, char ***argvp)
{
	PCL_PROFILE(MPI_Init);
	//	init mpi
	MPI_Init(argcp, argvp);
	SetErrHandler();
}

////////////////////////////////////////////////////////////////////////
void Abort(int errorcode)
{
	PCL_PROFILE(pclAbort);
	MPI_Abort(MPI_COMM_WORLD, 1);
}

////////////////////////////////////////////////////////////////////////
void Finalize()
{
	PCL_PROFILE(pclFinalize);
	MPI_Finalize();
}

////////////////////////////////////////////////////////////////////////
int ProcRank()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

////////////////////////////////////////////////////////////////////////
int NumProcs()
{
	int numProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
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
#if MPI_VERSION > 1
	MPI_Comm_create_errhandler( MPIErrorHandler, &newerr );
	MPI_Comm_set_errhandler( MPI_COMM_WORLD, newerr );
#else
	MPI_Errhandler_create( MPIErrorHandler, &newerr );
	MPI_Errhandler_set( MPI_COMM_WORLD, newerr );
#endif

}

}//	end of namespace
