//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d4

#include "mpi.h"
#include "pcl_base.h"
#include "common/log.h"
#include "pcl_profiling.h"
#include "pcl_errhandler.h"

namespace pcl
{

////////////////////////////////////////////////////////////////////////
void Init(int *argcp, char ***argvp)
{
	PCL_PROFILE(MPI_Init);
	//	init mpi
	MPI_Init(argcp, argvp);
	SetMPIErrHandler();
}

////////////////////////////////////////////////////////////////////////
void Finalize()
{
	PCL_PROFILE(pclFinalize);
	MPI_Finalize();
}

////////////////////////////////////////////////////////////////////////
int GetProcRank()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

int GetNumProcesses()
{
	int numProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
	return numProcesses;
}

}//	end of namespace
