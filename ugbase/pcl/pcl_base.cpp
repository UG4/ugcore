//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d4

#include "mpi.h"
#include "pcl_base.h"
#include "common/log.h"

namespace pcl
{

static int OUTPUT_PROC_RANK = 0;

////////////////////////////////////////////////////////////////////////
void Init(int argc, char* argv[])
{
//	init mpi
	MPI_Init(&argc, &argv);
}

////////////////////////////////////////////////////////////////////////
void Finalize()
{
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

////////////////////////////////////////////////////////////////////////
int GetOutputProcRank()
{
	return OUTPUT_PROC_RANK;
}

////////////////////////////////////////////////////////////////////////
void SetOutputProcRank(int rank)
{
	if(rank < -1 || rank >= GetNumProcesses()){
		UG_LOG("-- PCL-LOG: Can't change output process: invalid rank specified: " << rank << ".\n");
		return;
	}
	
	OUTPUT_PROC_RANK = rank;
}

}//	end of namespace
