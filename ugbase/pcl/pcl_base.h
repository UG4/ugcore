//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m06 d4

#ifndef __H__PCL__PCL_BASE__
#define __H__PCL__PCL_BASE__

namespace pcl
{

////////////////////////////////////////////////////////////////////////
///	call this method before any other pcl-operations.
//void Init(int argc, char* argv[]);
void Init(int *argcp, char ***argvp);

///	call this method right before quitting your application
void Finalize();


///	returns the number of processes
int GetNumProcesses();

///	returns the rank of the process
int GetProcRank();

///	returns the rank of the process for which the output shall be printed to the screen.
int GetOutputProcRank();

///	sets the rank of the process for which the output shall be printed.
void SetOutputProcRank(int rank);

///	returns true if the current process has to output data
inline bool IsOutputProc()
{
	int opr = GetOutputProcRank();
	return (GetProcRank() == opr) || (opr == -1);
}

}//	end of namespace
#endif
