// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 15.02.2011 (m,d,y)
 
#include <iostream>
#include <string>
#include "../registry.h"
#include "../ug_bridge.h"
#include "pcl/pcl.h"

namespace ug{
namespace bridge{

#ifdef UG_PARALLEL

template<typename T>
T ParallelMin(T t)
{
	pcl::ProcessCommunicator pc;
	return pc.allreduce(t, PCL_RO_MIN);
}

template<typename T>
T ParallelMax(T t)
{
	pcl::ProcessCommunicator pc;
	return pc.allreduce(t, PCL_RO_MAX);
}

template<typename T>
T ParallelSum(T t)
{
	pcl::ProcessCommunicator pc;
	return pc.allreduce(t, PCL_RO_SUM);
}

bool RegisterPCLInterface(Registry& reg, const char* parentGroup)
{
	std::string grpStr(parentGroup);
	grpStr.append("/pcl");

	reg.add_function("GetNumProcesses", &pcl::GetNumProcesses, grpStr.c_str(),
					"NumProcs", "", "Returns the number of active processes.");

	reg.add_function("GetProcessRank", &pcl::GetProcRank, grpStr.c_str(),
					"ProcRank", "", "Returns the rank of the current process.");

	reg.add_function("GetOutputProcessRank", &pcl::GetOutputProcRank, grpStr.c_str(),
					"ProcRank", "", "Returns the rank of the process whose output is logged.");

	reg.add_function("SetOutputProcessRank", &pcl::SetOutputProcRank, grpStr.c_str(),
					"", "ProcRank", "Sets the rank of the process whose output is logged.");

	reg.add_function("IsOutputProcess", &pcl::IsOutputProc, grpStr.c_str(),
					"", "", "Returns true if the current process is the output process.");

	reg.add_function("SynchronizeProcesses", &pcl::SynchronizeProcesses, grpStr.c_str(),
					"", "", "Waits until all active processes reached this point.");

	reg.add_function("ParallelMin", &ParallelMin<double>, grpStr.c_str(), "tmax", "t", "returns the maximum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelMax", &ParallelMax<double>, grpStr.c_str(), "tmin", "t", "returns the minimum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelSum", &ParallelSum<double>, grpStr.c_str(), "tsum", "t", "returns the sum of t over all processes. note: you have to assure that all processes call this function.");

	return true;
}

#else // UG_PARALLEL

///	Dummy method for serial compilation always returning 1
static int GetNumProcessesDUMMY()	{return 1;}

///	Dummy method for serial compilation always returning 0
static int GetProcRankDUMMY()				{return 0;}

///	Dummy method for serial compilation always returning 0
static int GetOutputProcRankDUMMY()		{return 0;}

///	Dummy method for serial compilation doing nothing
static void SetOutputProcRankDUMMY(int)			{}

///	Dummy method for serial compilation always returning true
static bool IsOutputProcDUMMY()			{return true;}

///	Dummy method for serial compilation doing nothing
static void SynchronizeProcessesDUMMY()			{}


template<typename T>
T ParallelMinDUMMY(T t)
{
	return t;
}

template<typename T>
T ParallelMaxDUMMY(T t)
{
	return t;
}

template<typename T>
T ParallelSumDUMMY(T t)
{
	return t;
}


bool RegisterPCLInterface(Registry& reg, const char* parentGroup)
{
	std::string grpStr(parentGroup);
	grpStr.append("/pcl");

	reg.add_function("GetNumProcesses", &GetNumProcessesDUMMY, grpStr.c_str(),
					"NumProcs", "", "Returns the number of active processes.");

	reg.add_function("GetProcessRank", &GetProcRankDUMMY, grpStr.c_str(),
					"ProcRank", "", "Returns the rank of the current process.");

	reg.add_function("GetOutputProcessRank", &GetOutputProcRankDUMMY, grpStr.c_str(),
					"ProcRank", "", "Returns the rank of the process whose output is logged.");

	reg.add_function("SetOutputProcessRank", &SetOutputProcRankDUMMY, grpStr.c_str(),
					"", "ProcRank", "Sets the rank of the process whose output is logged.");

	reg.add_function("IsOutputProcess", &IsOutputProcDUMMY, grpStr.c_str(),
					"", "", "Returns true if the current process is the output process.");

	reg.add_function("SynchronizeProcesses", &SynchronizeProcessesDUMMY, grpStr.c_str(),
					"", "", "Waits until all active processes reached this point.");


	reg.add_function("ParallelMinDUMMY", &ParallelMin<double>, grpStr.c_str(), "tmax", "t", "returns the maximum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelMaxDUMMY", &ParallelMax<double>, grpStr.c_str(), "tmin", "t", "returns the minimum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelSumDUMMY", &ParallelSum<double>, grpStr.c_str(), "tsum", "t", "returns the sum of t over all processes. note: you have to assure that all processes call this function.");

	return true;
}

#endif //UG_PARALLEL

}// end of namespace
}// end of namespace
