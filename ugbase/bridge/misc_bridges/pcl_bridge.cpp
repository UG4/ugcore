// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 15.02.2011 (m,d,y)
 
#include <iostream>
#include <string>
#include "registry/registry.h"
#include "bridge/bridge.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

using namespace std;

namespace ug{
namespace bridge{

/// \defgroup pcl_bridge PCL Bridge
/// \ingroup misc_bridge
/// \{

#ifdef UG_PARALLEL

static bool PclDebugBarrierEnabled()
{
#ifdef PCL_DEBUG_BARRIER_ENABLED
	return true;
#else
	return false;
#endif
}

static void PclDebugBarrierAll()
{
	PCL_DEBUG_BARRIER_ALL();
}

static bool PclAllProcsTrue(bool bTrue){
	return pcl::AllProcsTrue(bTrue);
}

template<typename T>
static T ParallelMin(T t)
{
	pcl::ProcessCommunicator pc;
	return pc.allreduce(t, PCL_RO_MIN);
}

template<typename T>
static T ParallelMax(T t)
{
	pcl::ProcessCommunicator pc;
	return pc.allreduce(t, PCL_RO_MAX);
}

template<typename T>
static T ParallelSum(T t)
{
	pcl::ProcessCommunicator pc;
	return pc.allreduce(t, PCL_RO_SUM);
}

void RegisterBridge_PCL(Registry& reg, string parentGroup)
{
	string grp(parentGroup);
	grp.append("/pcl");

	reg.add_function("PclDebugBarrierEnabled", &PclDebugBarrierEnabled, grp,
					"Enabled", "", "Returns the whether debug barriers are enabled.");

	reg.add_function("PclDebugBarrierAll", &PclDebugBarrierAll, grp,
					 "", "", "Synchronizes all parallel processes if the executable"
							 "has been compiled with PCL_DEBUG_BARRIER=ON");

	reg.add_function("NumProcs", &pcl::NumProcs, grp,
					"NumProcs", "", "Returns the number of active processes.");

	reg.add_function("ProcRank", &pcl::ProcRank, grp,
					"ProcRank", "", "Returns the rank of the current process.");

	reg.add_function("SynchronizeProcesses", &pcl::SynchronizeProcesses, grp,
					"", "", "Waits until all active processes reached this point.");

	reg.add_function("AllProcsTrue", &PclAllProcsTrue, grp,
					 "boolean", "boolean", "Returns true if all processes call the method with true.");

	reg.add_function("ParallelMin", &ParallelMin<double>, grp, "tmax", "t", "returns the maximum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelMax", &ParallelMax<double>, grp, "tmin", "t", "returns the minimum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelSum", &ParallelSum<double>, grp, "tsum", "t", "returns the sum of t over all processes. note: you have to assure that all processes call this function.");
}

#else // UG_PARALLEL

static bool PclDebugBarrierEnabledDUMMY()
{
	return false;
}

static bool PclDebugBarrierAllDUMMY()
{
	return false;
}

///	Dummy method for serial compilation always returning 1
static int NumProcsDUMMY()	{return 1;}

///	Dummy method for serial compilation always returning 0
static int ProcRankDUMMY()				{return 0;}

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

bool AllProcsTrueDUMMY(bool bTrue)
{
	return bTrue;
}

void RegisterBridge_PCL(Registry& reg, string parentGroup)
{
	string grp(parentGroup);
	grp.append("/PCL");

	reg.add_function("PclDebugBarrierEnabled", &PclDebugBarrierEnabledDUMMY, grp,
					"Enabled", "", "Returns the whether debug barriers are enabled.");

	reg.add_function("PclDebugBarrierAll", &PclDebugBarrierAllDUMMY, grp,
					 "", "", "Synchronizes all parallel processes if the executable"
							 "has been compiled with PCL_DEBUG_BARRIER=ON");

	reg.add_function("NumProcs", &NumProcsDUMMY, grp,
					"NumProcs", "", "Returns the number of active processes.");

	reg.add_function("ProcRank", &ProcRankDUMMY, grp,
					"ProcRank", "", "Returns the rank of the current process.");

	reg.add_function("SynchronizeProcesses", &SynchronizeProcessesDUMMY, grp,
					"", "", "Waits until all active processes reached this point.");

	reg.add_function("AllProcsTrue", &AllProcsTrueDUMMY, grp,
					 "boolean", "boolean", "Returns true if all processes call the method with true.");

	reg.add_function("ParallelMin", &ParallelMinDUMMY<double>, grp, "tmax", "t", "returns the maximum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelMax", &ParallelMaxDUMMY<double>, grp, "tmin", "t", "returns the minimum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelSum", &ParallelSumDUMMY<double>, grp, "tsum", "t", "returns the sum of t over all processes. note: you have to assure that all processes call this function.");
}

#endif //UG_PARALLEL

// end group pcl_bridge
/// \}

}// end of namespace
}// end of namespace
