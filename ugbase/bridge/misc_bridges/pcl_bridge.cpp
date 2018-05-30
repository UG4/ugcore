/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

template<typename T>
static vector<T> ParallelVecMin(const vector<T>& t)
{
	vector<T> tmp;
	pcl::ProcessCommunicator pc;
	pc.allreduce(t, tmp, PCL_RO_MIN);
	return tmp;
}

template<typename T>
static vector<T> ParallelVecMax(const vector<T>& t)
{
	vector<T> tmp;
	pcl::ProcessCommunicator pc;
	pc.allreduce(t, tmp, PCL_RO_MAX);
	return tmp;
}

template<typename T>
static vector<T> ParallelVecSum(const vector<T>& t)
{
	vector<T> tmp;
	pcl::ProcessCommunicator pc;
	pc.allreduce(t, tmp, PCL_RO_SUM);
	return tmp;
}


void RegisterBridge_PCL(Registry& reg, string parentGroup)
{
	string grp(parentGroup);
	grp.append("/pcl");

	reg.add_function("DisableMPIInit", &pcl::DisableMPIInit, grp, "", "",
	                 "Tells PCL to not call MPI_Init and MPI_Finalize.");

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

	reg.add_function("ParallelMin", &ParallelMin<double>, grp, "tmax", "t", "returns the minimum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelMax", &ParallelMax<double>, grp, "tmin", "t", "returns the maximum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelSum", &ParallelSum<double>, grp, "tsum", "t", "returns the sum of t over all processes. note: you have to assure that all processes call this function.");

	reg.add_function("ParallelVecMin", &ParallelVecMin<double>, grp, "tmax", "t", "returns the minimum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelVecMax", &ParallelVecMax<double>, grp, "tmin", "t", "returns the maximum of t over all processes. note: you have to assure that all processes call this function.");
	reg.add_function("ParallelVecSum", &ParallelVecSum<double>, grp, "tsum", "t", "returns the sum of t over all processes. note: you have to assure that all processes call this function.");
}

#else // UG_PARALLEL

void DisableMPIInitDUMMY ()
{}

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

	reg.add_function("DisableMPIInit", &DisableMPIInitDUMMY, grp, "", "",
	                 "Tells PCL to not call MPI_Init and MPI_Finalize.");

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
