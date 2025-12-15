/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Jan Friebertshäuser
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


#include "lib_disc/ordering_strategies/algorithms/downwindorder.h"

//ø #include "common/common.h"
#include <iostream>
#include <utility> // for pair
//ø #include <map> // for graph structure

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/dof_position_util.h"

namespace ug {
/**
 * For Debugging:
 */
DebugID LIB_DISC_ORDER("LIB_DISC_ORDER");

// Numbers vertices based on a adjacency matrix.
void NumeriereKnoten(const std::vector<std::vector<size_t> > &vvConnections,
		std::vector<bool> &vVisited, std::vector<size_t> & vAncestorsCount, std::vector<size_t> & vNewIndex, size_t & N, size_t v)
{
	vVisited[v] = true;
	vNewIndex[v] = N;
	N++;
	const std::vector<size_t> connections = vvConnections[v];
	for (auto AdjacencIter = connections.begin(); AdjacencIter != connections.end(); ++AdjacencIter)
	{
		vAncestorsCount[*AdjacencIter]--;
		if (vAncestorsCount[*AdjacencIter] == 0)
			NumeriereKnoten(vvConnections, vVisited, vAncestorsCount, vNewIndex, N, *AdjacencIter);
	}
}

// Calculates Downwind-Numbering for one DofDistribution only.
template <typename TDomain>
void OrderDownwindForDofDist(SmartPtr<DoFDistribution> dd, ConstSmartPtr<TDomain> domain,
		SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > spVelocity,
		number time, int si, number threshold)
{
	static constexpr int dim = TDomain::dim;
	const size_t num_ind = dd->num_indices();
	using pos_type = std::pair<MathVector<dim>, size_t>;
	using adjacency_type = std::vector<std::vector<size_t> >;

	//	get positions of indices
	typename std::vector<pos_type> vPositions;
	ExtractPositions(domain, dd, vPositions);

	// get adjacency vector of vectors
	adjacency_type vvConnections;
	dd->get_connections(vvConnections);

	// Check vector sizes match
	if (vvConnections.size() != num_ind)
		UG_THROW("OrderDownstreamForDofDist: "
				"Adjacency list of dimension " << num_ind << " expected, got "<< vvConnections.size());

	if (vPositions.size() != num_ind)
		UG_THROW("OrderDownstreamForDofDist: "
				"Position list of dimension " << num_ind << " expected, got "<< vPositions.size());

	// init helper structures
	std::vector<size_t> vNewIndex(num_ind, 0);
	std::vector<size_t> vAncestorsCount(num_ind, 0);
	std::vector<bool> vVisited(num_ind, false);

	// remove connections that are not in stream direction
	adjacency_type::iterator VertexIter;
	std::vector<size_t>::iterator AdjIter;
	std::vector<size_t> vAdjacency;

	MathVector<TDomain::dim> vVel1, vPos1, vPos2, vDir1_2;
	size_t i;
	for (VertexIter = vvConnections.begin(), i=0; VertexIter != vvConnections.end(); VertexIter++, i++)
	{
		UG_DLOG(LIB_DISC_ORDER, 2, "Filtering vertex " << i << " adjacency vector." <<std::endl);
		// count how many vertex were kept / removed per adjacency vector
		#ifdef UG_ENABLE_DEBUG_LOGS
		size_t initialcount = VertexIter->size();
		#endif
		size_t kept = 0, removed = 0;
		// get position and velocity of first trait
		vPos1 = vPositions.at(i).first;
		(*spVelocity)(vVel1, vPos1, time, si);
		if (VecLengthSq(vVel1) == 0 )
		{
			// if the velocity is zero at this trait it does not interfere with others
			// NOTE: otherwise this trait would be downwind-connected to all of it's neighbors
			// NOTE: VertexIter-> will access inner vector functions (*VertexIter) is the inner vector.
			removed = VertexIter->size();
			VertexIter->clear();
		}
		else {
			AdjIter = VertexIter->begin();
			while (AdjIter != VertexIter->end())
			{
				// get position of second trait
				vPos2 = vPositions.at(*AdjIter).first;

				// get difference vector as direction vector
				VecSubtract(vDir1_2, vPos2, vPos1);

				// compute angle between velocity and direction vector
				number anglex1_2 = VecAngle(vDir1_2, vVel1);

				// if angle is smaller then threshold continue else remove connection
				if (anglex1_2 <= threshold && i != *AdjIter)
				{
					vAncestorsCount.at(*AdjIter) += 1;
					++AdjIter;
					kept++;
				} else {
					AdjIter = VertexIter->erase(AdjIter);
					removed++;
				}
			}
		}
		UG_DLOG(LIB_DISC_ORDER, 2, "Kept: " << kept << ", removed: " << removed << " of " << initialcount
				<< " entries in adjacency matrix." << std::endl << std::endl);
	}
	// calculate downwindorder
	// Find vertexes without any ancestors and start NumeriereKnoten on them.
	size_t v,N;
	for (v=0, N=0; v < vvConnections.size(); v++)
	{
		if (vAncestorsCount[v] == 0 && !vVisited[v])
		{
			NumeriereKnoten(vvConnections, vVisited, vAncestorsCount, vNewIndex, N, v);
		}
	}

	// sanity check
	if (N < vvConnections.size()){
		size_t fails = 0;
		for (v=0; v < vvConnections.size(); v++) {
			if (!vVisited[v]) {
				UG_DLOG(LIB_DISC_ORDER, 2, v << "was not visited, has unresolved ancestors: " << vAncestorsCount[v] << std::endl);
				fails ++;
			}
		}
		UG_THROW("OrderDownwindForDist failed, " << fails << " traits unvisited." << std::endl);
	}

	//	reorder traits
	dd->permute_indices(vNewIndex);
}

// Calculates Downwind Numbering for all DofDistributions of one Domain.
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace,
		SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > spVelocity, number threshold)
{
	UG_LOG ("OrderDownwind: This function is obsolete and may cause problems. Avoid it! Alternatives: Ordering strategies in solvers etc.\n");
	
	// TODO: implement for variable time and subset
	int si = 0;
	std::vector<SmartPtr<DoFDistribution> > vDD = approxSpace.dof_distributions();
	UG_DLOG(LIB_DISC_ORDER, 2, "Starting DownwindOrdering." << std::endl);
	for(size_t i = 0; i < vDD.size(); ++i){
		number time = 0.0;
		UG_DLOG(LIB_DISC_ORDER, 2, "Ordering Domain Distribution " << i << "." << std::endl);
		OrderDownwindForDofDist<TDomain>(vDD[i], approxSpace.domain(), spVelocity, time, si, threshold);
	}
}

// Calculates Downwind Numbering for all DofDistributions of one Domain.
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace,
		SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > spVelocity)
{
	OrderDownwind<TDomain>(approxSpace, spVelocity, PI/4.0);
}

// Calculates Downwind Numbering for all DofDistributions of one Domain.
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace,
					const std::vector<number>& vVel)
{
	static constexpr int dim = TDomain::dim;
	if(vVel.size() != dim){
		UG_THROW("OrderDownstream: Velocity field of dimension " << dim << " expected, got "<< vVel.size());
	}

	OrderDownwind<TDomain>(approxSpace,  SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)) );
}

// Calculates Downwind Numbering for all DofDistributions of one Domain.
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace,
					const std::vector<number>& vVel, number threshold)
{
	static constexpr int dim = TDomain::dim;
	if(vVel.size() != dim){
		UG_THROW("OrderDownstream: Velocity field of dimension " << dim << " expected, got "<< vVel.size());
	}

	OrderDownwind<TDomain>(approxSpace,  SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)), threshold );
}


#ifdef UG_FOR_LUA

// Calculates Downwind Numbering for all DofDistributions of one Domain.
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace, const char* strVelocity)
{
	static constexpr int dim = TDomain::dim;

	SmartPtr<UserData<MathVector<dim>, dim> > spVelocity
	 = make_sp(new LuaUserData<MathVector<dim>, dim>(strVelocity));

	OrderDownwind<TDomain>(approxSpace, spVelocity);
}

// Calculates Downwind Numbering for all DofDistributions of one Domain.
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace, const char* strVelocity, number threshold)
{
	static constexpr int dim = TDomain::dim;

	SmartPtr<UserData<MathVector<dim>, dim> > spVelocity = make_sp(new LuaUserData<MathVector<dim>, dim>(strVelocity));

	OrderDownwind<TDomain>(approxSpace, spVelocity, threshold);
}

#ifdef UG_DIM_1
template void OrderDownwind<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, const char* strVelocity);
template void OrderDownwind<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, const char* strVelocity, number threshold);
#endif
#ifdef UG_DIM_2
template void OrderDownwind<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, const char* strVelocity);
template void OrderDownwind<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, const char* strVelocity, number threshold);
#endif
#ifdef UG_DIM_3
template void OrderDownwind<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, const char* strVelocity);
template void OrderDownwind<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, const char* strVelocity, number threshold);
#endif

#endif

#ifdef UG_DIM_1
template void OrderDownwind<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, SmartPtr<UserData<MathVector<Domain1d::dim>, Domain1d::dim> > spVelocity);
template void OrderDownwind<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, SmartPtr<UserData<MathVector<Domain1d::dim>, Domain1d::dim> > spVelocity, number threshold);
template void OrderDownwind<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, const std::vector<number>& vVel);
template void OrderDownwind<Domain1d>(ApproximationSpace<Domain1d>& approxSpace, const std::vector<number>& vVel, number threshold);
#endif

#ifdef UG_DIM_2
template void OrderDownwind<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, SmartPtr<UserData<MathVector<Domain2d::dim>, Domain2d::dim> > spVelocity);
template void OrderDownwind<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, SmartPtr<UserData<MathVector<Domain2d::dim>, Domain2d::dim> > spVelocity, number threshold);
template void OrderDownwind<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, const std::vector<number>& vVel);
template void OrderDownwind<Domain2d>(ApproximationSpace<Domain2d>& approxSpace, const std::vector<number>& vVel, number threshold);
#endif

#ifdef UG_DIM_3
template void OrderDownwind<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, SmartPtr<UserData<MathVector<Domain3d::dim>, Domain3d::dim> > spVelocity);
template void OrderDownwind<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, SmartPtr<UserData<MathVector<Domain3d::dim>, Domain3d::dim> > spVelocity, number threshold);
template void OrderDownwind<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, const std::vector<number>& vVel);
template void OrderDownwind<Domain3d>(ApproximationSpace<Domain3d>& approxSpace, const std::vector<number>& vVel, number threshold);
#endif

} // end namespace ug
