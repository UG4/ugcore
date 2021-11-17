/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Jan Friebertshäuser, Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS__DOWNWINDORDER__
#define __H__UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS__DOWNWINDORDER__

#include <vector>
#include <functional> // for lambda signature
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "bindings/lua/lua_user_data.h"
#include "common/log.h"

namespace ug{

/**
 * Numbers vertices based on a adjacency matrix.
 * Depth first traversal approach.
 * based on "NumeriereKnoten in chapter Downwind-Numbering in
 * <<Finite Volumen- und Mehrgitterverfahren für elliptische Randwertprobleme>>
 * by Jürgen Bey
 *
 * @param[in] vvConnections an adjacency vector of vectors of vertex indices
 * @param[in/out] vVisited a vector of booleans. vVisited[vertex] = true only if vertex was visited by NumeriereKnoten
 * @param[in] vAncestorsCount a vector of count of ancestor vertices per vertex.
 * @param[out] vNewIndex a vector of new indices per vector, usable for reordering
 * @param[in/out] N the highest Index up until now.
 * @param[in] v the vertex to inspect.
 */
void NumeriereKnoten(const std::vector<std::vector<size_t> > &vvConnections,
		std::vector<bool> &vVisited, std::vector<size_t> & vAncestorsCount, std::vector<size_t> & vNewIndex, size_t & N, size_t v);

/** Calculates Downwind-Numbering for one DofDistribution only.
 * uses NumeriereKnoten based on ideas of Jürgen Bey.
 *
 * @param[in] spDd DofDistribution to apply Downwind-Numbering on.
 * @param[in] spDomain Domain of the problem
 * @param[in] spVelocity Velocity field
 * @param[in] threshold Threshold for the angle between a connection and the Velocity field in radians.
 */
template <typename TDomain>
void OrderDownwindForDofDist(SmartPtr<DoFDistribution> spDd, ConstSmartPtr<TDomain> spDomain,
		SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > spVelocity, number threshold);

/**
 * Calculates Downwind Numbering for all DofDistributions of one Domain.
 * @param[in] approxSpace the domain.
 * @param[in] spVelocity the velocity field.
 */
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace,
					SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > spVelocity);

/**
 * Calculates Downwind Numbering for all DofDistributions of one Domain.
 * @param[in] approxSpace the domain.
 * @param[in] spVelocity the velocity field.
 * @param threshold threshold Threshold for the angle between a connection and the Velocity field in radians.
 */
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace,
					SmartPtr<UserData<MathVector<TDomain::dim>, TDomain::dim> > spVelocity, number threshold);

/**
 * Calculates Downwind Numbering for all DofDistributions of one Domain.
 * @param approxSpace the domain.
 * @param vVel a fixed velocity vector for a homogeneous velocity field.
 */
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace,
					const std::vector<number>& vVel);

/**
 * Calculates Downwind Numbering for all DofDistributions of one Domain.
 * @param approxSpace the domain.
 * @param vVel a fixed velocity vector for a homogeneous velocity field.
 * @param threshold threshold Threshold for the angle between a connection and the Velocity field in radians.
 */
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace,
					const std::vector<number>& vVel, number threshold);

#ifdef UG_FOR_LUA
/**
 * Calculates Downwind Numbering for all DofDistributions of one Domain.
 * @param approxSpace the domain.
 * @param strVelocity a lua callable to calculate the velocity field.
 */
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace, const char* strVelocity);

/**
 * Calculates Downwind Numbering for all DofDistributions of one Domain.
 * @param approxSpace the domain.
 * @param strVelocity a lua callable to calculate the velocity field.
 * @param threshold threshold Threshold for the angle between a connection and the Velocity field in radians.
 */
template <typename TDomain>
void OrderDownwind(ApproximationSpace<TDomain>& approxSpace, const char* strVelocity, number threshold);
#endif
} // end namespace ug

#endif
