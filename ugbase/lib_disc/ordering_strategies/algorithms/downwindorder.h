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
#include "lib_algebra/cpu_algebra_types.h"
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



/**
 * Numbers vertices based on a adjacency matrix. Depth first traversal approach.
 * Based on "NumeriereKnoten in chapter Downwind-Numbering in
 * <<Finite Volumen- und Mehrgitterverfahren für elliptische Randwertprobleme>>
 * by Jürgen Bey.
 *
 * @param[in] 		vvConnections	an adjacency vector of vectors of vertex indices
 * @param[in/out]	vVisited 		a vector of booleans. vVisited[vertex] = true,
 *									only if vertex was visited by NumberVertices
 * @param[in]		vAncestorsCount a vector of count of ancestor vertices per vertex.
 * @param[out]		vNewIndex 		a vector of new indices per vector, usable for reordering
 * @param[in/out]	N 				the highest Index up until now.
 * @param[in]		v 				the vertex to inspect.
 */
void NumberVertices(const std::vector<std::vector<size_t> >& vvConnections,
					std::vector<bool>&   vVisited,
					std::vector<size_t>& vAncestorsCount,
					std::vector<size_t>& vNewIndex,
					size_t & N, size_t v);

/** Calculates Downwind-Numbering for one DofDistribution only.
 * uses NumeriereKnoten based on ideas of Jürgen Bey.
 *
 * @param[in] spDd 		DofDistribution to apply Downwind-Numbering on.
 * @param[in] m 		stiffness matrix
 */
template <typename TDomain>
void OrderDownwindStiffForDofDist(SmartPtr<DoFDistribution> spDd, const CPUAlgebra::matrix_type& m);

/**
 * Calculates downwind numbering based on coupling in stiffness matrix m
 * for all DofDistributions of one Domain.
 *
 * @param[in] approxSpace 	the domain
 * @param[in] m 			stiffness matrix
 */
template <typename TDomain>
void OrderDownwindStiff(ApproximationSpace<TDomain>& approxSpace, const CPUAlgebra::matrix_type& m);

/*

template <typename TAlgebra, typename TDomain, typename O_t>
class DownwindStiffOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	/// Grid function type for the solution
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

	/// Position attachment type
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> Position_t;

	DownwindStiffOrdering(){ std::cout << "DownwindStiffOrdering: HALLO" << std::endl;}

	/// clone constructor
	DownwindStiffOrdering( const DownwindStiffOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new DownwindStiffOrdering<TAlgebra, TDomain, O_t>(*this));
	}

	void compute(){
		mat = NULL;
	}

	void check(){
		if(!is_permutation(o)){
			UG_THROW(name() << "::check: Not a permutation!");
		}
	}

	O_t& ordering(){
		return o;
	}

	void init(M_t* A, const V_t& V){
		#ifdef UG_ENABLE_DEBUG_LOGS
		//UG_LOG("Using " << name() << " in " << m_dir << " direction\n");
		#endif

		mat = A;
	}

	void init(M_t*){
		UG_THROW(name() << "::init: Cannot initialize smoother without a geometry. Specify the 2nd argument for init!");
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	virtual const char* name() const {return "DownwindStiffOrdering";}

private:
	O_t o;
	M_t* mat;

	//const char *m_dir;
	//std::vector<Position_t> m_vPositions;
};

*/



} // end namespace ug

#endif
