/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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


#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_NATIVE_CUTHILL_MCKEE_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_NATIVE_CUTHILL_MCKEE_ORDERING__

#include <vector>

#include "IOrderingAlgorithm.h"
#include "util.cpp"

//debug
#include "common/error.h"
#include "common/log.h"

namespace ug{

/// returns an array describing the needed index mapping for Cuthill-McKee ordering
/**
 * This function computes an index mapping that transforms an index-graph by a
 * Cuthill-McKee ordering. For each index, a vector of all adjacent indices
 * must be passed.
 *
 * There are two intended ways of this function to work, depending on the flag bPreserveConsec:
 * If set to true (default), the function expects vvNeighbor to contain adjacency information
 * sorted by associated geometric object, i.e.:
 * Geometric object 0
 *   DoF 0
 *   DoF 1
 *   DoF 2
 *   ...
 * Geometric object 1
 *   DoF 0
 *   DoF 1
 *   DoF 2
 *   ...
 * ...
 *
 * None of the first indices of any geometric object must be unconnected!
 * Only the first DoF index of any geometric object must be given neighbor information!
 * This is required to guarantee that DoF indices associated to a geometric object stay consecutive.
 *
 * If set to false, indices without any adjacency information given in vvNeighbour
 * will be sorted to the end. Nothing is preserved. This is applicable if vvNeighbour
 * contains the final non-zero entries of a sparse matrix to be re-ordered, e.g., as a
 * pre-processing step for ILU(-T) pre-conditioning.
 *
 *
 * On exit, the index field vNewIndex is filled with the index mapping:
 * newInd = vNewIndex[oldInd]
 *
 * \param[out]	vNewIndex		vector returning new index for old index
 * \param[in]	vvNeighbour		vector of adjacent indices for each index
 * \param[in]	bReverse		flag if "reverse Cuthill-McKee" is used
 * \param[in]   bPreserveConsec flag indicating whether ordering is done for DofDistribution
 * \returns		flag if ordering was successful
 */
void ComputeCuthillMcKeeOrder(std::vector<size_t>& vNewIndex,
                              std::vector<std::vector<size_t> >& vvNeighbour,
                              bool bReverse = true,
			       bool bPreserveConsec = true);


template <typename TAlgebra, typename O_t>
class NativeCuthillMcKeeOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	NativeCuthillMcKeeOrdering() : m_bReverse(false) {}

	/// clone constructor
	NativeCuthillMcKeeOrdering( const NativeCuthillMcKeeOrdering<TAlgebra, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new NativeCuthillMcKeeOrdering<TAlgebra, O_t>(*this));
	}

	void compute(){
		std::vector<std::vector<size_t> > neighbors;
		neighbors.resize(m->num_rows());

		for(size_t i=0; i<m->num_rows(); i++)
		{
			for(typename M_t::row_iterator i_it = m->begin_row(i); i_it != m->end_row(i); ++i_it){
				neighbors[i].push_back(i_it.index());
			}
		}

		ComputeCuthillMcKeeOrder(o, neighbors, m_bReverse, true);

		m = NULL;

		#ifdef UG_DEBUG
		check();
		#endif
	}

	void check(){
		UG_COND_THROW(!is_permutation(o), name() << "::check: Not a permutation!");
	}

	O_t& ordering(){
		return o;
	}

	void init(M_t* A, const V_t&){
		init(A);
	}

	void init(M_t* A){
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		m = A;
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	virtual const char* name() const {
		if(m_bReverse){
			return "ReverseNativeCuthillMcKeeOrdering (ug4 version)";
		}
		else{
			return "NativeCuthillMcKeeOrdering (ug4 version)";
		}
	}
private:
	O_t o;
	M_t* m;

	bool m_bReverse;
};


} // end namespace ug

#endif
