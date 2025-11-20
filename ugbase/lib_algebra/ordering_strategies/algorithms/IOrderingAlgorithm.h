/*
 * Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch
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

#ifndef UG_BASE_LIB_ALGEBRA_ORDERING_STRATEGIES_ALGORITHMS_IORDERINGALGORITHM_H
#define UG_BASE_LIB_ALGEBRA_ORDERING_STRATEGIES_ALGORITHMS_IORDERINGALGORITHM_H

#include "common/util/smart_pointer.h"

namespace ug{

/*
	Interface for ordering algorithms. O_t denotes the ordering container type.
	Ordering algorithms have to implement
		compute() - triggers computation of an ordering
	as well as
		ordering() - returns a reference to an ordering of type O_t
	O_t usually is a std::vector<size_t>
	G_t (not required here) denotes the underlying graph type, e.g.,
                                                        a boost graph
*/

template <typename TAlgebra, typename O_t=std::vector<size_t> >
class IOrderingAlgorithm{
public:
	using M_t = typename TAlgebra::matrix_type;
	using V_t = typename TAlgebra::vector_type;

	IOrderingAlgorithm()= default;
	virtual ~IOrderingAlgorithm()= default;

	virtual void compute() = 0;
	virtual void check() = 0;

	virtual O_t& ordering() = 0;

	//lib_algebra only
	virtual void init(M_t*) = 0;
	//lib_disc
	virtual void init(M_t*, const V_t&) = 0;
	//lib_algebra only, induced matrix
	virtual void init(M_t*, const std::vector<size_t>&) = 0;
	//lib_disc, induced matrix
	virtual void init(M_t*, const V_t&, const std::vector<size_t>&) = 0;

	virtual SmartPtr<IOrderingAlgorithm > clone() = 0;

	virtual const char* name() const = 0;
};

}

#endif
