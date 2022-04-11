/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_IORDERING_PREPROCESSOR__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_IORDERING_PREPROCESSOR__

#include "common/util/smart_pointer.h"

namespace ug{

/*
	Interface for preprocessor algorithms used for ordering algorithms.
	Ordering algorithms have to implement
		compute() - triggers preprocessing
	as well as
		init() - initializes preprocessing with a graph
	G_t denotes the underlying graph type, e.g., a boost graph
*/

template <typename TAlgebra, typename G_t>
class IOrderingPreprocessor{
public:
	IOrderingPreprocessor(){}
	virtual ~IOrderingPreprocessor(){}

	virtual void compute() = 0;

	//lib_algebra only
	virtual void init(const typename TAlgebra::vector_type&, G_t&) = 0;

	virtual SmartPtr<IOrderingPreprocessor<TAlgebra, G_t> > clone() = 0;

	virtual const char* name() const = 0;
};

} //namespace

#endif //guard
