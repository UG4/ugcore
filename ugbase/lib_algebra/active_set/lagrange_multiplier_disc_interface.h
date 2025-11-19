/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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

#ifndef LAGRANGE_MULTIPLIER_DISC_INTERFACE_H_
#define LAGRANGE_MULTIPLIER_DISC_INTERFACE_H_

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/common/multi_index.h"

namespace ug{

template <typename TDomain, typename TGridFunction>
class ILagrangeMultiplierDisc
{
	private:
	///	own type
		using this_type = ILagrangeMultiplierDisc<TDomain, TGridFunction>;

	public:
	///	Domain type
		using domain_type = TDomain;

	///	World dimension
		static constexpr int dim = TDomain::dim;

	public:
		ILagrangeMultiplierDisc(){};

	/// Virtual destructor
		virtual ~ILagrangeMultiplierDisc() {}

		virtual void lagrange_multiplier(TGridFunction& lagMult, const TGridFunction& u,
				std::vector<DoFIndex> vActiveSet, std::vector<int> vActiveSubsets) = 0;
};

} //end namespace ug

#endif