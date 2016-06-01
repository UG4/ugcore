/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include <cmath>
#include "subdivision_rules_piecewise_loop.h"

using namespace std;

namespace ug
{

SubdivRules_PLoop::
SubdivRules_PLoop()
{
//	precalculate betas
//	the number is quite arbitrary here.
	size_t numPrecals = 16;
	m_betas.resize(numPrecals);
	for(size_t i = 0; i < numPrecals; ++i)
		m_betas[i] = calculate_beta(i);
}

SubdivRules_PLoop::
SubdivRules_PLoop(const SubdivRules_PLoop& src)
{
//	since this method won't ever be executed it can stay empty
}

SubdivRules_PLoop& SubdivRules_PLoop::
operator=(const SubdivRules_PLoop& src)
{
//	since this method won't ever be executed it can stay empty
	return *this;
}


number SubdivRules_PLoop::
get_beta(size_t valency) const
{
	if(valency < m_betas.size())
		return m_betas[valency];
		
	return calculate_beta(valency);
}

number SubdivRules_PLoop::
calculate_beta(size_t valency) const
{
	if(valency == 6)
		return 0.0625;
		
	if(valency > 0){
		const number tmp = 0.375 + 0.25 * cos((2.0*PI)/(number)valency);
		return (0.625 - tmp*tmp)/(number)valency;
	}

	return 0;
}


}//	end of namespace
