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

#ifndef __H__UG__LIB_DISC__DOF_MANAGER__BOOST_INTERFACE__HPP
#define __H__UG__LIB_DISC__DOF_MANAGER__BOOST_INTERFACE__HPP


namespace ug{

template <typename G_t, typename O_t>
void boost_Cuthill_McKee_ordering(G_t &g, O_t &o, bool reverse);

template <typename G_t, typename O_t>
void boost_minimum_degree_ordering(G_t &G, O_t &O);


template <typename TDomain>
void OrderBoostMinimumDegree(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

template <typename TDomain>
void OrderBoostCuthillMcKee(SmartPtr<ApproximationSpace<TDomain> > approxSpace, bool reverse);


} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__BOOST_INTERFACE__HPP */
