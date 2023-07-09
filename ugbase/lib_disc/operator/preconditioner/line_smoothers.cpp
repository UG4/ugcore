/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#include "line_smoothers.h"
namespace ug{
template<int dim>
bool ComparePosDimYDir(const std::pair<MathVector<dim>, size_t> &p1,
					   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

template<>
bool ComparePosDimYDir<1>(const std::pair<MathVector<1>, size_t> &p1,
						  const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDimYDir<2>(const std::pair<MathVector<2>, size_t> &p1,
						  const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]) {
		return (p1.first[1] < p2.first[1]);
	}
	else {
		return (p1.first[0] < p2.first[0]);
	}
};

// Order for 3D
template<>
bool ComparePosDimYDir<3>(const std::pair<MathVector<3>, size_t> &p1,
						  const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[2]==p2.first[2]){
		if (p1.first[0]==p2.first[0]) {
			return (p1.first[1] < p2.first[1]);
		}
		else {
			return (p1.first[0] < p2.first[0]);
		}
	}
	else{
		return (p1.first[2] < p2.first[2]);
	}
};




template<int dim>
bool ComparePosDimZDir(const std::pair<MathVector<dim>, size_t> &p1,
					   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

// Order for 1D
template<>
bool ComparePosDimZDir<1>(const std::pair<MathVector<1>, size_t> &p1,
						  const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDimZDir<2>(const std::pair<MathVector<2>, size_t> &p1,
						  const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]) {
		return (p1.first[1] < p2.first[1]);
	}
	else {
		return (p1.first[0] < p2.first[0]);
	}
};

// Order for 3D
template<>
bool ComparePosDimZDir<3>(const std::pair<MathVector<3>, size_t> &p1,
						  const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]){
		if (p1.first[0]==p2.first[0]) {
			return (p1.first[2] < p2.first[2]);
		}
		else {
			return (p1.first[0] < p2.first[0]);
		}
	}
	else{
		return (p1.first[1] < p2.first[1]);
	}
};

} // namespace ug
