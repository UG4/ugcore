/*
 * Copyright (c) 2011-2022:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel, Lukas Larisch
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
#ifndef IG_UGBASE_LIB_DISC_ORDERING_STRATEGIES_ALGORITHMS_LEXORDER_COMPARATORS_H
#define IG_UGBASE_LIB_DISC_ORDERING_STRATEGIES_ALGORITHMS_LEXORDER_COMPARATORS_H

namespace ug {


// Order for 1D
template<int dim, size_t orderDim>
bool ComparePosDim(const std::pair<MathVector<dim>, size_t> &p1,
                   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

template<>
bool ComparePosDim<1, 0>(const std::pair<MathVector<1>, size_t> &p1,
                                const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]<p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDim<2,0>(const std::pair<MathVector<2>, size_t> &p1,
                      const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]) {
		return (p1.first[1] < p2.first[1]);
	}
	else {
		return (p1.first[0] < p2.first[0]);
	}
};
template<>
bool ComparePosDim<2,1>(const std::pair<MathVector<2>, size_t> &p1,
                      const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]) {
		return (p1.first[0] < p2.first[0]);
	}
	else {
		return (p1.first[1] < p2.first[1]);
	}
};

// Order for 3D
template<>
bool ComparePosDim<3,0>(const std::pair<MathVector<3>, size_t> &p1,
                      const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[2]==p2.first[2]){
		if (p1.first[1]==p2.first[1]) {
			return p1.first[0] < p2.first[0];
		}
		else {
			return (p1.first[1] < p2.first[1]);
		}
	}
	else{
		return (p1.first[2] < p2.first[2]);
	}
};
template<>
bool ComparePosDim<3,1>(const std::pair<MathVector<3>, size_t> &p1,
                      const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]){
		if (p1.first[2]==p2.first[2]) {
			return p1.first[1] < p2.first[1];
		}
		else {
			return (p1.first[2] < p2.first[2]);
		}
	}
	else{
		return (p1.first[0] < p2.first[0]);
	}
};
template<>
bool ComparePosDim<3,2>(const std::pair<MathVector<3>, size_t> &p1,
                      const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]){
		if (p1.first[0]==p2.first[0]) {
			return p1.first[2] < p2.first[2];
		}
		else {
			return (p1.first[0] < p2.first[0]);
		}
	}
	else{
		return (p1.first[1] < p2.first[1]);
	}
};


//Decreasing

// Order for 1D
template<int dim, size_t orderDim>
bool ComparePosDimDec(const std::pair<MathVector<dim>, size_t> &p1,
                   const std::pair<MathVector<dim>, size_t> &p2)
{return false;}

template<>
bool ComparePosDimDec<1, 0>(const std::pair<MathVector<1>, size_t> &p1,
                                const std::pair<MathVector<1>, size_t> &p2)
{
	return (p1.first[0]>p2.first[0]);
};

// Order for 2D
template<>
bool ComparePosDimDec<2,0>(const std::pair<MathVector<2>, size_t> &p1,
                      const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]) {
		return (p1.first[1] > p2.first[1]);
	}
	else {
		return (p1.first[0] > p2.first[0]);
	}
};
template<>
bool ComparePosDimDec<2,1>(const std::pair<MathVector<2>, size_t> &p1,
                      const std::pair<MathVector<2>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]) {
		return (p1.first[0] > p2.first[0]);
	}
	else {
		return (p1.first[1] > p2.first[1]);
	}
};

// Order for 3D
template<>
bool ComparePosDimDec<3,0>(const std::pair<MathVector<3>, size_t> &p1,
                      const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[2]==p2.first[2]){
		if (p1.first[1]==p2.first[1]) {
			return p1.first[0] > p2.first[0];
		}
		else {
			return (p1.first[1] > p2.first[1]);
		}
	}
	else{
		return (p1.first[2] > p2.first[2]);
	}
};
template<>
bool ComparePosDimDec<3,1>(const std::pair<MathVector<3>, size_t> &p1,
                      const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[0]==p2.first[0]){
		if (p1.first[2]==p2.first[2]) {
			return p1.first[1] > p2.first[1];
		}
		else {
			return (p1.first[2] > p2.first[2]);
		}
	}
	else{
		return (p1.first[0] > p2.first[0]);
	}
};
template<>
bool ComparePosDimDec<3,2>(const std::pair<MathVector<3>, size_t> &p1,
                      const std::pair<MathVector<3>, size_t> &p2)
{
	if (p1.first[1]==p2.first[1]){
		if (p1.first[0]==p2.first[0]) {
			return p1.first[2] > p2.first[2];
		}
		else {
			return (p1.first[0] > p2.first[0]);
		}
	}
	else{
		return (p1.first[1] > p2.first[1]);
	}
};


}

#endif