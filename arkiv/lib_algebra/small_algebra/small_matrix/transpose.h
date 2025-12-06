/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__UG__COMMON__TRANSPOSE_H_
#define __H__UG__COMMON__TRANSPOSE_H_


template<typename T>
class TE_TRANSPOSED
{
public:
	using value_type = typename T::value_type;
	TE_TRANSPOSED(const T&_t) : t(_t) {}
	inline size_t num_rows() const
	{
		return t.num_cols();
	}

	inline size_t num_cols() const
	{
		return t.num_rows();
	}

	const value_type &operator () (size_t r, size_t c) const
	{
		return t(c, r);
	}

	value_type &operator () (size_t r, size_t c)
	{
		return t(c, r);
	}
private:
	const T &t;
};

template<typename T>
inline TE_TRANSPOSED<T> transpose(const T &t)
{
	return TE_TRANSPOSED<T>(t);
}

inline const double &transpose(const double &t)
{
	return t;
}

inline double &transpose(double &t)
{
	return t;
}

#endif