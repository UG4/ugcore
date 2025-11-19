/*
 * Copyright (c) 2011-2013:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_ALGEBRA__STL_DEBUG__
#define __H__LIB_ALGEBRA__STL_DEBUG__

#include <vector>

#ifndef NDEBUG
#include "common/error.h"

namespace ug{

template<typename T, class Allocator = std::allocator<T> >
class stdvector : public std::vector<T, Allocator>
{
private:
	using super = std::vector<T, Allocator>;
public:
	using size_type = typename super::size_type;

private:

	inline void size_check(size_t i) const
	{
		UG_COND_THROW(i >= super::size(), "accessing element " << i << " but vector only has size " << super::size() << ".")
	}

public:
	explicit stdvector(const Allocator& a= Allocator()) : super(a) { }
	explicit stdvector(size_type n, const T& value = T(), const Allocator& a= Allocator()) : super(n, value, a) { }

	stdvector(const std::vector<T, Allocator> &x) : super(x) { }
	stdvector(const stdvector<T, Allocator> &x) : super(x) { }

	inline typename super::reference operator[] (size_t i)
	{
		size_check(i);
		return super::operator[](i);
	}
	inline typename super::const_reference operator[] (size_t i) const
	{
		size_check(i);
		return super::operator[](i);
	}
};

}
#else

#define stdvector std::vector

#endif


#endif