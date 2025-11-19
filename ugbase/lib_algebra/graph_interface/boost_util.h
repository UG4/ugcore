/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Author: Felix Salfelder, 2022
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

#ifndef UG_GRAPH_INTERFACE_BOOST_UTIL_H
#define UG_GRAPH_INTERFACE_BOOST_UTIL_H

#include <utility>
#include <boost/iterator/filter_iterator.hpp>

namespace ug{
namespace util{

namespace{
template<class G>
class noloop{
public:
	explicit noloop(const G& g) : _g(g) {}
	template<class E>
	bool operator()(E const&e) const{
		return boost::source(e, _g) != boost::target(e, _g);
	}
private:
	G const& _g;
};
}

template<class T, class G>
std::pair<boost::filter_iterator<noloop<G>, T>,
          boost::filter_iterator<noloop<G>, T>> omit_loops(std::pair<T, T> const& p, G const& g)
{
	noloop<G> P(g);
	using f = boost::filter_iterator<noloop<G>, T>;
	return std::make_pair(f(P, p.first, p.second), f(P, p.second, p.second));
}

} // util
} // ug

#endif
