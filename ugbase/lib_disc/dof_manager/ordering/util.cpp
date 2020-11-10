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
 
#ifndef __H__UG__LIB_DISC__DOF_MANAGER__UTIL__CPP
#define __H__UG__LIB_DISC__DOF_MANAGER__UTIL__CPP

#ifndef error
#define error() \
	std::cerr << "error " << __FILE__ << ":" << __LINE__ << ":" << __func__ << "\n"
#endif


#include "types.hpp"


/* vector<BOOL> is indead a vector of booleans, but vector<bool> isnt! */
class BOOL{
public:
	BOOL() : value_(bool()){}
	/* explicit */ BOOL(bool const& t): value_(t) {}
	// /* explicit */ operator bool&() { return value_; }
	/* explicit */ operator bool() const { return value_; }
private:
	char value_;
};


template <typename G_t>
void copy_graph(G_t &orig, Graph_t &copy){
	copy = Graph_t(boost::num_vertices(orig));

	typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd; 
	for(boost::tie(eIt, eEnd) = boost::edges(orig); eIt != eEnd; ++eIt){
		//if(!boost::edge(boost::source(*eIt, copy), boost::target(*eIt, copy), copy).second){
			boost::add_edge(boost::source(*eIt, orig), boost::target(*eIt, orig), copy);
		//}
	}
}


template <typename O_t>
bool is_permutation(O_t &o){
	std::vector<BOOL> container(o.size(), false);
	for(unsigned i = 0; i < o.size(); ++i){
		if(!container[o[i]]){
			container[o[i]] = true;
		}
		else{
			return false; //no doubles allowed
		}
	}
	
	for(unsigned i = 0; i < o.size(); ++i){
		if(!container[i]){
			return false;
		}
	}

	return true;
}



#endif
