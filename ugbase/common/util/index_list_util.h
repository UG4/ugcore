/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_index_list_util
#define __H__UG_index_list_util

#include <cctype>
#include <cstring>
#include <vector>
#include <sstream>

namespace ug {

/** Range strings are comma separated values and ranges, e.g.: "0,1,2-4,8,9-11"*/
template <class ind_t>
std::string IndexListToRangeString (const std::vector<ind_t>& inds)
{
	using namespace std;
	stringstream ss;

	for(size_t i = 0; i < inds.size(); ++i)
	{
		if(i > 0)
			ss << ",";

		const ind_t first = inds[i];
		ind_t cur = first;
		while((i + 1 < inds.size()) && inds[i+1] == cur + 1){
			++cur;
			++i;
		}

		if(first == cur)
			ss << first;
		else if(first + 1 == cur)
			ss << first << "," << cur;
		else
			ss << first << "-" << cur;
	}

	return ss.str();
}


/** Range strings are comma separated values and ranges, e.g.: "0,1,2-4,8,9-11"*/
template <class ind_t>
void RangeStringToIndexList (std::vector<ind_t>& indsOut, const char* rangeString)
{
	indsOut.clear();

	const char* cur = rangeString;
	
	bool readingRange = false;
	ind_t lastInd = 0;

	while(*cur) {
		ind_t ind = 0;
		while(*cur && *cur != ',' && *cur != '-'){
			if(isdigit(*cur))
				ind = ind * 10 + *cur - '0';
			++cur;
		}

		if(readingRange){
			readingRange = false;
			for(ind_t i = lastInd; i <= ind; ++i)
				indsOut.push_back (i);
		}
		else{
			switch(*cur){
				case 0: indsOut.push_back(ind); break;
				case ',': indsOut.push_back(ind); break;
				case '-': lastInd = ind; readingRange = true; break;
			}
		}

		if(*cur)
			++cur;		
	}
}

	
}//	end of namespace

#endif	//__H__UG_index_list_util
