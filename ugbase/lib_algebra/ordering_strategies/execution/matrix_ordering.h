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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_EXECUTION_MATRIX_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_EXECUTION_MATRIX_ORDERING__

#include <vector>

#include "../../../common/code_marker.h" //error()

namespace ug{

/*
	Orders a source matrix according to an ordering and stores it as target matrix.
*/

template <typename matrix_type, typename O_t=std::vector<size_t> >
void reorder_matrix(matrix_type& target, matrix_type& source, O_t& o){
	target.resize_and_clear(source.num_rows(), source.num_cols());

	for(size_t r = 0; r < source.num_rows(); ++r){
		size_t Pr = o[r];
		for(typename matrix_type::const_row_iterator it = source.begin_row(r); it != source.end_row(r); ++it){
			size_t Pc = o[it.index()];
			target(Pr, Pc) = it.value();
		}
	}
}

} //namespace

#endif //guard
