/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__MAP_ALGEBRA__SPARSEMATRIX_PRINT__
#define  __H__UG__MAP_ALGEBRA__SPARSEMATRIX_PRINT__

#include "mapsparsematrix.h"
#include "common/common.h"

namespace ug {

//!
//! print to console whole SparseMatrix
template<typename T>
void MapSparseMatrix<T>::print(const char * const text) const
{
	UG_LOG("================= MapSparseMatrix " << num_rows() << "x" << num_cols() << " =================\n");
	for(size_t i=0; i < num_rows(); i++)
		printrow(i);
}


//!
//! print the row row to the console
template<typename T>
void MapSparseMatrix<T>::printrow(size_t row) const
{
	UG_LOG("row " << row << ": ");
	for(const_row_iterator it=begin_row(row); it != end_row(row); ++it)
	{
		if(it.value() == 0.0) continue;
		UG_LOG(" ");
		UG_LOG("(" << it.index() << " -> " << it.value() << ")");
	}

	UG_LOG("\n");
}

template<typename T>
void MapSparseMatrix<T>::printtype() const
{
	std::cout << *this;
}

}
#endif
