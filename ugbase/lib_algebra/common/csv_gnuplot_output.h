/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

/*
 *      (based on connection_viewer_output.h by mrupp)
 */

#ifndef CSVGNUPLOTOUTPUT_H_
#define CSVGNUPLOTOUTPUT_H_

#include <iostream>
#include <fstream>
#include <vector>

namespace ug
{



// WriteVectorCSV (derived from ConnectionViewer output)
//--------------------------------------------------
/**
 * \brief writes to a file in somewhat SparseMatrix-market format (for connection viewer)
 * \param filename Filename to write matrix to
 * \param b Vector
 * \param positions Positions, there has to be one position for each i in (0, ..., max(A.num_rows(), A.num_cols())).
 * \param dimensions	Dimensions of Positions
 */
template<typename Vector_type, typename postype>
void WriteVectorCSV(const char *filename, const Vector_type &b, postype *positions, int dimensions)
{
	std::fstream file(filename, std::ios::out);
	size_t rows = b.size();

	// write status
	file << "# dim=" << dimensions<<  ", rows="<< rows <<std::endl;
	// write data
	if(dimensions == 1)
		for(size_t i=0; i < rows; i++)
				file << positions[i][0] << ", " << b[i] <<std::endl;

	else if(dimensions == 2)
		for(size_t i=0; i < rows; i++)
			file << positions[i][0] << ", " << positions[i][1] << ", " << b[i]<< std::endl;
	else
		for(size_t i=0; i < rows; i++)
		  file << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2] << ", " << b[i] << std::endl;

}




} // namespace ug
#endif
