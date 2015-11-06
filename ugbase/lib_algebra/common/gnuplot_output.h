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

#ifndef GNUPLOT_OUTPUT_H_
#define GNUPLOT_OUTPUT_H_

template<typename Vector_type, typename postype>
void WriteVectorGnuplot(std::string filename, const Vector_type &v,
		postype *positions, int dimensions, const Vector_type *compareVec=NULL)
{
	size_t N = A.num_rows();
	std::fstream f((filename+".sh").c_str(), std::ios::out);
	f << "#!/bin/bash\n"
			"cat > gnuplotTemporaryFile <<EOF\n"
			"set dgrid3d N,N\n"
			"set style data lines\n"
			"set pm3d\n"
			"splot \"-\" pal\n";
	if(dimensions == 1)
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " " << v[i] << "\n";
	else if(dimensions == 2)
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " " << positions[i][1] << " " << v[i] << "\n";
	else
		for(size_t i=0; i < N; i++)
			f << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << " " << v[i] << "\n";
	f <<	"e\n"
			"EOF\n"
			"cat gnuplotTemporaryFile | gnuplot -persist\n"
			"rm gnuplotTemporaryFile\n";
}


#endif /* GNUPLOT_OUTPUT_H_ */
