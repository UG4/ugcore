/*
 * connection_viewer_output.h
 *
 *  Created on: 03.08.2011
 *      Author: anaegel
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
#endif /* CSVOUTPUT_H_ */
