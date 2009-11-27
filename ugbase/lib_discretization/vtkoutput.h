/*
 * vtkoutput.h
 *
 *  Created on: 06.07.2009
 *      Author: andreasvogel
 */

#ifndef VTKOUTPUT_H_
#define VTKOUTPUT_H_

#include <cstdio>
#include "lib_algebra/lib_algebra.h"
#include "lib_grid/lib_grid.h"
#include "lib_discretization/lib_discretization.h"

namespace ug{

class VTKOutput{

	public:
		bool print(NumericalSolution& u, const char* filename, double Time = 0.0);
		bool print_subset(NumericalSolution& u, int subsetIndex, const char* filename, double Time = 0.0);

	private:
		bool write_prolog(FILE* file, double Time);
		bool write_piece_prolog(FILE* file);
		bool write_subset(SubsetHandler& sh, int subsetIndex, NumericalSolution& u, FILE* File);
		bool InitNumbers(SubsetHandler& sh, uint subIndex);
		bool write_points(FILE* File, SubsetHandler& sh, uint subsetIndex);
		bool write_elements(FILE* File,SubsetHandler& sh, uint subsetIndex);
		bool write_scalar(FILE* File, NumericalSolution& u, int comp, SubsetHandler& sh, uint subsetIndex);
		bool write_epilog(FILE* file);
		bool write_piece_epilog(FILE* file);

		template <class TElem>
		bool VTKOutput::write_elements_connectivity(FILE* File, SubsetHandler& sh, uint subsetIndex);
		template <class TElem>
		bool VTKOutput::write_elements_offsets(FILE* File, SubsetHandler& sh, uint subsetIndex, int& n);
		template <class TElem>
		bool VTKOutput::write_elements_types(FILE* File, SubsetHandler& sh, uint subsetIndex);


	private:
		struct {
			int noVertices;
			int noElements;
			int noConnections;
		} Numbers;

	protected:
		typedef ug::Attachment<uint> ADOFIndex;

	protected:
		ADOFIndex m_aDOFIndex;
		Grid::VertexAttachmentAccessor<ADOFIndex> m_aaDOFIndexVRT;
		uint m_numberOfDOF;


};



}


#endif /* VTKOUTPUT_H_ */
