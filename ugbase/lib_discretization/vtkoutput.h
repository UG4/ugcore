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
		bool print(NumericalSolution& u, int level, const char* filename, double Time = 0.0);
		bool print_subset(NumericalSolution& u, int level, int subsetIndex, const char* filename, double Time = 0.0);

	private:
		bool write_prolog(FILE* file, double Time);
		bool write_piece_prolog(FILE* file);
		bool write_subset(SubsetHandler& sh, int subsetIndex, NumericalSolution& u, int level, FILE* File);
		bool init_subset(SubsetHandler& sh, uint subIndex, int level);
		bool write_points(FILE* File, SubsetHandler& sh, uint subsetIndex);
		bool write_elements(FILE* File,SubsetHandler& sh, uint subsetIndex);
		bool write_scalar(FILE* File, NumericalSolution& u, int comp, SubsetHandler& sh, uint subsetIndex);
		bool write_epilog(FILE* file);
		bool write_piece_epilog(FILE* file);

		template <class TElem>
		bool VTKOutput::write_elements_connectivity(FILE* File, typename geometry_traits<TElem>::iterator iterBegin, typename geometry_traits<TElem>::iterator iterEnd, SubsetHandler& sh, uint subsetIndex);
		template <class TElem>
		bool VTKOutput::write_elements_offsets(FILE* File, typename geometry_traits<TElem>::iterator iterBegin, typename geometry_traits<TElem>::iterator iterEnd,SubsetHandler& sh, uint subsetIndex, int& n);
		template <class TElem>
		bool VTKOutput::write_elements_types(FILE* File, typename geometry_traits<TElem>::iterator iterBegin, typename geometry_traits<TElem>::iterator iterEnd, SubsetHandler& sh, uint subsetIndex);


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

	protected:
		int _level;
		geometry_traits<Vertex>::iterator _iterVRT, _iterBeginVRT, _iterEndVRT;
		geometry_traits<Triangle>::iterator _iterTRIANGLE, _iterBeginTRIANGLE, _iterEndTRIANGLE;
		geometry_traits<Quadrilateral>::iterator _iterQUADRILATERAL, _iterBeginQUADRILATERAL, _iterEndQUADRILATERAL;


};



}


#endif /* VTKOUTPUT_H_ */
