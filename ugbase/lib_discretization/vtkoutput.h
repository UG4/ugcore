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

template <int d>
class VTKOutput{

	public:
		bool print(NumericalSolution<d>& u, int level, const char* filename, double Time = 0.0);
		bool print_subset(NumericalSolution<d>& u, int level, int subsetIndex, const char* filename, double Time = 0.0);

	private:
		bool write_prolog(FILE* file, double Time);
		bool write_piece_prolog(FILE* file);
		bool write_subset(ISubsetHandler& sh, int subsetIndex, NumericalSolution<d>& u, int level, FILE* File);
		bool init_subset(ISubsetHandler& sh, uint subIndex, int level);
		bool write_points(FILE* File, ISubsetHandler& sh, uint subsetIndex, typename Domain<d>::position_attachment_type* aPos);
		bool write_elements(FILE* File,ISubsetHandler& sh, uint subsetIndex);
		bool write_scalar(FILE* File, NumericalSolution<d>& u, int comp, ISubsetHandler& sh, uint subsetIndex);
		bool write_epilog(FILE* file);
		bool write_piece_epilog(FILE* file);

		template <class TElem>
		bool write_elements_connectivity(FILE* File, typename geometry_traits<TElem>::iterator iterBegin, typename geometry_traits<TElem>::iterator iterEnd, ISubsetHandler& sh, uint subsetIndex);
		template <class TElem>
		bool write_elements_offsets(FILE* File, typename geometry_traits<TElem>::iterator iterBegin, typename geometry_traits<TElem>::iterator iterEnd,ISubsetHandler& sh, uint subsetIndex, int& n);
		template <class TElem>
		bool write_elements_types(FILE* File, typename geometry_traits<TElem>::iterator iterBegin, typename geometry_traits<TElem>::iterator iterEnd, ISubsetHandler& sh, uint subsetIndex);


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
