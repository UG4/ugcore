/*
 * vtkoutput.h
 *
 *  Created on: 06.07.2009
 *      Author: andreasvogel
 */

#ifndef VTKOUTPUT_H_
#define VTKOUTPUT_H_

// extern libraries
#include <cstdio>

// other ug modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/function_spaces/grid_function_space.h"

namespace ug{


/// ATTENTION: This class uses heavily the mark-function of the grid. Do not use any member function while having called begin_mark()
template <typename TDiscreteFunction>
class VTKOutput{
	public:
		typedef TDiscreteFunction discrete_function_type;

	public:
		bool print(discrete_function_type& u, const char* filename, double Time = 0.0);
		bool print_subset(discrete_function_type& u, int subsetIndex, const char* filename, double Time = 0.0);

	private:
		bool write_prolog(FILE* file, double Time);
		bool write_piece_prolog(FILE* file);
		bool write_subset(FILE* File, discrete_function_type& u,int subsetIndex, int dim);
		bool init_subset(discrete_function_type& u, int subsetIndex, int dim);
		bool write_points(FILE* File, discrete_function_type& u, int subsetIndex, int dim);
		bool write_elements(FILE* File,discrete_function_type& u, int subsetIndex, int dim);
		bool write_scalar(FILE* File, discrete_function_type& u, uint fct, int subsetIndex, int dim);
		bool write_epilog(FILE* file);
		bool write_piece_epilog(FILE* file);

		template <typename TElem>
		void count_elem_conn(discrete_function_type& u, int subsetIndex,
								typename geometry_traits<TElem>::iterator iterBegin,
								typename geometry_traits<TElem>::iterator iterEnd);

		template <typename TElem>
		bool write_points_elementwise(	FILE* File, discrete_function_type& u,
											typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd, int& n);

		template <typename TElem>
		bool write_elements_connectivity(	FILE* File,
											typename geometry_traits<TElem>::iterator iterBegin,
											typename geometry_traits<TElem>::iterator iterEnd);

		template <typename TElem>
		bool write_elements_offsets(	FILE* File,
										typename geometry_traits<TElem>::iterator iterBegin,
										typename geometry_traits<TElem>::iterator iterEnd, int& n);

		template <typename TElem>
		bool write_elements_types(	FILE* File,
									typename geometry_traits<TElem>::iterator iterBegin,
									typename geometry_traits<TElem>::iterator iterEnd);

		template <typename TElem>
		bool write_scalar_elementwise(	FILE* File,
										discrete_function_type& u, uint fct,
										typename geometry_traits<TElem>::iterator iterBegin,
										typename geometry_traits<TElem>::iterator iterEnd);


	private:
		struct {
			int noVertices;
			int noElements;
			int noConnections;
		} Numbers;

	protected:
		typedef ug::Attachment<int> ADOFIndex;

	protected:

		Grid* m_grid;

		ADOFIndex m_aDOFIndex;
		Grid::VertexAttachmentAccessor<ADOFIndex> m_aaDOFIndexVRT;
		uint m_numberOfDOF;

};



} // namespace ug

#include "vtkoutput_impl.h"

#endif /* VTKOUTPUT_H_ */
