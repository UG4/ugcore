/*
 * vtkoutput.h
 *
 *  Created on: 06.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__IO__VTKOUTPUT__
#define __H__LIB_DISCRETIZATION__IO__VTKOUTPUT__

// extern libraries
#include <vector>

// other ug modules
#include "lib_grid/lg_base.h"

namespace ug{


/// ATTENTION: This class uses heavily the mark-function of the grid. Do not use any member function while having called begin_mark()
template <typename TDiscreteFunction>
class VTKOutput{
	public:
		typedef TDiscreteFunction discrete_function_type;

	public:
		bool begin_timeseries(const char* filename, discrete_function_type& u);
		bool end_timeseries(const char* filename, discrete_function_type& u);

		bool print(const char* filename, discrete_function_type& u, size_t step = 0, number time = 0.0);
		bool print_subset(const char* filename, discrete_function_type& u, int si, size_t step = 0, number time = 0.0);

	private:
		bool write_prolog(FILE* file, double Time);
		bool write_piece_prolog(FILE* file);
		bool write_subset(FILE* File, discrete_function_type& u,int si, int dim);
		bool init_subset(discrete_function_type& u, int si, int dim);
		bool write_points(FILE* File, discrete_function_type& u, int si, int dim);
		bool write_elements(FILE* File,discrete_function_type& u, int si, int dim);
		bool write_scalar(FILE* File, discrete_function_type& u, uint fct, int si, int dim);
		bool write_epilog(FILE* file);
		bool write_piece_epilog(FILE* file);

		bool write_pvtu(discrete_function_type& u, const char* filename, int si, size_t step, number time);
		bool write_time_pvd(discrete_function_type& u, const char* filename);

		bool write_pvd(discrete_function_type& u, const char* filename, size_t timestep = 0, number time = 0.0);

		bool vtu_filename(char *nameOut, const char *nameIn, int rank, int si, size_t step);
		bool pvtu_filename(char *nameOut, const char *nameIn, int si, size_t step);
		bool pvd_filename(char *nameOut, const char *nameIn);
		bool pvd_time_filename(char *nameOut, const char *nameIn, size_t timestep);
		bool is_valid_filename(const char *nameIn);

		template <typename TElem>
		void count_elem_conn(discrete_function_type& u, int si,
								typename geometry_traits<TElem>::const_iterator iterBegin,
								typename geometry_traits<TElem>::const_iterator iterEnd);

		template <typename TElem>
		bool write_points_elementwise(	FILE* File, discrete_function_type& u,
											typename geometry_traits<TElem>::const_iterator iterBegin,
											typename geometry_traits<TElem>::const_iterator iterEnd, int& n);

		template <typename TElem>
		bool write_elements_connectivity(	FILE* File,
											typename geometry_traits<TElem>::const_iterator iterBegin,
											typename geometry_traits<TElem>::const_iterator iterEnd);

		template <typename TElem>
		bool write_elements_offsets(	FILE* File,
										typename geometry_traits<TElem>::const_iterator iterBegin,
										typename geometry_traits<TElem>::const_iterator iterEnd, int& n);

		template <typename TElem>
		bool write_elements_types(	FILE* File,
									typename geometry_traits<TElem>::const_iterator iterBegin,
									typename geometry_traits<TElem>::const_iterator iterEnd);

		template <typename TElem>
		bool write_scalar_elementwise(	FILE* File,
										discrete_function_type& u, uint fct,
										typename geometry_traits<TElem>::const_iterator iterBegin,
										typename geometry_traits<TElem>::const_iterator iterEnd, int si);


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

	protected:
		static const size_t NAMESIZE = 64*64;

		char m_seriesname[NAMESIZE];
		discrete_function_type* m_u;
		std::vector<number> m_timestep;
};



} // namespace ug

#include "vtkoutput_impl.h"

#endif /* __H__LIB_DISCRETIZATION__IO__VTKOUTPUT__ */
