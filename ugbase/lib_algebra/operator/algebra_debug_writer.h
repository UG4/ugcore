
#ifndef __H__LIB_ALGEBRA__OPERATOR__ALGEBRA_DEBUG_WRITER__
#define __H__LIB_ALGEBRA__OPERATOR__ALGEBRA_DEBUG_WRITER__

#include "debug_writer.h"
#include "lib_algebra/common/connection_viewer_output.h"
#include "common/util/file_util.h"

namespace ug{


/// Debug writer for connection viewer (based on algebraic information + vector positions only)
template <typename TAlgebra>
class AlgebraDebugWriter : public IDebugWriter<TAlgebra>
{
	public:
	///	type of matrix
		typedef TAlgebra algebra_type;

	///	type of vector
		typedef typename algebra_type::vector_type vector_type;

	///	type of matrix
		typedef typename algebra_type::matrix_type matrix_type;

	/// type of base
		typedef IDebugWriter<TAlgebra> base_type;
		using base_type::get_base_dir;

	public:
	///	Constructor
		AlgebraDebugWriter() : base_type() {}

	///	write vector
		virtual void write_vector(const vector_type& vec,
		                          const char* filename)
		{
			switch (base_type::current_dimension())
			{
				case 1: write_vector_dim<1>(vec, filename); break;
				case 2: write_vector_dim<2>(vec, filename); break;
				case 3: write_vector_dim<3>(vec, filename); break;
				default: UG_ASSERT(0, "Dimension not implemented.");
			}
		}

	///	write matrix
		virtual void write_matrix(const matrix_type& mat,
		                          const char* filename)
		{
		//	check name
			if( !FileTypeIs( filename, ".mat" ) ) {
				UG_THROW( "Only '.mat' format supported for matrices, but"
				          " filename is '" << filename << "'.");
			}
		// write to file
			switch (base_type::current_dimension())
			{
				case 1: write_matrix_dim<1>(mat, filename); break;
				case 2: write_matrix_dim<2>(mat, filename); break;
				case 3: write_matrix_dim<3>(mat, filename); break;
				default: UG_ASSERT(0, "Dimension not implemented.");
			}

		}
	private:
		/// auxiliary function for vectors
		template <int dim>
		void write_vector_dim(const vector_type& vec, const char* filename)
		{
			std::string name = get_base_dir() + "/" + filename;
			const std::vector<MathVector<dim> > &posvec = base_type::template get_positions<dim>();
			// check size
			if(vec.size() > posvec.size())
				UG_THROW("'AlgebraDebugWriter::write_vector':"
						 " Number of positions does not match.\n");

			// write connection viewer output to file
			ConnectionViewer::WriteVectorPar<vector_type, MathVector<dim> >(name, vec, &posvec[0], dim);

		}

		/// auxiliary function for matrices
		template <int dim>
		void write_matrix_dim(const matrix_type& mat,
				                          const char* filename)
		{
			std::string name = get_base_dir() + "/" + filename;
			const std::vector<MathVector<dim> > &posvec = base_type::template get_positions<dim>();
		// check size
			if(mat.num_rows() > posvec.size() || mat.num_cols() > posvec.size())
							UG_THROW("'AlgebraDebugWriter::write_matrix':"
									" Number of positions does not match.\n");
			// write to connection viewer
			ConnectionViewer::WriteMatrixPar<matrix_type, MathVector<dim> >
								( name, mat, &base_type::template get_positions<dim>()[0], dim);

		}

};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__ALGEBRA_DEBUG_WRITER__ */
