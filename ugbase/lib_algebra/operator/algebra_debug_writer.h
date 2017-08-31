/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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
									" Number of positions does not match: "<<   mat.num_rows() << ">" << posvec.size()<<"\n");
			// write to connection viewer
			ConnectionViewer::WriteMatrixPar<matrix_type, MathVector<dim> >
								( name, mat, &base_type::template get_positions<dim>()[0], dim);

		}

};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__ALGEBRA_DEBUG_WRITER__ */
