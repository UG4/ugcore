/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_overlap_writer
#define __H__UG_overlap_writer

#include <vector>

#include "common/error.h"
#include "lib_algebra/common/connection_viewer_output.h"
#include "lib_algebra/operator/debug_writer.h"
#include "parallelization_util.h"

namespace ug {

///	Writes overlapping matrices and vectors
/** Since local positions of overlapping matrices and vectors are not present
 * in the DebugWriter, we have to first communicate missing positions and then
 * use those positions to write the vectors and matrices.
 *
 * \sa CreateOverlap
 */
template <typename TAlgebra>
class OverlapWriter {
public:
	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;

	OverlapWriter () : m_dim (-1)	{}

	/** vector_t has to be a MathVector<dim> compatible type.
	 * \param nonOverlapVecSize		Size of the underlying non-overlapping vector.
	 * \param nonOverlapPositions	Array of length 'nonOverlapVecSize' containing
	 *								the positions of each vector entry in the underlying
	 *								non-overlapping vector.
	 * \param overlapVecSize		Size of the overlapping vector.
	 */
	template <typename vector_t>
	void init(const AlgebraLayouts& layouts,
	          size_t nonOverlapVecSize,
	          vector_t* nonOverlapPositions,
	          size_t overlapVecSize)
	{
		using namespace std;

		UG_COND_THROW(nonOverlapVecSize > overlapVecSize,
		              "nonOverlapVecSize > overlapVecSize. "
		              "This should not be the case!");

		m_dim = vector_t::Size;
		vector<MathVector<vector_t::Size> >& pos = get_pos (Int2Type<vector_t::Size>());

		pos.resize(overlapVecSize);
		for(size_t i = 0; i < nonOverlapVecSize; ++i){
			pos[i] = nonOverlapPositions[i];
		}

	//	copy positions from slave-overlap to master-overlap
		CopyValues(&pos, layouts.slave_overlap(), layouts.master_overlap(),
		           &layouts.comm());
	}

	/**
	 * \param dbgWriter			Used to extract position data. The size of the
	 *							position data is assumed to be the size of the
	 *							underlying non-overlapping vector.
	 * \param overlapVecSize	Size of the overlapping vector.
	 */
	inline
	void init(const AlgebraLayouts& layouts,
	          IVectorDebugWriter<vector_type>& dbgWriter,
	          size_t overlapVecSize)
	{
		m_dim = dbgWriter.get_dim();

		switch(m_dim) {
			case 1: init(layouts, dbgWriter.template get_positions<1>().size(),
			             &dbgWriter.template get_positions<1>().front(),
			             overlapVecSize);
					break;
			case 2: init(layouts, dbgWriter.template get_positions<2>().size(),
			             &dbgWriter.template get_positions<2>().front(),
			             overlapVecSize);
					break;
			case 3: init(layouts, dbgWriter.template get_positions<3>().size(),
			             &dbgWriter.template get_positions<3>().front(),
			             overlapVecSize);
					break;
			default: UG_THROW("Unsupported dimension: " << m_dim); break;
		}
	}


	/**	Writes a matrix or a vector to a connection viewer file.
	 * T has to be either of type TAlgebra::vector_type or TAlgebra::matrix_type.*/
	template <typename T>
	void write(const T& t, std::string name)
	{
		UG_COND_THROW (m_dim == -1, "Call 'OverlapWriter::init before calling "
		               "OverlapWriter::write");

		switch(m_dim) {
			case 1: write_dim_<1>(t, name); break;
			case 2: write_dim_<2>(t, name); break;
			case 3: write_dim_<3>(t, name); break;
		}
	}


private:
	template <int dim>
	void write_dim_(const vector_type& v, std::string name)
	{
		std::vector<MathVector<dim> >& pos = get_pos(Int2Type<dim>());
		ConnectionViewer::WriteVectorPar(name + ".vec", v, &pos.front(), m_dim);
	}

	template <int dim>
	void write_dim_(const matrix_type& A, std::string name)
	{
		std::vector<MathVector<dim> >& pos = get_pos(Int2Type<dim>());
		ConnectionViewer::WriteMatrixPar(name + ".mat", A, &pos.front(), m_dim);
	}

	std::vector<MathVector<1> >& get_pos(Int2Type<1>) {return m_pos1d;}
	std::vector<MathVector<2> >& get_pos(Int2Type<2>) {return m_pos2d;}
	std::vector<MathVector<3> >& get_pos(Int2Type<3>) {return m_pos3d;}

	//	MEMBER VARIABLES
	int m_dim;

	std::vector<MathVector<1> > m_pos1d;
	std::vector<MathVector<2> > m_pos2d;
	std::vector<MathVector<3> > m_pos3d;
};

}//	end of namespace

#endif